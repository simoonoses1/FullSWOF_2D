from __future__ import annotations

import contextlib
import csv
import io
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
import numpy as np
import rasterio
from affine import Affine
from rasterio import features
from rasterio.enums import Resampling
from rasterio.fill import fillnodata
from rasterio.transform import rowcol
from rasterio.windows import Window
from rasterio.windows import transform as window_transform

from .infiltration import InfiltrationModel
from .solver import OilSpillSolver, SolverConfig


DEFAULT_MAX_DEM_CELLS = 25_000_000
MIB = 1024 * 1024


@dataclass
class DemInfo:
    path: str
    width: int
    height: int
    dx: float
    dy: float
    source_width: int
    source_height: int
    source_dx: float
    source_dy: float
    target_resolution_m: float | None
    resampling_method: str | None
    interpolated_dem_path: str | None
    crs: str | None
    nodata: float | None
    z_min: float
    z_max: float
    invalid_total: int
    invalid_remaining: int


@dataclass
class SourceSite:
    name: str
    safe_name: str
    x: float
    y: float
    row: int
    col: int
    in_dem: bool
    layer: str | None = None
    fid: int | str | None = None


@dataclass
class SourceCurve:
    mode: str
    times: list[float]
    flows: list[float]
    volume_total_m3: float | None = None
    duration_s: float | None = None
    csv_path: str | None = None

    @classmethod
    def constant(cls, volume_total_m3: float, duration_s: float) -> "SourceCurve":
        if volume_total_m3 < 0.0:
            raise ValueError("volume_total_m3 must be >= 0")
        if duration_s <= 0.0:
            raise ValueError("duration_s must be > 0")
        q = volume_total_m3 / duration_s
        return cls(
            mode="constant",
            times=[0.0, float(duration_s)],
            flows=[float(q), float(q)],
            volume_total_m3=float(volume_total_m3),
            duration_s=float(duration_s),
        )

    @classmethod
    def from_csv(cls, path: str | Path) -> "SourceCurve":
        csv_path = Path(path)
        with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames:
                raise ValueError("source CSV has no header")
            fields = {name.strip().lower(): name for name in reader.fieldnames}
            time_col = fields.get("time_s")
            flow_col = fields.get("flow_m3_s")
            volume_col = fields.get("volume_m3")
            if time_col is None or (flow_col is None and volume_col is None):
                raise ValueError("source CSV must contain time_s and flow_m3_s or volume_m3")

            times: list[float] = []
            values: list[float] = []
            value_col = flow_col or volume_col
            assert value_col is not None
            for row in reader:
                times.append(float(row[time_col]))
                values.append(float(row[value_col]))

        _validate_increasing_times(times)
        if flow_col is not None:
            if any(v < 0.0 for v in values):
                raise ValueError("flow_m3_s values must be >= 0")
            return cls(mode="csv_flow", times=times, flows=values, csv_path=str(csv_path))

        if any(v < 0.0 for v in values):
            raise ValueError("volume_m3 values must be >= 0")
        if any(values[i + 1] < values[i] for i in range(len(values) - 1)):
            raise ValueError("volume_m3 must be cumulative and non-decreasing")
        interval_flows = [
            (values[i + 1] - values[i]) / (times[i + 1] - times[i])
            for i in range(len(times) - 1)
        ]
        return cls(
            mode="csv_volume",
            times=times,
            flows=interval_flows,
            volume_total_m3=values[-1] if values else 0.0,
            csv_path=str(csv_path),
        )

    def flow_at(self, t: float) -> float:
        t = float(t)
        if not self.times:
            return 0.0
        if self.mode == "constant":
            duration = self.duration_s if self.duration_s is not None else self.times[-1]
            return max(float(self.flows[0]), 0.0) if 0.0 <= t <= duration else 0.0
        if self.mode == "csv_flow":
            if t < self.times[0] or t > self.times[-1]:
                return 0.0
            return max(float(np.interp(t, self.times, self.flows)), 0.0)
        if self.mode == "csv_volume":
            if len(self.times) < 2:
                return 0.0
            for i, q in enumerate(self.flows):
                if self.times[i] <= t < self.times[i + 1]:
                    return max(float(q), 0.0)
            return 0.0
        raise ValueError(f"unsupported source curve mode: {self.mode}")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SourceCurve":
        return cls(
            mode=str(data["mode"]),
            times=[float(v) for v in data.get("times", [])],
            flows=[float(v) for v in data.get("flows", [])],
            volume_total_m3=_optional_float(data.get("volume_total_m3")),
            duration_s=_optional_float(data.get("duration_s")),
            csv_path=data.get("csv_path"),
        )


def _optional_float(value: Any) -> float | None:
    return None if value is None or value == "" else float(value)


def _validate_increasing_times(times: list[float]) -> None:
    if len(times) < 2:
        raise ValueError("at least two time rows are required")
    if any(t < 0.0 for t in times):
        raise ValueError("time_s values must be >= 0")
    if any(times[i + 1] <= times[i] for i in range(len(times) - 1)):
        raise ValueError("time_s values must be strictly increasing")


def source_curve_from_config(config: dict[str, Any]) -> SourceCurve:
    mode = config.get("mode", "constant")
    if mode == "constant":
        return SourceCurve.constant(
            float(config.get("volume_total_m3", 0.0)),
            float(config.get("duration_s", 1.0)),
        )
    if mode == "csv":
        return SourceCurve.from_csv(str(config["csv_path"]))
    if mode in {"csv_flow", "csv_volume"}:
        return SourceCurve.from_dict(config)
    raise ValueError(f"unsupported source curve mode: {mode}")


def _optional_positive_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    number = float(value)
    if number <= 0.0:
        raise ValueError("target_resolution_m must be > 0")
    return number


def _resampling_from_name(name: str | None) -> tuple[str, Resampling]:
    normalized = (name or "bilinear").strip().lower()
    methods = {
        "bilinear": Resampling.bilinear,
        "cubic": Resampling.cubic,
    }
    if normalized not in methods:
        raise ValueError("resampling_method must be 'bilinear' or 'cubic'")
    return normalized, methods[normalized]


def interpolated_dem_path(path: str | Path, target_resolution_m: float, resampling_method: str) -> Path:
    source = Path(path)
    method_name, _ = _resampling_from_name(resampling_method)
    resolution = f"{float(target_resolution_m):g}".replace(".", "p")
    return source.with_name(f"{source.stem}_interp_{resolution}m_{method_name}.tif")


def decimated_dem_path(path: str | Path, decimate: int) -> Path:
    source = Path(path)
    return source.with_name(f"{source.stem}_decim_{int(decimate)}x_average.tif")


def _source_cell_size(src: rasterio.io.DatasetReader) -> tuple[float, float]:
    return abs(float(src.transform.a)), abs(float(src.transform.e))


def _target_grid(
    src: rasterio.io.DatasetReader,
    decimate: int = 1,
    target_resolution_m: float | None = None,
) -> tuple[int, int, Any]:
    target = _optional_positive_float(target_resolution_m)
    if target is not None:
        src_dx, src_dy = _source_cell_size(src)
        width_m = src.width * src_dx
        height_m = src.height * src_dy
        out_cols = max(1, int(math.ceil(width_m / target)))
        out_rows = max(1, int(math.ceil(height_m / target)))
    elif decimate > 1:
        out_rows = max(1, src.height // decimate)
        out_cols = max(1, src.width // decimate)
    else:
        out_rows = src.height
        out_cols = src.width
    transform = src.transform * src.transform.scale(src.width / out_cols, src.height / out_rows)
    return out_rows, out_cols, transform


def _check_grid_size(rows: int, cols: int, max_cells: int = DEFAULT_MAX_DEM_CELLS) -> None:
    cells = int(rows) * int(cols)
    if cells > max_cells:
        raise ValueError(
            f"DEM grid would have {cells:,} cells after resampling ({cols}x{rows}), "
            f"above the safe limit of {max_cells:,}. Use a coarser target resolution, "
            "increase Decimate, or crop the DEM before running."
        )


def ensure_interpolated_dem(
    path: str | Path,
    target_resolution_m: float,
    resampling_method: str = "bilinear",
    max_cells: int = DEFAULT_MAX_DEM_CELLS,
) -> Path:
    target = _optional_positive_float(target_resolution_m)
    assert target is not None
    method_name, _ = _resampling_from_name(resampling_method)
    source_path = Path(path)
    out_path = interpolated_dem_path(source_path, target, method_name)
    if not shutil.which("gdalwarp"):
        raise RuntimeError("gdalwarp was not found in PATH; DEM interpolation to disk requires GDAL")

    with rasterio.open(source_path) as src:
        out_rows, out_cols, _ = _target_grid(src, target_resolution_m=target)
        _check_grid_size(out_rows, out_cols, max_cells=max_cells)

    if out_path.exists() and out_path.stat().st_mtime >= source_path.stat().st_mtime:
        return out_path

    tmp_path = out_path.with_name(f"{out_path.stem}.tmp{out_path.suffix}")
    if tmp_path.exists():
        tmp_path.unlink()
    cmd = [
        "gdalwarp",
        "-overwrite",
        "-multi",
        "-wm",
        "512",
        "-r",
        method_name,
        "-ts",
        str(out_cols),
        str(out_rows),
        "-of",
        "GTiff",
        "-co",
        "TILED=YES",
        "-co",
        "COMPRESS=DEFLATE",
        "-co",
        "BIGTIFF=IF_SAFER",
        str(source_path),
        str(tmp_path),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    if out_path.exists():
        out_path.unlink()
    tmp_path.replace(out_path)
    return out_path


def ensure_decimated_dem(
    path: str | Path,
    decimate: int,
    max_cells: int = DEFAULT_MAX_DEM_CELLS,
) -> Path:
    decimate = int(decimate)
    if decimate <= 1:
        return Path(path)
    source_path = Path(path)
    out_path = decimated_dem_path(source_path, decimate)
    if not shutil.which("gdalwarp"):
        raise RuntimeError("gdalwarp was not found in PATH; DEM decimation to disk requires GDAL")

    with rasterio.open(source_path) as src:
        out_rows, out_cols, _ = _target_grid(src, decimate=decimate)
        _check_grid_size(out_rows, out_cols, max_cells=max_cells)

    if out_path.exists() and out_path.stat().st_mtime >= source_path.stat().st_mtime:
        return out_path

    tmp_path = out_path.with_name(f"{out_path.stem}.tmp{out_path.suffix}")
    if tmp_path.exists():
        tmp_path.unlink()
    cmd = [
        "gdalwarp",
        "-overwrite",
        "-multi",
        "-wm",
        "512",
        "-r",
        "average",
        "-ts",
        str(out_cols),
        str(out_rows),
        "-of",
        "GTiff",
        "-co",
        "TILED=YES",
        "-co",
        "COMPRESS=DEFLATE",
        "-co",
        "BIGTIFF=IF_SAFER",
        str(source_path),
        str(tmp_path),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    if out_path.exists():
        out_path.unlink()
    tmp_path.replace(out_path)
    return out_path


def ensure_prepared_dem_file(
    path: str | Path,
    decimate: int = 1,
    target_resolution_m: float | None = None,
    resampling_method: str = "bilinear",
    max_cells: int = DEFAULT_MAX_DEM_CELLS,
) -> Path:
    target = _optional_positive_float(target_resolution_m)
    if target is not None:
        return ensure_interpolated_dem(path, target, resampling_method, max_cells=max_cells)
    if int(decimate) > 1:
        return ensure_decimated_dem(path, int(decimate), max_cells=max_cells)
    return Path(path)


def prepare_dem(
    path: str | Path,
    decimate: int = 1,
    target_resolution_m: float | None = None,
    resampling_method: str = "bilinear",
    max_cells: int = DEFAULT_MAX_DEM_CELLS,
) -> tuple[np.ndarray, dict[str, Any], float, float]:
    target = _optional_positive_float(target_resolution_m)
    effective_path = ensure_prepared_dem_file(
        path,
        decimate=decimate,
        target_resolution_m=target,
        resampling_method=resampling_method,
        max_cells=max_cells,
    )

    with rasterio.open(effective_path) as src:
        if target is not None or decimate > 1:
            dem = src.read(1)
            profile = src.profile
            transform = src.transform
        elif decimate <= 1:
            _check_grid_size(src.height, src.width, max_cells=max_cells)
            dem = src.read(1)
            profile = src.profile
            transform = src.transform

    dx = abs(transform.a)
    dy = abs(transform.e)
    return dem.astype(float), profile, dx, dy


def prepare_site_dem_window(
    path: str | Path,
    source_row: int,
    source_col: int,
    buffer_m: float,
    decimate: int = 1,
    target_resolution_m: float | None = None,
    resampling_method: str = "bilinear",
    max_cells: int = DEFAULT_MAX_DEM_CELLS,
) -> tuple[np.ndarray, dict[str, Any], float, float, dict[str, Any]]:
    effective_path = ensure_prepared_dem_file(
        path,
        decimate=decimate,
        target_resolution_m=target_resolution_m,
        resampling_method=resampling_method,
        max_cells=max_cells,
    )
    with rasterio.open(effective_path) as src:
        dx, dy = _source_cell_size(src)
        buffer_cells_y = max(1, int(math.ceil(float(buffer_m) / dy)))
        buffer_cells_x = max(1, int(math.ceil(float(buffer_m) / dx)))
        row = int(source_row)
        col = int(source_col)
        if not (0 <= row < src.height and 0 <= col < src.width):
            raise ValueError(f"source row/col ({row}, {col}) is outside prepared DEM grid {src.height}x{src.width}")

        row_start = max(0, row - buffer_cells_y)
        row_stop = min(src.height, row + buffer_cells_y + 1)
        col_start = max(0, col - buffer_cells_x)
        col_stop = min(src.width, col + buffer_cells_x + 1)
        crop_rows = row_stop - row_start
        crop_cols = col_stop - col_start
        _check_grid_size(crop_rows, crop_cols, max_cells=max_cells)

        window = Window(col_start, row_start, crop_cols, crop_rows)
        dem = src.read(1, window=window)
        transform = window_transform(window, src.transform)
        profile = src.profile.copy()
        profile.update(height=crop_rows, width=crop_cols, transform=transform)

    crop_info = {
        "enabled": True,
        "prepared_dem_path": str(effective_path),
        "buffer_m": float(buffer_m),
        "global_row": row,
        "global_col": col,
        "local_row": row - row_start,
        "local_col": col - col_start,
        "row_start": row_start,
        "row_stop": row_stop,
        "col_start": col_start,
        "col_stop": col_stop,
        "height": crop_rows,
        "width": crop_cols,
    }
    return dem.astype(float), profile, dx, dy, crop_info


def fill_nodata_local_average(arr: np.ndarray, invalid: np.ndarray, max_iters: int = 64) -> tuple[np.ndarray, np.ndarray]:
    work = arr.copy()
    rem = invalid.copy()
    for _ in range(max_iters):
        if not rem.any():
            break
        sums = np.zeros_like(work, dtype=float)
        count = np.zeros_like(work, dtype=int)
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]:
            src = work[max(0, dr): work.shape[0] + min(0, dr), max(0, dc): work.shape[1] + min(0, dc)]
            src_mask = rem[max(0, dr): work.shape[0] + min(0, dr), max(0, dc): work.shape[1] + min(0, dc)]
            dst_r0, dst_r1 = max(0, -dr), work.shape[0] - max(0, dr)
            dst_c0, dst_c1 = max(0, -dc), work.shape[1] - max(0, dc)
            valid = ~src_mask
            sums[dst_r0:dst_r1, dst_c0:dst_c1] += src * valid
            count[dst_r0:dst_r1, dst_c0:dst_c1] += valid.astype(int)
        fillable = rem & (count > 0)
        if not fillable.any():
            break
        work[fillable] = sums[fillable] / count[fillable]
        rem[fillable] = False
    return work, rem


def validate_dem(dem: np.ndarray, nodata: float | None) -> tuple[np.ndarray, dict[str, Any]]:
    work = dem.astype(float, copy=True)
    invalid = ~np.isfinite(work)
    if nodata is not None:
        invalid |= work == nodata
    invalid_total = int(invalid.sum())

    fillnodata_remaining = invalid.copy()
    if invalid_total and np.any(~invalid):
        try:
            filled = fillnodata(work, mask=(~invalid).astype("uint8"), max_search_distance=1000.0, smoothing_iterations=2)
            work = np.asarray(filled, dtype=float)
            fillnodata_remaining = ~np.isfinite(work)
            if nodata is not None:
                fillnodata_remaining |= work == nodata
        except Exception:
            fillnodata_remaining = invalid.copy()

    local_iters = max(64, min(max(work.shape) if work.ndim == 2 else 64, 512))
    clean, remaining = fill_nodata_local_average(work, fillnodata_remaining, max_iters=local_iters)

    fallback_filled = int(remaining.sum())
    if fallback_filled:
        valid = np.isfinite(clean) & ~remaining
        if not np.any(valid):
            raise ValueError("DEM contains no valid cells")
        fill_value = float(np.nanmedian(clean[valid]))
        clean[remaining] = fill_value
        remaining = ~np.isfinite(clean)

    stats = {
        "invalid_total": invalid_total,
        "invalid_remaining": int(remaining.sum()),
        "invalid_filled_with_fallback": fallback_filled,
        "z_min": float(np.nanmin(clean)),
        "z_max": float(np.nanmax(clean)),
    }
    return clean, stats


def inspect_dem(
    path: str | Path,
    decimate: int = 1,
    target_resolution_m: float | None = None,
    resampling_method: str = "bilinear",
    progress: Callable[[int, str], None] | None = None,
    max_cells: int = DEFAULT_MAX_DEM_CELLS,
) -> dict[str, Any]:
    def emit(percent: int, message: str) -> None:
        if progress is not None:
            progress(percent, message)

    emit(5, "Abriendo DEM y preparando grilla")
    target = _optional_positive_float(target_resolution_m)
    method_name = _resampling_from_name(resampling_method)[0] if target is not None else None
    with rasterio.open(path) as src:
        grid_rows, grid_cols, _ = _target_grid(src, decimate=decimate, target_resolution_m=target)
        _check_grid_size(grid_rows, grid_cols, max_cells=max_cells)
    emit(15, "Leyendo DEM" if target is None else f"Guardando DEM interpolado con {method_name}")
    dem, profile, dx, dy = prepare_dem(path, decimate, target, resampling_method, max_cells=max_cells)
    emit(60, "Validando nodata y calculando estadísticas")
    _, stats = validate_dem(dem, profile.get("nodata"))
    emit(85, "Leyendo metadatos del DEM original")
    with rasterio.open(path) as src:
        source_dx, source_dy = _source_cell_size(src)
        source_width = int(src.width)
        source_height = int(src.height)
    info = DemInfo(
        path=str(path),
        width=int(profile["width"]),
        height=int(profile["height"]),
        dx=float(dx),
        dy=float(dy),
        source_width=source_width,
        source_height=source_height,
        source_dx=float(source_dx),
        source_dy=float(source_dy),
        target_resolution_m=target,
        resampling_method=method_name,
        interpolated_dem_path=str(interpolated_dem_path(path, target, method_name)) if target is not None else None,
        crs=str(profile.get("crs")) if profile.get("crs") else None,
        nodata=_optional_float(profile.get("nodata")),
        z_min=stats["z_min"],
        z_max=stats["z_max"],
        invalid_total=stats["invalid_total"],
        invalid_remaining=stats["invalid_remaining"],
    )
    emit(100, "DEM listo")
    return asdict(info)


def sanitize_site_name(name: Any) -> str:
    text = str(name).strip() if name is not None else "site"
    text = re.sub(r'[<>:"/\\|?*\x00-\x1f]+', "_", text)
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"_+", "_", text)
    text = text.strip("._ ")
    return (text or "site")[:80]


def point_to_rowcol(x: float, y: float, transform: Any, height: int, width: int) -> tuple[int, int, bool]:
    row, col = rowcol(transform, x, y)
    row_i = int(row)
    col_i = int(col)
    return row_i, col_i, 0 <= row_i < height and 0 <= col_i < width


def ogr_available() -> bool:
    return bool(shutil.which("ogrinfo") and shutil.which("ogr2ogr"))


def list_vector_layers(vector_path: str | Path) -> list[dict[str, Any]]:
    if not shutil.which("ogrinfo"):
        raise RuntimeError("ogrinfo was not found in PATH")
    result = subprocess.run(
        ["ogrinfo", "-json", "-ro", "-summary", str(vector_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    data = json.loads(result.stdout)
    layers = []
    for layer in data.get("layers", []):
        fields_raw = layer.get("fields", [])
        fields = []
        for field in fields_raw:
            if isinstance(field, dict):
                fields.append(field.get("name"))
            else:
                fields.append(str(field))
        geom_types = []
        for geom_field in layer.get("geometryFields", []):
            if isinstance(geom_field, dict):
                geom_types.append(geom_field.get("type"))
        layers.append(
            {
                "name": layer.get("name"),
                "feature_count": layer.get("featureCount"),
                "fields": [f for f in fields if f],
                "geometry_types": [g for g in geom_types if g],
            }
        )
    return layers


def choose_default_site_layer(layers: list[dict[str, Any]], name_field: str = "Name") -> str | None:
    for layer in layers:
        fields = {str(f).lower() for f in layer.get("fields", [])}
        geom_types = " ".join(str(g).lower() for g in layer.get("geometry_types", []))
        if name_field.lower() in fields and "point" in geom_types:
            return layer.get("name")
    for layer in layers:
        fields = {str(f).lower() for f in layer.get("fields", [])}
        if name_field.lower() in fields:
            return layer.get("name")
    return layers[0].get("name") if layers else None


def load_source_sites(
    vector_path: str | Path,
    dem_path: str | Path,
    layer_name: str | None = None,
    name_field: str = "Name",
    decimate: int = 1,
    target_resolution_m: float | None = None,
    resampling_method: str = "bilinear",
) -> list[dict[str, Any]]:
    if not shutil.which("ogr2ogr"):
        raise RuntimeError("ogr2ogr was not found in PATH")
    if target_resolution_m is not None:
        _resampling_from_name(resampling_method)
    with rasterio.open(dem_path) as src:
        crs = src.crs
        height, width, transform = _target_grid(src, decimate=decimate, target_resolution_m=target_resolution_m)

    if layer_name is None:
        layer_name = choose_default_site_layer(list_vector_layers(vector_path), name_field=name_field)
    if not layer_name:
        raise ValueError("no vector layer was found")

    with tempfile.TemporaryDirectory(prefix="fullswof_sites_") as tmp:
        geojson_path = Path(tmp) / "sites.geojson"
        cmd = ["ogr2ogr", "-f", "GeoJSON", str(geojson_path), str(vector_path), str(layer_name)]
        if crs:
            cmd[3:3] = ["-t_srs", str(crs)]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        with geojson_path.open("r", encoding="utf-8") as f:
            data = json.load(f)

    sites: list[dict[str, Any]] = []
    for idx, feature in enumerate(data.get("features", [])):
        props = feature.get("properties") or {}
        resolved_name_field = _resolve_name_field(props, name_field)
        if resolved_name_field is None:
            available = ", ".join(str(k) for k in props.keys())
            raise ValueError(
                f"field {name_field!r} was not found in vector layer {layer_name!r}. "
                f"Available fields: {available}"
            )
        geom = feature.get("geometry") or {}
        coords = _representative_coordinates(geom)
        if coords is None:
            continue
        x, y = coords
        row, col, in_dem = point_to_rowcol(float(x), float(y), transform, height, width)
        name = str(props.get(resolved_name_field))
        site = SourceSite(
            name=name,
            safe_name=sanitize_site_name(name),
            x=float(x),
            y=float(y),
            row=row,
            col=col,
            in_dem=in_dem,
            layer=layer_name,
            fid=feature.get("id", idx),
        )
        sites.append(asdict(site))
    if not sites and data.get("features"):
        raise ValueError(
            f"layer {layer_name!r} has features, but no point/line/polygon geometry could be converted to source sites"
        )
    return sites


def _resolve_name_field(properties: dict[str, Any], requested: str) -> str | None:
    if requested in properties:
        return requested
    lower = {str(key).lower(): key for key in properties.keys()}
    for candidate in (requested, "Name", "NOMBRE", "nombre", "name"):
        key = lower.get(candidate.lower())
        if key is not None:
            return str(key)
    return None


def _representative_coordinates(geometry: dict[str, Any]) -> tuple[float, float] | None:
    geom_type = str(geometry.get("type", "")).lower()
    coords = geometry.get("coordinates")
    if geom_type == "point" and coords and len(coords) >= 2:
        return float(coords[0]), float(coords[1])
    if geom_type == "multipoint":
        return _mean_xy(coords)
    if geom_type == "linestring":
        return _mean_xy(coords)
    if geom_type == "multilinestring":
        return _mean_xy(_flatten_one_level(coords))
    if geom_type == "polygon":
        return _polygon_centroid(coords[0] if coords else None)
    if geom_type == "multipolygon":
        return _multipolygon_centroid(coords)
    if geom_type == "geometrycollection":
        for child in geometry.get("geometries", []):
            child_coords = _representative_coordinates(child)
            if child_coords is not None:
                return child_coords
    return None


def _flatten_one_level(items: Any) -> list[Any]:
    if not items:
        return []
    out: list[Any] = []
    for item in items:
        if isinstance(item, list):
            out.extend(item)
    return out


def _mean_xy(points: Any) -> tuple[float, float] | None:
    if not points:
        return None
    xy = []
    for point in points:
        if point and len(point) >= 2:
            xy.append((float(point[0]), float(point[1])))
    if not xy:
        return None
    arr = np.asarray(xy, dtype=float)
    return float(np.mean(arr[:, 0])), float(np.mean(arr[:, 1]))


def _polygon_centroid(ring: Any) -> tuple[float, float] | None:
    if not ring or len(ring) < 3:
        return _mean_xy(ring)
    points = [(float(p[0]), float(p[1])) for p in ring if p and len(p) >= 2]
    if len(points) < 3:
        return _mean_xy(points)
    if points[0] != points[-1]:
        points.append(points[0])

    twice_area = 0.0
    cx = 0.0
    cy = 0.0
    for (x0, y0), (x1, y1) in zip(points[:-1], points[1:]):
        cross = x0 * y1 - x1 * y0
        twice_area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(twice_area) < 1e-12:
        return _mean_xy(points)
    return cx / (3.0 * twice_area), cy / (3.0 * twice_area)


def _multipolygon_centroid(multipolygon: Any) -> tuple[float, float] | None:
    if not multipolygon:
        return None
    weighted: list[tuple[float, float, float]] = []
    fallback_points = []
    for polygon in multipolygon:
        ring = polygon[0] if polygon else None
        centroid = _polygon_centroid(ring)
        area = abs(_polygon_twice_area(ring)) if ring else 0.0
        if centroid is not None and area > 1e-12:
            weighted.append((centroid[0], centroid[1], area))
        elif centroid is not None:
            fallback_points.append(centroid)
    if weighted:
        total = sum(item[2] for item in weighted)
        return (
            sum(x * area for x, _, area in weighted) / total,
            sum(y * area for _, y, area in weighted) / total,
        )
    return _mean_xy(fallback_points)


def _polygon_twice_area(ring: Any) -> float:
    if not ring or len(ring) < 3:
        return 0.0
    points = [(float(p[0]), float(p[1])) for p in ring if p and len(p) >= 2]
    if len(points) < 3:
        return 0.0
    if points[0] != points[-1]:
        points.append(points[0])
    return sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(points[:-1], points[1:]))


def default_model_params() -> dict[str, Any]:
    return {
        "t_end": 1200.0,
        "cfl": 0.3,
        "friction_model": "manning",
        "manning_n": 0.03,
        "darcy_f": 0.02,
        "fluid_density": 900.0,
        "dynamic_viscosity": 0.005,
        "yield_stress": 0.0,
        "h_viscous_min": 1e-4,
        "k_oil": 5e-6,
        "capillary_suction": 0.0,
        "porosity_deficit": 0.2,
        "evap_rate": 0.0,
        "deg_rate": 0.0,
        "snapshot_interval": 10.0,
        "roi_buffer_m": 500.0,
        "video_fps": 2,
        "video_3d_palette": "oil_dark",
    }


def build_solver_config(
    params: dict[str, Any],
    nx: int,
    ny: int,
    dx: float,
    dy: float,
    source_row: int,
    source_col: int,
    curve: SourceCurve,
) -> SolverConfig:
    merged = {**default_model_params(), **params}
    return SolverConfig(
        nx=nx,
        ny=ny,
        dx=dx,
        dy=dy,
        t_end=float(merged["t_end"]),
        cfl=float(merged["cfl"]),
        friction_model=str(merged["friction_model"]),
        manning_n=float(merged["manning_n"]),
        darcy_f=float(merged["darcy_f"]),
        fluid_density=float(merged["fluid_density"]),
        dynamic_viscosity=float(merged["dynamic_viscosity"]),
        yield_stress=float(merged["yield_stress"]),
        h_viscous_min=float(merged["h_viscous_min"]),
        evaporation_rate=float(merged["evap_rate"]),
        degradation_rate=float(merged["deg_rate"]),
        point_source_row=int(source_row),
        point_source_col=int(source_col),
        point_source_flow_rate=curve.flow_at,
        roi_buffer_m=float(merged["roi_buffer_m"]),
    )


def build_infiltration(params: dict[str, Any]) -> InfiltrationModel:
    merged = {**default_model_params(), **params}
    return InfiltrationModel(
        saturated_hydraulic_conductivity=float(merged["k_oil"]),
        capillary_suction=float(merged["capillary_suction"]),
        porosity_deficit=float(merged["porosity_deficit"]),
    )


def run_site_case(case: dict[str, Any]) -> dict[str, Any]:
    started = time.time()
    out_dir = Path(case["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "run.log"
    log_buffer = io.StringIO()
    verbose = bool(case.get("verbose", False))
    stream_cm = contextlib.nullcontext() if verbose else contextlib.redirect_stdout(log_buffer)
    error_cm = contextlib.nullcontext() if verbose else contextlib.redirect_stderr(log_buffer)
    try:
        with stream_cm, error_cm:
            result = _run_site_case_inner(case, out_dir)
        result["status"] = "completed"
    except Exception as exc:
        result = {
            "status": "failed",
            "error": str(exc),
            "site": case.get("site", {}),
            "out_dir": str(out_dir),
        }
    finally:
        if not verbose:
            log_path.write_text(log_buffer.getvalue(), encoding="utf-8")

    finished = time.time()
    result["started_at"] = started
    result["finished_at"] = finished
    result["elapsed_s"] = finished - started
    result["log_path"] = str(log_path)
    return result


def _run_site_case_inner(case: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    target_resolution_m = _optional_positive_float(case.get("target_resolution_m"))
    resampling_method = str(case.get("resampling_method", "bilinear"))
    site_global = dict(case["site"])
    curve = SourceCurve.from_dict(case["source_curve"])
    params = case.get("model_params", {})
    merged_params = {**default_model_params(), **params}
    domain_buffer_m = float(merged_params.get("domain_buffer_m", merged_params["roi_buffer_m"]))

    dem, profile, dx, dy, crop_info = prepare_site_dem_window(
        case["dem_path"],
        int(site_global["row"]),
        int(site_global["col"]),
        buffer_m=domain_buffer_m,
        decimate=int(case.get("decimate", 1)),
        target_resolution_m=target_resolution_m,
        resampling_method=resampling_method,
    )
    dem, dem_stats = validate_dem(dem, profile.get("nodata"))
    if dem_stats["invalid_remaining"] > 0:
        raise ValueError(f"DEM contains {dem_stats['invalid_remaining']} unresolved invalid cells after filling")

    site = {
        **site_global,
        "global_row": crop_info["global_row"],
        "global_col": crop_info["global_col"],
        "row": crop_info["local_row"],
        "col": crop_info["local_col"],
    }
    ny, nx = dem.shape
    cfg = build_solver_config(params, nx, ny, dx, dy, int(site["row"]), int(site["col"]), curve)
    infiltration = build_infiltration(params)

    solver = OilSpillSolver(cfg, np.zeros_like(dem), z=dem, infiltration=infiltration)
    outputs, summary = solver.run(
        snapshot_interval=float(params.get("snapshot_interval", default_model_params()["snapshot_interval"])),
        verbose=bool(case.get("solver_verbose", False)),
    )

    output_options = {
        "arrays": True,
        "rasters": True,
        "maps": True,
        "diagnostic_report": True,
        "videos": True,
        "isochrones": True,
        "snapshots": True,
        **case.get("output_options", {}),
    }
    if output_options["arrays"]:
        np.save(out_dir / "final_h.npy", solver.h)
        np.save(out_dir / "final_infiltration.npy", solver.cum_infiltration)
        np.save(out_dir / "max_h.npy", solver.max_h)
    max_thickness_raster_path = None
    if output_options["rasters"]:
        max_thickness_raster_path = export_single_band_raster(
            solver.max_h,
            profile,
            out_dir / "max_thickness_m.tif",
        )
    snapshots_path = None
    if output_options["snapshots"]:
        snapshots_path = save_video_reexport_context(outputs, dem, profile, out_dir)
    snapshot_stats = _snapshot_storage_stats(out_dir, snapshots_path, outputs)

    domain_boundary_reached = wet_touches_boundary(solver.h)
    metadata = {
        **summary,
        **snapshot_stats,
        "dem_stats": dem_stats,
        "site": site,
        "dem_path": case["dem_path"],
        "vector_path": case.get("vector_path"),
        "source_curve": curve.to_dict(),
        "model_params": params,
        "output_options": output_options,
        "max_thickness_raster_path": max_thickness_raster_path,
        "snapshots_path": snapshots_path,
        "dem_resampling": {
            "decimate": int(case.get("decimate", 1)),
            "target_resolution_m": target_resolution_m,
            "resampling_method": _resampling_from_name(resampling_method)[0] if target_resolution_m is not None else None,
            "prepared_dem_path": crop_info["prepared_dem_path"],
            "interpolated_dem_path": (
                str(interpolated_dem_path(case["dem_path"], target_resolution_m, resampling_method))
                if target_resolution_m is not None
                else None
            ),
            "dx": float(dx),
            "dy": float(dy),
            "width": int(nx),
            "height": int(ny),
        },
        "domain_crop": crop_info,
        "domain_boundary_reached": domain_boundary_reached,
    }

    if output_options["maps"]:
        _save_basic_maps(solver, out_dir)

    isochrones_path = None
    isochrones_error = None
    if output_options["isochrones"]:
        interval_s = float(params.get("isochrone_interval_s", 60.0))
        try:
            isochrones_path = export_isochrones_gpkg(
                outputs,
                profile,
                out_dir,
                interval_s=interval_s,
                site=site,
                dem=dem,
                h_max=solver.max_h,
                model_params=merged_params,
                summary=summary,
                source_curve=curve.to_dict(),
            )
        except Exception as exc:
            isochrones_error = str(exc)

    if output_options["diagnostic_report"] or output_options["videos"]:
        from Examples.run_dem_methodology_test import export_mp4, export_mp4_3d, generate_diagnostic_report

        report_sim = type("ReportSim", (), {"profile": profile, "V_clipped_accum": 0.0, "z": dem, "cfg": solver.cfg})()
        if output_options["diagnostic_report"]:
            generate_diagnostic_report(outputs, report_sim, curve.flow_at, str(out_dir))
        if output_options["videos"]:
            fps = int(params.get("video_fps", default_model_params()["video_fps"]))
            video_3d_palette = str(params.get("video_3d_palette", default_model_params()["video_3d_palette"]))
            export_mp4(outputs, dem, profile, str(out_dir / "spill_2d.mp4"), fps=fps)
            export_mp4_3d(outputs, dem, profile, str(out_dir / "spill_3d.mp4"), fps=fps, spill_palette=video_3d_palette)

    metadata["isochrones_path"] = isochrones_path
    metadata["isochrones_error"] = isochrones_error
    with (out_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    return {
        "site": site,
        "out_dir": str(out_dir),
        "summary_path": str(out_dir / "summary.json"),
        "max_thickness_raster_path": max_thickness_raster_path,
        "snapshots_path": snapshots_path,
        **snapshot_stats,
        "isochrones_path": isochrones_path,
        "isochrones_error": isochrones_error,
        "domain_boundary_reached": domain_boundary_reached,
        "final_surface_mass": summary.get("final_surface_mass"),
        "mass_closure_error_pct": summary.get("mass_closure_error_pct"),
    }


def _save_basic_maps(solver: OilSpillSolver, out_dir: Path) -> None:
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    im0 = ax[0].imshow(solver.h, origin="lower", cmap="inferno")
    ax[0].set_title("Final oil thickness h [m]")
    fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)
    im1 = ax[1].imshow(solver.cum_infiltration, origin="lower", cmap="Blues")
    ax[1].set_title("Cumulative infiltration [m]")
    fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_dir / "diagnostic_maps.png", dpi=180)
    plt.close(fig)


def export_single_band_raster(
    array: np.ndarray,
    profile: dict[str, Any],
    path: str | Path,
    *,
    nodata: float | None = None,
) -> str:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    data = np.asarray(array, dtype=np.float32)
    raster_profile: dict[str, Any] = {
        "driver": "GTiff",
        "height": int(data.shape[0]),
        "width": int(data.shape[1]),
        "count": 1,
        "dtype": "float32",
        "transform": profile.get("transform"),
        "compress": "deflate",
        "predictor": 3,
        "zlevel": 4,
    }
    if nodata is not None:
        raster_profile["nodata"] = float(nodata)
    if profile.get("crs") is not None:
        raster_profile["crs"] = profile.get("crs")

    with rasterio.open(out_path, "w", **raster_profile) as dst:
        dst.write(data, 1)
    return str(out_path)


def wet_touches_boundary(h: np.ndarray, threshold: float = 1e-9) -> bool:
    arr = np.asarray(h)
    if arr.size == 0:
        return False
    return bool(
        np.any(arr[0, :] > threshold)
        or np.any(arr[-1, :] > threshold)
        or np.any(arr[:, 0] > threshold)
        or np.any(arr[:, -1] > threshold)
    )


def _path_size_bytes(path: str | Path | None) -> int:
    if path is None:
        return 0
    try:
        path_obj = Path(path)
        return int(path_obj.stat().st_size) if path_obj.exists() else 0
    except OSError:
        return 0


def _bytes_to_mib(value: int | float) -> float:
    return float(value) / float(MIB)


def _snapshot_uncompressed_bytes(outputs: list[dict[str, Any]]) -> int:
    h_bytes = 0
    for snap in outputs:
        h = snap.get("h")
        if h is not None:
            h_bytes += int(np.asarray(h).size) * np.dtype(np.float32).itemsize
    time_bytes = len(outputs) * np.dtype(np.float64).itemsize
    return int(h_bytes + time_bytes)


def _snapshot_storage_stats(out_dir: str | Path, snapshots_path: str | Path | None, outputs: list[dict[str, Any]]) -> dict[str, Any]:
    if snapshots_path is None:
        return {
            "snapshots_enabled": False,
            "snapshot_count": 0,
            "snapshots_archive_shape": None,
            "snapshots_size_bytes": 0,
            "snapshots_size_mb": 0.0,
            "snapshots_uncompressed_bytes": 0,
            "snapshots_uncompressed_mb": 0.0,
            "video_reexport_context_size_bytes": 0,
            "video_reexport_context_size_mb": 0.0,
        }

    out_path = Path(out_dir)
    snapshots_bytes = _path_size_bytes(snapshots_path)
    context_bytes = _path_size_bytes(out_path / "video_context.json")
    dem_bytes = _path_size_bytes(out_path / "domain_dem.npy")
    total_context_bytes = snapshots_bytes + context_bytes + dem_bytes
    uncompressed_bytes = _snapshot_uncompressed_bytes(outputs)
    first_h = next((snap.get("h") for snap in outputs if snap.get("h") is not None), None)
    if first_h is None:
        archive_shape = [len(outputs), 0, 0]
    else:
        h_shape = np.asarray(first_h).shape
        archive_shape = [len(outputs), *[int(value) for value in h_shape]]
    return {
        "snapshots_enabled": True,
        "snapshot_count": len(outputs),
        "snapshots_archive_shape": archive_shape,
        "snapshots_size_bytes": snapshots_bytes,
        "snapshots_size_mb": _bytes_to_mib(snapshots_bytes),
        "snapshots_uncompressed_bytes": uncompressed_bytes,
        "snapshots_uncompressed_mb": _bytes_to_mib(uncompressed_bytes),
        "video_reexport_context_size_bytes": total_context_bytes,
        "video_reexport_context_size_mb": _bytes_to_mib(total_context_bytes),
    }


def _profile_transform_to_json(profile: dict[str, Any]) -> list[float] | None:
    transform = profile.get("transform")
    if transform is None:
        return None
    return [float(transform.a), float(transform.b), float(transform.c), float(transform.d), float(transform.e), float(transform.f)]


def _profile_transform_from_json(values: list[float] | None) -> Any:
    if not values:
        return None
    return Affine(float(values[0]), float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]))


def save_video_reexport_context(
    outputs: list[dict[str, Any]],
    dem: np.ndarray,
    profile: dict[str, Any],
    out_dir: str | Path,
) -> str:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    snapshots_path = out_path / "snapshots.npz"
    dem_path = out_path / "domain_dem.npy"
    context_path = out_path / "video_context.json"

    times = np.asarray([float(snap.get("time", np.nan)) for snap in outputs], dtype=np.float64)
    h_stack = np.stack([np.asarray(snap["h"], dtype=np.float32) for snap in outputs], axis=0)
    np.savez_compressed(snapshots_path, time=times, h=h_stack)
    np.save(dem_path, np.asarray(dem, dtype=np.float32))

    context = {
        "snapshots_path": str(snapshots_path),
        "domain_dem_path": str(dem_path),
        "snapshot_count": int(h_stack.shape[0]),
        "snapshots_archive_shape": [int(value) for value in h_stack.shape],
        "snapshots_uncompressed_bytes": int(times.nbytes + h_stack.nbytes),
        "snapshots_size_bytes": _path_size_bytes(snapshots_path),
        "domain_dem_size_bytes": _path_size_bytes(dem_path),
        "profile": {
            "height": int(profile.get("height", dem.shape[0])),
            "width": int(profile.get("width", dem.shape[1])),
            "crs": str(profile.get("crs")) if profile.get("crs") is not None else None,
            "transform": _profile_transform_to_json(profile),
        },
    }
    context_path.write_text(json.dumps(context, indent=2), encoding="utf-8")
    context["video_context_size_bytes"] = _path_size_bytes(context_path)
    context["video_reexport_context_size_bytes"] = (
        int(context["snapshots_size_bytes"])
        + int(context["domain_dem_size_bytes"])
        + int(context["video_context_size_bytes"])
    )
    context_path.write_text(json.dumps(context, indent=2), encoding="utf-8")
    return str(snapshots_path)


def load_video_reexport_context(out_dir: str | Path, t_end_s: float | None = None) -> tuple[list[dict[str, Any]], np.ndarray, dict[str, Any]]:
    out_path = Path(out_dir)
    context_path = out_path / "video_context.json"
    if not context_path.exists():
        raise FileNotFoundError(f"video_context.json was not found in {out_path}")
    context = json.loads(context_path.read_text(encoding="utf-8"))

    snapshots_path = Path(context.get("snapshots_path") or out_path / "snapshots.npz")
    dem_path = Path(context.get("domain_dem_path") or out_path / "domain_dem.npy")
    if not snapshots_path.exists():
        raise FileNotFoundError(f"snapshots archive was not found: {snapshots_path}")
    if not dem_path.exists():
        raise FileNotFoundError(f"local DEM array was not found: {dem_path}")

    with np.load(snapshots_path) as data:
        times = np.asarray(data["time"], dtype=float)
        h_stack = np.asarray(data["h"], dtype=np.float32)

    if t_end_s is not None:
        mask = times <= float(t_end_s) + 1e-12
        if not np.any(mask):
            mask[0] = True
        times = times[mask]
        h_stack = h_stack[mask]

    outputs = [{"time": float(times[i]), "h": h_stack[i]} for i in range(len(times))]
    dem = np.load(dem_path)
    profile_data = context.get("profile", {})
    profile: dict[str, Any] = {
        "height": int(profile_data.get("height", dem.shape[0])),
        "width": int(profile_data.get("width", dem.shape[1])),
    }
    transform = _profile_transform_from_json(profile_data.get("transform"))
    if transform is not None:
        profile["transform"] = transform
    if profile_data.get("crs") is not None:
        profile["crs"] = profile_data.get("crs")
    return outputs, dem, profile


def _video_suffix(fps: int, t_end_s: float | None) -> str:
    end_label = "full" if t_end_s is None else f"to_{float(t_end_s):g}s".replace(".", "p")
    return f"fps{int(fps)}_{end_label}"


def _existing_video_or_gif(mp4_path: Path) -> str | None:
    if mp4_path.exists():
        return str(mp4_path)
    gif_path = mp4_path.with_suffix(".gif")
    if gif_path.exists():
        return str(gif_path)
    return None


def reexport_site_videos(
    out_dir: str | Path,
    fps: int = 2,
    t_end_s: float | None = None,
    include_2d: bool = True,
    include_3d: bool = True,
    video_3d_palette: str = "oil_dark",
) -> dict[str, Any]:
    fps = int(fps)
    if fps <= 0:
        raise ValueError("fps must be > 0")
    if t_end_s is not None and float(t_end_s) <= 0.0:
        raise ValueError("t_end_s must be > 0 when provided")

    outputs, dem, profile = load_video_reexport_context(out_dir, t_end_s=t_end_s)
    if not outputs:
        raise ValueError("no snapshots are available for video re-export")

    from Examples.run_dem_methodology_test import export_mp4, export_mp4_3d

    out_path = Path(out_dir)
    suffix = _video_suffix(fps, t_end_s)
    result: dict[str, Any] = {
        "status": "completed",
        "fps": fps,
        "t_end_s": t_end_s,
        "frames": len(outputs),
        "video_3d_palette": video_3d_palette,
    }
    if include_2d:
        path_2d = out_path / f"spill_2d_{suffix}.mp4"
        export_mp4(outputs, dem, profile, str(path_2d), fps=fps)
        result["video_2d_path"] = _existing_video_or_gif(path_2d)
    if include_3d:
        path_3d = out_path / f"spill_3d_{suffix}.mp4"
        export_mp4_3d(outputs, dem, profile, str(path_3d), fps=fps, spill_palette=video_3d_palette)
        result["video_3d_path"] = _existing_video_or_gif(path_3d)
    return result


def compute_arrival_map(outputs: list[dict[str, Any]], threshold: float = 1e-9) -> np.ndarray | None:
    ref = None
    for snap in outputs:
        h = snap.get("h")
        if h is not None:
            ref = np.asarray(h)
            break
    if ref is None:
        return None

    arrival = np.full(ref.shape, np.nan, dtype=float)
    for snap in outputs:
        h = snap.get("h")
        if h is None:
            continue
        h_arr = np.asarray(h, dtype=float)
        t = float(snap.get("time", np.nan))
        hit = np.isfinite(h_arr) & (h_arr > threshold) & ~np.isfinite(arrival)
        arrival[hit] = t
    return arrival


def _last_snapshot_array(outputs: list[dict[str, Any]], key: str, shape: tuple[int, ...]) -> np.ndarray | None:
    for snap in reversed(outputs):
        value = snap.get(key)
        if value is None:
            continue
        array = np.asarray(value, dtype=float)
        if array.shape == shape:
            return array
    return None


def _max_snapshot_array(outputs: list[dict[str, Any]], key: str, shape: tuple[int, ...]) -> np.ndarray | None:
    max_array: np.ndarray | None = None
    for snap in outputs:
        value = snap.get(key)
        if value is None:
            continue
        array = np.asarray(value, dtype=float)
        if array.shape != shape:
            continue
        if max_array is None:
            max_array = array.copy()
        else:
            np.maximum(max_array, array, out=max_array)
    return max_array


def _same_shape_array(value: np.ndarray | None, shape: tuple[int, ...]) -> np.ndarray | None:
    if value is None:
        return None
    array = np.asarray(value, dtype=float)
    return array if array.shape == shape else None


def _affine_cell_area(transform: Any) -> float:
    if transform is None:
        return 1.0
    try:
        return abs(float(transform.a) * float(transform.e) - float(transform.b) * float(transform.d))
    except AttributeError:
        return 1.0


def _zone_values(array: np.ndarray | None, mask: np.ndarray) -> np.ndarray:
    if array is None:
        return np.asarray([], dtype=float)
    values = np.asarray(array, dtype=float)[mask]
    return values[np.isfinite(values)]


def _zone_max(array: np.ndarray | None, mask: np.ndarray) -> float | None:
    values = _zone_values(array, mask)
    return float(np.max(values)) if values.size else None


def _zone_mean(array: np.ndarray | None, mask: np.ndarray) -> float | None:
    values = _zone_values(array, mask)
    return float(np.mean(values)) if values.size else None


def _zone_volume(array: np.ndarray | None, mask: np.ndarray, cell_area: float) -> float | None:
    values = _zone_values(array, mask)
    return float(np.sum(values) * cell_area) if values.size else None


def _isochrone_model_properties(
    model_params: dict[str, Any] | None,
    summary: dict[str, Any] | None,
    source_curve: dict[str, Any] | None,
    transform: Any,
) -> dict[str, Any]:
    params = {**default_model_params(), **(model_params or {})}
    props: dict[str, Any] = {
        "t_end_s": float(params["t_end"]),
        "cfl": float(params["cfl"]),
        "roi_buf_m": float(params["roi_buffer_m"]),
        "friction": str(params["friction_model"]),
        "manning_n": float(params["manning_n"]),
        "darcy_f": float(params["darcy_f"]),
        "rho_kg_m3": float(params["fluid_density"]),
        "mu_pa_s": float(params["dynamic_viscosity"]),
        "yield_pa": float(params["yield_stress"]),
        "h_visc_min": float(params["h_viscous_min"]),
        "k_oil_m_s": float(params["k_oil"]),
        "cap_suct": float(params["capillary_suction"]),
        "por_def": float(params["porosity_deficit"]),
        "evap_rate": float(params["evap_rate"]),
        "deg_rate": float(params["deg_rate"]),
        "snap_int_s": float(params["snapshot_interval"]),
        "dx_m": abs(float(transform.a)) if transform is not None else None,
        "dy_m": abs(float(transform.e)) if transform is not None else None,
        "cell_area": _affine_cell_area(transform),
    }
    if summary:
        props.update(
            {
                "max_cfl": summary.get("max_cfl"),
                "final_vol": summary.get("final_surface_mass"),
                "infil_vol": summary.get("infiltrated_mass"),
                "mass_err_pct": summary.get("mass_closure_error_pct"),
            }
        )
    if source_curve:
        props.update(
            {
                "src_mode": source_curve.get("mode"),
                "src_vol_m3": source_curve.get("volume_total_m3"),
                "src_dur_s": source_curve.get("duration_s"),
            }
        )
    return props


def _isochrone_zone_properties(
    mask: np.ndarray,
    *,
    h_max: np.ndarray | None,
    h_final: np.ndarray | None,
    z_f: np.ndarray | None,
    infiltration: np.ndarray | None,
    cell_area: float,
) -> dict[str, Any]:
    return {
        "area_m2": float(np.count_nonzero(mask) * cell_area),
        "h_max_m": _zone_max(h_max, mask),
        "h_final_max_m": _zone_max(h_final, mask),
        "h_final_mean_m": _zone_mean(h_final, mask),
        "z_f_m": _zone_max(z_f, mask),
        "z_f_mean_m": _zone_mean(z_f, mask),
        "inf_acc_m": _zone_mean(infiltration, mask),
        "inf_acc_max_m": _zone_max(infiltration, mask),
        "inf_acc_vol_m3": _zone_volume(infiltration, mask, cell_area),
    }


def export_isochrones_gpkg(
    outputs: list[dict[str, Any]],
    profile: dict[str, Any],
    out_dir: str | Path,
    interval_s: float = 60.0,
    threshold: float = 1e-9,
    site: dict[str, Any] | None = None,
    dem: np.ndarray | None = None,
    h_max: np.ndarray | None = None,
    model_params: dict[str, Any] | None = None,
    summary: dict[str, Any] | None = None,
    source_curve: dict[str, Any] | None = None,
) -> str | None:
    if interval_s <= 0.0:
        raise ValueError("interval_s must be > 0")
    if not shutil.which("ogr2ogr"):
        raise RuntimeError("ogr2ogr was not found in PATH")

    arrival = compute_arrival_map(outputs, threshold=threshold)
    if arrival is None or not np.isfinite(arrival).any():
        return None
    transform = profile.get("transform")
    if transform is None:
        raise ValueError("profile transform is required to export isochrone polygons")

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    nodata = -9999.0
    arrival_filled = np.where(np.isfinite(arrival), arrival, nodata).astype("float32")
    raster_path = out_path / "arrival_time_s.tif"
    gpkg_path = out_path / "isocronas.gpkg"
    geojson_path = out_path / "_isocronas_tmp.geojson"

    raster_profile = {
        "driver": "GTiff",
        "height": arrival_filled.shape[0],
        "width": arrival_filled.shape[1],
        "count": 1,
        "dtype": "float32",
        "nodata": nodata,
        "transform": profile.get("transform"),
    }
    if profile.get("crs") is not None:
        raster_profile["crs"] = profile.get("crs")

    with rasterio.open(raster_path, "w", **raster_profile) as dst:
        dst.write(arrival_filled, 1)

    finite = arrival[np.isfinite(arrival)]
    max_time = float(np.max(finite))
    if max_time <= 0.0:
        time_steps = [0.0]
    else:
        time_steps = [float(v) for v in np.arange(interval_s, max_time, interval_s)]
        time_steps.append(max_time)

    site = site or {}
    shape = arrival.shape
    cell_area = _affine_cell_area(transform)
    h_max_array = _same_shape_array(h_max, shape)
    if h_max_array is None:
        h_max_array = _max_snapshot_array(outputs, "h", shape)
    h_final = _last_snapshot_array(outputs, "h", shape)
    infiltration = _last_snapshot_array(outputs, "cumulative_infiltration", shape)
    dem_array = _same_shape_array(dem, shape)
    z_f = dem_array + h_final if dem_array is not None and h_final is not None else None
    model_properties = _isochrone_model_properties(model_params, summary, source_curve, transform)

    vector_features = []
    previous_time: float | None = None
    for time_s in time_steps:
        if previous_time is None:
            interval_mask = np.isfinite(arrival) & (arrival <= time_s)
            from_s = 0.0
        else:
            interval_mask = np.isfinite(arrival) & (arrival > previous_time) & (arrival <= time_s)
            from_s = float(previous_time)
        if not interval_mask.any():
            previous_time = time_s
            continue
        for geom, value in features.shapes(interval_mask.astype("uint8"), mask=interval_mask, transform=transform):
            if int(value) != 1:
                continue
            component_mask = (
                features.rasterize(
                    [(geom, 1)],
                    out_shape=arrival.shape,
                    transform=transform,
                    fill=0,
                    dtype="uint8",
                ).astype(bool)
                & interval_mask
            )
            if not component_mask.any():
                continue
            vector_features.append(
                {
                    "type": "Feature",
                    "properties": {
                        "time_s": float(time_s),
                        "from_s": from_s,
                        "to_s": float(time_s),
                        "kind": "interval",
                        "site": site.get("name"),
                        "row": site.get("row"),
                        "col": site.get("col"),
                        **_isochrone_zone_properties(
                            component_mask,
                            h_max=h_max_array,
                            h_final=h_final,
                            z_f=z_f,
                            infiltration=infiltration,
                            cell_area=cell_area,
                        ),
                        **model_properties,
                    },
                    "geometry": geom,
                }
            )
        previous_time = time_s

    if not vector_features:
        return None

    geojson = {"type": "FeatureCollection", "features": vector_features}
    geojson_path.write_text(json.dumps(geojson), encoding="utf-8")

    if gpkg_path.exists():
        gpkg_path.unlink()
    for suffix in ("-wal", "-shm"):
        sidecar = Path(str(gpkg_path) + suffix)
        if sidecar.exists():
            sidecar.unlink()

    cmd = ["ogr2ogr", "-f", "GPKG", str(gpkg_path), str(geojson_path), "-nln", "isocronas"]
    if profile.get("crs") is not None:
        cmd.extend(["-a_srs", str(profile.get("crs"))])
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    finally:
        if geojson_path.exists():
            geojson_path.unlink()
    return str(gpkg_path)


def build_cases(payload: dict[str, Any]) -> list[dict[str, Any]]:
    out_root = Path(payload["out_dir"])
    source_curve = source_curve_from_config(payload["source_curve"])
    cases = []
    used_names: dict[str, int] = {}
    decimate = int(payload.get("decimate", 1))
    target_resolution_m = _optional_positive_float(payload.get("target_resolution_m"))
    resampling_method = str(payload.get("resampling_method", "bilinear"))
    if target_resolution_m is not None:
        _resampling_from_name(resampling_method)
    with rasterio.open(payload["dem_path"]) as src:
        grid_height, grid_width, grid_transform = _target_grid(
            src,
            decimate=decimate,
            target_resolution_m=target_resolution_m,
        )
    for site in payload.get("sites", []):
        site = dict(site)
        if "x" in site and "y" in site:
            row, col, in_dem = point_to_rowcol(float(site["x"]), float(site["y"]), grid_transform, grid_height, grid_width)
            site.update({"row": row, "col": col, "in_dem": in_dem})
        if not site.get("in_dem", False):
            continue
        safe_name = sanitize_site_name(site.get("safe_name") or site.get("name"))
        count = used_names.get(safe_name, 0)
        used_names[safe_name] = count + 1
        folder = safe_name if count == 0 else f"{safe_name}_{count + 1}"
        cases.append(
            {
                "dem_path": payload["dem_path"],
                "vector_path": payload.get("vector_path"),
                "decimate": decimate,
                "target_resolution_m": target_resolution_m,
                "resampling_method": resampling_method,
                "site": site,
                "source_curve": source_curve.to_dict(),
                "model_params": payload.get("model_params", {}),
                "output_options": payload.get("output_options", {}),
                "out_dir": str(out_root / folder),
                "verbose": False,
                "solver_verbose": False,
            }
        )
    if not cases:
        raise ValueError("no selected sites inside the DEM")
    return cases


def write_runs_index(out_dir: str | Path, records: list[dict[str, Any]], payload: dict[str, Any]) -> str:
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=True)
    index_path = root / "runs_index.json"
    data = {
        "dem_path": payload.get("dem_path"),
        "vector_path": payload.get("vector_path"),
        "created_at": time.time(),
        "records": records,
    }
    with index_path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    return str(index_path)
