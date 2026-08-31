from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from Examples.run_dem_methodology_test import OIL_SPILL_CMAPS, _build_dem_cmap_norm, _build_oil_spill_cmap

import fullswof_oil.dem_workflow as dem_workflow
from fullswof_oil.dem_workflow import (
    SourceCurve,
    build_cases,
    compute_arrival_map,
    export_isochrones_gpkg,
    inspect_dem,
    interpolated_dem_path,
    load_source_sites,
    load_video_reexport_context,
    ogr_available,
    point_to_rowcol,
    prepare_dem,
    run_site_case,
    save_video_reexport_context,
    sanitize_site_name,
    validate_dem,
)


def test_constant_source_curve_from_volume_duration() -> None:
    curve = SourceCurve.constant(volume_total_m3=20.0, duration_s=10.0)

    assert curve.flow_at(0.0) == pytest.approx(2.0)
    assert curve.flow_at(5.0) == pytest.approx(2.0)
    assert curve.flow_at(11.0) == pytest.approx(0.0)


def test_csv_flow_source_curve(tmp_path) -> None:
    path = tmp_path / "flow.csv"
    path.write_text("time_s,flow_m3_s\n0,0\n10,2\n20,0\n", encoding="utf-8")

    curve = SourceCurve.from_csv(path)

    assert curve.mode == "csv_flow"
    assert curve.flow_at(5.0) == pytest.approx(1.0)
    assert curve.flow_at(20.1) == pytest.approx(0.0)


def test_csv_volume_source_curve_uses_interval_rates(tmp_path) -> None:
    path = tmp_path / "volume.csv"
    path.write_text("time_s,volume_m3\n0,0\n10,20\n20,30\n", encoding="utf-8")

    curve = SourceCurve.from_csv(path)

    assert curve.mode == "csv_volume"
    assert curve.flow_at(5.0) == pytest.approx(2.0)
    assert curve.flow_at(15.0) == pytest.approx(1.0)
    assert curve.flow_at(25.0) == pytest.approx(0.0)


def test_csv_source_curve_rejects_bad_time_order(tmp_path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text("time_s,flow_m3_s\n0,1\n0,2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="strictly increasing"):
        SourceCurve.from_csv(path)


def _assert_colormap_has_no_blue_dominant_colors(cmap) -> None:
    colors = cmap(np.linspace(0.0, 1.0, 256))[:, :3]
    blue_dominant = (colors[:, 2] > colors[:, 0]) & (colors[:, 2] > colors[:, 1])

    assert not np.any(blue_dominant)


def test_3d_video_palettes_have_no_blue_dominant_colors() -> None:
    dem_cmap, _ = _build_dem_cmap_norm(np.array([[0.0, 10.0], [20.0, 30.0]]))

    _assert_colormap_has_no_blue_dominant_colors(dem_cmap)
    for palette in OIL_SPILL_CMAPS:
        _assert_colormap_has_no_blue_dominant_colors(_build_oil_spill_cmap(palette))


def test_dem_palette_uses_dark_green_low_end() -> None:
    dem_cmap, _ = _build_dem_cmap_norm(np.array([[0.0, 10.0], [20.0, 30.0]]))
    low_color = np.asarray(dem_cmap(0.0))[:3]

    assert low_color[1] > low_color[0]
    assert low_color[1] > low_color[2]
    assert float(np.max(low_color)) < 0.35


def test_default_3d_video_palette_avoids_near_white_high_end() -> None:
    high_color = np.asarray(_build_oil_spill_cmap()(1.0))[:3]

    assert float(np.min(high_color)) < 0.75
    assert float(np.max(high_color)) < 0.95


def test_reexport_site_videos_passes_selected_3d_palette(monkeypatch, tmp_path) -> None:
    import Examples.run_dem_methodology_test as methodology

    captured: dict[str, str] = {}

    def fake_load_video_context(out_dir, t_end_s=None):
        outputs = [{"time": 0.0, "h": np.ones((2, 2), dtype=float)}]
        dem = np.zeros((2, 2), dtype=float)
        profile = {"height": 2, "width": 2, "transform": from_origin(0.0, 2.0, 1.0, 1.0)}
        return outputs, dem, profile

    def fake_export_mp4(outputs, dem, profile, out_path, **kwargs):
        Path(out_path).write_bytes(b"2d")

    def fake_export_mp4_3d(outputs, dem, profile, out_path, **kwargs):
        captured["spill_palette"] = kwargs["spill_palette"]
        Path(out_path).write_bytes(b"3d")

    monkeypatch.setattr(dem_workflow, "load_video_reexport_context", fake_load_video_context)
    monkeypatch.setattr(methodology, "export_mp4", fake_export_mp4)
    monkeypatch.setattr(methodology, "export_mp4_3d", fake_export_mp4_3d)

    result = dem_workflow.reexport_site_videos(
        tmp_path,
        fps=3,
        include_2d=False,
        include_3d=True,
        video_3d_palette="amber",
    )

    assert captured["spill_palette"] == "amber"
    assert result["video_3d_palette"] == "amber"
    assert result["video_3d_path"]


def test_sanitize_site_name_removes_invalid_path_chars() -> None:
    assert sanitize_site_name(' Pozo / A:*? " ') == "Pozo_A"
    assert sanitize_site_name("") == "site"


def test_point_to_rowcol_inside_and_outside_dem() -> None:
    transform = from_origin(100.0, 200.0, 10.0, 10.0)

    row, col, inside = point_to_rowcol(115.0, 185.0, transform, height=5, width=5)
    assert (row, col, inside) == (1, 1, True)

    row, col, inside = point_to_rowcol(10.0, 10.0, transform, height=5, width=5)
    assert inside is False


def _write_test_dem(path, width: int = 6, height: int = 6, pixel_size: float = 10.0) -> None:
    transform = from_origin(100.0, 200.0, pixel_size, pixel_size)
    data = np.zeros((height, width), dtype=np.float32)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=width,
        height=height,
        count=1,
        dtype="float32",
        transform=transform,
    ) as dst:
        dst.write(data, 1)


@pytest.mark.skipif(not shutil.which("gdalwarp"), reason="gdalwarp CLI not available")
def test_prepare_dem_resamples_to_target_resolution(tmp_path) -> None:
    dem_path = tmp_path / "dem_30m.tif"
    _write_test_dem(dem_path, width=4, height=3, pixel_size=30.0)

    dem, profile, dx, dy = prepare_dem(
        dem_path,
        target_resolution_m=10.0,
        resampling_method="bilinear",
    )

    assert dem.shape == (9, 12)
    assert profile["width"] == 12
    assert profile["height"] == 9
    assert dx == pytest.approx(10.0)
    assert dy == pytest.approx(10.0)
    assert interpolated_dem_path(dem_path, 10.0, "bilinear").exists()


def test_prepare_dem_rejects_unknown_resampling_method(tmp_path) -> None:
    dem_path = tmp_path / "dem.tif"
    _write_test_dem(dem_path)

    with pytest.raises(ValueError, match="bilinear"):
        prepare_dem(dem_path, target_resolution_m=5.0, resampling_method="lanczos")


@pytest.mark.skipif(not shutil.which("gdalwarp"), reason="gdalwarp CLI not available")
def test_prepare_dem_rejects_too_many_resampled_cells(tmp_path) -> None:
    dem_path = tmp_path / "dem_30m.tif"
    _write_test_dem(dem_path, width=10, height=10, pixel_size=30.0)

    with pytest.raises(ValueError, match="above the safe limit"):
        prepare_dem(dem_path, target_resolution_m=1.0, max_cells=1000)


@pytest.mark.skipif(not ogr_available(), reason="GDAL/OGR CLI not available")
def test_load_source_sites_from_geojson(tmp_path) -> None:
    dem_path = tmp_path / "dem.tif"
    _write_test_dem(dem_path)
    vector_path = tmp_path / "sites.geojson"
    vector_path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "name": "sites",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"Name": "Site A"},
                        "geometry": {"type": "Point", "coordinates": [115.0, 185.0]},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    sites = load_source_sites(vector_path, dem_path, layer_name="sites")

    assert sites[0]["name"] == "Site A"
    assert sites[0]["row"] == 1
    assert sites[0]["col"] == 1
    assert sites[0]["in_dem"] is True


@pytest.mark.skipif(not ogr_available(), reason="GDAL/OGR CLI not available")
def test_load_source_sites_uses_target_resolution_grid(tmp_path) -> None:
    dem_path = tmp_path / "dem_30m.tif"
    _write_test_dem(dem_path, width=4, height=4, pixel_size=30.0)
    vector_path = tmp_path / "sites.geojson"
    vector_path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "name": "sites",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"Name": "Site A"},
                        "geometry": {"type": "Point", "coordinates": [115.0, 185.0]},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    sites = load_source_sites(
        vector_path,
        dem_path,
        layer_name="sites",
        target_resolution_m=10.0,
        resampling_method="cubic",
    )

    assert sites[0]["row"] == 1
    assert sites[0]["col"] == 1
    assert sites[0]["in_dem"] is True


@pytest.mark.skipif(not ogr_available(), reason="GDAL/OGR CLI not available")
def test_load_source_sites_uses_polygon_centroid_and_nombre_fallback(tmp_path) -> None:
    dem_path = tmp_path / "dem.tif"
    _write_test_dem(dem_path)
    vector_path = tmp_path / "pads.geojson"
    vector_path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "name": "pads",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"NOMBRE": "Pad A"},
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [
                                [
                                    [110.0, 190.0],
                                    [120.0, 190.0],
                                    [120.0, 180.0],
                                    [110.0, 180.0],
                                    [110.0, 190.0],
                                ]
                            ],
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    sites = load_source_sites(vector_path, dem_path, layer_name="pads", name_field="Name")

    assert sites[0]["name"] == "Pad A"
    assert sites[0]["row"] == 1
    assert sites[0]["col"] == 1
    assert sites[0]["in_dem"] is True


def test_inspect_dem_reports_basic_stats(tmp_path) -> None:
    dem_path = tmp_path / "dem.tif"
    _write_test_dem(dem_path)

    info = inspect_dem(dem_path)

    assert info["width"] == 6
    assert info["height"] == 6
    assert info["crs"] is None


def test_inspect_dem_reports_progress(tmp_path) -> None:
    dem_path = tmp_path / "dem.tif"
    _write_test_dem(dem_path)
    updates = []

    inspect_dem(dem_path, progress=lambda percent, message: updates.append((percent, message)))

    assert updates[0][0] == 5
    assert updates[-1] == (100, "DEM listo")
    assert any("Validando" in message for _, message in updates)


def test_validate_dem_fills_large_nodata_hole() -> None:
    dem = np.arange(400, dtype=float).reshape(20, 20)
    dem[5:15, 5:15] = -9999.0

    clean, stats = validate_dem(dem, nodata=-9999.0)

    assert stats["invalid_total"] == 100
    assert stats["invalid_remaining"] == 0
    assert np.isfinite(clean).all()
    assert not np.any(clean == -9999.0)


def test_compute_arrival_map_uses_first_wet_time() -> None:
    outputs = [
        {"time": 0.0, "h": np.array([[0.0, 0.0], [0.0, 0.0]])},
        {"time": 5.0, "h": np.array([[1.0, 0.0], [0.0, 0.0]])},
        {"time": 10.0, "h": np.array([[2.0, 1.0], [0.0, 0.0]])},
    ]

    arrival = compute_arrival_map(outputs)

    assert arrival is not None
    assert arrival[0, 0] == pytest.approx(5.0)
    assert arrival[0, 1] == pytest.approx(10.0)
    assert np.isnan(arrival[1, 0])


@pytest.mark.skipif(not shutil.which("ogr2ogr"), reason="ogr2ogr CLI not available")
def test_export_isochrones_gpkg_from_synthetic_arrivals_as_polygons(tmp_path) -> None:
    outputs = []
    for step, time_s in enumerate([1.0, 2.0, 3.0, 4.0], start=1):
        h = np.zeros((5, 5), dtype=float)
        h[:, :step] = 1.0
        infiltration = np.zeros((5, 5), dtype=float)
        infiltration[:, :step] = 0.05 * step
        outputs.append({"time": time_s, "h": h, "cumulative_infiltration": infiltration})
    profile = {"transform": from_origin(0.0, 50.0, 10.0, 10.0), "crs": None}
    dem = np.full((5, 5), 100.0, dtype=float)
    model_params = {"friction_model": "viscous", "dynamic_viscosity": 0.02, "t_end": 4.0, "cfl": 0.25}
    summary = {"max_cfl": 0.24, "final_surface_mass": 500.0, "infiltrated_mass": 20.0}
    source_curve = {"mode": "constant", "volume_total_m3": 10.0, "duration_s": 2.0}

    path = export_isochrones_gpkg(
        outputs,
        profile,
        tmp_path,
        interval_s=1.0,
        dem=dem,
        model_params=model_params,
        summary=summary,
        source_curve=source_curve,
    )

    assert path is not None
    assert (tmp_path / "arrival_time_s.tif").exists()
    assert (tmp_path / "isocronas.gpkg").exists()
    if shutil.which("ogrinfo"):
        info = subprocess.run(
            ["ogrinfo", "-so", str(tmp_path / "isocronas.gpkg"), "isocronas"],
            check=True,
            capture_output=True,
            text=True,
        )
        assert "Polygon" in info.stdout
        assert "Feature Count: 4" in info.stdout
        assert "h_max_m: Real" in info.stdout
        assert "z_f_m: Real" in info.stdout
        assert "inf_acc_m: Real" in info.stdout
        assert "friction: String" in info.stdout
        assert "mu_pa_s: Real" in info.stdout
        features_text = subprocess.run(
            ["ogrinfo", "-al", "-geom=NO", str(tmp_path / "isocronas.gpkg"), "isocronas"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        assert 'kind (String) = interval' in features_text
        assert "from_s (Real) = 1" in features_text
        assert "to_s (Real) = 2" in features_text
        assert "friction (String) = viscous" in features_text
        assert "src_vol_m3 (Real) = 10" in features_text


def test_video_reexport_context_saves_and_filters_snapshots(tmp_path) -> None:
    dem = np.zeros((3, 3), dtype=float)
    profile = {"height": 3, "width": 3, "transform": from_origin(0.0, 30.0, 10.0, 10.0)}
    outputs = [
        {"time": 0.0, "h": np.zeros((3, 3), dtype=float)},
        {"time": 10.0, "h": np.ones((3, 3), dtype=float)},
        {"time": 20.0, "h": np.full((3, 3), 2.0, dtype=float)},
    ]

    path = save_video_reexport_context(outputs, dem, profile, tmp_path)
    loaded_outputs, loaded_dem, loaded_profile = load_video_reexport_context(tmp_path, t_end_s=10.0)
    context = json.loads((tmp_path / "video_context.json").read_text(encoding="utf-8"))

    assert (tmp_path / "video_context.json").exists()
    assert Path(path).exists()
    assert context["snapshot_count"] == 3
    assert context["snapshots_size_bytes"] > 0
    assert context["snapshots_uncompressed_bytes"] > 0
    assert len(loaded_outputs) == 2
    assert loaded_outputs[-1]["time"] == pytest.approx(10.0)
    assert loaded_dem.shape == (3, 3)
    assert loaded_profile["transform"].a == pytest.approx(10.0)


def test_run_site_case_light_integration_without_media(tmp_path) -> None:
    dem_path = tmp_path / "dem.tif"
    out_dir = tmp_path / "out"
    _write_test_dem(dem_path, width=8, height=8)
    site = {
        "name": "Site A",
        "safe_name": "Site_A",
        "x": 135.0,
        "y": 165.0,
        "row": 3,
        "col": 3,
        "in_dem": True,
    }
    payload = {
        "dem_path": str(dem_path),
        "out_dir": str(out_dir),
        "decimate": 1,
        "sites": [site],
        "source_curve": {"mode": "constant", "volume_total_m3": 0.01, "duration_s": 0.1},
        "model_params": {"t_end": 0.2, "snapshot_interval": 0.1, "cfl": 0.3},
        "output_options": {
            "arrays": True,
            "maps": False,
            "diagnostic_report": False,
            "videos": False,
            "isochrones": False,
            "snapshots": False,
        },
    }
    case = build_cases(payload)[0]

    result = run_site_case(case)
    summary = json.loads((out_dir / "Site_A" / "summary.json").read_text(encoding="utf-8"))

    assert result["status"] == "completed"
    assert (out_dir / "Site_A" / "summary.json").exists()
    assert (out_dir / "Site_A" / "final_h.npy").exists()
    assert (out_dir / "Site_A" / "max_h.npy").exists()
    assert result["max_thickness_raster_path"] == str(out_dir / "Site_A" / "max_thickness_m.tif")
    with rasterio.open(out_dir / "Site_A" / "max_thickness_m.tif") as src:
        max_thickness = src.read(1)
        assert src.width == 8
        assert src.height == 8
        assert src.transform == from_origin(100.0, 200.0, 10.0, 10.0)
        assert max_thickness.dtype == np.float32
        final_max = float(np.max(np.load(out_dir / "Site_A" / "final_h.npy")))
        assert float(np.max(max_thickness)) + 1e-9 >= final_max
    assert result["snapshots_path"] is None
    assert result["snapshots_enabled"] is False
    assert result["snapshots_size_bytes"] == 0
    assert summary["snapshots_enabled"] is False
    assert not (out_dir / "Site_A" / "snapshots.npz").exists()


def test_run_site_case_uses_local_domain_crop(tmp_path) -> None:
    dem_path = tmp_path / "dem.tif"
    out_dir = tmp_path / "out"
    _write_test_dem(dem_path, width=9, height=9)
    site = {
        "name": "Site A",
        "safe_name": "Site_A",
        "x": 145.0,
        "y": 155.0,
        "row": 4,
        "col": 4,
        "in_dem": True,
    }
    payload = {
        "dem_path": str(dem_path),
        "out_dir": str(out_dir),
        "decimate": 1,
        "sites": [site],
        "source_curve": {"mode": "constant", "volume_total_m3": 0.01, "duration_s": 0.1},
        "model_params": {"t_end": 0.2, "snapshot_interval": 0.1, "cfl": 0.3, "roi_buffer_m": 10.0},
        "output_options": {
            "arrays": True,
            "maps": False,
            "diagnostic_report": False,
            "videos": False,
            "isochrones": False,
        },
    }

    result = run_site_case(build_cases(payload)[0])
    final_h = np.load(out_dir / "Site_A" / "final_h.npy")
    summary = json.loads((out_dir / "Site_A" / "summary.json").read_text(encoding="utf-8"))

    assert result["status"] == "completed"
    assert final_h.shape == (3, 3)
    assert summary["site"]["global_row"] == 4
    assert summary["site"]["row"] == 1
    assert summary["domain_crop"]["row_start"] == 3
