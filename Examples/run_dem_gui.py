from __future__ import annotations

import argparse
import json
import os
import threading
import time
import uuid
import webbrowser
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import numpy as np
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fullswof_oil.dem_workflow import (
    build_cases,
    default_model_params,
    inspect_dem,
    list_vector_layers,
    load_source_sites,
    ogr_available,
    reexport_site_videos,
    run_site_case,
    source_curve_from_config,
    write_runs_index,
)


RUNS: dict[str, dict[str, Any]] = {}
RUNS_LOCK = threading.Lock()
DEM_INFO_JOBS: dict[str, dict[str, Any]] = {}
DEM_INFO_LOCK = threading.Lock()
VIDEO_JOBS: dict[str, dict[str, Any]] = {}
VIDEO_JOBS_LOCK = threading.Lock()
DIALOG_LOCK = threading.Lock()


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _dialog_select_file(title: str, filetypes: list[tuple[str, str]]) -> str:
    import tkinter as tk
    from tkinter import filedialog

    with DIALOG_LOCK:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        try:
            path = filedialog.askopenfilename(title=title, filetypes=filetypes)
        finally:
            root.destroy()
    return path or ""


def _dialog_select_dir(title: str) -> str:
    import tkinter as tk
    from tkinter import filedialog

    with DIALOG_LOCK:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        try:
            path = filedialog.askdirectory(title=title)
        finally:
            root.destroy()
    return path or ""


def _set_run(job_id: str, **updates: Any) -> None:
    with RUNS_LOCK:
        RUNS.setdefault(job_id, {}).update(updates)


def _set_dem_info_job(job_id: str, **updates: Any) -> None:
    with DEM_INFO_LOCK:
        DEM_INFO_JOBS.setdefault(job_id, {}).update(updates)


def _set_video_job(job_id: str, **updates: Any) -> None:
    with VIDEO_JOBS_LOCK:
        VIDEO_JOBS.setdefault(job_id, {}).update(updates)


def _run_video_job(job_id: str, payload: dict[str, Any]) -> None:
    try:
        _set_video_job(job_id, status="running", message="Reexportando videos", updated_at=time.time())
        t_end_raw = payload.get("t_end_s")
        t_end_s = None if t_end_raw in (None, "") else float(t_end_raw)
        result = reexport_site_videos(
            payload["out_dir"],
            fps=int(payload.get("fps", 2)),
            t_end_s=t_end_s,
            include_2d=bool(payload.get("video_2d", True)),
            include_3d=bool(payload.get("video_3d", True)),
            video_3d_palette=str(payload.get("video_3d_palette", "oil_dark")),
        )
        _set_video_job(job_id, **result, message="Videos listos", updated_at=time.time())
    except Exception as exc:
        _set_video_job(job_id, status="failed", message=str(exc), error=str(exc), updated_at=time.time())


def _run_dem_info_job(job_id: str, payload: dict[str, Any]) -> None:
    def progress(percent: int, message: str) -> None:
        _set_dem_info_job(job_id, status="running", progress=percent, message=message, updated_at=time.time())

    try:
        progress(0, "En cola")
        info = inspect_dem(
            payload["path"],
            int(payload.get("decimate", 1)),
            payload.get("target_resolution_m"),
            payload.get("resampling_method", "bilinear"),
            progress=progress,
        )
        _set_dem_info_job(
            job_id,
            status="completed",
            progress=100,
            message="DEM listo",
            info=info,
            updated_at=time.time(),
        )
    except Exception as exc:
        _set_dem_info_job(
            job_id,
            status="failed",
            progress=100,
            message=str(exc),
            error=str(exc),
            updated_at=time.time(),
        )


def _update_record(job_id: str, index: int, **updates: Any) -> None:
    with RUNS_LOCK:
        job = RUNS.setdefault(job_id, {})
        records = job.setdefault("records", [])
        records[index].update(updates)
        job["completed"] = sum(1 for record in records if record.get("status") in {"completed", "failed"})


def _run_job(job_id: str, payload: dict[str, Any]) -> None:
    try:
        cases = build_cases(payload)
        max_workers = max(1, int(payload.get("max_workers", 1)))
        initial_records = [
            {
                "status": "queued",
                "message": "En cola",
                "site": case.get("site", {}),
                "out_dir": case.get("out_dir"),
                "error": "",
            }
            for case in cases
        ]
        _set_run(job_id, status="running", total=len(cases), completed=0, records=initial_records, message="Running")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            case_iter = iter(enumerate(cases))
            futures: dict[Any, tuple[int, dict[str, Any]]] = {}

            def submit_next() -> bool:
                try:
                    index, case = next(case_iter)
                except StopIteration:
                    return False
                _update_record(job_id, index, status="running", message="Procesando")
                futures[executor.submit(run_site_case, case)] = (index, case)
                return True

            for _ in range(min(max_workers, len(cases))):
                submit_next()

            while futures:
                done, _ = wait(futures.keys(), return_when=FIRST_COMPLETED)
                for future in done:
                    index, case = futures.pop(future)
                    try:
                        record = future.result()
                    except Exception as exc:
                        record = {
                            "status": "failed",
                            "error": str(exc),
                            "site": case.get("site", {}),
                            "out_dir": case.get("out_dir"),
                        }
                    _update_record(job_id, index, **record, message=record.get("status", "done"))
                    submit_next()
        with RUNS_LOCK:
            records = list(RUNS[job_id].get("records", []))
        index_path = write_runs_index(payload["out_dir"], records, payload)
        _set_run(job_id, status="completed", message="Completed", index_path=index_path)
    except Exception as exc:
        _set_run(job_id, status="failed", message=str(exc), error=str(exc))


class DemGuiHandler(BaseHTTPRequestHandler):
    def log_message(self, fmt: str, *args: Any) -> None:
        print(f"[GUI] {self.address_string()} - {fmt % args}")

    def do_GET(self) -> None:
        if self.path in {"/", "/index.html"}:
            self._send_html(INDEX_HTML)
            return
        if self.path == "/api/defaults":
            self._send_json({"model_params": default_model_params(), "ogr_available": ogr_available()})
            return
        self.send_error(HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:
        try:
            payload = self._read_json()
            if self.path == "/api/select-file":
                self._api_select_file(payload)
            elif self.path == "/api/select-dir":
                self._send_json({"path": _dialog_select_dir("Seleccione carpeta de salida")})
            elif self.path == "/api/dem-info":
                self._send_json(
                    inspect_dem(
                        payload["path"],
                        int(payload.get("decimate", 1)),
                        payload.get("target_resolution_m"),
                        payload.get("resampling_method", "bilinear"),
                    )
                )
            elif self.path == "/api/start-dem-info":
                job_id = uuid.uuid4().hex
                with DEM_INFO_LOCK:
                    DEM_INFO_JOBS[job_id] = {
                        "status": "queued",
                        "progress": 0,
                        "message": "En cola",
                        "created_at": time.time(),
                        "updated_at": time.time(),
                    }
                thread = threading.Thread(target=_run_dem_info_job, args=(job_id, payload), daemon=True)
                thread.start()
                self._send_json({"job_id": job_id})
            elif self.path == "/api/dem-info-status":
                job_id = payload["job_id"]
                with DEM_INFO_LOCK:
                    status = DEM_INFO_JOBS.get(job_id)
                if status is None:
                    self._send_json({"status": "missing", "message": "Unknown DEM job"}, status=HTTPStatus.NOT_FOUND)
                else:
                    self._send_json(status)
            elif self.path == "/api/vector-layers":
                self._send_json({"layers": list_vector_layers(payload["path"])})
            elif self.path == "/api/sites":
                sites = load_source_sites(
                    payload["vector_path"],
                    payload["dem_path"],
                    payload.get("layer_name") or None,
                    payload.get("name_field", "Name"),
                    int(payload.get("decimate", 1)),
                    payload.get("target_resolution_m"),
                    payload.get("resampling_method", "bilinear"),
                )
                self._send_json({"sites": sites})
            elif self.path == "/api/validate-curve":
                curve = source_curve_from_config(payload)
                self._send_json({"curve": curve.to_dict()})
            elif self.path == "/api/start-run":
                job_id = uuid.uuid4().hex
                with RUNS_LOCK:
                    RUNS[job_id] = {
                        "status": "queued",
                        "message": "Queued",
                        "created_at": time.time(),
                        "total": 0,
                        "completed": 0,
                        "records": [],
                    }
                thread = threading.Thread(target=_run_job, args=(job_id, payload), daemon=True)
                thread.start()
                self._send_json({"job_id": job_id})
            elif self.path == "/api/status":
                job_id = payload["job_id"]
                with RUNS_LOCK:
                    status = RUNS.get(job_id)
                if status is None:
                    self._send_json({"status": "missing", "message": "Unknown job"}, status=HTTPStatus.NOT_FOUND)
                else:
                    self._send_json(status)
            elif self.path == "/api/start-video-export":
                job_id = uuid.uuid4().hex
                with VIDEO_JOBS_LOCK:
                    VIDEO_JOBS[job_id] = {
                        "status": "queued",
                        "message": "En cola",
                        "created_at": time.time(),
                        "updated_at": time.time(),
                    }
                thread = threading.Thread(target=_run_video_job, args=(job_id, payload), daemon=True)
                thread.start()
                self._send_json({"job_id": job_id})
            elif self.path == "/api/video-export-status":
                job_id = payload["job_id"]
                with VIDEO_JOBS_LOCK:
                    status = VIDEO_JOBS.get(job_id)
                if status is None:
                    self._send_json({"status": "missing", "message": "Unknown video job"}, status=HTTPStatus.NOT_FOUND)
                else:
                    self._send_json(status)
            else:
                self.send_error(HTTPStatus.NOT_FOUND)
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)

    def _api_select_file(self, payload: dict[str, Any]) -> None:
        kind = payload.get("kind", "")
        if kind == "dem":
            path = _dialog_select_file("Seleccione DEM GeoTIFF", [("GeoTIFF", "*.tif *.tiff"), ("Todos", "*.*")])
        elif kind == "vector":
            path = _dialog_select_file("Seleccione puntos GPKG/SHP", [("Vector", "*.gpkg *.shp"), ("Todos", "*.*")])
        elif kind == "csv":
            path = _dialog_select_file("Seleccione curva CSV", [("CSV", "*.csv"), ("Todos", "*.*")])
        else:
            path = _dialog_select_file("Seleccione archivo", [("Todos", "*.*")])
        self._send_json({"path": path})

    def _read_json(self) -> dict[str, Any]:
        size = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(size).decode("utf-8") if size else "{}"
        return json.loads(raw or "{}")

    def _send_json(self, data: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(data, ensure_ascii=False, default=_json_default).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_html(self, html: str) -> None:
        body = html.encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def run_server(host: str, port: int, open_browser: bool) -> None:
    server = ThreadingHTTPServer((host, port), DemGuiHandler)
    url = f"http://{host}:{server.server_port}"
    print(f"FullSWOF DEM GUI listening at {url}")
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("Stopping GUI server")
    finally:
        server.server_close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Local web GUI for FullSWOF DEM multi-site runs")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    return parser


INDEX_HTML = r"""<!doctype html>
<html lang="es">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>FullSWOF DEM Multi-Sitio</title>
  <style>
    :root {
      color-scheme: light;
      font-family: "Segoe UI", Arial, sans-serif;
      --bg: #f6f7f9;
      --panel: #ffffff;
      --line: #cdd3da;
      --text: #1f2933;
      --muted: #5b6773;
      --accent: #0f766e;
      --accent-dark: #115e59;
      --warn: #9a3412;
      --bad: #b91c1c;
      --good: #047857;
    }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--text);
    }
    header {
      background: #263238;
      color: white;
      padding: 14px 22px;
      display: flex;
      align-items: baseline;
      gap: 14px;
    }
    header h1 {
      font-size: 18px;
      margin: 0;
      font-weight: 650;
    }
    header span {
      color: #d8dee4;
      font-size: 13px;
    }
    main {
      padding: 18px;
      display: grid;
      grid-template-columns: minmax(360px, 480px) minmax(520px, 1fr);
      gap: 16px;
    }
    section {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 14px;
      margin-bottom: 14px;
    }
    h2 {
      font-size: 15px;
      margin: 0 0 12px;
    }
    label {
      display: block;
      font-size: 12px;
      color: var(--muted);
      margin-bottom: 4px;
    }
    input, select {
      box-sizing: border-box;
      width: 100%;
      min-height: 32px;
      border: 1px solid var(--line);
      border-radius: 4px;
      padding: 6px 8px;
      font: inherit;
      background: white;
    }
    input[type="checkbox"] {
      width: auto;
      min-height: 0;
    }
    button {
      border: 1px solid #0f766e;
      background: var(--accent);
      color: white;
      border-radius: 4px;
      min-height: 32px;
      padding: 6px 10px;
      font-weight: 600;
      cursor: pointer;
    }
    button.secondary {
      background: white;
      color: var(--accent-dark);
    }
    button:disabled {
      opacity: 0.55;
      cursor: not-allowed;
    }
    .row {
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px;
      align-items: end;
      margin-bottom: 10px;
    }
    .grid2 {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px;
    }
    .grid3 {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 10px;
    }
    .muted {
      color: var(--muted);
      font-size: 12px;
    }
    .status {
      font-size: 13px;
      white-space: pre-wrap;
      border-left: 3px solid var(--line);
      padding-left: 8px;
      min-height: 20px;
    }
    .site-table {
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
    }
    .site-table th, .site-table td {
      border-bottom: 1px solid #e3e7eb;
      padding: 6px;
      text-align: left;
    }
    .site-table th {
      background: #f2f5f7;
      position: sticky;
      top: 0;
    }
    .table-wrap {
      max-height: 360px;
      overflow: auto;
      border: 1px solid #e3e7eb;
      border-radius: 4px;
    }
    .ok { color: var(--good); font-weight: 650; }
    .bad { color: var(--bad); font-weight: 650; }
    .warn { color: var(--warn); font-weight: 650; }
    @media (max-width: 980px) {
      main { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
<header>
  <h1>FullSWOF DEM Multi-Sitio</h1>
  <span>GUI web local para corridas independientes por punto fuente</span>
</header>
<main>
  <div>
    <section>
      <h2>Archivos</h2>
      <div class="row">
        <div><label>DEM GeoTIFF</label><input id="demPath" placeholder="Seleccione o pegue ruta .tif"></div>
        <button onclick="selectFile('dem','demPath')">Seleccionar</button>
      </div>
      <div class="row">
        <div><label>Puntos GPKG/SHP</label><input id="vectorPath" placeholder="Seleccione o pegue ruta .gpkg/.shp"></div>
        <button onclick="selectFile('vector','vectorPath')">Seleccionar</button>
      </div>
      <div class="row">
        <div><label>Carpeta de salida</label><input id="outDir" placeholder="Seleccione carpeta raíz de salida"></div>
        <button onclick="selectDir()">Seleccionar</button>
      </div>
      <div class="grid2">
        <div><label>Decimate</label><input id="decimate" type="number" min="1" step="1" value="1"></div>
        <div><label>Campo identificador</label><input id="nameField" value="Name"></div>
        <div><label>Resolución objetivo (m)</label><input id="targetResolution" type="number" min="0.001" step="1" placeholder="Ej. 10"></div>
        <div><label>Interpolación DEM</label><select id="resamplingMethod"><option value="bilinear">bilinear</option><option value="cubic">cubic</option></select></div>
      </div>
      <div style="margin-top:10px; display:flex; gap:8px;">
        <button id="loadDemButton" class="secondary" onclick="loadDem()">Leer DEM</button>
        <button class="secondary" onclick="loadLayers()">Leer capas</button>
        <button class="secondary" onclick="loadSites()">Cargar sitios</button>
      </div>
      <p id="fileStatus" class="status"></p>
    </section>

    <section>
      <h2>Curva fuente</h2>
      <div class="grid2">
        <div>
          <label>Modo</label>
          <select id="curveMode" onchange="toggleCurveMode()">
            <option value="constant">Volumen + duración</option>
            <option value="csv">CSV</option>
          </select>
        </div>
        <div><label>Duración (s)</label><input id="durationS" type="number" step="1" value="100"></div>
      </div>
      <div id="constantCurve" style="margin-top:10px;">
        <label>Volumen total (m3)</label><input id="volumeTotal" type="number" step="0.001" value="50">
      </div>
      <div id="csvCurve" style="margin-top:10px; display:none;">
        <div class="row">
          <div><label>CSV time_s + flow_m3_s o volume_m3</label><input id="csvPath" placeholder="Ruta CSV"></div>
          <button onclick="selectFile('csv','csvPath')">Seleccionar</button>
        </div>
      </div>
    </section>

    <section>
      <h2>Paralelo y salidas</h2>
      <div class="grid2">
        <div><label>Max workers</label><input id="maxWorkers" type="number" min="1" step="1" value="1"></div>
        <div><label>FPS videos</label><input id="videoFps" type="number" min="1" step="1" value="2"></div>
        <div><label>Intervalo isocronas (s)</label><input id="isochroneInterval" type="number" min="0.001" step="1" value="60"></div>
        <div><label>Paleta derrame 3D</label><select id="video3dPalette" onchange="syncReexportPalette()"><option value="oil_dark">Petróleo oscuro</option><option value="oil_fire">Fuego</option><option value="amber">Ámbar</option></select></div>
      </div>
      <div class="grid2" style="margin-top:10px;">
        <label><input id="maps" type="checkbox" checked> Mapas</label>
        <label><input id="diagnosticReport" type="checkbox" checked> Reporte</label>
        <label><input id="videos" type="checkbox" checked> Videos 2D/3D</label>
        <label><input id="arrays" type="checkbox" checked> Arrays NPY</label>
        <label><input id="rasters" type="checkbox" checked> Rasters GeoTIFF</label>
        <label><input id="isochrones" type="checkbox" checked> Isocronas por intervalo GPKG</label>
        <label><input id="snapshots" type="checkbox" checked> Snapshots para reexportar video</label>
      </div>
    </section>
  </div>

  <div>
    <section>
      <h2>Parámetros del modelo</h2>
      <div class="grid3">
        <div><label>t_end (s)</label><input id="tEnd" type="number" step="1" value="1200"></div>
        <div><label>CFL</label><input id="cfl" type="number" step="0.01" value="0.3"></div>
        <div><label>ROI buffer (m)</label><input id="roiBuffer" type="number" step="1" value="500"></div>
        <div><label>Fricción</label><select id="frictionModel"><option>manning</option><option>none</option><option>darcy-weisbach</option><option>viscous</option></select></div>
        <div><label>Manning n</label><input id="manningN" type="number" step="0.001" value="0.03"></div>
        <div><label>Darcy f</label><input id="darcyF" type="number" step="0.001" value="0.02"></div>
        <div><label>Densidad (kg/m3)</label><input id="fluidDensity" type="number" step="1" value="900"></div>
        <div><label>Viscosidad (Pa*s)</label><input id="dynamicViscosity" type="number" step="0.001" value="0.005"></div>
        <div><label>h visc min (m)</label><input id="hViscousMin" type="number" step="0.0001" value="0.0001"></div>
        <div><label>K oil</label><input id="kOil" type="number" step="0.000001" value="0.000005"></div>
        <div><label>Capillary suction</label><input id="capillarySuction" type="number" step="0.001" value="0"></div>
        <div><label>Porosity deficit</label><input id="porosityDeficit" type="number" step="0.01" value="0.2"></div>
        <div><label>Evap rate</label><input id="evapRate" type="number" step="0.000001" value="0"></div>
        <div><label>Deg rate</label><input id="degRate" type="number" step="0.000001" value="0"></div>
        <div><label>Snapshot interval</label><input id="snapshotInterval" type="number" step="1" value="10"></div>
      </div>
    </section>

    <section>
      <h2>Sitios</h2>
      <div class="grid2" style="margin-bottom:8px;">
        <div><label>Capa</label><select id="layerName"></select></div>
        <div><label>Filtro Name</label><input id="siteFilter" placeholder="Filtrar..." oninput="renderSites()"></div>
      </div>
      <div style="margin-bottom:8px; display:flex; gap:8px;">
        <button class="secondary" onclick="setAllSites(true)">Seleccionar visibles</button>
        <button class="secondary" onclick="setAllSites(false)">Quitar visibles</button>
      </div>
      <div class="table-wrap"><table class="site-table" id="sitesTable"></table></div>
      <p id="siteStatus" class="muted"></p>
    </section>

    <section>
      <h2>Ejecución</h2>
      <button id="runButton" onclick="startRun()">Correr sitios seleccionados</button>
      <p id="runStatus" class="status"></p>
      <div class="grid2" style="margin:10px 0;">
        <div><label>FPS reexport</label><input id="reexportFps" type="number" min="1" step="1" value="2"></div>
        <div><label>Cortar video hasta t (s)</label><input id="videoCutoffS" type="number" min="0.001" step="1" placeholder="Vacio = completo"></div>
        <div><label>Paleta derrame 3D reexport</label><select id="reexport3dPalette" onchange="reexportPaletteTouched = true"><option value="oil_dark">Petróleo oscuro</option><option value="oil_fire">Fuego</option><option value="amber">Ámbar</option></select></div>
        <label><input id="reexport2d" type="checkbox" checked> Video 2D</label>
        <label><input id="reexport3d" type="checkbox" checked> Video 3D</label>
      </div>
      <p id="videoExportStatus" class="status"></p>
      <div class="table-wrap" style="max-height:240px;"><table class="site-table" id="resultsTable"></table></div>
    </section>
  </div>
</main>
<script>
let layers = [];
let sites = [];
let selected = new Set();
let currentJob = null;
let pollTimer = null;
let currentDemJob = null;
let demPollTimer = null;
let demFetchFailures = 0;
let currentVideoJob = null;
let videoPollTimer = null;
let reexportPaletteTouched = false;

async function api(path, payload = {}) {
  let res;
  try {
    res = await fetch(path, {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(payload)
    });
  } catch (err) {
    throw new Error(`No se pudo conectar con el backend local (${err.message}). Revise que run_dem_gui.py siga abierto y recargue la pagina.`);
  }
  const text = await res.text();
  let data = {};
  if (text) {
    try {
      data = JSON.parse(text);
    } catch (_err) {
      throw new Error(`El backend respondio ${res.status}, pero no envio JSON. Recargue la pagina o reinicie run_dem_gui.py.`);
    }
  }
  if (!res.ok) throw new Error(data.error || data.message || 'Error');
  return data;
}

function val(id) { return document.getElementById(id).value; }
function num(id) { return Number(val(id)); }
function optNum(id) {
  const text = val(id).trim();
  return text === '' ? null : Number(text);
}
function checked(id) { return document.getElementById(id).checked; }
function setStatus(id, text, cls='') {
  const el = document.getElementById(id);
  el.className = 'status ' + cls;
  el.textContent = text;
}

function paletteLabel(value) {
  return {
    oil_dark: 'Petróleo oscuro',
    oil_fire: 'Fuego',
    amber: 'Ámbar'
  }[value] || value;
}

function syncReexportPalette() {
  if (!reexportPaletteTouched) {
    document.getElementById('reexport3dPalette').value = val('video3dPalette');
  }
}

async function selectFile(kind, target) {
  const data = await api('/api/select-file', {kind});
  if (data.path) document.getElementById(target).value = data.path;
}

async function selectDir() {
  const data = await api('/api/select-dir', {});
  if (data.path) document.getElementById('outDir').value = data.path;
}

async function loadDem() {
  try {
    if (!val('demPath').trim()) {
      setStatus('fileStatus', 'Seleccione o pegue la ruta del DEM antes de leerlo.', 'bad');
      return;
    }
    if (demPollTimer) clearInterval(demPollTimer);
    document.getElementById('loadDemButton').disabled = true;
    demFetchFailures = 0;
    setStatus('fileStatus', '0% | En cola', 'warn');
    let data;
    try {
      data = await api('/api/start-dem-info', demGridPayload({path: val('demPath')}));
    } catch (err) {
      setStatus('fileStatus', `${err.message}\nIntentando lectura directa sin barra de avance...`, 'warn');
      await loadDemDirectFallback();
      return;
    }
    currentDemJob = data.job_id;
    demPollTimer = setInterval(pollDemInfo, 500);
    await pollDemInfo();
  } catch (err) {
    document.getElementById('loadDemButton').disabled = false;
    setStatus('fileStatus', err.message, 'bad');
  }
}

async function loadDemDirectFallback() {
  try {
    const info = await api('/api/dem-info', demGridPayload({path: val('demPath')}));
    document.getElementById('loadDemButton').disabled = false;
    setStatus('fileStatus', `DEM listo\n${formatDemInfo(info)}`, 'ok');
  } catch (err) {
    document.getElementById('loadDemButton').disabled = false;
    setStatus('fileStatus', err.message, 'bad');
  }
}

async function pollDemInfo() {
  if (!currentDemJob) return;
  try {
    const data = await api('/api/dem-info-status', {job_id: currentDemJob});
    demFetchFailures = 0;
    const progress = Number(data.progress || 0);
    const message = data.message || '';
    const cls = data.status === 'failed' ? 'bad' : (data.status === 'completed' ? 'ok' : 'warn');
    if (data.status === 'completed' && data.info) {
      if (demPollTimer) clearInterval(demPollTimer);
      document.getElementById('loadDemButton').disabled = false;
      setStatus('fileStatus', `${progress}% | ${message}\n${formatDemInfo(data.info)}`, cls);
    } else if (data.status === 'failed') {
      if (demPollTimer) clearInterval(demPollTimer);
      document.getElementById('loadDemButton').disabled = false;
      setStatus('fileStatus', `${progress}% | ${message}`, cls);
    } else {
      setStatus('fileStatus', `${progress}% | ${message}`, cls);
    }
  } catch (err) {
    demFetchFailures += 1;
    if (demFetchFailures <= 6) {
      setStatus('fileStatus', `Reconectando con backend local (${demFetchFailures}/6)...\n${err.message}`, 'warn');
      return;
    }
    if (demPollTimer) clearInterval(demPollTimer);
    document.getElementById('loadDemButton').disabled = false;
    setStatus('fileStatus', `${err.message}\nSi el proceso de Python se cerro, pruebe una resolucion objetivo mas gruesa o recorte el DEM.`, 'bad');
  }
}

function formatDemInfo(info) {
  const source = `${info.source_width}x${info.source_height} @ ${info.source_dx}x${info.source_dy}`;
  const method = info.resampling_method ? ` | ${info.resampling_method} objetivo=${info.target_resolution_m}` : '';
  const interpolated = info.interpolated_dem_path ? `\nInterpolado: ${info.interpolated_dem_path}` : '';
  return `DEM ${info.width}x${info.height} | dx=${info.dx} dy=${info.dy} | fuente=${source}${method} | CRS=${info.crs || 'sin CRS'} | z=${info.z_min.toFixed(2)}..${info.z_max.toFixed(2)}${interpolated}`;
}

async function loadLayers() {
  try {
    const data = await api('/api/vector-layers', {path: val('vectorPath')});
    layers = data.layers || [];
    const select = document.getElementById('layerName');
    select.innerHTML = layers.map(l => `<option value="${escapeHtml(l.name)}">${escapeHtml(l.name)} (${(l.fields||[]).join(', ')})</option>`).join('');
    setStatus('fileStatus', `${layers.length} capa(s) vectoriales detectadas`, 'ok');
  } catch (err) {
    setStatus('fileStatus', err.message, 'bad');
  }
}

async function loadSites() {
  try {
    const data = await api('/api/sites', {
      dem_path: val('demPath'),
      vector_path: val('vectorPath'),
      layer_name: val('layerName'),
      name_field: val('nameField'),
      decimate: num('decimate'),
      target_resolution_m: optNum('targetResolution'),
      resampling_method: val('resamplingMethod')
    });
    sites = data.sites || [];
    selected = new Set(sites.filter(s => s.in_dem).map((_, i) => String(i)));
    renderSites();
    setStatus('fileStatus', `${sites.length} sitio(s) cargados`, 'ok');
  } catch (err) {
    setStatus('fileStatus', err.message, 'bad');
  }
}

function renderSites() {
  const filter = val('siteFilter').toLowerCase();
  const rows = sites.map((s, i) => ({s, i})).filter(x => String(x.s.name).toLowerCase().includes(filter));
  document.getElementById('sitesTable').innerHTML =
    '<tr><th></th><th>Name</th><th>Row</th><th>Col</th><th>Estado</th></tr>' +
    rows.map(({s, i}) => `<tr>
      <td><input type="checkbox" ${selected.has(String(i)) ? 'checked' : ''} ${s.in_dem ? '' : 'disabled'} onchange="toggleSite(${i}, this.checked)"></td>
      <td>${escapeHtml(s.name)}</td><td>${s.row}</td><td>${s.col}</td>
      <td class="${s.in_dem ? 'ok' : 'bad'}">${s.in_dem ? 'dentro DEM' : 'fuera DEM'}</td>
    </tr>`).join('');
  document.getElementById('siteStatus').textContent = `${selected.size} sitio(s) seleccionados`;
}

function toggleSite(i, isChecked) {
  if (isChecked) selected.add(String(i)); else selected.delete(String(i));
  renderSites();
}

function setAllSites(isChecked) {
  const filter = val('siteFilter').toLowerCase();
  sites.forEach((s, i) => {
    if (s.in_dem && String(s.name).toLowerCase().includes(filter)) {
      if (isChecked) selected.add(String(i)); else selected.delete(String(i));
    }
  });
  renderSites();
}

function toggleCurveMode() {
  const csv = val('curveMode') === 'csv';
  document.getElementById('constantCurve').style.display = csv ? 'none' : 'block';
  document.getElementById('csvCurve').style.display = csv ? 'block' : 'none';
}

function sourceCurvePayload() {
  if (val('curveMode') === 'csv') return {mode: 'csv', csv_path: val('csvPath')};
  return {mode: 'constant', volume_total_m3: num('volumeTotal'), duration_s: num('durationS')};
}

function demGridPayload(extra = {}) {
  return {
    decimate: num('decimate'),
    target_resolution_m: optNum('targetResolution'),
    resampling_method: val('resamplingMethod'),
    ...extra
  };
}

function modelParams() {
  return {
    t_end: num('tEnd'), cfl: num('cfl'), roi_buffer_m: num('roiBuffer'),
    friction_model: val('frictionModel'), manning_n: num('manningN'), darcy_f: num('darcyF'),
    fluid_density: num('fluidDensity'), dynamic_viscosity: num('dynamicViscosity'),
    yield_stress: 0, h_viscous_min: num('hViscousMin'),
    k_oil: num('kOil'), capillary_suction: num('capillarySuction'), porosity_deficit: num('porosityDeficit'),
    evap_rate: num('evapRate'), deg_rate: num('degRate'), snapshot_interval: num('snapshotInterval'),
    video_fps: num('videoFps'), video_3d_palette: val('video3dPalette'), isochrone_interval_s: num('isochroneInterval')
  };
}

async function startRun() {
  try {
    document.getElementById('runButton').disabled = true;
    const selectedSites = [...selected].map(i => sites[Number(i)]);
    const payload = {
      dem_path: val('demPath'),
      vector_path: val('vectorPath'),
      out_dir: val('outDir'),
      decimate: num('decimate'),
      target_resolution_m: optNum('targetResolution'),
      resampling_method: val('resamplingMethod'),
      max_workers: num('maxWorkers'),
      sites: selectedSites,
      source_curve: sourceCurvePayload(),
      model_params: modelParams(),
      output_options: {
        arrays: checked('arrays'), maps: checked('maps'),
        diagnostic_report: checked('diagnosticReport'), videos: checked('videos'),
        rasters: checked('rasters'), isochrones: checked('isochrones'), snapshots: checked('snapshots')
      }
    };
    const data = await api('/api/start-run', payload);
    currentJob = data.job_id;
    setStatus('runStatus', `Job ${currentJob} iniciado`);
    pollTimer = setInterval(pollStatus, 1500);
    await pollStatus();
  } catch (err) {
    document.getElementById('runButton').disabled = false;
    setStatus('runStatus', err.message, 'bad');
  }
}

async function pollStatus() {
  if (!currentJob) return;
  const data = await api('/api/status', {job_id: currentJob});
  const total = data.total || 0;
  const completed = data.completed || 0;
  const runClass = data.status === 'failed' ? 'bad' : (data.status === 'completed' ? 'ok' : 'warn');
  const records = data.records || [];
  const snapshotMb = totalSnapshotsMb(records);
  const snapshotText = snapshotMb > 0 ? `\nSnapshots guardados: ${snapshotMb.toFixed(2)} MB` : '';
  setStatus('runStatus', `${data.status}: ${data.message || ''}\n${completed}/${total} completadas${snapshotText}${data.index_path ? '\nIndex: ' + data.index_path : ''}`, runClass);
  renderResults(records);
  if (data.status === 'completed' || data.status === 'failed') {
    clearInterval(pollTimer);
    document.getElementById('runButton').disabled = false;
  }
}

function renderResults(records) {
  document.getElementById('resultsTable').innerHTML =
    '<tr><th>Sitio</th><th>Estado</th><th>Mensaje</th><th>Tiempo (s)</th><th>Salida</th><th>Espesor max raster</th><th>Isocronas por intervalo</th><th>Snapshots</th><th>Video</th><th>Log</th><th>Error/Aviso</th></tr>' +
    records.map(r => `<tr>
      <td>${escapeHtml((r.site && r.site.name) || '')}</td>
      <td class="${statusClass(r.status)}">${escapeHtml(r.status || '')}</td>
      <td>${escapeHtml(r.message || '')}</td>
      <td>${formatSeconds(r.elapsed_s)}</td>
      <td>${escapeHtml(r.out_dir || '')}</td>
      <td>${escapeHtml(r.max_thickness_raster_path || '')}</td>
      <td>${escapeHtml(r.isochrones_path || '')}</td>
      <td>${formatMb(r.snapshots_size_mb)}</td>
      <td>${videoButtonHtml(r)}</td>
      <td>${escapeHtml(r.log_path || '')}</td>
      <td>${escapeHtml(r.error || r.isochrones_error || boundaryWarning(r) || '')}</td>
    </tr>`).join('');
}

function videoButtonHtml(record) {
  if (record.status !== 'completed' || !record.snapshots_path) return '';
  return `<button class="secondary" onclick="startVideoExport('${escapeJsArg(record.out_dir || '')}')">Reexportar</button>`;
}

async function startVideoExport(outDir) {
  try {
    if (!outDir) throw new Error('No hay carpeta de salida para reexportar.');
    if (videoPollTimer) clearInterval(videoPollTimer);
    const palette = val('reexport3dPalette');
    setStatus('videoExportStatus', `En cola\nPaleta derrame 3D: ${paletteLabel(palette)}`, 'warn');
    const data = await api('/api/start-video-export', {
      out_dir: outDir,
      fps: num('reexportFps'),
      t_end_s: optNum('videoCutoffS'),
      video_2d: checked('reexport2d'),
      video_3d: checked('reexport3d'),
      video_3d_palette: palette
    });
    currentVideoJob = data.job_id;
    videoPollTimer = setInterval(pollVideoExport, 1000);
    await pollVideoExport();
  } catch (err) {
    setStatus('videoExportStatus', err.message, 'bad');
  }
}

async function pollVideoExport() {
  if (!currentVideoJob) return;
  try {
    const data = await api('/api/video-export-status', {job_id: currentVideoJob});
    const cls = data.status === 'failed' ? 'bad' : (data.status === 'completed' ? 'ok' : 'warn');
    const paths = [
      data.video_2d_path ? `2D: ${data.video_2d_path}` : '',
      data.video_3d_path ? `3D: ${data.video_3d_path}` : ''
    ].filter(Boolean).join('\n');
    const palette = data.video_3d_palette ? `\nPaleta derrame 3D: ${paletteLabel(data.video_3d_palette)}` : '';
    setStatus('videoExportStatus', `${data.status}: ${data.message || ''}${palette}${paths ? '\n' + paths : ''}`, cls);
    if (data.status === 'completed' || data.status === 'failed') clearInterval(videoPollTimer);
  } catch (err) {
    if (videoPollTimer) clearInterval(videoPollTimer);
    setStatus('videoExportStatus', err.message, 'bad');
  }
}

function boundaryWarning(record) {
  return record.domain_boundary_reached ? 'El derrame toco el borde del recorte; aumente ROI buffer (m).' : '';
}

function formatSeconds(value) {
  const number = Number(value);
  return Number.isFinite(number) ? number.toFixed(1) : '';
}

function formatMb(value) {
  const number = Number(value);
  return Number.isFinite(number) && number > 0 ? `${number.toFixed(2)} MB` : '';
}

function totalSnapshotsMb(records) {
  return records.reduce((sum, record) => {
    const value = Number(record.snapshots_size_mb);
    return Number.isFinite(value) ? sum + value : sum;
  }, 0);
}

function statusClass(status) {
  if (status === 'completed') return 'ok';
  if (status === 'failed') return 'bad';
  return 'warn';
}

function escapeHtml(text) {
  return String(text ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'}[c]));
}

function escapeJsArg(text) {
  return String(text ?? '').replace(/\\/g, '\\\\').replace(/'/g, "\\'").replace(/\n/g, ' ');
}

toggleCurveMode();
syncReexportPalette();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_server(args.host, args.port, open_browser=not args.no_browser)
