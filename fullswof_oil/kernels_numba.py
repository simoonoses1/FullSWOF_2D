from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange

    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover - optional dependency
    NUMBA_AVAILABLE = False


EPS = 1e-10


if NUMBA_AVAILABLE:

    @njit(parallel=True, fastmath=False)
    def max_wave_speed_numba(h: np.ndarray, hu: np.ndarray, hv: np.ndarray, g: float, eps: float = EPS) -> float:
        ny, nx = h.shape
        row_max = np.zeros(ny, dtype=np.float64)
        for i in prange(ny):
            local_max = 0.0
            for j in range(nx):
                hij = h[i, j]
                u = hu[i, j] / hij if hij > eps else 0.0
                v = hv[i, j] / hij if hij > eps else 0.0
                c = np.sqrt(g * hij) if hij > 0.0 else 0.0
                s = max(abs(u) + c, abs(v) + c)
                if s > local_max:
                    local_max = s
            row_max[i] = local_max
        out = 0.0
        for i in range(ny):
            if row_max[i] > out:
                out = row_max[i]
        return out


    @njit(parallel=True, fastmath=False)
    def conservative_update_numba(
        h: np.ndarray,
        hu: np.ndarray,
        hv: np.ndarray,
        z: np.ndarray,
        dt: float,
        dx: float,
        dy: float,
        g: float,
        eps: float = EPS,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ny, nx = h.shape
        hp = np.zeros((ny + 2, nx + 2), dtype=np.float64)
        hup = np.zeros((ny + 2, nx + 2), dtype=np.float64)
        hvp = np.zeros((ny + 2, nx + 2), dtype=np.float64)
        zp = np.zeros((ny + 2, nx + 2), dtype=np.float64)

        for i in prange(ny):
            for j in range(nx):
                hp[i + 1, j + 1] = h[i, j]
                hup[i + 1, j + 1] = hu[i, j]
                hvp[i + 1, j + 1] = hv[i, j]
                zp[i + 1, j + 1] = z[i, j]

        for i in prange(ny):
            hp[i + 1, 0] = h[i, 0]
            hp[i + 1, nx + 1] = h[i, nx - 1]
            hup[i + 1, 0] = -hu[i, 0]
            hup[i + 1, nx + 1] = -hu[i, nx - 1]
            hvp[i + 1, 0] = hv[i, 0]
            hvp[i + 1, nx + 1] = hv[i, nx - 1]
            zp[i + 1, 0] = z[i, 0]
            zp[i + 1, nx + 1] = z[i, nx - 1]

        for j in prange(nx + 2):
            hp[0, j] = hp[1, j]
            hp[ny + 1, j] = hp[ny, j]
            hup[0, j] = hup[1, j]
            hup[ny + 1, j] = hup[ny, j]
            hvp[0, j] = -hvp[1, j]
            hvp[ny + 1, j] = -hvp[ny, j]
            zp[0, j] = zp[1, j]
            zp[ny + 1, j] = zp[ny, j]

        fx_h = np.zeros((ny, nx + 1), dtype=np.float64)
        fx_hu = np.zeros((ny, nx + 1), dtype=np.float64)
        fx_hv = np.zeros((ny, nx + 1), dtype=np.float64)
        sx = np.zeros((ny, nx + 1), dtype=np.float64)

        for i in prange(ny):
            ii = i + 1
            for j in range(nx + 1):
                h_l = hp[ii, j]
                h_r = hp[ii, j + 1]
                hu_l = hup[ii, j]
                hu_r = hup[ii, j + 1]
                hv_l = hvp[ii, j]
                hv_r = hvp[ii, j + 1]
                z_l = zp[ii, j]
                z_r = zp[ii, j + 1]

                dz = z_r - z_l
                h_l_hr = max(0.0, h_l - max(0.0, dz))
                h_r_hr = max(0.0, h_r - max(0.0, -dz))

                u_l = hu_l / h_l if h_l > eps else 0.0
                u_r = hu_r / h_r if h_r > eps else 0.0
                v_l = hv_l / h_l if h_l > eps else 0.0
                v_r = hv_r / h_r if h_r > eps else 0.0

                hl_hul = h_l_hr * u_l
                hl_hvl = h_l_hr * v_l
                hr_hur = h_r_hr * u_r
                hr_hvr = h_r_hr * v_r

                fl_h = hl_hul
                fl_hu = hl_hul * u_l + 0.5 * g * h_l_hr * h_l_hr
                fl_hv = hl_hul * v_l

                fr_h = hr_hur
                fr_hu = hr_hur * u_r + 0.5 * g * h_r_hr * h_r_hr
                fr_hv = hr_hur * v_r

                c_l = np.sqrt(g * h_l_hr) if h_l_hr > 0.0 else 0.0
                c_r = np.sqrt(g * h_r_hr) if h_r_hr > 0.0 else 0.0
                smax = max(abs(u_l) + c_l, abs(u_r) + c_r)

                fx_h[i, j] = 0.5 * (fl_h + fr_h) - 0.5 * smax * (h_r_hr - h_l_hr)
                fx_hu[i, j] = 0.5 * (fl_hu + fr_hu) - 0.5 * smax * (hr_hur - hl_hul)
                fx_hv[i, j] = 0.5 * (fl_hv + fr_hv) - 0.5 * smax * (hr_hvr - hl_hvl)

                sx[i, j] = 0.5 * g * (
                    (h_l_hr - h_l) * (h_l_hr + h_l)
                    + (h_r - h_r_hr) * (h_r + h_r_hr)
                    + (h_l + h_r) * dz
                )

        fy_h = np.zeros((ny + 1, nx), dtype=np.float64)
        fy_hu = np.zeros((ny + 1, nx), dtype=np.float64)
        fy_hv = np.zeros((ny + 1, nx), dtype=np.float64)
        sy = np.zeros((ny + 1, nx), dtype=np.float64)

        for i in prange(ny + 1):
            for j in range(nx):
                jj = j + 1
                h_b = hp[i, jj]
                h_t = hp[i + 1, jj]
                hu_b = hup[i, jj]
                hu_t = hup[i + 1, jj]
                hv_b = hvp[i, jj]
                hv_t = hvp[i + 1, jj]
                z_b = zp[i, jj]
                z_t = zp[i + 1, jj]

                dz = z_t - z_b
                h_b_hr = max(0.0, h_b - max(0.0, dz))
                h_t_hr = max(0.0, h_t - max(0.0, -dz))

                u_b = hu_b / h_b if h_b > eps else 0.0
                u_t = hu_t / h_t if h_t > eps else 0.0
                v_b = hv_b / h_b if h_b > eps else 0.0
                v_t = hv_t / h_t if h_t > eps else 0.0

                hb_hub = h_b_hr * u_b
                hb_hvb = h_b_hr * v_b
                ht_hut = h_t_hr * u_t
                ht_hvt = h_t_hr * v_t

                fl_h = hb_hvb
                fl_hu = hb_hvb * u_b
                fl_hv = hb_hvb * v_b + 0.5 * g * h_b_hr * h_b_hr

                fr_h = ht_hvt
                fr_hu = ht_hvt * u_t
                fr_hv = ht_hvt * v_t + 0.5 * g * h_t_hr * h_t_hr

                c_b = np.sqrt(g * h_b_hr) if h_b_hr > 0.0 else 0.0
                c_t = np.sqrt(g * h_t_hr) if h_t_hr > 0.0 else 0.0
                smax = max(abs(v_b) + c_b, abs(v_t) + c_t)

                fy_h[i, j] = 0.5 * (fl_h + fr_h) - 0.5 * smax * (h_t_hr - h_b_hr)
                fy_hu[i, j] = 0.5 * (fl_hu + fr_hu) - 0.5 * smax * (ht_hut - hb_hub)
                fy_hv[i, j] = 0.5 * (fl_hv + fr_hv) - 0.5 * smax * (ht_hvt - hb_hvb)

                sy[i, j] = 0.5 * g * (
                    (h_b_hr - h_b) * (h_b_hr + h_b)
                    + (h_t - h_t_hr) * (h_t + h_t_hr)
                    + (h_b + h_t) * dz
                )

        h_new = h.copy()
        hu_new = hu.copy()
        hv_new = hv.copy()
        cdx = dt / dx
        cdy = dt / dy
        for i in prange(ny):
            for j in range(nx):
                h_new[i, j] -= cdx * (fx_h[i, j + 1] - fx_h[i, j]) + cdy * (fy_h[i + 1, j] - fy_h[i, j])
                hu_new[i, j] -= cdx * (fx_hu[i, j + 1] - fx_hu[i, j]) + cdy * (fy_hu[i + 1, j] - fy_hu[i, j])
                hv_new[i, j] -= cdx * (fx_hv[i, j + 1] - fx_hv[i, j]) + cdy * (fy_hv[i + 1, j] - fy_hv[i, j])

                hu_new[i, j] -= cdx * (sx[i, j + 1] - sx[i, j])
                hv_new[i, j] -= cdy * (sy[i + 1, j] - sy[i, j])

                if h_new[i, j] < 0.0:
                    h_new[i, j] = 0.0
                if h_new[i, j] < eps:
                    hu_new[i, j] = 0.0
                    hv_new[i, j] = 0.0

        return h_new, hu_new, hv_new

else:

    def max_wave_speed_numba(h: np.ndarray, hu: np.ndarray, hv: np.ndarray, g: float, eps: float = EPS) -> float:  # pragma: no cover
        raise RuntimeError("numba is not available")

    def conservative_update_numba(
        h: np.ndarray,
        hu: np.ndarray,
        hv: np.ndarray,
        z: np.ndarray,
        dt: float,
        dx: float,
        dy: float,
        g: float,
        eps: float = EPS,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # pragma: no cover
        raise RuntimeError("numba is not available")
