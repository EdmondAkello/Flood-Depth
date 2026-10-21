# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Optional accelerators with pure NumPy / GDAL fallbacks.

FloodDepth never requires anything beyond NumPy and GDAL. When SciPy is
importable (it ships with most QGIS installs), a few kernels use it:

    label           scipy.ndimage.label          (else gdal.Polygonize + Rasterize)
    box_sum         scipy.ndimage.uniform_filter (else cumulative sums)
    nearest index   scipy.ndimage.distance_transform_edt (else jump flooding)
    FFT             scipy.fft with all cores     (else numpy.fft)

Both paths are tested against each other (tests/test_accel.py). Set the
environment variable FLOODDEPTH_NO_SCIPY=1 to force the fallbacks.

Threads: the uncertainty ensemble runs its members in a thread pool. GDAL,
NumPy and SciPy release the GIL in the heavy calls, so members overlap on
multi-core machines with no extra process (nothing is spawned).
"""

import os

import numpy as np

_SCIPY = None


def scipy_ok():
    """True when SciPy (ndimage + fft) is usable and not disabled."""
    global _SCIPY
    if os.environ.get("FLOODDEPTH_NO_SCIPY", "").strip() not in ("", "0"):
        return False
    if _SCIPY is None:
        try:
            import scipy.fft  # noqa: F401
            import scipy.ndimage  # noqa: F401
            _SCIPY = True
        except Exception:  # missing, or a broken install
            _SCIPY = False
    return _SCIPY


def describe():
    if scipy_ok():
        import scipy
        return f"SciPy {scipy.__version__} (optional accelerator) used"
    if os.environ.get("FLOODDEPTH_NO_SCIPY", "").strip() not in ("", "0"):
        return "SciPy disabled (FLOODDEPTH_NO_SCIPY): NumPy / GDAL paths used"
    return "SciPy not available: NumPy / GDAL paths used"


def workers(n_tasks=None, cap=8):
    """Thread count for parallel work: the CPU count (minus one, kept free for
    QGIS), at least 1, at most `cap` and `n_tasks`. FLOODDEPTH_THREADS overrides."""
    env = os.environ.get("FLOODDEPTH_THREADS", "").strip()
    if env.isdigit() and int(env) > 0:
        n = int(env)
    else:
        n = max(1, (os.cpu_count() or 1) - 1)
    n = min(n, cap)
    if n_tasks is not None:
        n = min(n, max(1, int(n_tasks)))
    return n


# ------------------------------------------------------------------ FFT

def _fast_len(n):
    """Smallest 2-3-5-smooth integer >= n (FFT-friendly)."""
    if scipy_ok():
        from scipy.fft import next_fast_len
        return next_fast_len(int(n), real=True)
    m = int(n)
    while True:
        k = m
        for p in (2, 3, 5):
            while k % p == 0:
                k //= p
        if k == 1:
            return m
        m += 1


def fast_shape(shape):
    return tuple(_fast_len(s) for s in shape)


def rfft2(a):
    if scipy_ok():
        import scipy.fft as sf
        return sf.rfft2(a, workers=-1)
    return np.fft.rfft2(a)


def irfft2(a, s):
    if scipy_ok():
        import scipy.fft as sf
        return sf.irfft2(a, s=s, workers=-1)
    return np.fft.irfft2(a, s=s)


# ------------------------------------------------------------------ ndimage

def box_sum(a, r):
    """Sum over a (2r+1)^2 window with zero padding, or None (use the fallback)."""
    if not scipy_ok():
        return None
    from scipy.ndimage import uniform_filter
    k = 2 * int(r) + 1
    return uniform_filter(np.asarray(a, dtype="float64"), size=k, mode="constant", cval=0.0) * float(k * k)


def label(mask, connectivity=8):
    """(labels int32, count) or None (use the fallback)."""
    if not scipy_ok():
        return None
    from scipy.ndimage import label as _label
    st = np.ones((3, 3), bool) if connectivity == 8 else np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)
    lab, n = _label(np.asarray(mask, bool), structure=st)
    return lab.astype("int32"), int(n)


def nearest_index(seeds, cell_x=1.0, cell_y=1.0):
    """(rows, cols, distance in map units) of the nearest seed for every cell,
    exact Euclidean (anisotropic cells), or None (use the fallback)."""
    if not scipy_ok():
        return None
    from scipy.ndimage import distance_transform_edt
    seeds = np.asarray(seeds, bool)
    if not seeds.any():
        return None
    d, (ri, ci) = distance_transform_edt(~seeds, sampling=(float(cell_y), float(cell_x)), return_indices=True)
    return ri.astype("int32"), ci.astype("int32"), d.astype("float32")
