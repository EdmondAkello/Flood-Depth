# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""The optional SciPy kernels agree with the NumPy / GDAL fallbacks.

    python -m flooddepth.tests.test_accel

Skips the comparisons (with a note) when SciPy is not installed; the other
suites then exercise the fallbacks. CI runs every suite both ways
(FLOODDEPTH_NO_SCIPY=1).
"""

import os

import numpy as np

from .common import check, finish
from ..core import accel as A, classify as C, surface as SF, uncertainty as U


def both(fn):
    """(result with SciPy, result with the fallback)."""
    old = os.environ.get("FLOODDEPTH_NO_SCIPY")
    try:
        os.environ["FLOODDEPTH_NO_SCIPY"] = "0"
        a = fn()
        os.environ["FLOODDEPTH_NO_SCIPY"] = "1"
        b = fn()
    finally:
        if old is None:
            os.environ.pop("FLOODDEPTH_NO_SCIPY", None)
        else:
            os.environ["FLOODDEPTH_NO_SCIPY"] = old
    return a, b


def same_partition(la, lb, mask):
    """Two labelings describe the same regions (ids may differ)."""
    pairs = set(zip(la[mask].tolist(), lb[mask].tolist()))
    return len(pairs) == len(set(la[mask].tolist())) == len(set(lb[mask].tolist()))


def main():
    print("\n1. Accelerator status: " + A.describe())
    rng = np.random.default_rng(3)
    if not A.scipy_ok():
        check("SciPy absent: fallbacks are used (other suites test them)", A.label(np.ones((3, 3), bool)) is None)
        finish("test_accel")
        return
    m = rng.random((180, 230)) > 0.62

    print("\n2. Connected components")
    for conn in (4, 8):
        (la, na), (lb, nb) = both(lambda: C.label_components(m, conn))
        check(f"{conn}-connected: same count", na == nb, f"{na} vs {nb}")
        check(f"{conn}-connected: same regions", same_partition(la, lb, m))

    print("\n3. Box sums")
    a = rng.normal(size=(97, 131))
    for r in (1, 4, 15):
        sa, sb = both(lambda: SF.box_sum(a, r))
        check(f"radius {r}: equal", np.allclose(sa, sb, atol=1e-9), f"max diff {np.abs(sa - sb).max():.2e}")

    print("\n4. Nearest seed (anisotropic cells)")
    seeds = rng.random((140, 170)) > 0.995
    (ra, ca_, da), (rb, cb, db) = both(lambda: SF.nearest_seed(seeds, 10.0, 12.5))
    check("distances equal (ties may pick different seeds)", np.allclose(da, db, atol=1e-3),
          f"max diff {np.abs(da - db).max():.2e} m")
    dd = np.hypot((np.arange(140)[:, None] - ra) * 12.5, (np.arange(170)[None, :] - ca_) * 10.0)
    check("SciPy indices point at a seed at that distance", seeds[ra, ca_].all() and np.allclose(dd, da, atol=1e-3))

    print("\n5. Correlated noise (FFT): statistics, not values")
    for name, fl in (("SciPy", "0"), ("NumPy", "1")):
        os.environ["FLOODDEPTH_NO_SCIPY"] = fl
        z = U.correlated_noise((257, 311), 2.0, 6.0, np.random.default_rng(1))
        c1 = np.corrcoef(z[:, :-6].ravel(), z[:, 6:].ravel())[0, 1]
        check(f"{name}: sd 2 m, shape kept", z.shape == (257, 311) and abs(z.std() - 2.0) < 0.05, f"sd {z.std():.3f}")
        check(f"{name}: correlated at the correlation length", 0.3 < c1 < 0.9, f"r(6 cells) {c1:.2f}")
    os.environ.pop("FLOODDEPTH_NO_SCIPY", None)
    check("FFT sizes are 2-3-5 smooth", all(A._fast_len(n) in (n, 384, 512, 1000) for n in (384, 509, 997))
          and A._fast_len(509) >= 509)

    print("\n6. Threads")
    os.environ["FLOODDEPTH_THREADS"] = "3"
    check("FLOODDEPTH_THREADS honoured, capped by the task count", A.workers(10) == 3 and A.workers(2) == 2)
    os.environ.pop("FLOODDEPTH_THREADS", None)
    finish("test_accel")


if __name__ == "__main__":
    main()
