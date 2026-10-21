# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Band recognition from file names and metadata (no files needed).

    python -m flooddepth.tests.test_sensors
"""

from .common import check, finish
from ..core import sensors as S


def info(path, count=1, alpha=None, descs=None, dtype="Byte", md=None, bmd=None):
    bands = []
    for i in range(1, count + 1):
        bands.append({"index": i, "description": (descs or [""] * count)[i - 1], "dtype": dtype,
                      "colorinterp": "", "metadata": (bmd or [{}] * count)[i - 1], "alpha": i == alpha})
    return {"path": path, "count": count, "alpha": alpha, "bands": bands, "metadata": md or {}, "grid": None}


def ident(paths, infos, overrides=None):
    lut = dict(zip(paths, infos))
    return S.identify(paths, lambda p: lut[p], overrides)


def test_names():
    print("\n1. File names")
    b = "2024-05-05-00_00_2024-05-05-23_59_Sentinel-2_L2A_{}_(Raw).tiff"
    ps = [b.format(x) for x in ("B03", "B08", "B11", "B8A", "B02")]
    sc = ident(ps, [info(p, 2, alpha=2) for p in ps])
    check("Browser Raw: green = B03", sc.sources["green"][0].endswith("B03_(Raw).tiff"))
    check("Browser Raw: NIR = B08, not B8A", "B08" in sc.sources["nir"][0])
    check("Browser Raw: SWIR1 = B11", "B11" in sc.sources["swir1"][0])
    check("kind = bands, family s2, date read", sc.kind == "bands" and sc.family == "s2" and sc.date == "2024-05-05")
    ps = ["T36MZE_20240505T074611_B03_10m.jp2", "T36MZE_20240505T074611_B8A_20m.jp2",
          "T36MZE_20240505T074611_B11_20m.jp2"]
    sc = ident(ps, [info(p, dtype="UInt16") for p in ps])
    check("SAFE: B8A used as NIR when B08 is absent", "B8A" in sc.sources["nir"][0])
    ps = ["LC09_L2SP_169060_20240506_20240507_02_T1_SR_B3.TIF", "LC09_L2SP_169060_20240506_20240507_02_T1_SR_B5.TIF",
          "LC09_L2SP_169060_20240506_20240507_02_T1_SR_B6.TIF"]
    sc = ident(ps, [info(p, dtype="UInt16") for p in ps])
    check("Landsat 9: B3 green, B5 NIR, B6 SWIR1", "SR_B3" in sc.sources["green"][0] and "SR_B5" in sc.sources["nir"][0]
          and "SR_B6" in sc.sources["swir1"][0])
    ps = ["LE07_L2SP_169060_20100506_SR_B2.TIF", "LE07_L2SP_169060_20100506_SR_B4.TIF"]
    sc = ident(ps, [info(p, dtype="UInt16") for p in ps])
    check("Landsat 7: B2 green, B4 NIR", "SR_B2" in sc.sources["green"][0] and "SR_B4" in sc.sources["nir"][0])


def test_composites():
    print("\n2. Rendered products")
    p = "2024-05-05-00_00_2024-05-05-23_59_Sentinel-2_L2A_False_color.tiff"
    sc = ident([p], [info(p, 3)])
    check("False color: NIR = channel 1, green = channel 3", sc.sources["nir"] == (p, 1) and sc.sources["green"] == (p, 3))
    p = "2024-05-05-00_00_2024-05-05-23_59_Sentinel-2_L2A_True_color.tiff"
    sc = ident([p], [info(p, 3)])
    check("True color: no index bands, kind rgb", sc.kind == "rgb" and not sc.has("nir"))
    p = "2024-05-06-00_00_2024-05-06-23_59_Sentinel-1_AWS-IW-VVVH_Enhanced_visualization.tiff"
    sc = ident([p], [info(p, 3, dtype="UInt16")])
    check("rendered Sentinel-1 refused (recorded)", sc.rendered_sar is not None and sc.kind == "unknown")
    p = "S1A_IW_GRDH_20240506_VV.tif"
    sc = ident([p], [info(p, 1, dtype="Float32")])
    check("Sentinel-1 VV file", sc.kind == "sar" and "vv" in sc.sources)


def test_metadata():
    print("\n3. Band metadata and overrides")
    p = "stack.tif"
    sc = ident([p], [info(p, 4, descs=["Blue", "Green", "Red", "NIR"])])
    check("descriptions: green = 2, NIR = 4", sc.sources["green"] == (p, 2) and sc.sources["nir"] == (p, 4))
    sc = ident([p], [info(p, 3, bmd=[{"WAVELENGTH": "560"}, {"WAVELENGTH": "0.865", "WAVELENGTH_UNITS": "micrometers"},
                                     {"WAVELENGTH": "1610"}])])
    check("wavelengths: green, NIR (um), SWIR1", sc.sources["green"] == (p, 1) and sc.sources["nir"] == (p, 2)
          and sc.sources["swir1"] == (p, 3))
    p2 = "S2B_MSIL2A_20240505_stack.tif"
    sc = ident([p2], [info(p2, 12, dtype="UInt16")])
    check("12-band file named Sentinel-2: B03 = 3, B08 = 8, B11 = 11", sc.sources["green"] == (p2, 3)
          and sc.sources["nir"] == (p2, 8) and sc.sources["swir1"] == (p2, 11))
    p3 = "unnamed.tif"
    sc = ident([p3], [info(p3, 6)])
    check("unnamed 6-band stack: nothing guessed", sc.kind == "unknown" and not sc.sources)
    sc = ident([p3], [info(p3, 6)], {"green": 2, "nir": 4})
    check("overrides win", sc.sources["green"] == (p3, 2) and sc.kind == "bands")
    try:
        ident([p3], [info(p3, 6)], {"nir": 9})
        check("override beyond the band count refused", False)
    except ValueError:
        check("override beyond the band count refused", True)
    sc = ident([p], [info(p, 4, descs=["Blue", "Green", "Red", "NIR"], dtype="UInt16",
                          md={"BOA_ADD_OFFSET_B3": "-1000"})])
    check("BOA_ADD_OFFSET metadata applied", sc.offset["green"] == -1000.0)


def test_scale():
    print("\n4. Reflectance scale")
    check("Byte -> /255", abs(S.reflectance_scale("Byte", 200) - 1 / 255) < 1e-12)
    check("UInt16 0-10000 -> /10000", S.reflectance_scale("UInt16", 6000) == 1e-4)
    check("UInt16 Browser 0-65535 -> /65535", abs(S.reflectance_scale("UInt16", 25000) - 1 / 65535) < 1e-12)
    check("Float reflectance -> 1", S.reflectance_scale("Float32", 0.9) == 1.0)


def test_scene_scale():
    """A x 65535 export whose blue band is dark (dry soil) must not be scaled
    band by band: blue alone looks like x 10000 and clear land became 'cloud'
    (Juba, Dec 2023: 61 % of a clear scene masked)."""
    print("\n5. Scene-wide reflectance scale (bright-cloud test)")
    import os
    import tempfile
    import numpy as np
    from osgeo import gdal
    from ..core import imagery as I, raster as R
    d = tempfile.mkdtemp()
    rng = np.random.default_rng(1)
    vals = {"B02": 2200, "B03": 4200, "B04": 3000, "B08": 20000, "B11": 14000}
    paths = []
    for b, v in vals.items():
        p = os.path.join(d, f"2023-12-03-00_00_2023-12-03-23_59_Sentinel-2_L2A_{b}_(Raw).tiff")
        ds = gdal.GetDriverByName("GTiff").Create(p, 60, 60, 1, gdal.GDT_UInt16)
        ds.SetGeoTransform((500000, 10, 0, 100000, 0, -10))
        sr = R.srs_from_epsg(32638)
        ds.SetProjection(sr.ExportToWkt())
        ds.GetRasterBand(1).WriteArray((v * rng.uniform(0.5, 1.6, (60, 60))).astype("uint16"))
        ds = None
        paths.append(p)
    sc = I.scene_from_paths(paths)
    g = R.grid_of(R.open_dataset(paths[0]))
    cl = I.bright_cloud(sc, g, 0.25, [])
    check("blue in a x 65535 export is not cloud", cl is not None and cl.mean() == 0.0, f"{cl.mean():.2f}")
    check("scene p99 taken from the brightest band", I.scene_p99(sc) > 12000)


def test_sar_uint16():
    """Browser 'Raw' Sentinel-1 in UINT16 is linear backscatter x 65535; it was
    read as dB, so nothing looked dark and no water was found (Nyando, 6 May 2024)."""
    print("\n6. Sentinel-1 UINT16 export (linear x 65535)")
    import os
    import tempfile
    import numpy as np
    from osgeo import gdal
    from ..core import imagery as I, raster as R
    d = tempfile.mkdtemp()
    rng = np.random.default_rng(2)
    water = np.zeros((80, 80), bool)
    water[20:50, 10:40] = True
    paths = []
    for pol, (w, l) in (("VV", (0.01, 0.15)), ("VH", (0.002, 0.03))):
        lin = np.where(water, w, l) * rng.gamma(4.0, 0.25, water.shape)
        p = os.path.join(d, f"2024-05-06-00_00_2024-05-06-23_59_Sentinel-1_AWS-IW-VVVH_{pol}_(Raw).tiff")
        ds = gdal.GetDriverByName("GTiff").Create(p, 80, 80, 2, gdal.GDT_UInt16)
        ds.SetGeoTransform((34.5, 0.0001, 0, -0.2, 0, -0.0001))
        ds.SetProjection(R.srs_from_epsg(4326).ExportToWkt())
        ds.GetRasterBand(1).WriteArray(np.clip(lin * 65535, 0, 65535).astype("uint16"))
        ds.GetRasterBand(2).WriteArray(np.ones((80, 80), "uint16"))
        ds = None
        paths.append(p)
    r = I.classify_image(paths)
    t = r["info"]["threshold"]["threshold"]
    share = r["info"]["water_share"]
    check("threshold in the SAR water range (dB)", -26 < t < -10, f"{t:.1f} dB")
    check("the water block is found", 0.1 < share < 0.2, f"{100 * share:.1f} % (true 14 %)")


if __name__ == "__main__":
    test_names()
    test_composites()
    test_metadata()
    test_scale()
    test_scene_scale()
    test_sar_uint16()
    finish("test_sensors")
