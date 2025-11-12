"""
Utility functions for fetch-arcticdem.
"""

import logging
from pathlib import Path

import numpy as np
import rioxarray as rxr  # noqa: F401
import xarray as xr
from osgeo import gdal


def setup_logging(output_dir, log_filename="arcticdem.log"):
    """
    Configure logging to both file and console.

    Parameters
    ----------
    output_dir : str or Path
        Directory where log file will be written
    log_filename : str
        Name of the log file

    Returns
    -------
    logger : logging.Logger
        Configured logger instance
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / log_filename

    logger = logging.getLogger("fetch_arcticdem")
    logger.setLevel(logging.INFO)

    if logger.handlers:
        logger.handlers.clear()

    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def filter_assets(assets):
    """
    Filter out None values from STAC assets dictionary.

    Parameters
    ----------
    assets : dict
        Dictionary of STAC assets

    Returns
    -------
    dict
        Filtered assets dictionary
    """
    return {k: v for k, v in assets.items() if v is not None}


def clean_json_types(value):
    """
    Recursively convert numpy types to native Python types for JSON serialization.

    Parameters
    ----------
    value : any
        Value to convert

    Returns
    -------
    any
        Converted value
    """
    if isinstance(value, (np.float32, np.float64)):
        return float(value)
    if isinstance(value, (np.int32, np.int64)):
        return int(value)
    if isinstance(value, list):
        return [clean_json_types(v) for v in value]
    if isinstance(value, dict):
        return {k: clean_json_types(v) for k, v in value.items()}
    if value is None:
        return None
    return value


def intersection_ratio(item_geom, aoi_geom, aoi_area):
    """
    Calculate the ratio of intersection area to AOI area.

    Parameters
    ----------
    item_geom : shapely.geometry
        Geometry of the STAC item
    aoi_geom : shapely.geometry
        Geometry of the area of interest
    aoi_area : float
        Area of the AOI

    Returns
    -------
    float
        Ratio of intersection area to AOI area
    """
    intersection = item_geom.intersection(aoi_geom)
    if intersection.is_empty:
        return 0.0
    return intersection.area / aoi_area


def gdaldem_hillshade_in_memory(da, z=1.0, s=1.0, subcommand="hillshade"):
    """
    Generate hillshade from DEM DataArray using GDAL in memory.

    Parameters
    ----------
    da : xr.DataArray
        Input DEM
    z : float
        Vertical exaggeration factor
    s : float
        Scale factor
    subcommand : str
        GDAL DEM processing subcommand

    Returns
    -------
    xr.DataArray
        Hillshade raster
    """
    arr = np.array(da, dtype=np.float32)
    ny, nx = arr.shape
    mem_drv = gdal.GetDriverByName("MEM")
    ds = mem_drv.Create("", nx, ny, 1, gdal.GDT_Float32)
    affine = da.rio.transform()
    geotransform = (affine.c, affine.a, affine.b, affine.f, affine.d, affine.e)
    ds.SetGeoTransform(geotransform)
    if hasattr(da.rio, "crs") and da.rio.crs is not None:
        ds.SetProjection(da.rio.crs.to_wkt())
    ds.GetRasterBand(1).WriteArray(arr)
    if hasattr(da.rio, "nodata") and da.rio.nodata is not None:
        ds.GetRasterBand(1).SetNoDataValue(float(da.rio.nodata))
    out_ds = gdal.DEMProcessing(
        "",
        ds,
        subcommand,
        format="MEM",
        zFactor=z,
        scale=s,
        multiDirectional=True,
        computeEdges=True,
    )
    result = out_ds.ReadAsArray()
    out_ds = None
    ds = None
    out_da = xr.DataArray(result, dims=("y", "x"), coords={"y": da.y.values, "x": da.x.values})
    out_da.rio.write_crs(da.rio.crs, inplace=True)
    out_da.rio.write_transform(affine, inplace=True)
    return out_da
