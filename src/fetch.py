"""
Fetch ArcticDEM DEMs from STAC and export as cloud-optimized GeoTIFFs.
"""

from pathlib import Path

import geopandas as gpd
import odc.stac
import pandas as pd
import pystac_client
import rasterio
import rioxarray as rxr  # noqa: F401
import stac_geoparquet.arrow

from .utils import filter_assets, gdaldem_hillshade_in_memory, intersection_ratio, setup_logging


def to_geopandas(collection):
    """
    Convert STAC ItemCollection to GeoDataFrame.

    Parameters
    ----------
    collection : pystac.ItemCollection
        STAC items to convert

    Returns
    -------
    gpd.GeoDataFrame
        Converted items
    """
    if len(collection) == 0:
        raise ValueError("ItemCollection is empty, cannot convert to GeoDataFrame")

    record_batch_reader = stac_geoparquet.arrow.parse_stac_items_to_arrow(collection)
    gf = gpd.GeoDataFrame.from_arrow(record_batch_reader)
    gf["assets"] = gf["assets"].apply(filter_assets)
    gf["dayofyear"] = gf["datetime"].dt.dayofyear

    return gf


def load_aoi_geometry(aoi_file=None, bounds=None):
    """
    Load AOI geometry from file or bounds.

    Parameters
    ----------
    aoi_file : str or Path, optional
        Path to AOI file (any format GeoPandas can read)
    bounds : list of float, optional
        WGS84 bounds [minx, miny, maxx, maxy]

    Returns
    -------
    tuple
        (GeoDataFrame, geometry, area)
    """
    if aoi_file is not None:
        gf = gpd.read_file(aoi_file)
        gf["geometry"] = gf["geometry"].envelope
        geom = gf.geometry.iloc[0]
        area = geom.area
        return gf, geom, area
    elif bounds is not None:
        from shapely.geometry import box

        gf = gpd.GeoDataFrame({"geometry": [box(*bounds)]}, crs="EPSG:4326")
        geom = gf.geometry.iloc[0]
        area = geom.area
        return gf, geom, area
    else:
        raise ValueError("Either aoi_file or bounds must be provided")


def fetch_arcticdem_stack(
    aoi_file=None,
    bounds=None,
    output_dir=".",
    stac_url="https://stac.pgc.umn.edu/api/v1/",
    collection="arcticdem-strips-s2s041-2m",
    date_range=None,
    min_valid_fraction=0.0,
    intersection_threshold=0.8,
    resolution=2,
    generate_hillshade=False,
    out_crs=None,
):
    """
    Fetch ArcticDEM DEMs from STAC and save as COGs.

    Parameters
    ----------
    aoi_file : str or Path, optional
        Path to AOI file
    bounds : list of float, optional
        WGS84 bounds [minx, miny, maxx, maxy]
    output_dir : str or Path
        Output directory
    stac_url : str
        STAC API URL
    collection : str
        STAC collection ID
    date_range : list of str, optional
        [min_date, max_date] in YYYY-MM-DD format
    min_valid_fraction : float
        Minimum valid pixel fraction (0.0-1.0)
    intersection_threshold : float
        Minimum intersection ratio with AOI (0.0-1.0)
    resolution : float
        Output resolution in meters
    generate_hillshade : bool
        Whether to generate hillshade
    out_crs : str, optional
        Output CRS as EPSG code (e.g., 'EPSG:32606')

    Returns
    -------
    None
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(output_dir, "fetch_arcticdem.log")
    logger.info("Starting ArcticDEM fetch process")
    logger.info(f"Output directory: {output_dir}")

    gf_aoi, aoi_geom, aoi_area = load_aoi_geometry(aoi_file, bounds)
    logger.info(f"Loaded AOI with area: {aoi_area:.2f} square units")

    catalog = pystac_client.Client.open(stac_url)
    search = catalog.search(collections=[collection], intersects=gf_aoi.geometry.envelope.item())
    items = search.item_collection()
    logger.info(f"Found {len(items)} items from STAC search")

    gf_stac = to_geopandas(items)
    gf_stac = gf_stac.drop("proj:geometry", axis=1, errors="ignore")

    gf_stac_filtered = gf_stac[
        gf_stac.geometry.apply(
            lambda item_geom: intersection_ratio(item_geom, aoi_geom, aoi_area)
            >= intersection_threshold
        )
    ]
    logger.info(
        f"Filtered by intersection ratio >= {intersection_threshold}: {len(gf_stac_filtered)} items"
    )

    if date_range is not None:
        min_date = pd.Timestamp(date_range[0])
        max_date = pd.Timestamp(date_range[1])
        if gf_stac_filtered["datetime"].dt.tz is not None:
            min_date = min_date.tz_localize("UTC")
            max_date = max_date.tz_localize("UTC")
        gf_stac_filtered = gf_stac_filtered[
            (gf_stac_filtered["datetime"] >= min_date) & (gf_stac_filtered["datetime"] <= max_date)
        ]
        logger.info(
            f"Filtered by date range {date_range[0]} to {date_range[1]}: {len(gf_stac_filtered)} items"
        )

    filtered_ids = set(gf_stac_filtered["id"])
    selected_items = [item for item in items if item.id in filtered_ids]

    if out_crs is None:
        target_crs = gf_stac_filtered["proj:code"].iloc[0]
    else:
        target_crs = out_crs

    logger.info(f"Target CRS: {target_crs}")

    gf_aoi_proj = gf_aoi.to_crs(target_crs)
    load_bbox = gf_aoi_proj.total_bounds

    bands_to_load = ["dem", "mask"]
    dask_chunks = {"x": 2048, "y": 2048}
    stac_cfg = {
        collection: {
            "assets": {
                "dem": {"nodata": -9999},
                "mask": {"nodata": 0},
            },
            "aliases": {
                "dem": "dem",
                "mask": "mask",
            },
        },
        "*": {"warnings": "ignore"},
    }

    processed_count = 0

    for idx, item in enumerate(selected_items):
        logger.info(f"Processing scene {idx + 1}/{len(selected_items)}: {item.id}")
        logger.info(f"Date: {item.datetime}")

        ds = odc.stac.load(
            [item],
            bands=bands_to_load,
            crs=target_crs,
            x=(load_bbox[0], load_bbox[2]),
            y=(load_bbox[1], load_bbox[3]),
            resolution=resolution,
            chunks=dask_chunks,
            groupby="time",
            stac_cfg=stac_cfg,
        )

        da_dem = ds["dem"].isel(time=0).squeeze()
        da_mask = ds["mask"].isel(time=0).squeeze()
        valid_mask = da_mask == 0
        nodata_val = (
            ds.dem.rio.nodata if hasattr(ds.dem, "rio") and ds.dem.rio.nodata is not None else -9999
        )
        da_dem_masked = da_dem.where(valid_mask, nodata_val)
        da_dem_masked.rio.write_nodata(nodata_val, inplace=True)

        valid_pixel_count = int(valid_mask.sum().compute())
        total_pixel_count = da_mask.size
        valid_fraction = float(valid_pixel_count) / total_pixel_count
        pixel_area_m2 = abs(da_dem.rio.resolution()[0] * da_dem.rio.resolution()[1])
        valid_area_km2 = valid_pixel_count * pixel_area_m2 / 1e6

        logger.info(f"Valid pixels: {valid_pixel_count}/{total_pixel_count} ({valid_fraction:.3f})")
        logger.info(f"Valid area: {valid_area_km2:.2f} km2")

        if valid_fraction < min_valid_fraction:
            logger.info(
                f"Skipping scene due to valid fraction {valid_fraction:.3f} < {min_valid_fraction}"
            )
            continue

        meta_data = {
            "datetime": str(item.datetime),
            "gsd": item.properties.get("gsd"),
            "instruments": item.properties.get("instruments"),
            "avg_convergence_angle": item.properties.get("pgc:avg_convergence_angle"),
            "avg_expected_height_accuracy_m": item.properties.get(
                "pgc:avg_expected_height_accuracy"
            ),
            "avg_sun_elevs": item.properties.get("pgc:avg_sun_elevs"),
            "image_ids": item.properties.get("pgc:image_ids"),
            "pairname": item.properties.get("pgc:pairname"),
            "coreg_rmse_m": item.properties.get("pgc:rmse"),
            "valid_pixel_count": valid_pixel_count,
            "total_pixel_count": total_pixel_count,
            "valid_fraction": valid_fraction,
            "valid_area_km2": valid_area_km2,
        }

        dem_filename = f"{item.id}_dem.tif"
        dem_out_path = output_dir / dem_filename

        with rasterio.open(
            dem_out_path,
            "w",
            driver="COG",
            height=da_dem_masked.shape[0],
            width=da_dem_masked.shape[1],
            count=1,
            dtype=da_dem_masked.dtype,
            crs=da_dem_masked.rio.crs,
            transform=da_dem_masked.rio.transform(),
            nodata=nodata_val,
            compress="DEFLATE",
        ) as dst:
            dst.write(da_dem_masked.values, 1)
            dst.update_tags(**{k: str(v) for k, v in meta_data.items()})
            dst.set_band_description(1, "elevation")

        logger.info(f"Wrote DEM COG: {dem_out_path}")

        if generate_hillshade:
            hill_da = gdaldem_hillshade_in_memory(da_dem_masked, z=1.0, s=1.0)
            hill_filename = f"{item.id}_hillshade.tif"
            hill_out_path = output_dir / hill_filename

            with rasterio.open(
                hill_out_path,
                "w",
                driver="COG",
                height=hill_da.shape[0],
                width=hill_da.shape[1],
                count=1,
                dtype="uint8",
                crs=hill_da.rio.crs,
                transform=hill_da.rio.transform(),
                nodata=0,
                compress="DEFLATE",
            ) as dst:
                dst.write(hill_da.values.astype("uint8"), 1)
                dst.update_tags(**{k: str(v) for k, v in meta_data.items()})
                dst.set_band_description(1, "multidirectional_hillshade")

            logger.info(f"Wrote hillshade COG: {hill_out_path}")

        processed_count += 1

    logger.info(f"Fetch complete. Processed {processed_count} scenes.")
