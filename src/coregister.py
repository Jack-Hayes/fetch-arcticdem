"""
Coregister ArcticDEM DEMs using xdem.
"""

import logging
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import rioxarray as rxr
import xarray as xr
import xdem
from osgeo import gdal

from .utils import gdaldem_hillshade_in_memory, setup_logging

# Enable GDAL exceptions to avoid FutureWarning
gdal.UseExceptions()


def build_coreg_pipeline(steps):
    """
    Build xdem coregistration pipeline from list of step names.

    Parameters
    ----------
    steps : list of str
        List of coregistration steps: 'VerticalShift', 'ICP', 'NuthKaab',
        'AffineCoreg', 'DhMinimize'

    Returns
    -------
    xdem.coreg pipeline
        Combined coregistration pipeline
    """
    step_map = {
        "VerticalShift": xdem.coreg.VerticalShift,
        "ICP": xdem.coreg.ICP,
        "NuthKaab": xdem.coreg.NuthKaab,
        "AffineCoreg": xdem.coreg.AffineCoreg,
        "DhMinimize": xdem.coreg.DhMinimize,
    }

    if not steps:
        raise ValueError("At least one coregistration step must be provided")

    pipeline = step_map[steps[0]]()
    for step_name in steps[1:]:
        pipeline = pipeline + step_map[step_name]()

    return pipeline


def select_reference_dem(dem_stack, ref_index=None, ref_date=None):
    """
    Select reference DEM from stack.

    Parameters
    ----------
    dem_stack : xr.DataArray
        Stack of DEMs with time dimension
    ref_index : int, optional
        Index of reference DEM
    ref_date : str, optional
        Date of reference DEM (YYYY-MM-DD)

    Returns
    -------
    tuple
        (reference_index, reference_DataArray)
    """
    if ref_index is not None:
        return ref_index, dem_stack.isel(time=ref_index)

    if ref_date is not None:
        ref_timestamp = pd.Timestamp(ref_date)
        time_diffs = [
            abs((pd.Timestamp(t.values) - ref_timestamp).total_seconds()) for t in dem_stack.time
        ]
        ref_index = int(np.argmin(time_diffs))
        return ref_index, dem_stack.isel(time=ref_index)

    valid_counts = []
    for i in range(len(dem_stack.time)):
        da = dem_stack.isel(time=i)
        valid_count = (~np.isnan(da.values)).sum()
        valid_counts.append(valid_count)

    max_valid = max(valid_counts)
    max_indices = [i for i, v in enumerate(valid_counts) if v == max_valid]

    if len(max_indices) == 1:
        ref_index = max_indices[0]
    else:
        times = [pd.Timestamp(dem_stack.time.values[i]) for i in max_indices]
        mean_time = pd.Timestamp(np.mean([t.value for t in times]))
        time_diffs = [abs((t - mean_time).total_seconds()) for t in times]
        ref_index = max_indices[np.argmin(time_diffs)]

    return ref_index, dem_stack.isel(time=ref_index)


def coregister_arcticdem_stack(
    input_dir,
    output_dir,
    ref_index=None,
    ref_date=None,
    coreg_steps=None,
    slope_min=2,
    slope_max=30,
    resolution=2,
    generate_hillshade=True,
    generate_slope_files=False,
):
    """
    Coregister ArcticDEM stack using xdem.

    Parameters
    ----------
    input_dir : str or Path
        Directory containing input DEM COGs
    output_dir : str or Path
        Output directory for coregistered DEMs
    ref_index : int, optional
        Index of reference DEM
    ref_date : str, optional
        Date of reference DEM (YYYY-MM-DD)
    coreg_steps : list of str, optional
        Coregistration steps, defaults to ['VerticalShift', 'ICP', 'NuthKaab']
    slope_min : float
        Minimum slope threshold for inlier mask (degrees)
    slope_max : float
        Maximum slope threshold for inlier mask (degrees)
    resolution : float
        Output resolution in meters
    generate_hillshade : bool
        Whether to generate hillshade for aligned DEMs
    generate_slope_files : bool
        Whether to output slope rasters

    Returns
    -------
    None
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(output_dir, "coregister_arcticdem.log")
    logger.info("Starting ArcticDEM coregistration process")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")

    # Configure xdem logging
    xdem_logger = logging.getLogger("xdem")
    xdem_logger.setLevel(logging.INFO)
    for handler in logger.handlers:
        xdem_logger.addHandler(handler)

    dem_files = sorted(glob(str(input_dir / "*_dem.tif")))

    if not dem_files:
        raise ValueError(f"No DEM files found in {input_dir}")

    logger.info(f"Found {len(dem_files)} DEM files")

    dems = []
    times = []

    for f in dem_files:
        da = rxr.open_rasterio(f, masked=True).squeeze("band", drop=True)
        da_utm = da.rio.reproject(da.rio.crs, resolution=resolution)

        with rasterio.open(f) as src:
            tags = src.tags()
            datetime_str = tags.get("datetime", None)

        if datetime_str:
            t = pd.to_datetime(datetime_str)
        else:
            t = pd.Timestamp.now()

        dems.append(da_utm)
        times.append(t)

    dem_stack = xr.concat(dems, dim=pd.Index(times, name="time"))
    dem_stack = dem_stack.sortby("time")

    logger.info(f"Loaded DEM stack with {len(dem_stack.time)} timestamps")

    ref_idx, ref_da = select_reference_dem(dem_stack, ref_index, ref_date)
    logger.info(
        f"Selected reference DEM at index {ref_idx}, date {pd.to_datetime(ref_da.time.values)}"
    )

    ref_dem = xdem.DEM.from_xarray(ref_da)

    logger.info("Generating inlier mask from reference DEM")
    slope_raster = xdem.terrain.slope(ref_dem)

    inlier_mask_array = (slope_raster.data.filled(np.nan) > slope_min) & (
        slope_raster.data.filled(np.nan) < slope_max
    )
    inlier_mask_array &= ~np.isnan(ref_dem.data.filled(np.nan))

    logger.info(f"Using {np.count_nonzero(inlier_mask_array)} pixels as stable ground")

    if generate_slope_files:
        slope_output = output_dir / "reference_slope.tif"
        slope_da = xr.DataArray(
            slope_raster.data,
            dims=("y", "x"),
            coords={"y": ref_da.y.values, "x": ref_da.x.values},
        )
        slope_da.rio.write_crs(ref_da.rio.crs, inplace=True)
        slope_da.rio.write_transform(ref_da.rio.transform(), inplace=True)
        slope_da.rio.to_raster(slope_output, driver="COG", compress="DEFLATE")
        logger.info(f"Wrote slope raster: {slope_output}")

    if coreg_steps is None:
        coreg_steps = ["VerticalShift", "ICP", "NuthKaab"]

    logger.info(f"Coregistration pipeline: {' + '.join(coreg_steps)}")

    # Open coregistration results file
    coreg_results_path = output_dir / "coregistration_results.txt"
    with open(coreg_results_path, "w") as coreg_file:
        coreg_file.write("ArcticDEM Coregistration Results\n")
        coreg_file.write("=" * 80 + "\n\n")
        coreg_file.write(
            f"Reference DEM: Index {ref_idx}, Date {pd.to_datetime(ref_da.time.values)}\n"
        )
        coreg_file.write(f"Coregistration pipeline: {' + '.join(coreg_steps)}\n")
        coreg_file.write(f"Slope thresholds: {slope_min}-{slope_max} degrees\n")
        coreg_file.write(f"Stable ground pixels: {np.count_nonzero(inlier_mask_array)}\n")
        coreg_file.write("\n" + "=" * 80 + "\n\n")

        for i in range(len(dem_stack.time)):
            tba_da = dem_stack.isel(time=i)
            original_attrs = tba_da.attrs.copy()
            datetime_str = str(pd.to_datetime(tba_da.time.values))

            if i == ref_idx:
                logger.info(
                    f"Processing DEM {i + 1}/{len(dem_stack.time)} (reference, no coregistration)"
                )

                tba_da.attrs["coreg_ref_time"] = datetime_str
                tba_da.attrs["datetime"] = datetime_str

                out_name = f"aligned_dem_{i:02d}_{pd.to_datetime(tba_da.time.values).strftime('%Y%m%dT%H%M%S')}.tif"
                out_path = output_dir / out_name
                tba_da.rio.to_raster(out_path, driver="COG", compress="DEFLATE")
                logger.info(f"Wrote aligned DEM (reference): {out_path}")

                # Write to results file
                coreg_file.write(f"DEM {i + 1}/{len(dem_stack.time)}: {datetime_str} (REFERENCE)\n")
                coreg_file.write("  No transformation applied (reference DEM)\n\n")

                if generate_hillshade:
                    hill_da = gdaldem_hillshade_in_memory(tba_da, z=1.0, s=1.0)
                    hill_name = f"aligned_hillshade_{i:02d}_{pd.to_datetime(tba_da.time.values).strftime('%Y%m%dT%H%M%S')}.tif"
                    hill_path = output_dir / hill_name

                    with rasterio.open(
                        hill_path,
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
                        dst.set_band_description(1, "multidirectional_hillshade")

                    logger.info(f"Wrote hillshade: {hill_path}")

                continue

            tba_dem = xdem.DEM.from_xarray(tba_da)

            logger.info(f"Processing DEM {i + 1}/{len(dem_stack.time)} ({datetime_str})")

            pipeline_i = build_coreg_pipeline(coreg_steps)

            logger.info("Fitting coregistration pipeline...")
            pipeline_i.fit(
                reference_elev=ref_dem, to_be_aligned_elev=tba_dem, inlier_mask=inlier_mask_array
            )
            logger.info("Coregistration fit complete")

            # Extract transformation information
            coreg_file.write(f"DEM {i + 1}/{len(dem_stack.time)}: {datetime_str}\n")

            # Get metadata from pipeline
            meta = pipeline_i.meta

            # Log affine transformation parameters
            if "outputs" in meta:
                outputs = meta["outputs"]

                # Extract shifts
                shift_x = outputs.get("shift_x", 0.0)
                shift_y = outputs.get("shift_y", 0.0)
                shift_z = outputs.get("shift_z", 0.0)

                coreg_file.write("  Shifts (m):\n")
                coreg_file.write(f"    Easting (X):  {shift_x:>10.3f}\n")
                coreg_file.write(f"    Northing (Y): {shift_y:>10.3f}\n")
                coreg_file.write(f"    Vertical (Z): {shift_z:>10.3f}\n")

                logger.info(
                    f"  Shifts - X: {shift_x:.3f} m, Y: {shift_y:.3f} m, Z: {shift_z:.3f} m"
                )

                # If affine matrix is available, extract rotations
                if "matrix" in outputs:
                    matrix = outputs["matrix"]
                    coreg_file.write("  Affine transformation matrix:\n")
                    for row in matrix:
                        coreg_file.write(f"    {row}\n")

                    # Extract rotations using xdem utility
                    try:
                        rotations = xdem.coreg.AffineCoreg.to_rotations(matrix)
                        coreg_file.write("  Rotations (degrees):\n")
                        coreg_file.write(f"    X-axis: {np.degrees(rotations[0]):>10.6f}\n")
                        coreg_file.write(f"    Y-axis: {np.degrees(rotations[1]):>10.6f}\n")
                        coreg_file.write(f"    Z-axis: {np.degrees(rotations[2]):>10.6f}\n")
                        logger.info(
                            f"  Rotations - X: {np.degrees(rotations[0]):.6f}°, "
                            f"Y: {np.degrees(rotations[1]):.6f}°, "
                            f"Z: {np.degrees(rotations[2]):.6f}°"
                        )
                    except Exception:
                        pass

                # Log iteration info if available
                if "last_iteration" in outputs:
                    last_iter = outputs["last_iteration"]
                    coreg_file.write(f"  Iterations completed: {last_iter}\n")
                    logger.info(f"  Iterations: {last_iter}")

                if "all_tolerances" in outputs:
                    tolerances = outputs["all_tolerances"]
                    coreg_file.write(f"  Final tolerance: {tolerances[-1]:.6f}\n")

            coreg_file.write("\n")

            aligned_dem_obj = pipeline_i.apply(elev=tba_dem)

            aligned_da = aligned_dem_obj.to_xarray()
            if "band" in aligned_da.dims:
                aligned_da = aligned_da.squeeze("band", drop=True)
            aligned_da = aligned_da.assign_coords(time=tba_da.time)

            aligned_da.attrs.update(original_attrs)
            aligned_da.attrs["datetime"] = datetime_str
            aligned_da.attrs["coreg_ref_time"] = str(pd.to_datetime(ref_da.time.values))

            out_name = f"aligned_dem_{i:02d}_{pd.to_datetime(tba_da.time.values).strftime('%Y%m%dT%H%M%S')}.tif"
            out_path = output_dir / out_name
            aligned_da.rio.to_raster(out_path, driver="COG", compress="DEFLATE")
            logger.info(f"Wrote aligned DEM: {out_path}")

            if generate_hillshade:
                hill_da = gdaldem_hillshade_in_memory(aligned_da, z=1.0, s=1.0)
                hill_name = f"aligned_hillshade_{i:02d}_{pd.to_datetime(tba_da.time.values).strftime('%Y%m%dT%H%M%S')}.tif"
                hill_path = output_dir / hill_name

                with rasterio.open(
                    hill_path,
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
                    dst.set_band_description(1, "multidirectional_hillshade")

                logger.info(f"Wrote hillshade: {hill_path}")

    logger.info(f"Coregistration results written to: {coreg_results_path}")
    logger.info("Coregistration complete")
