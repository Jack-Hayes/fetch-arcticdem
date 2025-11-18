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


def _calculate_nmad(dem1_values, dem2_values):
    """
    Calculate NMAD between two DEM arrays.

    Parameters
    ----------
    dem1_values : np.ndarray
        First DEM array
    dem2_values : np.ndarray
        Second DEM array

    Returns
    -------
    float
        NMAD value, or np.nan if insufficient valid pixels
    """
    dh = dem2_values - dem1_values
    valid_mask = ~np.isnan(dh)

    if np.count_nonzero(valid_mask) < 100:
        return np.nan

    return xdem.spatialstats.nmad(dh[valid_mask])


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
        "Deramp": xdem.coreg.Deramp,
    }

    if not steps:
        raise ValueError("At least one coregistration step must be provided")

    pipeline = step_map[steps[0]]()
    for step_name in steps[1:]:
        pipeline = pipeline + step_map[step_name]()

    return pipeline


def select_reference_dem(dem_stack, ref_index=None, ref_date=None, subsample=10):
    """
    Select reference DEM from stack.

    If no ref_index or ref_date provided, selects DEM with lowest median NMAD
    to all other DEMs in the stack.

    Parameters
    ----------
    dem_stack : xr.DataArray
        Stack of DEMs with time dimension
    ref_index : int, optional
        Index of reference DEM
    ref_date : str, optional
        Date of reference DEM (YYYY-MM-DD)
    subsample : int
        Subsampling factor for NMAD calculation (default: 10)

    Returns
    -------
    tuple
        (reference_index, reference_DataArray)
    """
    if ref_index is not None:
        return ref_index, dem_stack.isel(time=ref_index)

    if ref_date is not None:
        ref_timestamp = pd.Timestamp(ref_date)
        time_diffs = np.abs(
            [(pd.Timestamp(t.values) - ref_timestamp).total_seconds() for t in dem_stack.time]
        )
        ref_index = int(np.argmin(time_diffs))
        return ref_index, dem_stack.isel(time=ref_index)

    # Automatic selection based on lowest median NMAD
    n_dems = len(dem_stack.time)

    if n_dems == 1:
        return 0, dem_stack.isel(time=0)

    # Subsample for efficiency
    dem_arrays = dem_stack.values[:, ::subsample, ::subsample]

    # Calculate all pairwise differences at once
    # Shape: (n_dems, n_dems, height, width)
    dh_matrix = dem_arrays[:, np.newaxis, :, :] - dem_arrays[np.newaxis, :, :, :]

    # Calculate NMAD for each pair efficiently
    median_nmads = np.full(n_dems, np.inf)

    for i in range(n_dems):
        # Get all differences for DEM i (excluding self-comparison)
        dh_i = np.concatenate([dh_matrix[i, :i], dh_matrix[i, i + 1 :]], axis=0)

        # Calculate valid pixels across all comparisons
        valid_mask = ~np.isnan(dh_i)

        # Calculate NMAD for each comparison
        nmads_i = []
        for j in range(n_dems - 1):
            dh_j = dh_i[j]
            valid_j = valid_mask[j]

            if np.count_nonzero(valid_j) >= 100:
                nmad_j = xdem.spatialstats.nmad(dh_j[valid_j])
                nmads_i.append(nmad_j)

        if nmads_i:
            median_nmads[i] = np.median(nmads_i)

    # Select DEM with lowest median NMAD
    ref_index = int(np.argmin(median_nmads))

    return ref_index, dem_stack.isel(time=ref_index)


def create_inlier_mask_from_dh(
    dem_stack, ref_idx, slope_min=2, slope_max=30, nmad_threshold=1.0, logger=None
):
    """
    Create inlier mask based on elevation differences between first and last DEM.
    Falls back to slope-based mask if insufficient valid pixels.

    Parameters
    ----------
    dem_stack : xr.DataArray
        Stack of DEMs with time dimension
    ref_idx : int
        Index of reference DEM
    slope_min : float
        Minimum slope threshold for fallback mask (degrees)
    slope_max : float
        Maximum slope threshold for fallback mask (degrees)
    nmad_threshold : float
        NMAD threshold multiplier for inlier selection
    logger : logging.Logger, optional
        Logger instance

    Returns
    -------
    np.ndarray
        Boolean inlier mask
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    ref_da = dem_stack.isel(time=ref_idx)
    ref_dem = xdem.DEM.from_xarray(ref_da)

    # Get first and last DEMs
    first_dem = dem_stack.isel(time=0)
    last_dem = dem_stack.isel(time=-1)

    # Calculate NMAD using helper function
    nmad = _calculate_nmad(first_dem.values, last_dem.values)

    # Check if we have valid NMAD calculation
    if np.isnan(nmad):
        logger.warning(
            "Insufficient valid pixels for NMAD-based inlier mask. "
            f"Falling back to slope-based mask ({slope_min}-{slope_max} degrees)."
        )
        slope_raster = xdem.terrain.slope(ref_dem)
        inlier_mask = (slope_raster.data.filled(np.nan) > slope_min) & (
            slope_raster.data.filled(np.nan) < slope_max
        )
        inlier_mask &= ~np.isnan(ref_dem.data.filled(np.nan))
        logger.info(f"Using slope-based inlier mask with {np.count_nonzero(inlier_mask)} pixels")
        return inlier_mask

    # Calculate elevation difference for mask creation
    dh = last_dem.values - first_dem.values
    valid_mask = ~np.isnan(dh)
    valid_count = np.count_nonzero(valid_mask)

    logger.info(f"Valid pixels for dh calculation: {valid_count}")
    logger.info(f"Elevation difference NMAD: {nmad:.3f} m")
    logger.info(f"Using NMAD threshold: {nmad_threshold} * NMAD = {nmad_threshold * nmad:.3f} m")

    # Create NMAD-based inlier mask
    inlier_mask = valid_mask & (np.abs(dh) < (nmad_threshold * nmad))
    inlier_pixels = np.count_nonzero(inlier_mask)

    logger.info(f"NMAD-based inlier mask contains {inlier_pixels} pixels")

    # Check if NMAD mask has sufficient pixels
    min_pixels_required = 100
    if inlier_pixels < min_pixels_required:
        logger.warning(
            f"NMAD-based mask has too few pixels ({inlier_pixels} < {min_pixels_required}). "
            f"Falling back to slope-based mask ({slope_min}-{slope_max} degrees)."
        )
        slope_raster = xdem.terrain.slope(ref_dem)
        inlier_mask = (slope_raster.data.filled(np.nan) > slope_min) & (
            slope_raster.data.filled(np.nan) < slope_max
        )
        inlier_mask &= ~np.isnan(ref_dem.data.filled(np.nan))
        logger.info(f"Using slope-based inlier mask with {np.count_nonzero(inlier_mask)} pixels")

    return inlier_mask


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
        Minimum slope threshold for fallback inlier mask (degrees)
    slope_max : float
        Maximum slope threshold for fallback inlier mask (degrees)
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

    # Create reference DEM object for coregistration
    ref_dem = xdem.DEM.from_xarray(ref_da)

    logger.info("Generating inlier mask from elevation differences")
    inlier_mask_array = create_inlier_mask_from_dh(
        dem_stack, ref_idx, slope_min, slope_max, nmad_threshold=1, logger=logger
    )

    logger.info(f"Using {np.count_nonzero(inlier_mask_array)} pixels as stable ground")

    if generate_slope_files:
        logger.info("Generating slope raster for reference DEM")
        slope_raster = xdem.terrain.slope(ref_dem)
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

            logger.info("Coregistration fit complete")

            # Extract transformation information
            coreg_file.write(f"DEM {i + 1}/{len(dem_stack.time)}: {datetime_str}\n")

            # Check if the pipeline is fully affine before trying to get a cumulative matrix
            is_affine_pipeline = all(step._is_affine for step in pipeline_i.pipeline)

            if is_affine_pipeline:
                # Get cumulative transformation from the final matrix
                final_matrix = pipeline_i.to_matrix()
                cumulative_shift_x = final_matrix[0, 3]
                cumulative_shift_y = final_matrix[1, 3]
                cumulative_shift_z = final_matrix[2, 3]

                # Write cumulative transformation (summary)
                coreg_file.write("  Cumulative Transformation:\n")
                coreg_file.write("    Shifts (m):\n")
                coreg_file.write(f"      Easting (X):  {cumulative_shift_x:>10.3f}\n")
                coreg_file.write(f"      Northing (Y): {cumulative_shift_y:>10.3f}\n")
                coreg_file.write(f"      Vertical (Z): {cumulative_shift_z:>10.3f}\n")
                coreg_file.write("    Transformation Matrix:\n")
                for row in final_matrix:
                    coreg_file.write(f"      {row}\n")

                # Extract rotations from final matrix
                try:
                    rotations = xdem.coreg.AffineCoreg.to_rotations(final_matrix)
                    coreg_file.write("    Rotations (degrees):\n")
                    coreg_file.write(f"      X-axis: {np.degrees(rotations[0]):>10.6f}\n")
                    coreg_file.write(f"      Y-axis: {np.degrees(rotations[1]):>10.6f}\n")
                    coreg_file.write(f"      Z-axis: {np.degrees(rotations[2]):>10.6f}\n")
                except Exception:
                    pass

                coreg_file.write("\n")

                # Log to console
                logger.info(
                    f"  Cumulative shifts - X: {cumulative_shift_x:.3f} m, "
                    f"Y: {cumulative_shift_y:.3f} m, Z: {cumulative_shift_z:.3f} m"
                )
            else:
                logger.info(
                    "  Pipeline contains non-affine steps. Cumulative matrix not calculated."
                )
                coreg_file.write(
                    "  Pipeline contains non-affine steps. Cumulative matrix not applicable.\n\n"
                )

            # Write detailed per-step transformations
            coreg_file.write("  Individual Step Transformations:\n")
            for step_idx, step in enumerate(pipeline_i.pipeline):
                step_name = type(step).__name__
                step_meta = step.meta

                coreg_file.write(f"\n    Step {step_idx + 1}: {step_name}\n")

                # Extract step-specific outputs
                if "outputs" in step_meta and "affine" in step_meta["outputs"]:
                    affine_outputs = step_meta["outputs"]["affine"]

                    # Write shifts if available
                    if (
                        "shift_x" in affine_outputs
                        or "shift_y" in affine_outputs
                        or "shift_z" in affine_outputs
                    ):
                        coreg_file.write("      Shifts (m):\n")
                        if "shift_x" in affine_outputs:
                            coreg_file.write(
                                f"        X: {float(affine_outputs['shift_x']):>10.6f}\n"
                            )
                        if "shift_y" in affine_outputs:
                            coreg_file.write(
                                f"        Y: {float(affine_outputs['shift_y']):>10.6f}\n"
                            )
                        if "shift_z" in affine_outputs:
                            coreg_file.write(
                                f"        Z: {float(affine_outputs['shift_z']):>10.6f}\n"
                            )

                    # Write matrix if available
                    if "matrix" in affine_outputs:
                        coreg_file.write("      Transformation Matrix:\n")
                        matrix = affine_outputs["matrix"]
                        for row in matrix:
                            coreg_file.write(f"        {row}\n")

                    # Write centroid if available (ICP)
                    if "centroid" in affine_outputs:
                        centroid = affine_outputs["centroid"]
                        coreg_file.write(
                            f"      Centroid: ({float(centroid[0]):.2f}, {float(centroid[1]):.2f}, {float(centroid[2]):.2f})\n"
                        )

                # Write iteration info if available
                if "outputs" in step_meta:
                    outputs = step_meta["outputs"]
                    if "random" in outputs and "subsample_final" in outputs["random"]:
                        subsample = outputs["random"]["subsample_final"]
                        coreg_file.write(f"      Subsample size: {subsample}\n")

            coreg_file.write("\n" + "=" * 80 + "\n\n")

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
