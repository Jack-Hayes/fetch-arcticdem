"""
Wrapper script to fetch and coregister ArcticDEM DEMs.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.coregister import coregister_arcticdem_stack
from src.fetch import fetch_arcticdem_stack


def main():
    parser = argparse.ArgumentParser(
        description="Fetch and coregister ArcticDEM DEMs (wrapper script)"
    )

    aoi_group = parser.add_mutually_exclusive_group(required=True)
    aoi_group.add_argument("--aoi_file", type=str, help="Path to AOI file")
    aoi_group.add_argument(
        "--bounds",
        type=float,
        nargs=4,
        metavar=("MINX", "MINY", "MAXX", "MAXY"),
        help="WGS84 bounds",
    )

    parser.add_argument(
        "--fetch_output_dir", type=str, required=True, help="Output directory for fetched DEMs"
    )

    parser.add_argument(
        "--coreg_output_dir", type=str, required=True, help="Output directory for coregistered DEMs"
    )

    parser.add_argument("--stac_url", type=str, default="https://stac.pgc.umn.edu/api/v1/")

    parser.add_argument("--collection", type=str, default="arcticdem-strips-s2s041-2m")

    parser.add_argument("--date_range", type=str, nargs=2, metavar=("START_DATE", "END_DATE"))

    parser.add_argument("--min_valid_fraction", type=float, default=0.0)

    parser.add_argument("--intersection_threshold", type=float, default=0.8)

    parser.add_argument("--resolution", type=float, default=2)

    parser.add_argument("--out_crs", type=str, default=None)

    ref_group = parser.add_mutually_exclusive_group()
    ref_group.add_argument("--ref_index", type=int)
    ref_group.add_argument("--ref_date", type=str)

    parser.add_argument(
        "--coreg_steps",
        type=str,
        nargs="+",
        choices=["VerticalShift", "ICP", "NuthKaab", "AffineCoreg", "DhMinimize"],
        default=None,
    )

    parser.add_argument("--slope_min", type=float, default=2)

    parser.add_argument("--slope_max", type=float, default=30)

    parser.add_argument("--generate_slope_files", action="store_true")

    args = parser.parse_args()

    print("=" * 60)
    print("STEP 1: FETCHING ARCTICDEM DEMS")
    print("=" * 60)

    fetch_arcticdem_stack(
        aoi_file=args.aoi_file,
        bounds=args.bounds,
        output_dir=args.fetch_output_dir,
        stac_url=args.stac_url,
        collection=args.collection,
        date_range=args.date_range,
        min_valid_fraction=args.min_valid_fraction,
        intersection_threshold=args.intersection_threshold,
        resolution=args.resolution,
        generate_hillshade=False,
        out_crs=args.out_crs,
    )

    print("\n" + "=" * 60)
    print("STEP 2: COREGISTERING DEMS")
    print("=" * 60)

    coregister_arcticdem_stack(
        input_dir=args.fetch_output_dir,
        output_dir=args.coreg_output_dir,
        ref_index=args.ref_index,
        ref_date=args.ref_date,
        coreg_steps=args.coreg_steps,
        slope_min=args.slope_min,
        slope_max=args.slope_max,
        resolution=args.resolution,
        generate_hillshade=True,
        generate_slope_files=args.generate_slope_files,
    )

    print("\n" + "=" * 60)
    print("FETCH AND COREGISTRATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
