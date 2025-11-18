"""
CLI script to coregister ArcticDEM DEMs using xdem.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.coregister import coregister_arcticdem_stack


def main():
    parser = argparse.ArgumentParser(description="Coregister ArcticDEM DEMs using xdem")

    parser.add_argument(
        "--input_dir", type=str, required=True, help="Directory containing input DEM COGs"
    )

    parser.add_argument(
        "--output_dir", type=str, required=True, help="Output directory for coregistered DEMs"
    )

    ref_group = parser.add_mutually_exclusive_group()
    ref_group.add_argument("--ref_index", type=int, help="Index of reference DEM in stack")
    ref_group.add_argument("--ref_date", type=str, help="Date of reference DEM (YYYY-MM-DD)")

    parser.add_argument(
        "--coreg_steps",
        type=str,
        nargs="+",
        choices=["VerticalShift", "ICP", "NuthKaab", "AffineCoreg", "DhMinimize", "Deramp"],
        default=None,
        help="Coregistration steps (default: VerticalShift ICP NuthKaab)",
    )

    parser.add_argument(
        "--slope_min",
        type=float,
        default=2,
        help="Minimum slope threshold for inlier mask in degrees (default: 2)",
    )

    parser.add_argument(
        "--slope_max",
        type=float,
        default=30,
        help="Maximum slope threshold for inlier mask in degrees (default: 30)",
    )

    parser.add_argument(
        "--resolution", type=float, default=2, help="Output resolution in meters (default: 2)"
    )

    parser.add_argument(
        "--generate_hillshade",
        action="store_true",
        default=True,
        help="Generate hillshade for aligned DEMs (default: True)",
    )

    parser.add_argument(
        "--no_hillshade",
        dest="generate_hillshade",
        action="store_false",
        help="Do not generate hillshade",
    )

    parser.add_argument(
        "--generate_slope_files",
        action="store_true",
        help="Output slope rasters used for coregistration",
    )

    args = parser.parse_args()

    coregister_arcticdem_stack(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        ref_index=args.ref_index,
        ref_date=args.ref_date,
        coreg_steps=args.coreg_steps,
        slope_min=args.slope_min,
        slope_max=args.slope_max,
        resolution=args.resolution,
        generate_hillshade=args.generate_hillshade,
        generate_slope_files=args.generate_slope_files,
    )


if __name__ == "__main__":
    main()
