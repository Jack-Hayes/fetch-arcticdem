"""
CLI script to fetch ArcticDEM DEMs from STAC.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.fetch import fetch_arcticdem_stack


def main():
    parser = argparse.ArgumentParser(
        description="Fetch ArcticDEM DEMs from STAC and save as cloud-optimized GeoTIFFs"
    )

    aoi_group = parser.add_mutually_exclusive_group(required=True)
    aoi_group.add_argument(
        "--aoi_file", type=str, help="Path to AOI file (any format GeoPandas can read)"
    )
    aoi_group.add_argument(
        "--bounds",
        type=float,
        nargs=4,
        metavar=("MINX", "MINY", "MAXX", "MAXY"),
        help="WGS84 bounds: minx miny maxx maxy",
    )

    parser.add_argument(
        "--output_dir", type=str, required=True, help="Output directory for results"
    )

    parser.add_argument(
        "--stac_url",
        type=str,
        default="https://stac.pgc.umn.edu/api/v1/",
        help="STAC API URL (default: https://stac.pgc.umn.edu/api/v1/)",
    )

    parser.add_argument(
        "--collection",
        type=str,
        default="arcticdem-strips-s2s041-2m",
        help="STAC collection ID (default: arcticdem-strips-s2s041-2m)",
    )

    parser.add_argument(
        "--date_range",
        type=str,
        nargs=2,
        metavar=("START_DATE", "END_DATE"),
        help="Date range as YYYY-MM-DD YYYY-MM-DD (optional)",
    )

    parser.add_argument(
        "--min_valid_fraction",
        type=float,
        default=0.5,
        help="Minimum valid pixel fraction required (0.0-1.0, default: 0.0)",
    )

    parser.add_argument(
        "--intersection_threshold",
        type=float,
        default=0.8,
        help="Minimum intersection ratio with AOI (0.0-1.0, default: 0.8)",
    )

    parser.add_argument(
        "--resolution", type=float, default=2, help="Output resolution in meters (default: 2)"
    )

    parser.add_argument(
        "--generate_hillshade", action="store_true", help="Generate hillshade rasters"
    )

    parser.add_argument(
        "--out_crs",
        type=str,
        default=None,
        help="Output CRS as EPSG code (e.g., EPSG:32606). Default is native CRS.",
    )

    args = parser.parse_args()

    fetch_arcticdem_stack(
        aoi_file=args.aoi_file,
        bounds=args.bounds,
        output_dir=args.output_dir,
        stac_url=args.stac_url,
        collection=args.collection,
        date_range=args.date_range,
        min_valid_fraction=args.min_valid_fraction,
        intersection_threshold=args.intersection_threshold,
        resolution=args.resolution,
        generate_hillshade=args.generate_hillshade,
        out_crs=args.out_crs,
    )


if __name__ == "__main__":
    main()
