#!/usr/bin/env python3

import argparse
from spatialdata_io import xenium
import pandas
import anndata

def main():
    pandas.options.mode.string_storage = "python"
    anndata.settings.allow_write_nullable_strings = True
    parser = argparse.ArgumentParser(
        description="Convert 10x Xenium output directory to a SpatialData Zarr store."
    )

    parser.add_argument(
        "xenium_dir",
        help="Path to Xenium output directory"
    )

    parser.add_argument(
        "--output_zarr",
        default="./data.zarr/",
        help="Path to output Zarr file"
    )

    args = parser.parse_args()

    print(f"Reading Xenium directory: {args.xenium_dir}")
    sdata = xenium(args.xenium_dir)

    print(f"Writing SpatialData Zarr: {args.output_zarr}")
    sdata.write(args.output_zarr, overwrite=True)

    print("Done.")


if __name__ == "__main__":
    main()