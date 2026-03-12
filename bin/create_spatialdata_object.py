#!/usr/bin/env python3

import argparse
import logging
import re
from pathlib import Path

import anndata
import ome_types
import pandas
import tifffile
from dask_image.imread import imread
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity
from spatialdata_io import xenium


def load_morphology_focus_v4(
    morphology_focus_dir: Path,
    image_models_kwargs: dict = {},
) -> Image2DModel:
    """Manually load v4 morphology focus images (ch0000_dapi.ome.tif naming scheme).

    This is a workaround for a bug in spatialdata-io where v4-style filenames are
    not detected when the experiment.xenium metadata reports a version < 2.0.0.
    """
    first_tiff_path = morphology_focus_dir / "ch0000_dapi.ome.tif"
    if not first_tiff_path.exists():
        raise FileNotFoundError(
            f"Expected v4 DAPI file not found at {first_tiff_path}. "
            "Check that your morphology_focus directory contains ch0000_dapi.ome.tif."
        )

    # Read channel names from OME-XML metadata embedded in the TIFF
    ome = ome_types.from_xml(tifffile.tiffcomment(first_tiff_path), validate=False)
    ome_channels = ome.images[0].pixels.channels
    channels = []
    for ome_ch in ome_channels:
        if ome_ch.name is None:
            raise ValueError(f"Found a channel without a name in {first_tiff_path}")
        match = re.match(r"Channel:(\d+)", ome_ch.id)
        if match is None:
            raise ValueError(
                f"Expected OME channel ID of the form 'Channel:<index>', found: {ome_ch.id}"
            )
        channels.append((int(match.group(1)), ome_ch.name))
    channel_names = [name for _, name in sorted(channels)]
    logging.info(f"Detected morphology focus channels: {channel_names}")

    # Suppress the tifffile warning about multi-file pyramids
    class IgnoreMultiFilePyramid(logging.Filter):
        def filter(self, record: logging.LogRecord) -> bool:
            return "OME series cannot read multi-file pyramids" not in record.getMessage()

    logger = tifffile.logger()
    logger.addFilter(IgnoreMultiFilePyramid())
    image = imread(first_tiff_path)
    logger.removeFilter(IgnoreMultiFilePyramid())

    return Image2DModel.parse(
        image,
        dims=("c", "y", "x"),
        transformations={"global": Identity()},
        c_coords=channel_names,
        rgb=None,
        **image_models_kwargs,
    )


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

    xenium_dir = Path(args.xenium_dir)

    # Use multiscale pyramids and explicit chunk sizes to avoid large-chunk warnings
    # and potential compression errors when writing.
    image_models_kwargs = {
        "scale_factors": [2, 2, 2, 2],
        "chunks": (1, 4096, 4096),
    }

    # Detect whether this dataset uses v4-style morphology focus filenames
    # (ch0000_dapi.ome.tif). These are not detected by spatialdata-io when the
    # experiment.xenium metadata reports a version < 2.0.0, requiring a manual load.
    morphology_focus_dir = xenium_dir / "morphology_focus"
    use_workaround = (morphology_focus_dir / "ch0000_dapi.ome.tif").exists()

    print(f"Reading Xenium directory: {xenium_dir}")
    if use_workaround:
        print("Detected v4-style morphology focus files — using manual load workaround.")
        sdata = xenium(xenium_dir, morphology_focus=False, image_models_kwargs=image_models_kwargs)
        sdata.images["morphology_focus"] = load_morphology_focus_v4(
            morphology_focus_dir, image_models_kwargs=image_models_kwargs
        )
    else:
        sdata = xenium(xenium_dir, morphology_focus=True, image_models_kwargs=image_models_kwargs)

    print(f"Writing SpatialData Zarr: {args.output_zarr}")
    sdata.write(args.output_zarr, overwrite=True)

    print("Done.")


if __name__ == "__main__":
    main()
