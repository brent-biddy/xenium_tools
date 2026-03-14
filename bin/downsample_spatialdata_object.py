#!/usr/bin/env python3

import argparse
import numpy as np
from pathlib import Path

import spatialdata
from spatialdata import SpatialData


def get_point_coordinate_columns(points_df):
    candidates = [
        ("x_location", "y_location"),
        ("x", "y"),
        ("global_x", "global_y"),
    ]
    for x_col, y_col in candidates:
        if x_col in points_df.columns and y_col in points_df.columns:
            return {"x": x_col, "y": y_col}

    raise ValueError(
        "Could not determine point coordinate columns. "
        f"Available columns: {list(points_df.columns)}"
    )


def get_point_transformations(points):
    attrs = getattr(points, "attrs", {})
    transform = attrs.get("transform", {})
    if transform:
        return dict(transform)
    return None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Spatially-aware downsampling of a Xenium SpatialData Zarr store."
    )
    parser.add_argument("--zarr_file", required=True, help="Path to the input SpatialData Zarr store")
    parser.add_argument("--bin_size", type=float, default=50,
                        help="Side length of spatial bins in microns (default: 50)")
    parser.add_argument("--fraction", type=float, default=0.1,
                        help="Fraction of cells to sample per bin (default: 0.1)")
    return parser.parse_args()


def spatial_downsample(sdata, bin_size=50, fraction=0.1):

    table = sdata.tables["table"]
    coords = table.obsm["spatial"]

    x = coords[:, 0]
    y = coords[:, 1]

    bin_x = np.floor(x / bin_size).astype(int)
    bin_y = np.floor(y / bin_size).astype(int)
    bin_ids = [f"{bx}_{by}" for bx, by in zip(bin_x, bin_y)]

    cell_ids = np.array(table.obs_names)

    bins = {}
    for cell_id, bin_id in zip(cell_ids, bin_ids):
        bins.setdefault(bin_id, []).append(cell_id)

    rng = np.random.default_rng()
    sampled_cells = []
    for cells_in_bin in bins.values():
        n_target = max(1, int(np.floor(len(cells_in_bin) * fraction)))
        if len(cells_in_bin) <= n_target:
            sampled_cells.extend(cells_in_bin)
        else:
            sampled_cells.extend(rng.choice(cells_in_bin, n_target, replace=False).tolist())

    sampled_set = set(sampled_cells)

    print(f"Original cell count: {len(cell_ids)}")
    print(f"Downsampled cell count: {len(sampled_cells)}")
    print(f"Total bins: {len(bins)}")

    # Filter table
    filtered_table = table[table.obs_names.isin(sampled_set)].copy()

    # Filter shapes — preserve coordinate transformation attrs
    filtered_shapes = {}
    for shape_name, shapes in sdata.shapes.items():
        filtered = shapes.loc[shapes.index.isin(sampled_set)].copy()
        filtered.attrs = shapes.attrs
        filtered_shapes[shape_name] = filtered

    # Filter transcripts to those assigned to sampled cells
    filtered_points = {}
    for point_name, points in sdata.points.items():
        points_df = points.compute()
        print(f"Processing points element: {point_name}")
        print(f"Point columns: {list(points_df.columns)}")
        print(f"Point attrs keys: {list(getattr(points, 'attrs', {}).keys())}")
        if "cell_id" in points_df.columns:
            points_df = points_df[points_df["cell_id"].isin(sampled_set)]
        coordinates = get_point_coordinate_columns(points_df)
        transformations = get_point_transformations(points)
        parse_kwargs = {
            "coordinates": coordinates,
        }
        if transformations is not None:
            parse_kwargs["transformations"] = transformations
        filtered = spatialdata.models.PointsModel.parse(points_df, **parse_kwargs)
        filtered_points[point_name] = filtered

    downsampled_sdata = SpatialData(
        images=sdata.images,
        shapes=filtered_shapes,
        points=filtered_points,
        tables={"table": filtered_table},
    )

    return downsampled_sdata


def main():
    args = parse_args()

    print(f"Loading SpatialData from: {args.zarr_file}")
    sdata = spatialdata.read_zarr(args.zarr_file)

    downsampled_sdata = spatial_downsample(sdata, bin_size=args.bin_size, fraction=args.fraction)

    output_path = "data_downsampled.zarr"
    print(f"Writing downsampled SpatialData to: {output_path}")
    downsampled_sdata.write(output_path, overwrite=True)

    print("Done.")


if __name__ == "__main__":
    main()
