#!/usr/bin/env python3
"""Per-nucleus measurements from a label volume and multi-channel image.

Computes volume, centroid, and per-channel mean/integrated intensity for
every labelled nucleus. Writes measurements.tsv.

Usage:
    python src/02_measure.py --image data/synthetic_blastocyst.ome.tif --labels outputs/labels.tif
"""
import argparse
import numpy as np
import pandas as pd
import tifffile
from pathlib import Path
from skimage.measure import regionprops_table


def load_multichannel(path):
    """Load a multi-channel 3D image as (C, Z, Y, X)."""
    img = tifffile.imread(str(path))
    if img.ndim == 4:
        if img.shape[0] <= 10:
            return img  # CZYX
        # ZCYX -> CZYX
        return np.moveaxis(img, 1, 0)
    elif img.ndim == 3:
        return img[np.newaxis]  # ZYX -> 1ZYX
    else:
        raise ValueError(f"unexpected image shape {img.shape}")


def main():
    parser = argparse.ArgumentParser(description="per-nucleus measurements")
    parser.add_argument("--image", required=True, help="multi-channel image")
    parser.add_argument("--labels", required=True, help="label volume")
    parser.add_argument("--channel-names", nargs="+",
                        default=["DAPI", "marker_A", "marker_B", "marker_C"],
                        help="channel names")
    parser.add_argument("--output-dir", default="outputs")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)

    img = load_multichannel(args.image)
    labels = tifffile.imread(str(args.labels))
    n_channels = img.shape[0]
    ch_names = args.channel_names[:n_channels]

    # basic morphology from the label volume
    props = regionprops_table(labels, properties=("label", "area", "centroid"))
    df = pd.DataFrame(props)
    df.rename(columns={
        "area": "volume_voxels",
        "centroid-0": "centroid_z",
        "centroid-1": "centroid_y",
        "centroid-2": "centroid_x",
    }, inplace=True)

    # per-channel intensity
    for c, name in enumerate(ch_names):
        ch_props = regionprops_table(labels, intensity_image=img[c],
                                     properties=("label", "intensity_mean"))
        ch_df = pd.DataFrame(ch_props)
        df[f"{name}_mean"] = ch_df["intensity_mean"]
        df[f"{name}_integrated"] = df[f"{name}_mean"] * df["volume_voxels"]

    out_path = out_dir / "measurements.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path}  ({len(df)} nuclei, {n_channels} channels)")


if __name__ == "__main__":
    main()
