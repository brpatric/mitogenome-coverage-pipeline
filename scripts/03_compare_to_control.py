#!/usr/bin/env python3
"""
03_compare_to_control.py

Purpose
-------
Validates the padding-and-folding correction by comparing the corrected
(folded) depth profile against the unpadded control mapping, across the
genome interior -- positions K+1 .. N-K, which padding does not touch and
where the unpadded control is free of the origin artifact.

Produces the validation result reported in the Results:
  - percent difference in interior mean depth
  - minimum interior depth, and its position, in each profile

Since pipeline version 4, 02_fold_and_plot.py runs this same comparison
automatically when prelim_filtered_depth.txt is present, and reports
identical values. This standalone script is kept for depth profiles
that were folded with an earlier version of 02_fold_and_plot.py.

How to run
----------
Like 02_fold_and_plot.py, this runs directly from Windows Python IDLE
(Run > Run Module, or F5) or from WSL with  python3 03_compare_to_control.py

Put these three files in the SAME FOLDER as this script:
    <prefix>_folded_depth.txt     from 02_fold_and_plot.py
    prelim_filtered_depth.txt     from run_pipeline.sh (Phase 7)
    padding_metadata.txt          from run_pipeline.sh (coverage_run/)

Output
------
    <prefix>_control_comparison.txt   (also printed to screen)
"""

import glob
import hashlib
import os
import platform
from datetime import datetime

import matplotlib
import numpy as np
import pandas as pd

__version__ = "4.4"   # coverage pipeline version

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONTROL_FILENAME = "prelim_filtered_depth.txt"
METADATA_FILENAME = "padding_metadata.txt"


def read_metadata(path):
    meta = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                meta[parts[0]] = parts[1]
    return meta


def md5_of_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance_lines(script_name, input_paths, meta):
    """Provenance block appended to each reported output file: what
    produced it, when, with which library versions, and from which inputs
    (MD5 checksums). The values above the block are unaffected."""
    lines = [
        "",
        "--- Provenance ---",
        f"Produced by       : {script_name} (pipeline version {__version__})",
        f"Date              : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Python            : {platform.python_version()}",
        f"NumPy / pandas / Matplotlib : {np.__version__} / {pd.__version__} / {matplotlib.__version__}",
    ]
    if "assembly_md5" in meta:
        lines.append(f"Assembly (MD5)    : {meta.get('assembly_file', 'assembly')} {meta['assembly_md5']}")
    lines.append(f"Padding built by  : pipeline version {meta.get('pipeline_version', 'earlier than 4.1 (not recorded)')}")
    for p in input_paths:
        lines.append(f"Input (MD5)       : {os.path.basename(p)} {md5_of_file(p)}")
    return lines


def pick_folded_file():
    hits = sorted(glob.glob(os.path.join(SCRIPT_DIR, "*_folded_depth.txt")))
    if not hits:
        return None
    if len(hits) == 1:
        return hits[0]
    print("More than one *_folded_depth.txt found:")
    for i, h in enumerate(hits, 1):
        print(f"  {i}. {os.path.basename(h)}")
    choice = input("Enter the number of the one to use: ").strip()
    return hits[int(choice) - 1]


def main():
    control_path = os.path.join(SCRIPT_DIR, CONTROL_FILENAME)
    metadata_path = os.path.join(SCRIPT_DIR, METADATA_FILENAME)
    folded_path = pick_folded_file()

    missing = []
    if folded_path is None:
        missing.append("<prefix>_folded_depth.txt")
    for p in (control_path, metadata_path):
        if not os.path.isfile(p):
            missing.append(os.path.basename(p))
    if missing:
        print(f"ERROR: missing from {SCRIPT_DIR}:")
        for m in missing:
            print(f"  - {m}")
        return

    meta = read_metadata(metadata_path)
    n = int(meta["original_length_N"])
    k = int(meta["padding_length_K"])

    folded = pd.read_csv(folded_path, sep="\t")["Depth"].to_numpy()
    control = pd.read_csv(control_path, sep="\t", header=None,
                          names=["Ref", "Pos", "Depth"])["Depth"].to_numpy()

    for name, arr in (("folded", folded), ("control", control)):
        if len(arr) != n:
            print(f"ERROR: {name} profile has {len(arr)} positions, "
                  f"but metadata says N = {n}. Files are from different runs?")
            return
    if n - 2 * k < 1:
        print(f"ERROR: genome (N = {n}) too short for an interior with K = {k}.")
        return

    # interior = 1-based positions K+1 .. N-K  ->  0-based slice [k : n-k]
    lo, hi = k + 1, n - k
    f_int = folded[k:n - k]
    c_int = control[k:n - k]

    f_mean = float(np.mean(f_int))
    c_mean = float(np.mean(c_int))
    pct_diff = 100 * (f_mean - c_mean) / c_mean

    f_min = int(np.min(f_int))
    c_min = int(np.min(c_int))
    f_min_pos = np.flatnonzero(f_int == f_min) + lo
    c_min_pos = np.flatnonzero(c_int == c_min) + lo
    same_min = (f_min == c_min) and np.array_equal(f_min_pos, c_min_pos)

    def fmt_pos(arr):
        shown = ", ".join(str(p) for p in arr[:10])
        return shown + (f" (+{len(arr) - 10} more)" if len(arr) > 10 else "")

    prefix = os.path.basename(folded_path)[:-len("_folded_depth.txt")]
    lines = [
        f"Folded profile : {os.path.basename(folded_path)}",
        f"Control profile: {CONTROL_FILENAME}",
        f"N = {n} bp, K = {k} bp",
        f"Interior compared: positions {lo}-{hi} ({hi - lo + 1} bp)",
        "",
        f"Interior mean depth, corrected : {f_mean:.4f}",
        f"Interior mean depth, control   : {c_mean:.4f}",
        f"Difference (corrected - control): {pct_diff:+.3f}%",
        "",
        f"Interior minimum, corrected : {f_min} at position(s) {fmt_pos(f_min_pos)}",
        f"Interior minimum, control   : {c_min} at position(s) {fmt_pos(c_min_pos)}",
        f"Same minimum at same position(s): {'YES' if same_min else 'NO'}",
    ]
    out_path = os.path.join(SCRIPT_DIR, f"{prefix}_control_comparison.txt")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
        fh.write("\n".join(provenance_lines("03_compare_to_control.py",
                                             [folded_path, metadata_path, control_path],
                                             meta)) + "\n")
    print("\n".join(lines))
    print(f"\nWritten to: {out_path}")


if __name__ == "__main__":
    main()
