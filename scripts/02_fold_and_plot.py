#!/usr/bin/env python3
"""
02_fold_and_plot.py

Purpose
-------
Takes the raw per-base depth from the PADDED reference (length N+K) and
"folds" it back onto the true circular genome (length N):

    folded[pos] = padded_depth[pos]                      for pos in (K, N]
    folded[pos] = padded_depth[pos] + padded_depth[N+pos] for pos in [1, K]

This adds back, onto positions 1..K, the depth from reads/pairs that were
placed in the duplicated tail because they truly span the biological
origin. The result is an artifact-free per-base depth profile of the real
genome, which is what gets reported and plotted.

If the unpadded control depth file (prelim_filtered_depth.txt) is present,
the script also validates the correction: it compares the corrected
profile with the control across the genome interior -- positions
K+1 .. N-K, which padding does not touch and where the control is free of
the origin artifact. If the file is absent, this step is skipped with a
note and everything else runs as normal.

How to run
----------
This version is written to be run directly from the Windows Python IDLE
shell (Run > Run Module, or F5) -- no command-line arguments needed.

Before running, put these files in the SAME FOLDER as this script:
    padded_depth.txt            Step 6 of the protocol (Phase 5 of run_pipeline.sh)
    padding_metadata.txt        Step 3 of the protocol (Phase 2 of run_pipeline.sh)
    prelim_filtered_depth.txt   Step 7 of the protocol (Phase 6 of run_pipeline.sh);
                                optional, needed only for the control comparison
When you run it, the program will ask you to type an output prefix (e.g.
MySpecies_sample01) -- this is used to name all the output files.

It only needs pandas, matplotlib, and numpy -- no bioinformatics tools --
so it also runs fine inside WSL if you'd rather use that instead.

Outputs (written next to this script, named using your output prefix)
-----------------------------------------------------------------------
    <prefix>_folded_depth.txt        per-base depth table for the true genome (Pos, Depth)
    <prefix>_coverage_stats.txt      mean/sd/min/max/breadth computed from folded depth
    <prefix>_linear_coverage.png     linear coverage plot
    <prefix>_circular_coverage.png   circular/polar coverage plot
    <prefix>_control_comparison.txt  interior comparison with the unpadded control
                                     (only if prelim_filtered_depth.txt is present)

03_compare_to_control.py produces the same comparison as a standalone
step, for depth profiles folded with an earlier version of this script.

Provenance: the statistics and comparison files end with a provenance block
recording the pipeline version, date, Python and library versions, and MD5
checksums of the input files (and, via padding_metadata.txt, of the
assembly). The values above the block are unchanged by it.
"""

import hashlib
import os
import math
import platform
from datetime import datetime

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

__version__ = "4.1.1"   # coverage pipeline version

# ============================= EDIT THESE ===================================
PLOT_TITLE = "Per-base Read Coverage"     # set to the full title you want shown, e.g. "Per-base Read Coverage: Genus species"
LINEAR_COLOR = "#648FFF"
MEAN_LINE_COLOR = "#FE6100"
# =============================================================================

# Input files are expected in the same folder as this script.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEPTH_FILENAME = "padded_depth.txt"
METADATA_FILENAME = "padding_metadata.txt"
CONTROL_FILENAME = "prelim_filtered_depth.txt"


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


def fold_depth(padded_depth, n, k):
    """padded_depth: 1D numpy array, 0-indexed array where index i
    corresponds to 1-based genomic position i+1 on the padded reference.
    Returns folded depth array of length n (index i -> position i+1)."""
    if len(padded_depth) != n + k:
        raise ValueError(
            f"padded_depth has {len(padded_depth)} rows but metadata says "
            f"N+K = {n + k}. Did you point this at the right depth file?"
        )
    folded = padded_depth[:n].copy()
    # positions 1..K (0-indexed 0..K-1) get the padded tail added back
    folded[:k] = folded[:k] + padded_depth[n:n + k]
    return folded


def compare_to_control(folded, control, n, k, folded_label):
    """Compare corrected and control profiles across the interior
    (1-based positions K+1 .. N-K). Returns the report as a list of lines,
    or None (with a printed note) if the comparison cannot be made.
    The comparison lines are identical to those of 03_compare_to_control.py."""
    if len(control) != n:
        print(f"NOTE: {CONTROL_FILENAME} has {len(control)} positions, but "
              f"metadata says N = {n}. Files are from different runs? "
              f"Control comparison skipped.")
        return None
    if n - 2 * k < 1:
        print(f"NOTE: genome (N = {n}) too short for an interior with "
              f"K = {k}. Control comparison skipped.")
        return None

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

    return [
        f"Folded profile : {folded_label}",
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


def main():
    depth_path = os.path.join(SCRIPT_DIR, DEPTH_FILENAME)
    metadata_path = os.path.join(SCRIPT_DIR, METADATA_FILENAME)
    control_path = os.path.join(SCRIPT_DIR, CONTROL_FILENAME)

    missing = [p for p in (depth_path, metadata_path) if not os.path.isfile(p)]
    if missing:
        print("ERROR: could not find the following required file(s) in this")
        print(f"script's folder ({SCRIPT_DIR}):")
        for p in missing:
            print(f"  - {os.path.basename(p)}")
        print("Copy padded_depth.txt and padding_metadata.txt into this folder")
        print("(from Steps 6 and 3 of the protocol) and run this again.")
        return

    out_prefix = input(
        "Enter an output prefix for this sample (e.g. MySpecies_sample01): "
    ).strip()
    if not out_prefix:
        out_prefix = "coverage_output"
        print(f"No prefix entered -- using default: '{out_prefix}'")

    meta = read_metadata(metadata_path)
    n = int(meta["original_length_N"])
    k = int(meta["padding_length_K"])
    header = meta.get("original_header", "assembly")

    data = pd.read_csv(depth_path, sep="\t", header=None,
                        names=["Ref", "Pos", "Depth"])
    padded_depth = data["Depth"].to_numpy()

    folded = fold_depth(padded_depth, n, k)
    positions = np.arange(1, n + 1)

    # ---- write folded depth table ----
    folded_df = pd.DataFrame({"Pos": positions, "Depth": folded})
    folded_name = f"{out_prefix}_folded_depth.txt"
    folded_out_path = os.path.join(SCRIPT_DIR, folded_name)
    folded_df.to_csv(folded_out_path, sep="\t", index=False)

    # ---- stats ----
    mean = float(np.mean(folded))
    sd = float(np.std(folded))
    dmin = int(np.min(folded))
    dmax = int(np.max(folded))
    breadth = int(np.sum(folded > 0))

    stats_out_path = os.path.join(SCRIPT_DIR, f"{out_prefix}_coverage_stats.txt")
    with open(stats_out_path, "w") as fh:
        fh.write(f"Reference: {header}\n")
        fh.write(f"Genome length (N): {n} bp\n")
        fh.write(f"Padding used (K): {k} bp\n")
        fh.write(f"Mean Depth: {mean:.4f}\n")
        fh.write(f"Std Dev: {sd:.4f}\n")
        fh.write(f"Min: {dmin}\n")
        fh.write(f"Max: {dmax}\n")
        fh.write(f"Breadth: {breadth} bp ({100 * breadth / n:.2f}% of genome)\n")
        fh.write("\n".join(provenance_lines("02_fold_and_plot.py",
                                             [depth_path, metadata_path], meta)) + "\n")

    print(f"Mean Depth: {mean:.2f}")
    print(f"Std Dev: {sd:.2f}")
    print(f"Min: {dmin}")
    print(f"Max: {dmax}")
    print(f"Breadth: {breadth} bp / {n} bp")
    print(f"Folded depth table -> {folded_out_path}")
    print(f"Stats -> {stats_out_path}")

    # =========================== LINEAR PLOT ================================
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.fill_between(positions, folded, color=LINEAR_COLOR, alpha=0.3)
    ax.plot(positions, folded, color=LINEAR_COLOR, lw=0.8)
    ax.axhline(y=mean, color=MEAN_LINE_COLOR, linestyle="--", lw=1.2,
               label=f"Mean: {mean:.1f}×")
    ax.set_title(PLOT_TITLE, fontsize=16)
    ax.set_xlabel("Genomic Position (bp)", fontsize=11)
    ax.set_ylabel("Read Depth", fontsize=11)
    ax.legend(fontsize=10)
    ax.set_xlim(positions.min(), positions.max())
    plt.tight_layout()
    linear_out_path = os.path.join(SCRIPT_DIR, f"{out_prefix}_linear_coverage.png")
    plt.savefig(linear_out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Linear plot -> {linear_out_path}")

    # ========================== CIRCULAR PLOT ================================
    theta = 2 * math.pi * positions / n
    fig2 = plt.figure(figsize=(7, 7))
    ax2 = fig2.add_subplot(111, projection="polar")
    ax2.plot(theta, folded, color=LINEAR_COLOR, lw=0.8)
    ax2.fill(theta, folded, color=LINEAR_COLOR, alpha=0.3)
    # mean depth reference circle
    mean_circle = np.full_like(theta, mean)
    ax2.plot(theta, mean_circle, color=MEAN_LINE_COLOR, linestyle="--", lw=1.2,
              label=f"Mean: {mean:.1f}×")
    ax2.set_theta_zero_location("N")
    ax2.set_theta_direction(-1)
    ax2.set_title(f"{PLOT_TITLE} (circular)", fontsize=12, pad=20)
    ax2.set_xlabel("")
    ax2.legend(loc="upper right", bbox_to_anchor=(1.15, 1.1), fontsize=9)
    plt.tight_layout()
    circular_out_path = os.path.join(SCRIPT_DIR, f"{out_prefix}_circular_coverage.png")
    plt.savefig(circular_out_path, dpi=300, bbox_inches="tight")
    plt.close(fig2)
    print(f"Circular plot -> {circular_out_path}")

    # ===================== CONTROL COMPARISON (VALIDATION) ===================
    print("")
    if not os.path.isfile(control_path):
        print(f"NOTE: {CONTROL_FILENAME} not found next to this script, so the")
        print("control comparison was skipped. To run it, copy that file here")
        print("(from Phase 6 of run_pipeline.sh / Step 7 of the protocol) and")
        print("run this script again.")
        return

    control = pd.read_csv(control_path, sep="\t", header=None,
                          names=["Ref", "Pos", "Depth"])["Depth"].to_numpy()
    lines = compare_to_control(folded, control, n, k, folded_name)
    if lines is None:
        return
    comparison_out_path = os.path.join(SCRIPT_DIR, f"{out_prefix}_control_comparison.txt")
    with open(comparison_out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
        fh.write("\n".join(provenance_lines("02_fold_and_plot.py",
                                             [depth_path, metadata_path, control_path],
                                             meta)) + "\n")
    print("\n".join(lines))
    print(f"Control comparison -> {comparison_out_path}")


if __name__ == "__main__":
    main()
