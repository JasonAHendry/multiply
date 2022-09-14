import os
import sys
import pandas as pd
import numpy as np

from multiply.util.printing import print_header, print_footer
from multiply.util.dirs import produce_dir
from .algorithms import PrimerDimerAlgorithm


def align(primer_csv):
    """
    Run a pairwise alignment between all primers in `primer_csv` to identify
    potential primer dimers

    """

    # Parameters
    SAVE_TOP = None

    # PARSE CLI
    t0 = print_header("MULTIPLY: Perform pairwise alignment to identify possible primer dimers")
    input_dir = os.path.dirname(primer_csv)
    output_dir = produce_dir(input_dir, "align")
    print("Parsing inputs...")
    print(f"  Primer CSV: {primer_csv}")
    print(f"  Output directory: {output_dir}")
    print("Done.\n")

    # LOAD DATA
    print("Loading primers...")
    primer_df = pd.read_csv(primer_csv)
    primer_df.reset_index(inplace=True, drop=True)  # precaution, to ensure ordering
    n_primers = primer_df.shape[0]
    print(f"  Found {n_primers} primers.")
    print("Done.\n")

    # SET MODEL
    model = PrimerDimerAlgorithm()
    model.load_parameters()

    # COMPUTE PAIRWISE
    # NB: Essential to keep track of ordering here.
    print("Computing pairwise alignments...")
    alignments = []
    pairwise_scores = np.zeros((n_primers, n_primers))
    fmt_str = "  {:<4} {:<14} {:4>}/{:<4}"
    print("  {:<4} {:<14} {:<}".format("#", "Primer", "Completed"))
    for i in range(n_primers):

        # Extract first primer sequecne
        primer1_seq, primer1_name = primer_df.iloc[i][["seq", "primer_name"]]

        for j in range(i, n_primers):

            # Extract second primer sequence
            primer2_seq, primer2_name = primer_df.iloc[j][["seq", "primer_name"]]

            # Align
            model.set_primers(primer1_seq, primer2_seq, primer1_name, primer2_name)
            model.align()

            # Save
            alignments.append(model.get_primer_alignment())
            pairwise_scores[i, j] = model.score
            pairwise_scores[j, i] = model.score

            # Print
            sys.stdout.write("\r")
            sys.stdout.flush()
            sys.stdout.write(fmt_str.format(i + 1, primer1_name, j + 1, n_primers))
        sys.stdout.write("\r")
        sys.stdout.flush()
    print("\nDone.\n")

    # SAVE AS CSV
    print("Saving outputs...")
    matrix_output_csv = f"{output_dir}/matrix.pairwise_scores.csv"
    pairwise_df = pd.DataFrame(
        pairwise_scores,
        index=primer_df["primer_name"],
        columns=primer_df["primer_name"],
    )
    pairwise_df.to_csv(matrix_output_csv)
    print(f"  Pairwise interaction matrix: {matrix_output_csv}")

    # ADDITIONAL OUTPUTS FOR HIGH-SCORING ALIGNMENTS
    if SAVE_TOP is None:
        SAVE_TOP = len(alignments)
    alignments.sort()
    alignment_df = pd.DataFrame([a for a in alignments[:SAVE_TOP]])
    alignment_df.insert(0, "rank", range(SAVE_TOP))
    alignment_df.to_csv(f"{output_dir}/table.alignment_scores.csv", index=False)
    with open(f"{output_dir}/alignment_diagrams.txt", "w") as fn:
        for ix, a in enumerate(alignments[:SAVE_TOP]):
            fn.write(f"Alignment Rank: {ix:05d}\n")
            fn.write(f"{a.alignment}\n\n")
    print("Done.\n")

    print_footer(t0)
