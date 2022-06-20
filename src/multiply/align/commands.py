import click
import os
import pandas as pd
import numpy as np

from multiply.util.dirs import produce_dir
from .algorithms import PrimerDimerAlgorithm


# ================================================================================
# Main function wrapped for Click CLI
#
# ================================================================================


@click.command(short_help="Search for primer dimers.")
@click.option(
    "-p",
    "--primer_csv",
    type=click.Path(exists=True),
    required=True,
    help="Path to candidate primer CSV file (e.g. `table.candidate_primers.csv`).",
)
def align(primer_csv):
    """
    Run a pairwise alignment between all primers in `primer_csv` to identify
    potential primer dimers

    """
    main(primer_csv)


# ================================================================================
# Main function, unwrapped
#
# ================================================================================


def main(primer_csv):
    # PARSE CLI
    input_dir = os.path.dirname(primer_csv)
    output_dir = produce_dir(input_dir, "align")

    # LOAD DATA
    primer_df = pd.read_csv(primer_csv)
    primer_df.reset_index(inplace=True, drop=True) # precaution, to ensure ordering
    n_primers = primer_df.shape[0]

    # SET MODEL
    model = PrimerDimerAlgorithm()
    model.load_parameters()

    # COMPUTE PAIRWISE
    # - Note how essential ordering is here
    print("Computing alignments...")
    alignments = []
    pairwise_scores = np.zeros((n_primers, n_primers))
    for i in range(n_primers):

        # Extract first primer sequecne
        primer1_seq, primer1_name = primer_df.iloc[i][["seq", "primer_name"]]

        for j in range(i, n_primers):

            # Extract second primer sequence
            primer2_seq, primer2_name = primer_df.iloc[j][["seq", "primer_name"]]

            # Align
            #print(primer1_name, primer2_name)
            model.set_primers(primer1_seq, primer2_seq, primer1_name, primer2_name)
            model.align()

            # Save
            alignments.append(model.get_primer_alignment())
            pairwise_scores[i, j] = model.score
            pairwise_scores[j, i] = model.score

    # SAVE AS NPY
    np.save(f"{output_dir}/matrix.pairwise_scores.npy", pairwise_scores)

    # SAVE AS CSV
    pairwise_df = pd.DataFrame(
        pairwise_scores, 
        index=primer_df["primer_name"],
        columns=primer_df["primer_name"]
    )
    pairwise_df.to_csv(f"{output_dir}/matrix.pairwise_scores.csv")

    # ADDITIONAL OUTPUTS FOR HIGH-SCORING ALIGNMENTS
    # Want to actually iterate and print these alignments...
    SAVE_TOP = 200
    alignments.sort()
    alignment_df = pd.DataFrame([a for a in alignments[:SAVE_TOP]])
    alignment_df.to_csv(f"{output_dir}/table.alignment_scores.csv", index=False)
    with open(f"{output_dir}/alignment_diagrams.txt", "w") as fn:
        for ix, a in enumerate(alignments[:SAVE_TOP]):
            fn.write(f"Alignment Index: {ix:05d}\n")
            fn.write(f"{a.alignment}\n\n")

    # OPTIONALLY -- visualise matrix
    print("Done.")
