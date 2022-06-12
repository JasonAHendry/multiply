import click
import os
import pandas as pd
import numpy as np

from multiply.util.dirs import produce_dir
from .algorithms import PrimerDimerAlgorithm


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
    Run a pairwise alignment algorithm between all candidate primers,
    scoring the likelihood of them producing a dimer

    """
    # PARSE CLI
    input_dir = os.path.dirname(primer_csv)
    output_dir = produce_dir(input_dir, "align")

    # LOAD DATA
    primer_df = pd.read_csv(primer_csv)
    n_primers = primer_df.shape[0]

    # SET MODEL
    model = PrimerDimerAlgorithm()

    # COMPUTE PAIRWISE
    # - Note how essential ordering is here
    alignments = []
    pairwise_scores = np.zeros((n_primers, n_primers))
    for i in range(n_primers):
        for j in range(i, n_primers):
            # Extract sequences
            primer1_seq = primer_df.iloc[i]["seq"]
            primer1_name = primer_df.iloc[i]["primer_name"]
            primer2_seq = primer_df.iloc[j]["seq"]
            primer2_name = primer_df.iloc[j]["primer_name"]

            # Align
            model.set_primers(primer1_seq, primer2_seq, primer1_name, primer2_name)
            model.align()

            # Save
            alignments.append(model.get_primer_alignment())
            pairwise_scores[i, j] = model.score
            pairwise_scores[j, i] = model.score

    # ADDITIONAL OUTPUTS FOR HIGH-SCORING ALIGNMENTS
    alignments.sort()
    alignment_df = pd.DataFrame([a for a in alignments[:100]])
    alignment_df.to_csv(f"{output_dir}/table.alignment_scores.csv", index=False)
    np.save(f"{output_dir}/matrix.pairwise_scores.npy", pairwise_scores)

    # OPTIONALLY -- visualise matrix
    print("Done.")
