import os
import click

import pandas as pd

from multiply.util.dirs import produce_dir
from multiply.util.io import write_primers_to_bed
from multiply.download.collection import genome_collection
from .bedtools import bedtools_intersect


# ================================================================================
# Main function wrapped for Click CLI
#
# ================================================================================


@click.command(short_help="Check candidate primers for SNPs.")
@click.option(
    "-p",
    "--primer_csv",
    type=click.Path(exists=True),
    required=True,
    help="Path to candidate primer CSV file (e.g. `table.candidate_primers.csv`).",
)
@click.option(
    "-g",
    "--genome_name",
    type=click.Choice(genome_collection),
    required=True,
    help="Name of genome to download.",
)
def snpcheck(primer_csv, genome_name):
    """
    Check primers in `primer_csv` for the presence of any SNPs, as defined
    by the `genome_name`

    Variation data for a given `genome_name` should be found in the
    `genomes/variation` folder.

    """
    main(primer_csv, genome_name)


# ================================================================================
# Main function, unwrapped
#
# ================================================================================


def main(primer_csv, genome_name):
    # PARSE CLI
    input_dir = os.path.dirname(primer_csv)
    output_dir = produce_dir(input_dir, "snpcheck")
    genome = genome_collection[genome_name]

    # LOAD DATA
    primer_df = pd.read_csv(primer_csv)

    # WRITE PRIMERS TO BED
    bed_path = f"{output_dir}/candidate_primers.bed"
    write_primers_to_bed(primer_df, bed_path)

    # CHECK FOR VARIATION DATA
    # ITERATE OVER VARIATION FILES
    # RUN BEDTOOLS INTERSECT -C
    # CONSTRUCT OUTPUT DATAFRAME
