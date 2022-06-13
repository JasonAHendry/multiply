import os
import click
import pandas as pd

from .blast import make_blast_db, run_blast
from multiply.download.collection import genome_collection
from multiply.util.dirs import produce_dir
from multiply.util.io import write_fasta_from_dict


@click.command(short_help="BLAST candidate primers.")
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
def blast(primer_csv, genome_name):
    """
    Search for off-target binding sites and amplicons using blast

    """

    # PARSE CLI
    input_dir = os.path.dirname(primer_csv)
    output_dir = produce_dir(input_dir, "blast")
    genome = genome_collection[genome_name]

    # LOAD DATA
    primer_df = pd.read_csv(primer_csv)
    n_primers = primer_df.shape[0]

    # WRITE PRIMER FASTA (for blast input)
    primer_dt = dict(zip(primer_df["primer_name"], primer_df["seq"]))
    primer_fasta = f"{output_dir}/cadidate_primer.fasta"
    write_fasta_from_dict(
        input_dt=primer_dt,
        output_fasta=primer_fasta
    )

    # CHECK FOR / CREATE BLAST DATABASE (will probably want this in download)
    print("  Making blast database...")
    blast_db = make_blast_db(genome.fasta_path)

    # RUN BLAST
    print("  Runnning blast...")
    blast_output=f"{output_dir}/table.blast_results.output"
    run_blast(
        db_path=blast_db,
        input_fasta=primer_fasta,
        output_blast=blast_output,
        output_fmt=6
    )

    print("Done.")
    # Look for off-target alignments more intelligently
    # - (1) Find probable binding sites (3') end
    # - (2) See if any of those produce amplicon
    # - (3) See if any of those interact with *target* amplicons


    # CREATE SUMMARY DATAFRAME(S) OF OUTPUT
