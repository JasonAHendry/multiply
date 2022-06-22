import os
import click
import json
import pandas as pd

from .runner import BlastRunner
from .annotator import BlastResultsAnnotator
from .offtarget import AmpliconFinder
from multiply.download.collection import genome_collection
from multiply.util.printing import print_header, print_footer
from multiply.util.dirs import produce_dir
from multiply.util.io import write_fasta_from_dict


# ================================================================================
# Main function wrapped for Click CLI
#
# ================================================================================


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
    BLAST primers in `primer_csv` against reference genome given by `genome_name` to
    find off-target binding sites and amplicons

    Visit `settings/blast` for blast settings.

    """
    main(primer_csv, genome_name)


# ================================================================================
# Main function, unwrapped
#
# ================================================================================


def main(primer_csv, genome_name):
    # PARSE CLI
    t0 = print_header()
    input_dir = os.path.dirname(primer_csv)
    print("Parsing inputs...")
    output_dir = produce_dir(input_dir, "blast")
    genome = genome_collection[genome_name]
    print(f"  Primer CSV: {primer_csv}")
    print(f"  Output directory: {output_dir}")
    print(f"Done.\n")

    # LOAD DATA
    print("Loading data...")
    primer_df = pd.read_csv(primer_csv)
    print(f"  Found {primer_df.shape[0]} primers.")

    # LOAD PARAMATERS
    params = json.load(open("settings/blast/parameters.json", "r"))
    param_str = [f"{k}={v}" for k, v in params.items()]
    print(f"  BLAST parameters: {', '.join(param_str)}")
    print("Done.\n")

    # WRITE PRIMER FASTA, for BLAST input
    primer_dt = dict(zip(primer_df["primer_name"], primer_df["seq"]))
    primer_fasta = f"{output_dir}/cadidate_primer.fasta"
    write_fasta_from_dict(input_dt=primer_dt, output_fasta=primer_fasta)

    # RUN BLAST
    print("Running BLAST...")
    blast_df = (
        BlastRunner(input_fasta=primer_fasta, reference_fasta=genome.fasta_path)
        .create_database()
        .run(
            word_size=params["word_size"],
            output_archive=f"{output_dir}/blast.candidate_primers.archive",
        )
        .reformat_output_as_table(
            output_table=f"{output_dir}/blast.candidate_primers.table"
        )
        .get_dataframe()
    )
    print("Done.\n")

    # ANNOTATE RESULTS
    print("Annotating results...")
    annotator = BlastResultsAnnotator(blast_df)
    annotator.build_annotation_dict(
        length_threshold=params["alignment_length_threshold"],
        evalue_threshold=params["evalue_threshold"],
    )
    annotator.add_annotations()
    annotator.summarise_by_primer(
        f"{output_dir}/table.blast.candidate_primers.summary.csv"
    )
    print(f"  Found {blast_df.shape[0]} BLAST hits.")
    print("Done.\n")

    # PREDICT AMPLICONS
    print("Predicting off-target amplicons...")
    bound_df = annotator.get_predicted_bound()
    print(f"  Found {bound_df.shape[0]} putative primer binding sites.")
    amplicon_finder = AmpliconFinder(bound_df)
    amplicon_finder.find_amplicons()
    amplicon_df = amplicon_finder.amplicon_df
    print(f"  Yeilding {amplicon_df.shape[0]} potential amplicons.")
    amplicon_df.to_csv(f"{output_dir}/table.predicted_amplicons.csv")

    amplicon_finder.create_ontarget_dataframe(primer_df)
    amplicon_finder.create_offtarget_dataframe()

    amplicon_finder.ontarget_df.to_csv(f"{output_dir}/table.ontarget_amplicons.csv")
    amplicon_finder.offtarget_df.to_csv(f"{output_dir}/table.offtarget_amplicons.csv")
    print("Done.\n")

    # PLOT
    print_footer(t0)