import os
import click
import json
import pandas as pd

from .runner import BlastRunner
from .annotator import BlastResultsAnnotator, ANNOTATIONS, PrimerBlastRecord
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
    #n_primers = primer_df.shape[0]

    # LOAD PARAMATERS
    params = json.load(open("settings/blast/parameters.json", "r"))
    print(params)

    # WRITE PRIMER FASTA, for BLAST input
    primer_dt = dict(zip(primer_df["primer_name"], primer_df["seq"]))
    primer_fasta = f"{output_dir}/cadidate_primer.fasta"
    write_fasta_from_dict(input_dt=primer_dt, output_fasta=primer_fasta)

    # RUN BLAST
    blast_df = (
        BlastRunner(input_fasta=primer_fasta, reference_fasta=genome.fasta_path)
        .create_database()
        .run(
            word_size=params["word_size"],
            output_archive=f"{output_dir}/blast.candidate_primers.archive")
        .reformat_output_as_table(
            output_table=f"{output_dir}/blast.candidate_primers.table"
        )
        .get_dataframe()
    )

    # ANNOTATE RESULTS
    annotator = BlastResultsAnnotator(blast_df)
    annotator.build_annotation_dict(
        length_threshold=params["alignment_length_threshold"],
        evalue_threshold=params["evalue_threshold"]
    )
    annotator.add_annotations()
    annotator.summarise_by_primer(f"{output_dir}/table.blast.candidate_primers.summary.csv")



    # # Insert additional columns
    # for annot_name, annot_func in ANNOTATIONS.items():
    #     blast_df.insert(
    #         blast_df.shape[1], annot_name, blast_df.apply(func=annot_func, axis=1)
    #     )

    # # Produce summary of alignments for each primer
    # primer_blast_records = [
    #     PrimerBlastRecord(
    #         primer_name=qseqid,
    #         primer_pair_name=qseqid[:-2],
    #         target_name=qseqid.split("_")[0],
    #         total_alignments=qseqid_df.shape[0],
    #         **qseqid_df[ANNOTATIONS].sum().to_dict(),
    #     )
    #     for qseqid, qseqid_df in blast_df.groupby("qseqid")
    # ]
    # blast_primer_df = (
    #     pd.DataFrame(primer_blast_records)
    #     .sort_values("predicted_bound", ascending=False)
    #     .to_csv(f"{output_dir}/table.blast.candidate_primers.summary.csv")
    # )



    # Look for off-target alignments more intelligently
    # - (1) Find probable binding sites (3') end
    # - (2) See if any of those produce amplicon(s)
    # - (3) See if any of those interact with *target* amplicons

    # CREATE SUMMARY DATAFRAME(S) OF OUTPUT

    # PLOT / VISUALISE
    # - What would *actually* be a nice visualisation?
