import click
import os
import pandas as pd
from multiply.view.plot import (
    SequencePlotter,
    GffPlotter,
    PrimerPlotter,
    CombinedPlotter,
)
from multiply.download.collection import genome_collection
from multiply.download.gff import load_gff
from multiply.util.io import load_fasta_as_dict
from multiply.util.dirs import produce_dir


@click.command(short_help="Visualise candidate primer locations.")
@click.option(
    "-i",
    "--input_dir",
    type=click.Path(exists=True),
    required=True,
    help="Input directory containing: `table.candidate_primers.csv`, `table.targets_overview.csv`, and `targets_sequence.fasta`.",
)
@click.option(
    "-g",
    "--genome_name",
    type=click.Choice(genome_collection),
    required=True,
    help="Name of genome to download.",
)
def view(input_dir, genome_name):
    """
    Visualise all candidate primers and amplicons for each target

    """
    # DEFINE AND CHECK FILE PATHS
    primer_path = f"{input_dir}/table.candidate_primers.csv"
    targets_path = f"{input_dir}/table.targets_overview.csv"
    fasta_path = f"{input_dir}/targets_sequence.fasta"
    for path in [primer_path, targets_path, fasta_path]:
        if not os.path.exists(path):
            raise FileNotFoundError(
                "Can't find {path}. Ensure it exists in input directory."
            )

    # Create output directory
    output_dir = produce_dir(input_dir, "view")

    # LOAD
    # Multiplex data
    primer_df = pd.read_csv(primer_path)
    targets_df = pd.read_csv(targets_path)
    seqs = load_fasta_as_dict(fasta_path)

    # Gff
    genome = genome_collection[genome_name]
    gff_df = load_gff(genome.gff_raw_download)

    # ITERATE OVER TARGETS, PLOT
    for target_id, row in targets_df.groupby("ID"):

        # Extract information
        target_primer_df = primer_df.query("target_id == @target_id")
        target_seq = [
            seq for header, seq in seqs.items() if header.startswith(f"ID={target_id}")
        ][0]

        # Prepare plotters
        seq_plotter = SequencePlotter(target_seq)
        gff_plotter = GffPlotter(
            gff=gff_df,
            chrom=str(row["chrom"].values[0]),
            start=int(row["pad_start"]),
            end=int(row["pad_end"]),
        )
        primer_plotter = PrimerPlotter(target_primer_df)

        # Compose
        comb_plotter = CombinedPlotter(
            sequence_plotter=seq_plotter,
            gff_plotter=gff_plotter,
            primer_plotter=primer_plotter,
        )

        # Plot
        comb_plotter.plot(
            start=int(row["pad_start"]), 
            end=int(row["pad_end"]),
            output_path=f"{input_dir}/view/{target_id}.pdf"
        )
