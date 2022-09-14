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
from multiply.util.printing import print_header, print_footer


def view(result_dir, genome_name):
    """
    View binding locations of candidate primers

    """
    # DEFINE AND CHECK FILE PATHS
    t0 = print_header()
    print("Parsing inputs...")
    primer_path = f"{result_dir}/table.candidate_primers.csv"
    targets_path = f"{result_dir}/table.targets_overview.csv"
    fasta_path = f"{result_dir}/targets_sequence.fasta"
    for path in [primer_path, targets_path, fasta_path]:
        if not os.path.exists(path):
            raise FileNotFoundError(
                "Can't find {path}. Ensure it exists in input directory."
            )
    output_dir = produce_dir(result_dir, "view")
    print(f"  Primer CSV: {primer_path}")
    print(f"  Targets CSV: {targets_path}")
    print(f"  Targets FASTA: {fasta_path}")
    print(f"  Genome: {genome_name}")
    print("Done.\n")

    # LOAD
    # Multiplex data
    primer_df = pd.read_csv(primer_path)
    targets_df = pd.read_csv(targets_path)
    seqs = load_fasta_as_dict(fasta_path)

    # Gff
    print("Loading GFF for gene body plots...")
    genome = genome_collection[genome_name]
    gff_df = load_gff(genome.gff_raw_download)
    print("Done.\n")

    # ITERATE OVER TARGETS, PLOT
    print("Plotting primer locations for each target...")
    for target_id, row in targets_df.groupby("ID"):

        # Extract information
        print(f"  {target_id}")
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
            title=f"{target_id} | {row['name'].values[0]}",
            output_path=f"{result_dir}/view/{target_id}.pdf",
        )
    print("Done.\n")

    print(f"Plots can be found in directory: {output_dir}\n")

    print_footer(t0)
