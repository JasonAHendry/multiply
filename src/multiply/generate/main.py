import pandas as pd

from multiply.util.printing import print_header, print_footer, print_parameters
from multiply.util.parsing import parse_parameters
from multiply.util.dirs import create_output_directory, produce_dir
from multiply.util.io import load_bed_as_dataframe
from multiply.download.collection import genome_collection
from multiply.generate.targets import Target, TargetSet
from multiply.generate.primer3 import Primer3Runner
from multiply.generate.primers import load_primer_pairs_from_primer3_output


def generate(design):
    """
    Generate a pool of candidate primers for a given `design` using
    primer3

    Visit the `settings/primer3` directory to observe or change primer3
    settings.

    """
    # PARSE CLI
    t0 = print_header("MULTIPLY: Generate candidate primers for all targets using Primer3")
    params = parse_parameters(design)
    genome = genome_collection[params["genome"]]
    params["output_dir"] = create_output_directory(params)
    print_parameters(design, params)

    # EXTRACT GENE INFORMATION
    print("Preparing targets...")
    genes = []
    if params["from_genes"]:
        gene_df = pd.read_csv(genome.gff_path)
        target_ids = params["target_ids"]
        gene_df.query("ID in @target_ids", inplace=True)
        genes = [Target.from_series(row) for _, row in gene_df.iterrows()]

        # Pretty ugly, would be nice to encapsulate
        if params["has_names"]:
            for gene in genes:
                gene.name = params["target_id_to_name"][gene.ID]
    print(f"  Found {len(genes)} gene(s).")

    # EXTRACT REGION INFORMATION
    regions = []
    if params["from_regions"]:
        region_df = load_bed_as_dataframe(params["region_bed"])
        regions = [Target.from_series(row) for _, row in region_df.iterrows()]
    print(f"  Found {len(regions)} region(s).")

    # MERGE
    print("  Merging genes and regions...")
    targets = genes + regions
    target_set = (
        TargetSet(targets)
        .check_size_compatible(params["max_size_bp"])
        .calc_pads()
        .extract_seqs(genome.fasta_path, include_pads=True)
        .to_csv(f"{params['output_dir']}/table.targets_overview.csv")
        .to_fasta(f"{params['output_dir']}/targets_sequence.fasta")
    )
    print("Done.\n")

    # RUN PRIMER3
    print("Running primer3...")
    primer3_output_dir = produce_dir(params["output_dir"], "primer3")
    primer3_runner = Primer3Runner()

    # Storage
    primer_pair_dt = {target.ID: [] for target in target_set.targets}

    # Iterate over settings
    for primer3_setting in params["primer3_settings"]:
        print(f"  Generating primers using {primer3_setting} settings...")

        primer3_runner.load_primer3_settings(primer3_setting)
        primer3_runner.set_amplicon_size_ranges(
            min_size_bp=params["min_size_bp"], max_size_bp=params["max_size_bp"]
        )

        # Iterate over targets
        for target in target_set.targets:
            primer3_runner.set_target(
                ID=target.ID,
                seq=target.seq,
                start=target.start,
                pad_start=target.pad_start,
                length=target.length,
            )
            primer3_runner.run(output_dir=primer3_output_dir)

            # Store
            primer_pairs = load_primer_pairs_from_primer3_output(
                primer3_runner.output_path, add_target=target
            )
            primer_pair_dt[target.ID].extend(primer_pairs)
    print("Done.\n")

    # REDUCE TO UNIQUE PAIRS
    print("Primer pairs discovered by primer3:")
    print(f"  {'Target':<15} {'Total':<10} {'Unique':<10}")
    for target_id, all_primer_pairs in primer_pair_dt.items():

        # Reduce to unique primer pairs, and give names
        # this definitely can be cleaned
        uniq_primer_pairs = set(all_primer_pairs)
        for ix, pair in enumerate(uniq_primer_pairs):
            pair.give_primers_names(primer_code=params["primer_code"], primer_ix=ix)

        # Print summary
        print(
            f"  {target_id:<15} {len(all_primer_pairs):<10} {len(uniq_primer_pairs):<10}"
        )

        # Store
        primer_pair_dt[target_id] = uniq_primer_pairs
    print("Done.\n")

    # Add tails, if provided
    if params["include_tails"]:
        print("Adding tails to discovered primers...")
        for target_id, primer_pairs in primer_pair_dt.items():
            for primer_pair in primer_pairs:
                primer_pair.F.add_tail(params["F_tail"])
                primer_pair.R.add_tail(params["R_tail"])
        print("Done.\n")

    # WRITE
    print("Writing output table...")
    primer_df = pd.DataFrame(
        [
            pair.get_primer_as_dict(direction)
            for _, primer_pairs in primer_pair_dt.items()
            for pair in primer_pairs
            for direction in ["F", "R"]
        ]
    )
    primer_df = primer_df[
        [
            "target_id",
            "target_name",
            "pair_name",
            "primer_name",
            "direction",
            "seq",
            "length",
            "tm",
            "gc",
            "chrom",
            "start",
            "product_bp",
            "pair_penalty",
        ]
    ]
    output_csv = f"{params['output_dir']}/table.candidate_primers.csv"
    primer_df.to_csv(output_csv, index=False)
    print(f"  to: {output_csv}")
    print("Done.\n")

    print_footer(t0)
