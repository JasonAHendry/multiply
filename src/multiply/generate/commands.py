import click
import pandas as pd
#from multiply.cli import design_path_option
from multiply.util.dirs import create_output_directory, produce_dir
from multiply.util.parsing import parse_parameters
from multiply.util.io import load_bed_as_dataframe
from multiply.download.collection import genome_collection
from multiply.generate.targets import Target, TargetSet
from multiply.generate.primer3 import Primer3Runner
from multiply.generate.primers import Primer, PrimerPair, load_primer_pairs_from_primer3_output


@click.command(short_help="Generate candidate primers.")
@click.option(
        "-d",
        "--design",
        type=str,
        required=True,
        help="Path to MULTIPLY design file (e.g. 'designs/pf-default.ini').",
    )
#@design_path_option
def generate(design):
    # PARSE CLI
    params = parse_parameters(design)
    genome = genome_collection[params["genome"]]

    # CREATE OUTPUT DIRECTORY
    params["output_dir"] = create_output_directory(params)

    # EXTRACT GENE INFORMATION
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
        

    # EXTRACT REGION INFORMATION
    regions = []
    if params["from_regions"]:
        region_df = load_bed_as_dataframe(params["region_bed"])
        # region_df = pd.read_csv(
        #     params["region_bed"], sep="\t", names=["seqname", "start", "end", "ID"]
        # )
        regions = [Target.from_series(row) for _, row in region_df.iterrows()]

    # MERGE
    targets = genes + regions
    target_set = (
        TargetSet(targets)
        .check_size_compatible(params["max_size_bp"])
        .calc_pads()
        .extract_seqs(genome.fasta_path, include_pads=True)
        .to_csv(f"{params['output_dir']}/table.targets_overview.csv")
        .to_fasta(f"{params['output_dir']}/targets_sequence.fasta")
    )

    # RUN PRIMER3
    primer3_output_dir = produce_dir(params["output_dir"], "primer3")
    primer3_runner = Primer3Runner()
    
    # Storage
    primer_pair_dt = {target.ID: [] for target in target_set.targets}

    # Iterate over settings
    for primer3_setting in params["primer3_settings"]:

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

            # UP TO HERE I THINK EVERYTHING IS GOOD

            # Store
            primer_pairs = load_primer_pairs_from_primer3_output(
                primer3_runner.output_path,
                add_target=target
            )

            primer_pair_dt[target.ID].extend(primer_pairs)

    # REDUCE TO UNIQUE PAIRS
    print("Primer pairs discovered by primer3:")
    print(f"{'Target':<15} {'Total':<10} {'Unique':<10}")
    for target_id, all_primer_pairs in primer_pair_dt.items():

        # Reduce to unique primer pairs, and give names
        # this definitely can be cleaned
        uniq_primer_pairs = set(all_primer_pairs)
        for ix, pair in enumerate(uniq_primer_pairs):
            pair.give_primers_names(
                primer_code=params["primer_code"],
                primer_ix=ix
            )

        # Print summary
        print(f"{target_id:<15} {len(all_primer_pairs):<10} {len(uniq_primer_pairs):<10}")

        # Store
        primer_pair_dt[target_id] = uniq_primer_pairs
    print("Done.")
    print("")

    # Add tails, if provided
    if params["include_tails"]:
        # TODO: incapsulate below
        for target_id, primer_pairs in primer_pair_dt.items():
            for primer_pair in primer_pairs:
                primer_pair.F.add_tail(params["F_tail"])
                primer_pair.R.add_tail(params["R_tail"])

    # WRITE
    primer_df = pd.DataFrame([
        pair.get_primer_as_dict(direction)
        for _, primer_pairs in primer_pair_dt.items()
        for pair in primer_pairs
        for direction in ["F", "R"]
    ])
    primer_df = primer_df[
        ["target_id", "target_name", "pair_name", "primer_name", "direction", 
        "seq", "length", "tm", "gc", "chrom", "start", "product_bp", "pair_penalty"]
    ]
    # if it has a proper start, it needs a chromosome
    # probably want to insert a candidate index column
    primer_df.to_csv(f"{params['output_dir']}/table.candidate_primers.csv", index=False)


