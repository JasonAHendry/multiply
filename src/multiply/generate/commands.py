import click
import pandas as pd
#from multiply.cli import design_path_option
from multiply.util.dirs import create_output_directory, produce_dir
from multiply.util.parsing import parse_parameters
from multiply.download.collection import genome_collection
from multiply.generate.targets import Target, TargetSet
from multiply.generate.primer3 import Primer3Runner


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
        gene_df = pd.read_csv(genome.standard_gff_path)
        target_ids = params["target_ids"]
        gene_df.query("ID in @target_ids", inplace=True)
        genes = [Target.from_series(row) for _, row in gene_df.iterrows()]

    # EXTRACT REGION INFORMATION
    regions = []
    if params["from_regions"]:
        region_df = pd.read_csv(
            params["region_bed"], sep="\t", names=["seqname", "start", "end", "ID"]
        )
        regions = [Target.from_series(row) for _, row in region_df.iterrows()]

    # MERGE
    targets = genes + regions
    target_set = (
        TargetSet(targets)
        .check_size_compatible(params["max_size_bp"])
        .calc_pads()
        .extract_seqs(genome.fasta_path, include_pads=True)
        .to_csv(f"{params['output_dir']}/target_genes.summary.csv")
    )

    # RUN PRIMER3
    primer3_output_dir = produce_dir(params["output_dir"], "primer3")
    primer3_runner = Primer3Runner()
    primer3_runner.load_primer3_settings("relaxed")
    primer3_runner.set_amplicon_size_ranges(
        min_size_bp=params["min_size_bp"], max_size_bp=params["max_size_bp"]
    )

    for target in target_set.targets:
        primer3_runner.set_target(
            ID=target.ID,
            seq=target.seq,
            start=target.start,
            pad_start=target.pad_start,
            length=target.length,
        )

        primer3_runner.run(output_dir=primer3_output_dir)
