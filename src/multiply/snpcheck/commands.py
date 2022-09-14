import click
from multiply.download.collection import genome_collection


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
    from .main import snpcheck

    snpcheck(primer_csv, genome_name)
