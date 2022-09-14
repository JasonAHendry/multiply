import click
from multiply.download.collection import genome_collection


@click.command(short_help="Download genome information.")
@click.option("--available", is_flag=True, help="List genomes available for download.")
@click.option("--all", is_flag=True, help="Download all genomes available.")
@click.option(
    "-g",
    "--genome_name",
    type=click.Choice(genome_collection),
    default=None,
    help="Name of genome to download.",
)
def download(available, all, genome_name):
    """
    Download genome information for a given `genome_name`

    This includes genome sequences (.fasta) and information
    about gene locations (.gff). Diversity information
    must be prepared separately.

    """
    from .main import download

    download(available, all, genome_name)
