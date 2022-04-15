import click
from multiply.download.collection import genome_collection
from multiply.download.gff import gff_standardisation_functions
from multiply.download.downloaders import GenomeDownloader


# @click.option("--all")
@click.command(short_help="Download genome information.")
@click.option("--available", is_flag=True, help="List genomes available for download.")
@click.option(
    "-g",
    "--genome_name",
    default=None,
    help="Genome to download; must exist in collection.",
)
def download(available, genome_name):
    """
    Utilities to download genome information for MULTIPLY

    This includes genome sequences (.fasta) and information
    about gene locations (.gff). Diversity information
    must be prepared separately.

    """
    # Print available
    if available:
        # TODO:
        # - This probably goes into a function on collection.py
        # logger.info(f"Found {len(genome_collection)} genomes in collection.")
        # logger.info("")
        print(f"Found {len(genome_collection)} genomes in collection.")
        print("")
        print(f"{'Name':25} {'Source':15}")
        for _, genome in genome_collection.items():
            print(f"{genome.name:25} {genome.source:15}")

    # Prepare downloader
    downloader = GenomeDownloader()

    # Download
    if genome_name is not None:
        genome = genome_collection[genome_name]
        downloader.set_genome(genome)
        downloader.download_fasta()
        downloader.download_gff()
        downloader.standardise_gff(gff_standardisation_functions[genome.source])
        downloader.close_logging()
