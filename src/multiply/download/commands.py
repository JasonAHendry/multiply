import click
from multiply.download.collection import genome_collection
from multiply.download.gff import gff_standardisation_functions
from multiply.download.downloaders import GenomeDownloader


# ================================================================================
# Main function wrapped for Click CLI
#
# ================================================================================


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
    main(available, all, genome_name)


# ================================================================================
# Main function, unwrapped
#
# ================================================================================
    

def main(available, all, genome_name):
    if available:
        genome_collection.display()
        return

    # Prepare downloader
    downloader = GenomeDownloader()

    if all:
        for genome_name, genome in genome_collection.items():
            downloader.set_genome(genome)
            downloader.download_fasta()
            downloader.download_gff()
            downloader.standardise_gff(gff_standardisation_functions[genome.source])
            downloader.close_logging()
    elif genome_name is not None:
        genome = genome_collection[genome_name]
        downloader.set_genome(genome)
        downloader.download_fasta()
        downloader.download_gff()
        downloader.standardise_gff(gff_standardisation_functions[genome.source])
        downloader.close_logging()
    else:
        print("Must specify options. Type 'multiply download --help' for details.")

