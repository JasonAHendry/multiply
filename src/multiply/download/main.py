from multiply.download.collection import genome_collection
from multiply.download.gff import gff_standardisation_functions
from multiply.download.fasta import unmask_fasta_info
from multiply.download.downloaders import GenomeDownloader


def download(available, all, genome_name):
    """
    Download genome information for a given `genome_name`

    This includes genome sequences (.fasta) and information
    about gene locations (.gff). Diversity information
    must be prepared separately.

    """
    if available:
        genome_collection.display()
        return

    # Prepare downloader
    downloader = GenomeDownloader()

    if all:
        for genome_name, genome in genome_collection.items():
            downloader.set_genome(genome)
            downloader.download_fasta(unmask_fasta_info[genome.source])
            downloader.download_gff()
            downloader.standardise_gff(gff_standardisation_functions[genome.source])
            downloader.close_logging()
    elif genome_name is not None:
        genome = genome_collection[genome_name]
        downloader.set_genome(genome)
        downloader.download_fasta(unmask_fasta_info[genome.source])
        downloader.download_gff()
        downloader.standardise_gff(gff_standardisation_functions[genome.source])
        downloader.close_logging()
    else:
        print("Must specify options. Type 'multiply download --help' for details.")
