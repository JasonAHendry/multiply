import configparser
from multiply.download.genomes import PlasmoDBFactory, EnsemblGenomesFactory

INI_PATH = "genomes/collection.ini"  # need to check this exists
factories = [PlasmoDBFactory, EnsemblGenomesFactory]


def gen_genome_collection(ini_path, factories):
    """
    Generate the a genome collection, based on
    the current genome collectionn `ini_path` and
    genome `factories`

    """
    # Load genome collection from source
    config = configparser.ConfigParser()
    config.read(ini_path)

    # Create collection
    genome_collection = {}
    for section in config.sections():

        # Extract parameters
        kwargs = dict(config.items(section))
        source = kwargs.pop("source")

        # Get suitable factory
        (factory,) = [factory() for factory in factories if factory.source == source]

        # Create genome
        genome = factory.create_genome(name=section, **kwargs)

        # Store
        genome_collection[genome.name] = genome

    return genome_collection


genome_collection = gen_genome_collection(ini_path=INI_PATH, factories=factories)
