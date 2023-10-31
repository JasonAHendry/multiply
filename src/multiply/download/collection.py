import os
import configparser
from multiply.util.exceptions import GenomeCollectionError
from multiply.util.definitions import ROOT_DIR
from multiply.download.genomes import PlasmoDBFactory, EnsemblGenomesFactory, RefSeqGenomesFactory




INI_PATH = f"{ROOT_DIR}/genomes/collection.ini"
FACTORIES = [PlasmoDBFactory, EnsemblGenomesFactory, RefSeqGenomesFactory]

# ================================================================================
# Define a collection of Genome objects
#
# ================================================================================


class GenomeCollection(dict):
    def __init__(self, ini_path, factories):
        """
        Hold a collection of Genome objects as a dictionary,
        defined from an initiation file located at `ini_path`

        params
            ini_path: str
                Path to .ini file containing information
                about defined genomes.
            factories: list[GenomeFactory]
                List of GenomeFactory instances, used to
                create individual genomes based on
                their source.
        
        """

        # Ensure points to valid file, then set
        if not os.path.isfile(ini_path):
            raise FileNotFoundError(f"No genome collection found at {ini_path}. Check file path is correct.")
        self._ini_path = ini_path

        # Read config object
        self._config = configparser.ConfigParser()
        self._config.read(ini_path)

        # Set the factories used to build Genome objects
        self._factories = factories
        self._factory_names = [f.source for f in self._factories]

    def populate(self):
        """
        Populate the genome collection

        """

        for section in self._config.sections():

            # Extract parameters
            kwargs = dict(self._config.items(section))
            source = kwargs.pop("source")

            # Get suitable factory
            if not source in self._factory_names:
                raise GenomeCollectionError(f"No support to download genomes from source '{source}'. "
                                           f" Supported sources: {', '.join(self._factory_names)}."
                                           )

            (factory,) = [factory() for factory in self._factories if factory.source == source]

            # Create genome
            genome = factory.create_genome(name=section, **kwargs)

            # Store
            self[genome.name] = genome

    def is_downloaded(self, genome_name):
        """
        Check if a given genome is downloaded
        
        """

        has_fasta = os.path.isfile(self[genome_name].fasta_path)
        has_gff = os.path.isfile(self[genome_name].gff_path)

        return has_fasta and has_gff

    def display(self):
        """
        Display the genome collection
        
        """

        print("SUMMARY")
        print(f"  Collection located at: {self._ini_path}")
        print(f"  No. genomes in collection: {len(self)}")
        print("")
        print("BREAKDOWN")
        print(f"  {'Name':25} {'Source':25} {'Downloaded [both FASTA and GFF]':10}")
        for genome_name, genome in self.items():
            print(f"  {genome.name:25} {genome.source:25} {self.is_downloaded(genome_name)!s:10}")


# ================================================================================
# Initialise the genome collection
# - Not sure I like this here, as it gets run on every import
# ================================================================================

genome_collection = GenomeCollection(INI_PATH, FACTORIES)
genome_collection.populate()