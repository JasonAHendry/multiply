from abc import ABC, abstractmethod
from dataclasses import dataclass


# ================================================================================
# Define a Genome
#
# ================================================================================


@dataclass(frozen=True)
class Genome:
    """Encapsulate all information for a given reference genome"""

    name: str
    source: str

    fasta_url: str
    fasta_raw_download: str
    fasta_path: str

    gff_url: str = ""
    gff_raw_download: str = ""
    gff_path: str = ""

    variation_dir: str = ""


# ================================================================================
# Create Genomes from different sources
#
# ================================================================================


class GenomeFactory(ABC):
    """Base class for creating genomes from different sources"""

    output_dir = "genomes/information"

    @abstractmethod
    def create_genome(self):
        pass


class PlasmoDBFactory(GenomeFactory):
    """
    Create, download, and postprocess genome information
    from plasmodb

    """

    source = "plasmodb"
    source_url = "https://plasmodb.org/common/downloads"
    release = 52

    def create_genome(self, name, genus, species, strain, include_variation=None):
        """Create a genome object from plasmodb"""

        # Process input information
        lineage = f"{genus.capitalize()[0]}{species.lower()}"

        # Get URL for species, strain
        data_url = f"{self.source_url}/release-{self.release}/{lineage}{strain}"

        # Prepare FASTA information
        fasta_fn = f"PlasmoDB-{self.release}_{lineage}{strain}_Genome.fasta"
        fasta_url = f"{data_url}/fasta/data/{fasta_fn}"
        fasta_raw_download = f"{self.output_dir}/{name}/{fasta_fn}"
        fasta_path = fasta_raw_download  # already decompressed

        # Prepare GFF information
        gff_fn = f"PlasmoDB-{self.release}_{lineage}{strain}.gff"
        gff_url = f"{data_url}/gff/data/{gff_fn}"
        gff_raw_download = f"{self.output_dir}/{name}/{gff_fn}"
        gff_path = gff_raw_download.replace(".gff", ".csv")

        # Create Genome
        genome = Genome(
            name=name,
            source=self.source,
            fasta_url=fasta_url,
            fasta_raw_download=fasta_raw_download,
            fasta_path=fasta_path,
            gff_url=gff_url,
            gff_raw_download=gff_raw_download,
            gff_path=gff_path,
            variation_dir=include_variation if include_variation is not None else "",
        )

        return genome


class EnsemblGenomesFactory(GenomeFactory):

    source = "ensemblgenomes"
    release = 52

    def create_genome(
        self, name, clade, genus, species, assembly, include_variation=None
    ):
        """Create a genome object from EnsemblGenomes"""

        # Define source URL
        source_url = "http://ftp.ensemblgenomes.org/pub/"
        source_url += f"{clade}/release-{self.release}"

        # Process input information
        lineage = f"{genus}_{species}".lower()

        # Prepare FASTA information
        fasta_fn = f"{lineage.capitalize()}.{assembly}.dna.toplevel.fa.gz"
        fasta_url = f"{source_url}/fasta/{lineage}/dna/{fasta_fn}"
        fasta_raw_download = f"{self.output_dir}/{name}/{fasta_fn}"
        fasta_path = fasta_raw_download.replace(".gz", "")

        # Prepare GFF information
        gff_fn = f"{lineage.capitalize()}.{assembly}.{self.release}.gff3.gz"
        gff_url = f"{source_url}/gff3/{lineage}/{gff_fn}"
        gff_raw_download = f"{self.output_dir}/{name}/{gff_fn}"
        gff_path = gff_raw_download.replace(".gff3.gz", ".csv")

        # Create Genome
        genome = Genome(
            name=name,
            source=self.source,
            fasta_url=fasta_url,
            fasta_raw_download=fasta_raw_download,
            fasta_path=fasta_path,
            gff_url=gff_url,
            gff_raw_download=gff_raw_download,
            gff_path=gff_path,
            variation_dir=include_variation if include_variation is not None else "",
        )

        return genome
