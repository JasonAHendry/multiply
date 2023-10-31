import os
import urllib.request
import pandas as pd
from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import lru_cache
from multiply.util.definitions import ROOT_DIR
from multiply.util.dirs import produce_dir
from multiply.util.exceptions import GenomeCollectionError


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

    include_variation: str = ""


# ================================================================================
# Create Genomes from different sources
#
# ================================================================================


class GenomeFactory(ABC):
    """
    Base class for creating genomes from different sources

    """

    output_dir = f"{ROOT_DIR}/genomes/information"

    @abstractmethod
    def create_genome(self):
        pass


class PlasmoDBFactory(GenomeFactory):
    """
    Create Genome objects from PlasmoDB

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
            include_variation=include_variation
            if include_variation is not None
            else "",
        )

        return genome


class EnsemblGenomesFactory(GenomeFactory):
    """
    Create Genome objects from EnsemblGenomes

    """

    source = "ensemblgenomes"
    release = 56

    clades = ["plants", "metazoa", "protists", "fungi", "bacteria"]

    def create_genome(
        self, name, clade, genus, species, assembly, include_variation=None
    ):
        """Create a genome object from EnsemblGenomes"""

        if not clade in self.clades:
            raise GenomeCollectionError(
                f"Provided clade '{clade}' not in clade list for {self.source} downloads:\n{', '.join(self.clades)}."
            )

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
            include_variation=include_variation
            if include_variation is not None
            else "",
        )

        return genome


class RefSeqGenomesFactory(GenomeFactory):
    """
    Create Genome objects from RefSeq Genome database

    https://www.ncbi.nlm.nih.gov/genome/

    """

    source = "refseq"
    refseq_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"

    clades = [
        "archea",
        "bacteria",
        "fungi",
        "invertebrate",
        "metagenomes",
        "mitochondrion",
        "plant",
        "plasmid",
        "protozoa",
        "unknown",
        "vertebrate_mammalian",
        "vertebrate_other",
        "viral",
    ]  # could be an Enum

    @staticmethod
    def exists_locally(file_path):
        return os.path.isfile(file_path)

    def _download_assembly_summary(
        self, species_url: str, genome_dir: str, check_already_exists: bool = True
    ) -> str:
        """
        Given a URL to a RefSeq species directory, download the `assembly_summary.txt`
        file

        Return the local path to the downloaded file

        """

        # Define local path to assembly, check if it exists
        assembly_txt = f"{genome_dir}/assembly_summary.txt"
        if check_already_exists and self.exists_locally(assembly_txt):
            return assembly_txt

        # Download
        assembly_url = f"{species_url}/assembly_summary.txt"
        try:
            urllib.request.urlretrieve(url=assembly_url, filename=assembly_txt)
        except urllib.error.HTTPError:
            raise GenomeCollectionError(
                f"Cannot find assembly information for {os.path.basename(genome_dir)} on RefSeq."
            )

        return assembly_txt

    @staticmethod
    @lru_cache(maxsize=None)
    def _extract_assembly_url(assembly_txt: str, assembly_name: str) -> str:
        """
        Extract the URL of a given RefSeq assembly from the `assembly_summary.txt`
        file

        TODO: This slows down the API --> Could I cache?

        """

        # Read the assembly text file
        assembly_df = pd.read_csv(assembly_txt, skiprows=1, sep="\t")
        assembly_df.rename(
            {"#assembly_accession": "assembly_accession"}, axis=1, inplace=True
        )

        # Query for target accession
        assert "assembly_accession" in assembly_df.columns
        assert "ftp_path" in assembly_df.columns
        assembly_df.query("assembly_accession == @assembly_name", inplace=True)

        if assembly_df.shape[0] != 1:
            raise GenomeCollectionError(
                f"Tried to find single unique URL for assembly {assembly_name}, "
                + f"but found {assembly_df.shape[0]} matches. Please confirm this assembly exists in the {assembly_txt} file."
            )

        return assembly_df.squeeze()["ftp_path"]

    def create_genome(
        self, name, clade, genus, species, assembly, include_variation=None
    ):
        """Create a genome object from RefSeq Genomes"""

        if not clade in self.clades:
            raise GenomeCollectionError(
                f"Provided clade '{clade}' not in clade list for {self.source} downloads:\n{', '.join(self.clades)}."
            )

        # Define source URL
        genome_dir = produce_dir(
            self.output_dir,
            f"{genus.lower().capitalize()}{species.lower().capitalize()}",
        )

        species_url = (
            f"{self.refseq_url}/{clade}/{genus.lower().capitalize()}_{species.lower()}"
        )
        assembly_txt = self._download_assembly_summary(species_url, genome_dir)
        source_url = self._extract_assembly_url(assembly_txt, assembly)
        assembly_full = os.path.basename(
            source_url
        )  # NB: Required because RefSeq modifies assembly name in some cases!

        # Prepare FASTA information
        fasta_fn = f"{assembly_full}_genomic.fna.gz"
        fasta_url = f"{source_url}/{fasta_fn}"
        fasta_raw_download = f"{self.output_dir}/{name}/{fasta_fn}"
        fasta_path = fasta_raw_download.replace(".gz", "")

        # Prepare GFF information
        gff_fn = f"{assembly_full}_genomic.gff.gz"
        gff_url = f"{source_url}/{gff_fn}"
        gff_raw_download = f"{self.output_dir}/{name}/{gff_fn}"
        gff_path = gff_raw_download.replace(".gff.gz", ".csv")

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
            include_variation=include_variation
            if include_variation is not None
            else "",
        )

        return genome
