import os
import gzip
import shutil
import logging
import urllib.request
from multiply.util.definitions import ROOT_DIR
from multiply.download.gff import load_gff
from multiply.download.fasta import convert_fasta_to_all_uppercase


# Prepare logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)


class GenomeDownloader:
    def __init__(self):
        """
        Download the FASTA and GFF files associated
        with a given Genome instance
        
        """
        
        self.genome = None
        self.log_handler = None

    def set_genome(self, genome):
        """
        Select genome to download
        
        """

        # Assign genome
        self.genome = genome

        # Prepare logging
        log_path = f"{ROOT_DIR}/genomes/information/{genome.name}/{genome.name}.log"
        self.produce_dir(log_path)

        # File Handler
        self.log_handler = logging.FileHandler(log_path)
        formatter = logging.Formatter("[%(asctime)s] %(message)s")
        self.log_handler.setFormatter(formatter)

        logger.addHandler(self.log_handler)

        logger.info(f"Name: {self.genome.name}")
        logger.info(f"Source: {self.genome.source}")
        logger.info("")

    @staticmethod
    def exists_locally(file_path):
        return os.path.isfile(file_path)

    @staticmethod
    def produce_dir(file_path):
        file_dir = os.path.dirname(file_path)
        if not os.path.isdir(file_dir):
            os.makedirs(file_dir)

    @staticmethod
    def decompress_file(input_file_path, decompressed_file_path):
        if not input_file_path.endswith(".gz"):
            return
    
        with gzip.open(input_file_path, "rb") as fi:
            with open(decompressed_file_path, "wb") as fo:
                shutil.copyfileobj(fi, fo)


    def _download(self, file_url, file_raw_path, decompress=False, file_decompressed_path=None):
        """
        Download a source file from `file_url` to a local path `file_raw_path`, and then
        optionally `decompress` to `file_decompressed_path`

        """

        logger.info(f"  Source URL: {file_url}")
        logger.info(f"  Destination path: {file_raw_path}")

        # Otherwise, download from the URL
        logger.info("  Downloading...")

        self.produce_dir(file_raw_path)
        urllib.request.urlretrieve(
            url=file_url, filename=file_raw_path
        )

        if decompress and file_decompressed_path is not None:
            logger.info("  Decompressing...")
            self.decompress_file(
                input_file_path=file_raw_path,
                decompressed_file_path=file_decompressed_path
            )

        logger.info("  Done.")
        logger.info("")
        
    def download_fasta(self, unmask=False):
        """
        Download .fasta information associated with Genome

        """

        logger.info("Downloading .fasta information.")
        if self.exists_locally(f"{ROOT_DIR}/{self.genome.fasta_path}"):
            logger.info("  Already downloaded.")
            logger.info("  Skipping.")
            logger.info("")
            return
        
        self._download(
            file_url=self.genome.fasta_url,
            file_raw_path=f"{ROOT_DIR}/{self.genome.fasta_raw_download}",
            decompress=True,
            file_decompressed_path=f"{ROOT_DIR}/{self.genome.fasta_path}"
        )
        
        if unmask:
            logger.info("This .fasta file requires soft-clipping to be unmasked.")
            logger.info("  Unmasking...")
            convert_fasta_to_all_uppercase(f"{ROOT_DIR}/{self.genome.fasta_path}")
            logger.info("  Done.")
            logger.info("")
        
    def download_gff(self):
        """
        Download .fasta information associated with Genome

        """

        logger.info("Downloading .gff information.")
        if self.exists_locally(f"{ROOT_DIR}/{self.genome.gff_raw_download}"):
            logger.info("  Already downloaded.")
            logger.info("  Skipping.")
            logger.info("")
            return
        
        self._download(
            file_url=self.genome.gff_url,
            file_raw_path=f"{ROOT_DIR}/{self.genome.gff_raw_download}",
            decompress=False
        )

    def standardise_gff(self, standardise_fn):
        """
        Process .gff information into a standardised format
        using a passed function `standardise_fn`
        
        """

        logger.info("Standardising .gff information.")
        logger.info(f"  Raw .gff: {ROOT_DIR}/{self.genome.gff_raw_download}")
        logger.info(f"  Standardised .gff: {ROOT_DIR}/{self.genome.gff_path}")

        # Skip if already downloaded
        if self.exists_locally(f"{ROOT_DIR}/{self.genome.gff_path}"):
            logger.info("  Already standardised.")
            logger.info("  Skipping.")
            logger.info("")
            return

        # Otherwise standardise
        logger.info(f"  Loading downloaded .gff...")
        gff = load_gff(f"{ROOT_DIR}/{self.genome.gff_raw_download}")
        logger.info(f"  Found {gff.shape[0]} entries in .gff.")
        logger.info(f"  Standardising...")
        standard_gff = standardise_fn(gff)
        logger.info(f"  {standard_gff.shape[0]} entries remain.")
        logger.info(f"  Example IDs: {', '.join(standard_gff.sample(6)['ID'])}")
        logger.info(f"  Writing...")
        standard_gff.to_csv(f"{ROOT_DIR}/{self.genome.gff_path}", index=False)
        logger.info("  Done.")
        logger.info("")

    def close_logging(self):
        """
        Close the logger specific for the genome download
        
        """
        
        logger.info("Log complete.")
        logger.info("")
        logger.removeHandler(self.log_handler)
