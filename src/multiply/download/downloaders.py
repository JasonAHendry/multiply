import os
import logging
import urllib.request
from multiply.download.gff import load_gff


# Prepare logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler()
logger.addHandler(stream_handler)


class GenomeDownloader:
    def __init__(self):
        self.genome = None
        self.log_handler = None

    def set_genome(self, genome):
        """Set genome to download"""

        # Assign genome
        self.genome = genome

        # Prepare logging
        log_path = f"genomes/information/{genome.name}/{genome.name}.log"
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

    def download_fasta(self):
        """Download .fasta information"""

        logger.info("Downloading .fasta information.")
        logger.info(f"  Source URL: {self.genome.fasta_url}")
        logger.info(f"  Destination path: {self.genome.fasta_path}")

        # Skip if already downloaded
        if self.exists_locally(self.genome.fasta_path):
            logger.info("  Already downloaded.")
            logger.info("  Skipping.")
            logger.info("")

            return

        # Otherwise, download from the URL
        logger.info("  Downloading...")

        self.produce_dir(self.genome.fasta_path)
        urllib.request.urlretrieve(
            url=self.genome.fasta_url, filename=self.genome.fasta_path
        )

        logger.info("  Done.")
        logger.info("")

    def download_gff(self):
        """Download .gff information"""

        logger.info("Downloading .gff information.")
        logger.info(f"  Source URL: {self.genome.gff_url}")
        logger.info(f"  Destination path: {self.genome.gff_path}")

        # Skip if already downloaded
        if self.exists_locally(self.genome.gff_path):
            logger.info("  Already downloaded.")
            logger.info("  Skipping.")
            logger.info("")

            return

        # Otherwise, download from the URL
        logger.info("  Downloading...")

        self.produce_dir(self.genome.gff_path)
        urllib.request.urlretrieve(
            url=self.genome.gff_url, filename=self.genome.gff_path
        )

        logger.info("  Done.")
        logger.info("")

    def standardise_gff(self, standardise_fn):
        """Process .gff information into a standardised format"""

        logger.info("Standardising .gff information.")
        logger.info(f"  Raw .gff: {self.genome.gff_path}")
        logger.info(f"  Standardised .gff: {self.genome.standard_gff_path}")

        # Skip if already downloaded
        if self.exists_locally(self.genome.standard_gff_path):
            logger.info("  Already standardised.")
            logger.info("  Skipping.")
            logger.info("")

            return

        # Otherwise standardise
        logger.info(f"  Loading downloaded .gff...")
        gff = load_gff(self.genome.gff_path)
        logger.info(f"  Found {gff.shape[0]} entries in .gff.")
        logger.info(f"  Standardising...")
        standard_gff = standardise_fn(gff)
        logger.info(f"  {standard_gff.shape[0]} entries remain.")
        logger.info(f"  Example IDs: {', '.join(standard_gff.sample(6)['ID'])}")
        logger.info(f"  Writing...")
        standard_gff.to_csv(self.genome.standard_gff_path, index=False)
        logger.info("  Done.")
        logger.info("")

    def close_logging(self):

        logger.info("Log complete.")
        logger.info("")
        logger.removeHandler(self.log_handler)
