import re
import os
import subprocess
import pandas as pd


# CONSTANTS
BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen"
WORD_SIZE = 7


class BlastRunner:

    BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen"

    def __init__(self, input_fasta, reference_fasta):
        """
        Run BLAST


        """
        self.input_fasta = input_fasta
        self.reference_fasta = reference_fasta

    def create_database(self):
        """
        Check if a `blast` database has already been generated for
        the `reference_fasta`, if not, create one.

        """

        # Define database name
        if not self.reference_fasta.endswith(".fasta"):
            "`reference_fasta` must be a .fasta file, ending with `.fasta`."
        self.db_path = re.sub(".fasta$", "", self.reference_fasta)

        # Check if database already exists
        db_suffixes = [".nhr", ".nin", ".nsq"]
        if all([os.path.isfile(f"{self.db_path}{suff}") for suff in db_suffixes]):
            print(f"BLAST database '{self.db_path}' already exists.")
            return self

        # Create database
        cmd = "makeblastdb"
        cmd += f" -in  {self.reference_fasta}"
        cmd += " -dbtype nucl -parse_seqids"
        cmd += f" -out {self.db_path}"

        subprocess.run(cmd, check=True, shell=True)

        return self

    def run(self, output_archive):
        """
        Run blast, writing a BLAST archive to `output_archive`

        Note that we run with output format 11 `-outfmt 11` to
        produce the archive; from this format you can convert to
        all other formats.

        """

        # Define command
        cmd = "blastn"
        cmd += f" -db {self.db_path}"
        cmd += f" -query {self.input_fasta}"
        cmd += f" -word_size {WORD_SIZE}"
        cmd += f" -outfmt 11"
        cmd += f" -out {output_archive}"

        # Run
        subprocess.run(cmd, check=True, shell=True)

        # Save as instance variable, for reformattings
        self.output_archive = output_archive

        return self

    def reformat_output_as_table(self, output_table):
        """
        Reformat the output from `self.run()`

        """

        # Define command
        cmd = "blast_formatter"
        cmd += f" -archive {self.output_archive}"
        cmd += f" -outfmt '6 {self.BLAST_COLS}'"
        cmd += f" -out {output_table}"

        # Run
        subprocess.run(cmd, check=True, shell=True)

        # Save
        self.output_table = output_table

        return self

    def _convert_to_dataframe(self):
        """
        Convert from tabular output to a pandas dataframe, and save

        """
        # Load as a dataframe
        self.blast_df = pd.read_csv(
            self.output_table, sep="\t", names=self.BLAST_COLS.split(" ")
        )

    def get_dataframe(self):

        self._convert_to_dataframe()

        return self.blast_df
