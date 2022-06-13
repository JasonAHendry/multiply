import re
import subprocess


# Constants -- move to .json
BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen"
WORD_SIZE = 7

def make_blast_db(source_fasta):
    """
    Make a `blast` database from a `source_fasta` file
    
    """
    
    assert source_fasta.endswith(".fasta")
    output_db = re.sub(".fasta$", "", source_fasta)
    
    cmd = "makeblastdb"
    cmd += f" -in  {source_fasta}" 
    cmd += " -dbtype nucl -parse_seqids"
    cmd += f" -out {output_db}"
    
    subprocess.run(cmd, check=True, shell=True)
    
    return output_db

def run_blast(db_path, input_fasta, output_blast, output_fmt=None):
    """
    Run `blast` on an `input_fasta` and blast database `db_path`
    
    - Want to allow for different output formats
        `blast_formatter` should be used for this
    - Is it possible to run blast only one time, and produce multiply outputs?
    
    """

    cmd = "blastn"
    cmd += f" -db {db_path}"
    cmd += f" -query {input_fasta}"
    cmd += f" -word_size {WORD_SIZE}"
    cmd += f" -out {output_blast}"
    if output_fmt == 6:
        cmd += f" -outfmt '6 {BLAST_COLS}'"
    
    subprocess.run(cmd, check=True, shell=True)