import subprocess


def convert_fasta_to_all_uppercase(fasta_path: str, dry_run=False) -> None:
    """ 
    Convert sequences within a FASTA file to all upper case
    
    Lowercase is an indication of soft-masking in some cases. For example,
    see FAQ 23 here:
    
    https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/genome/doc/ftpfaq/#alignmask
    
    Where they advise the following:
    
    awk '{if(/^[^>]/)$0=toupper($0);print $0}' genomic.fna > genomic.unmasked.fna
    
    """
    
    fasta_unmasked_path = f"{fasta_path}.unmasked"    
    
    cmd_unmask = "awk '{if(/^[^>]/)$0=toupper($0);print $0}'"
    cmd_unmask += f" {fasta_path} > {fasta_unmasked_path}"
    cmd_mv = f"mv -f {fasta_unmasked_path} {fasta_path}"
    cmd = " && ".join([cmd_unmask, cmd_mv])
    
    if dry_run:
        print(cmd)
        return

    subprocess.run(cmd, shell=True, check=True)


unmask_fasta_info = {
    "plasmodb": False,
    "ensemblgenomes": False,
    "refseq": True,
    "vectorbase": False
}

