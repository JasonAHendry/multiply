import subprocess
import pandas as pd
from dataclasses import dataclass


def write_fasta_from_bed(bed_path, reference_fasta_path, output_path, verbose=False):
    """
    Write a .fasta file from a .bed file, using bedtools

    params
        bed_path : str
            Path to a .bed file.
        reference_fasta_path : str
            Path to reference genome .fasta.
        output_path : str
            Path to write .fasta file.
        verbose : bool
            Print to stdout?

    returns
        None

    """

    cmd = "bedtools getfasta -name -fi %s -bed %s" % (reference_fasta_path, bed_path)
    cmd += " > %s" % output_path

    if verbose:
        print("  Running bedtools...")
        print("  %s" % cmd)

    subprocess.run(cmd, shell=True, check=True)

    return None


def load_bed_as_dataframe(bed_path):
    """ Load a .bed file into a dataframe """
    
    @dataclass
    class BedRecord:
        seqname: str
        start: int
        end: int
        ID: str
        name: str=""
    records = []

    with open(bed_path, "r") as bed:
        for l in bed:
            if l.startswith("#"):
                continue
            fields = l.strip().split("\t")
            record = BedRecord(seqname=fields[0], start=int(fields[1]), end=int(fields[2]), ID=fields[3])
            if len(fields) == 5:
                record.name = fields[4]
            records.append(record)
    
    return pd.DataFrame(records)


def targets_to_bed(targets, bed_path, include_pads=True):
    """Write a set of targets to a bed file"""

    # Open bed file at `bed_path`
    with open(bed_path, "w") as bed:

        # Write header
        bed.write("# MUTLIPLY: Targets .bed file\n")
        if include_pads:
            bed.write("# Note that pads have been included.\n")

        # Iterate over targets
        for target in targets:

            # Extract start and end, depending on pad inclusion
            s = target.pad_start if include_pads else target.start
            e = target.pad_end if include_pads else target.end

            # Write to bed
            text = f"{target.chrom}\t{s}\t{e}\t{target.ID}\t{target.name}\n"
            bed.write(text)
