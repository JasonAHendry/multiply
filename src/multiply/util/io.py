import subprocess
import pandas as pd
from dataclasses import dataclass
from multiply.util.exceptions import BEDFormattingError


# ================================================================================
# Loading from ad writing to BED files
#
# ================================================================================


def load_bed_as_dataframe(bed_path):
    """Load a .bed file into a dataframe"""

    @dataclass
    class BedRecord:
        seqname: str
        start: int
        end: int
        ID: str
        name: str = ""

    records = []

    with open(bed_path, "r") as bed:
        for l in bed:
            if l.startswith("#"):
                continue
            fields = l.strip().split("\t")
            
            if not len(fields) == 4:
                raise BEDFormattingError(f"In bed file: '{bed_path}' ...\n" 
                                         f"...incorrect formatting of this line: '{repr(l)}'.")

            record = BedRecord(
                seqname=fields[0],
                start=int(fields[1]),
                end=int(fields[2]),
                ID=fields[3],
            )
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


def write_primers_to_bed(primer_df, output_path):
    """
    Write a primers data frame to a `.bed` file

    params
        primer_df: pandas DataFrame, shape(n_columns, n_primers)
            A pandas dataframe, typically created from the
            `multiply generate` command; where each row
            contains information about a specific primer.
        output_path: str
            Path to output `.bed` file.

    """

    with open(output_path, "w") as bed:
        for _, row in primer_df.iterrows():

            start = row["start"]
            if row["direction"] == "F":
                start = row["start"]
                end = start + row["length"]
            else:
                end = row["start"]
                start = end - row["length"]

            bed.write(f"{row['chrom']}\t{start}\t{end}\t{row['primer_name']}\n")


def write_amplicons_to_bed(primer_df, output_path):
    """
    Use primers data to write a bed file of amplicons
    
    params
        primer_df: pandas DataFrame, shape(n_columns, n_primers)
            A pandas dataframe, typically created from the
            `multiply generate` command; where each row
            contains information about a specific primer.
        output_path: str
            Path to output `.bed` file.
    
    """

    # Create BED file
    with open(output_path, "w") as bed:
        for target, target_df in primer_df.groupby("target_id"):
            
            # Get primer information for this target
            F_info = target_df.query("direction == 'F'").squeeze()
            R_info = target_df.query("direction == 'R'").squeeze()
            
            # Create a BED record
            record = f"{F_info['chrom']}"
            record += f"\t{F_info['start']}\t{R_info['start']}"
            record += f"\t{target}\n"
            
            # Write
            bed.write(record)


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


# ================================================================================
# Handling FASTA files
#
# ================================================================================


def load_fasta_as_dict(fasta_path):
    """
    Load a `.fasta` file as a dictionary

    """

    dt = {}
    with open(fasta_path, "r") as fasta:

        for line in fasta:

            # Extract header and sequence
            if line.startswith(">"):
                header = line[1:].rstrip()
            seq = fasta.readline().rstrip()

            # Ensure unique
            if header in dt:
                raise ValueError(
                    f"Header {header} found more than once in {fasta_path}."
                )

            # Add
            dt[header] = seq

    return dt


def write_fasta_from_dict(input_dt, output_fasta):
    """
    Write a `.fasta` file to `output_fasta` from an input dictionary
    `input_dt`

    """
    with open(output_fasta, "w") as fasta:
        for header, seq in input_dt.items():
            fasta.write(f">{header}\n")
            fasta.write(f"{seq}\n")
