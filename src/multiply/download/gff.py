import gzip
import pandas as pd
from dataclasses import dataclass


def load_gff(gff_path):
    """Load a gene feature format (.gff) file into a pandas DataFrame"""

    # Define fields
    @dataclass
    class gffEntry:
        seqname: str
        source: str
        feature: str
        start: int
        end: int
        score: str
        strand: str
        frame: str
        attribute: str

    # Prepare to load, handle compression
    if gff_path.endswith(".gz"):
        binary_gff = True
        open_gff = gzip.open
    else:
        binary_gff = False
        open_gff = open

    # Open gff
    entries = []
    with open_gff(gff_path) as gff:

        # Iterate over rows
        for line in gff:
            # Decode as necessary
            if binary_gff:
                line = line.decode()

            # Skip if info
            if line.startswith("#"):
                continue

            # Extract gff fields
            fields = line.strip().split("\t")
            entry = gffEntry(*fields)

            # Store
            entries.append(entry)

    return pd.DataFrame(entries)


def standardise_PlasmoDB_gff(gff_df, restrict_to=["protein_coding_gene"]):
    """Processing entails adding an ID and Name column"""

    # Query for relevant features
    standard_df = gff_df.query("feature in @restrict_to")

    # Extract ID and Name fields
    IDs = []
    Names = []
    for attributes in standard_df["attribute"]:

        # Split attributes
        fields = attributes.split(";")

        # Extract fields
        ID = [f.split("=")[1] for f in fields if f.startswith("ID=")]
        Name = [f.split("=")[1] for f in fields if f.startswith("Name=")]

        # Store
        IDs.append(ID[0])
        Names.append(Name[0] if Name else None)

    # Add to DataFrame
    standard_df.insert(9, "ID", IDs)
    standard_df.insert(10, "name", Names)

    return standard_df


def standardise_EnsemblGenomes_gff(gff_df, restrict_to=["gene"]):
    """Processing entails adding an ID and Name column"""

    # Query for relevant features
    standard_df = gff_df.query("feature in @restrict_to")

    # Extract ID and Name fields
    IDs = []
    Names = []
    for attributes in standard_df["attribute"]:

        # Split attributes
        fields = attributes.split(";")

        # Extract fields
        ID = [f.split("=")[1] for f in fields if f.startswith("ID=")]
        Name = [f.split("=")[1] for f in fields if f.startswith("Name=")]

        # Store
        IDs.append(ID[0].split(":")[1])
        Names.append(Name[0] if Name else None)

    # Add to DataFrame
    standard_df.insert(9, "ID", IDs)
    standard_df.insert(10, "name", Names)

    return standard_df


# Prepare .gff standardisation
gff_standardisation_functions = {
    "plasmodb": standardise_PlasmoDB_gff,
    "ensemblgenomes": standardise_EnsemblGenomes_gff,
}
