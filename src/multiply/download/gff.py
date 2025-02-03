import gzip
import pandas as pd
from dataclasses import dataclass


# ================================================================================
# Loading and adding columns to Genome Feature Format files
#
# ================================================================================


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
    
    # Coerce data types
    gff_df = pd.DataFrame(entries)
    gff_df["start"] = gff_df["start"].astype("int")
    gff_df["end"] = gff_df["end"].astype("int")

    return gff_df


def add_gff_attributes(input_df, field_names=["Parent", "ID", "Name"]):
    """ Add specific attributes as columns to .gff """
    
    # Dictionary to store new columns
    dt = {field_name: [] for field_name in field_names}
    
    # Iterate over attributes
    for attributes in input_df["attribute"]:
        
        # Extract fields
        fields = attributes.split(";")
        
        # Iterate over names to add
        for field_name in field_names:
            
            # Extract value
            value = [f.split("=")[1] for f in fields if f.startswith(f"{field_name}=")]
            
            # Add to dictionary, if exists
            dt[field_name].append(value[0] if value else None)
    
    # Add to data frame
    df = pd.DataFrame(dt)
    
    # Need to reset input index
    input_df.reset_index(inplace=True, drop=True)
            
    return pd.concat([input_df, df], axis=1)


# ================================================================================
# Functions to specifically 'standardise' GFF files from different sources
# 
# Unfortunately there are differences in feature keywords and interpretations
# that need to be remedied for seemless downstream usage
#
# ================================================================================


def standardise_PlasmoDB_gff(gff_df):
    """
    Standardise GFF dataframe download from PlasmoDB
    
    For later releases of PlasmoDB this is really quite cumbersome,
    as for genes with multiple exons, there is no single record that 
    indicates the start and end of the ORF / CDS.
    
    It's slightly annoying that we lose the name / attribute columns
    
    """
    
    standard_df = gff_df.query("feature in 'CDS'")
    standard_df = add_gff_attributes(
        input_df=standard_df, 
        field_names=["ID", "Parent", "Name"]
    )
    standard_df.rename({"Name": "name"}, axis=1, inplace=True)

    # Define fields
    @dataclass
    class StandardisedGffEntry:
        seqname: str
        source: str
        feature: str
        start: int
        end: int
        score: str
        strand: str
        frame: str
        attribute: str
        ID: str
        name: str
    
    # Iterate over CDS parents
    records = []
    for _, gdf in standard_df.groupby("Parent"):

        # Extract start + stop across all CDS
        start = gdf["start"].min()
        end = gdf["end"].max()

        # Populate (largely from first row)
        kwarg_columns = ["seqname", "source", "feature", "score", "strand", "attribute", "name"]
        first_row = gdf.iloc[0]
        record = StandardisedGffEntry(
            start=start, 
            end=end,
            frame=None,
            ID=first_row["Parent"].split(".")[0],
            **first_row[kwarg_columns].to_dict()
        )

        # Store
        records.append(record)
        
    return pd.DataFrame(records)


def standardise_EnsemblGenomes_gff(gff_df, restrict_to=["gene"]):
    """
    Standardise GFF dataframe download from EnsemblGenomes
    
    """

    # Query for relevant features
    standard_df = gff_df.query("feature in @restrict_to")
    standard_df = add_gff_attributes(
        input_df=standard_df,
        field_names=["gene_id", "Name"]
    )
    standard_df.rename({"gene_id":"ID", "Name":"name"}, axis=1, inplace=True)

    return standard_df


def standardise_RefSeqGenomes_gff(gff_df, restrict_to=["gene"], source_only=["RefSeq", "BestRefSeq"]):
    """
    Standardise GFF dataframe download from RefSeq Genome database
    
    """
    
    # Query for relevant features
    standard_df = gff_df.query("feature in @restrict_to and source in @source_only")
    standard_df = add_gff_attributes(
        input_df=standard_df,
        field_names=["Name", "ID"]
    )
    standard_df.rename({"Name": "name"}, axis=1, inplace=True)
    standard_df["ID"] = [s[5:] for s in standard_df["ID"]]  # drop prefix 'gene-'
    
    return standard_df

def standardise_VectorBaseGenomes_gff(gff_df, restrict_to=["protein_coding_gene"], source_only=["VEuPathDB"]):
    """
    Standardise GFF dataframe download from RefSeq Genome database
    
    """
    
    # Query for relevant features
    standard_df = gff_df.query("feature in @restrict_to and source in @source_only")
    standard_df = add_gff_attributes(
        input_df=standard_df,
        field_names=["Name", "ID"]
    )
    standard_df.rename({"Name": "name"}, axis=1, inplace=True)
        
    return standard_df

# Prepare .gff standardisation
gff_standardisation_functions = {
    "plasmodb": standardise_PlasmoDB_gff,
    "vectorbase": standardise_VectorBaseGenomes_gff,
    "ensemblgenomes": standardise_EnsemblGenomes_gff,
    "refseq": standardise_RefSeqGenomes_gff,
}
