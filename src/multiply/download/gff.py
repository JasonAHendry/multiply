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
    df.columns = [c.lower() for c in df.columns]
    
    # Need to reset input index
    input_df.reset_index(inplace=True, drop=True)
            
    return pd.concat([input_df, df], axis=1)


def standardise_PlasmoDB_gff_old(gff_df, restrict_to=["protein_coding_gene"]):
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


def standardise_PlasmoDB_gff(gff_df):
    """
    Processing entails adding an ID and Name column
    
    For later releases of PlasmoDB wthis is really quite cumbersome,
    as for genes with multiple exons, there is no single record that 
    indicates the start and end of the ORF / CDS.
    
    It's slightly annoying that we lose the name / attribute columns
    
    """
    
    standard_df = gff_df.query("feature in 'CDS'")
    standard_df = add_gff_attributes(
        input_df=standard_df, 
        field_names=["ID", "Parent", "Name"]
    )

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
    for g, gdf in standard_df.groupby("parent"):

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
            ID=first_row["parent"].split(".")[0],
            #name=None,
            **first_row[kwarg_columns].to_dict()
        )

        records.append(record)
        
    return pd.DataFrame(records)


def standardise_EnsemblGenomes_gff(gff_df, restrict_to=["gene"]):
    """Processing entails adding an ID and Name column"""

    # Query for relevant features
    standard_df = gff_df.query("feature in @restrict_to")
    standard_df = add_gff_attributes(
        input_df=standard_df,
        field_names=["ID", "Name"]
    )

    # # Extract ID and Name fields
    # IDs = []
    # Names = []
    # for attributes in standard_df["attribute"]:

    #     # Split attributes
    #     fields = attributes.split(";")

    #     # Extract fields
    #     ID = [f.split("=")[1] for f in fields if f.startswith("ID=")]
    #     Name = [f.split("=")[1] for f in fields if f.startswith("Name=")]

    #     # Store
    #     IDs.append(ID[0].split(":")[1])
    #     Names.append(Name[0] if Name else None)

    # # Add to DataFrame
    # standard_df.insert(9, "ID", IDs)
    # standard_df.insert(10, "name", Names)

    return standard_df


# Prepare .gff standardisation
gff_standardisation_functions = {
    "plasmodb": standardise_PlasmoDB_gff,
    "ensemblgenomes": standardise_EnsemblGenomes_gff,
}
