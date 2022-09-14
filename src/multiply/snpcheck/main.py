import os
import json

import pandas as pd
from functools import partial, reduce

from multiply.util.printing import print_header, print_footer
from multiply.util.dirs import produce_dir
from multiply.util.io import write_primers_to_bed
from multiply.download.collection import genome_collection
from .bedtools import bedtools_intersect


def snpcheck(primer_csv, genome_name):
    """
    Check primers in `primer_csv` for the presence of any SNPs, as defined
    by the `genome_name`

    Variation data for a given `genome_name` should be found in the
    `genomes/variation` folder.

    """
    # PARSE CLI
    t0 = print_header("MULTIPLY: Identify SNPs in primers")
    input_dir = os.path.dirname(primer_csv)
    output_dir = produce_dir(input_dir, "snpcheck")
    genome = genome_collection[genome_name]

    # LOAD DATA
    primer_df = pd.read_csv(primer_csv)

    # WRITE PRIMERS TO BED
    bed_path = f"{output_dir}/candidate_primers.bed"
    write_primers_to_bed(primer_df, bed_path)

    # GATHER VARIATION DATA
    print("Gathering variation data...")
    if not genome.include_variation:
        print(f"No variation data is available for {genome_name}. Exiting.")
        return
    print(f"  Found data at: {genome.include_variation}")
    variation_dt = json.load(open(genome.include_variation, "r"))
    print(f"  Includes: {', '.join(variation_dt)}")
    print("Done.")
    print("")

    # ITERATE OVER VARIATION FILES
    # - Might want to also output information about overlaps
    dfs = []
    for pop, pop_fn in variation_dt.items():

        # Parse
        print(f"Looking for intersections with: {pop}")
        print(f"  Variation file: {pop_fn}")
        print(f"  Intersecting...")

        # Intersect
        pop_output_path = f"{output_dir}/primers.snp_counts.{pop}.bed"
        bedtools_intersect(
            a=bed_path, b=pop_fn, flags=["-c"], output_path=pop_output_path
        )
        print(f"  Output: {pop_output_path}")

        # Load as dataframe
        df = pd.read_csv(
            pop_output_path,
            names=["chrom", "primer_start", "primer_end", "primer_name", pop],
            sep="\t",
        )

        # Summary
        print(f"  Total No. primers with SNPs: {(df[pop] > 0).sum()}")
        print(f"  Total No. SNPs found in primers: {df[pop].sum()}")

        # Store
        dfs.append(df)
        print("Done.")
        print("")

    # MERGE RESULT, WRITE SUMMARY
    print("Merging and writing summary CSV...")
    # Path
    output_path = f"{output_dir}/table.candidate_primers.snp_counts.csv"
    print(f"  to: {output_path}")

    # Merge & write
    merge_func = partial(
        pd.merge, on=["chrom", "primer_start", "primer_end", "primer_name"]
    )
    merged_df = reduce(merge_func, dfs)
    merged_df.to_csv(output_path, index=False)
    print("Done.")
    print("")
    print_footer(t0)
