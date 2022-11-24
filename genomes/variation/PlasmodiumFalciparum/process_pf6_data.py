# Small script to process biallelic SNPs from Pf6
# into only common SNPs, partitioned by geography
#
# Uses a version of Pf6 downloaded locally on my machine
#
# Subpopulations in Pf6: CAF, EAF, ESEA, OCE, SAM, SAS, WAF, WSEA, GLOBAL
#
# J.Hendry, 2022/06/20

import click
import subprocess
import pandas as pd
import numpy as np


# PARAMETERS
N_CHROM = 14
CHROMS = np.arange(1, N_CHROM + 1)
VCF_DIR = "../../../../pf-snpcheck/data/pf6/vcfs"
VCF_FILE_TEMP = "Pf_60_public_Pf3D7_%.02d_v3.high_quality_biallelic_snps.vcf.gz"
METADATA_PATH = "../../../../pf-snpcheck/data/pf6/metadata/Pf_60_candidate_public_metadata_20170919_cleaned.csv"
MIN_AF = 0.05


@click.command()
@click.option(
    "-t",
    "--target_population",
    required=True,
    type=str,
    help="Target population from Pf6.",
)
def main(target_population):
    # Introduction
    print("Processing Pf6 data for MULTIPLY")
    print(f"  VCF source dir.: {VCF_DIR}")
    print(f"  Minimum allele frequency: {MIN_AF}")
    print("")

    # Load metadata
    print("Loading metadata...")
    metadata = pd.read_csv(METADATA_PATH)
    print(f"  From: {METADATA_PATH}")
    print(f"  Total no. samples {metadata.shape[0]}")
    print(f"Done.")
    print("")

    # Extract samples per population
    pop_samples_dt = {
        pop: pop_df["sample"].tolist() for pop, pop_df in metadata.groupby("population")
    }
    pop_samples_dt["GLOBAL"] = metadata["sample"].tolist()
    print(f"Target populations available: {', '.join(pop_samples_dt)}")

    # Iterate over populations
    for population, samples in pop_samples_dt.items():

        if population != target_population:
            continue

        pop_output_file = f"pf6.filtered_snps.{population}.bed"

        print(f"Population: {population}")
        print(f"  No. samples: {len(samples)}")
        print(f"  Samples: {', '.join(samples[:3])}...")

        for chrom in CHROMS:
            print(f"  Processing chromosome {chrom}...")

            cmd = "bcftools view"
            cmd += f" -s {','.join(samples)}"
            cmd += f" -q {MIN_AF}"
            cmd += f" {VCF_DIR}/{VCF_FILE_TEMP % chrom}"
            cmd += f" | bcftools query"
            cmd += " -f '%CHROM\\t%POS\\t%END\\t%ID\\n'"
            cmd += f" - >> {pop_output_file}"

            subprocess.run(cmd, check=True, shell=True)

        print("Done.")
        print("")


if __name__ == "__main__":
    main()
