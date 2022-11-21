import click
import pandas as pd
from dataclasses import dataclass, field
from multiply.download.collection import genome_collection


# --------------------------------------------------------------------------------
# REGIONS -> BED file
#
# --------------------------------------------------------------------------------


@dataclass(order=True)
class TargetRegion:
    chrom: str
    start: int
    end: int
    name: str = field(default=None, compare=False)


def write_target_regions_to_bed(bed_path, regions):
    """Write a list of target regions to a BED file"""
    with open(bed_path, "w") as bed:
        for region in regions:
            bed.write(f"{region.chrom}\t{region.start}\t{region.end}\t{region.name}\n")


# --------------------------------------------------------------------------------
# MAIN
#
# --------------------------------------------------------------------------------


@click.command(short_help="Define multiplex PCR pools for contiguous tiling.")
@click.option(
    "-t", "--target_id", help="Target ID of gene to tile around.", required=True
)
@click.option("-m", "--min_amp_size", default=2000, help="Minimum amplicon size in bp.")
@click.option("-M", "--max_amp_size", default=3000, help="Maximum amplicon size in bp.")
@click.option("-p", "--pool_size", default=10, help="Number of amplicons per pool.")
@click.option(
    "-g",
    "--genome_name",
    type=click.Choice(genome_collection),
    default="PlasmodiumFalciparum",
    help="Name of genome to download.",
)
def create_pools(target_id, min_amp_size, max_amp_size, pool_size, genome_name):
    """
    Define two multiplex PCR pools for contiguous tiling around a given
    target gene, output as BED files

    """
    # Input parameters
    print("Input parameters")
    print(f"  Target ID: {target_id}")
    print(f"  Minimnum amplicon size (bp): {min_amp_size}")
    print(f"  Maximum amplicon size (bp): {max_amp_size}")
    print(f"  Number amplicons per pool: {pool_size}")
    print(f"  Target genome: {genome_name}")
    print("Done.\n")

    # Define genome
    genome = genome_collection[genome_name]

    # Define central region
    print("Defining central region...")
    gff = pd.read_csv(genome.gff_path)
    target_info = gff.query("ID == @target_id").squeeze()
    target_chrom = target_info["seqname"]
    target_size = int(target_info["end"]) - int(target_info["start"])
    if target_size > max_amp_size:
        raise ValueError(
            f"Script does not support target size of {target_size}bp"
            f"larger than maximum amplicon size of {max_amp_size}bp."
        )
    amp_excess = min_amp_size - target_size
    central_start = target_info["start"] - int(amp_excess / 2)
    central_end = target_info["end"] + int(amp_excess / 2)

    # Create region
    central = [
        TargetRegion(
            chrom=target_chrom,
            start=central_start,
            end=central_end,
            name=f"{target_id}",
        )
    ]

    # Define upstream and downstream regions, satisfying pool size
    downstream = [
        TargetRegion(
            chrom=target_chrom,
            start=central_end + (i - 1) * min_amp_size,
            end=central_end + i * min_amp_size,
            name=f"dwnstrm{i:02d}",
        )
        for i in range(1, pool_size)
    ]
    upstream = [
        TargetRegion(
            chrom=target_chrom,
            start=central_start - i * min_amp_size,
            end=central_start - (i - 1) * min_amp_size,
            name=f"upstrm{i:02d}",
        )
        for i in range(1, pool_size + 1)
    ]
    regions = upstream + central + downstream
    regions.sort()

    # Split into pools
    poolA = regions[::2]
    poolB = regions[1::2]
    print(f"Created {len(regions)}.")
    print(f"Split into two pools:")
    print(f"  Pool A: {len(poolA)}")
    print(f"  Pool B: {len(poolB)}")

    # Write as BED files
    print("Writing to BED files...")
    write_target_regions_to_bed(bed_path=f"tiling.{target_id}.poolA.bed", regions=poolA)
    write_target_regions_to_bed(bed_path=f"tiling.{target_id}.poolB.bed", regions=poolB)
    print("Done.\n")


if __name__ == "__main__":
    create_pools()
