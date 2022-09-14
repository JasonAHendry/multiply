import click
from multiply.download.collection import genome_collection


@click.command(short_help="Visualise candidate primer locations.")
@click.option(
    "-r",
    "--result_dir",
    type=click.Path(exists=True),
    required=True,
    help="Result directory containing: `table.candidate_primers.csv`, `table.targets_overview.csv`, and `targets_sequence.fasta`.",
)
@click.option(
    "-g",
    "--genome_name",
    type=click.Choice(genome_collection),
    required=True,
    help="Name of genome to download.",
)
def view(result_dir, genome_name):
    """
    View binding locations of candidate primers

    """
    from .main import view

    view(result_dir, genome_name)
