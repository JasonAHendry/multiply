import click
from multiply.download.collection import genome_collection


@click.command(short_help="BLAST candidate primers.")
@click.option(
    "-p",
    "--primer_csv",
    type=click.Path(exists=True),
    required=True,
    help="Path to candidate primer CSV file (e.g. `table.candidate_primers.csv`).",
)
@click.option(
    "-g",
    "--genome_name",
    type=click.Choice(genome_collection),
    required=True,
    help="Name of genome to download.",
)
def blast(primer_csv, genome_name):
    """
    BLAST primers in `primer_csv` against reference genome given by `genome_name` to
    find off-target binding sites and amplicons

    Visit `settings/blast` for blast settings.

    """
    from .main import blast

    blast(primer_csv, genome_name)
