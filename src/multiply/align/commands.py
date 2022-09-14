import click


@click.command(short_help="Search for primer dimers.")
@click.option(
    "-p",
    "--primer_csv",
    type=click.Path(exists=True),
    required=True,
    help="Path to candidate primer CSV file (e.g. `table.candidate_primers.csv`).",
)
def align(primer_csv):
    """
    Run a pairwise alignment between all primers in `primer_csv` to identify
    potential primer dimers

    """
    from .main import align

    align(primer_csv)
