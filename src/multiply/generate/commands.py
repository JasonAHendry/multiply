import click


@click.command(short_help="Generate candidate primers.")
@click.option(
    "-d",
    "--design",
    type=str,
    required=True,
    help="Path to MULTIPLY design file (e.g. 'designs/pf-default.ini').",
)
def generate(design):
    """
    Generate a pool of candidate primers for a given `design` using
    primer3

    Visit the `settings/primer3` directory to observe or change primer3
    settings.

    """
    from .main import generate

    generate(design)
