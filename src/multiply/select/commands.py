import click
from .selectors import selector_collection


@click.command(short_help="Select optimal multiplex(es).")
@click.option(
    "-r",
    "--result_dir",
    type=click.Path(exists=True),
    required=True,
    help="Path to results directory for multiplex design (e.g. `results/2022-06-11_pf-default`).",
)
@click.option(
    "-a",
    "--algorithm",
    type=click.Choice(selector_collection),
    default="Greedy",
    help="Search algorithm for optimal multiplex. Note that `BruteForce` is exceedingly slow for large multiplexes.",
)
def select(result_dir, algorithm):
    """
    Select optimal multiplex(es) from a set of candidate primers

    It is assumed that the following commands have been run:
    `multiply generate`, `multiply align`, `multiply blast`.

    Information about primer quality, primer dimers, and off-target
    binding sites are fed into a cost function, which is then minimised.

    """
    from .main import select

    select(result_dir, algorithm)
