import click

# import logging
from multiply.download.commands import download
from multiply.generate.commands import generate


# ================================================================
# Decorators for commonly used options
#
# ================================================================


def design_path_option(fn):
    fn = click.option(
        "-d",
        "--design",
        type=str,
        help="Path to MULTIPLY design file (e.g. 'designs/pf-default.ini').",
    )(fn)
    return fn


# ================================================================
# Entry point for all commands
#
# ================================================================


# TODO: I think this is where I set logging verbosity
# - Then we want that verbosity level to pass to the sub-modules, probably
@click.group()
def cli():
    """
    MULTIPLY: Design multiplex PCRs in silico

    """
    # # Prepare logger
    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.INFO)
    # stream_handler = logging.StreamHandler()
    # logger.addHandler(stream_handler)
    pass


cli.add_command(download)
cli.add_command(generate)


# ================================================================
# Complete pipeline
#
# ================================================================


# @cli.command(short_help="Run complete pipeline.")
# @design_path_option
# def pipeline(design_path):
#     """
#     Run the complete multiplex PCR primer design process

#     """
#     pass


# # ================================================================
# # Individual commands
# #
# # ================================================================


# @cli.command(short_help="Map to P.f. reference.")
# @design_path_option
# def generate(design_path):
#     """
#     Map .fastq files to the P. falciparum reference genome

#     """
#     from multiply import generate

#     generate.main(design_path)

if __name__ == "__main__":
    cli()
