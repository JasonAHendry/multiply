import click
import logging

from multiply.util.parsing import parse_parameters

import multiply.download.commands as download
import multiply.generate.commands as generate
import multiply.view.commands as view
import multiply.snpcheck.commands as snpcheck
import multiply.align.commands as align
import multiply.blast.commands as blast
import multiply.select.commands as select



# ================================================================
# Entry point for all sub-commands
#
# ================================================================


# TODO: I think this is where I set logging verbosity
# - Then we want that verbosity level to pass to the sub-modules, probably
@click.group()
def cli():
    """
    Multiplex PCR design, in silico

    """
    # Prepare logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    logger.addHandler(stream_handler)
    
    pass


# ================================================================
# Individual sub-commands
#
# ================================================================


cli.add_command(download.download)
cli.add_command(generate.generate)
cli.add_command(view.view)
cli.add_command(snpcheck.snpcheck)
cli.add_command(align.align)
cli.add_command(blast.blast)
cli.add_command(select.select)



# ================================================================
# Run as a pipeline
#
# ================================================================


@click.command(short_help="Run the full multiply pipeline.")
@click.option(
        "-d",
        "--design",
        type=str,
        required=True,
        help="Path to MULTIPLY design file (e.g. 'designs/pf-default.ini').",
    )
def pipeline(design):
    """
    Run the complete `multiply` pipeline
    
    """
    # PARSE INPUT PARAMETERS
    params = parse_parameters(design)
    primer_csv = f"{params['output_dir']}/table.candidate_primers.csv"

    # Pipeline
    generate.main(design=design)
    view.main(result_dir=params['output_dir'], genome_name=params['genome'])
    snpcheck.main(primer_csv=primer_csv, genome_name=params['genome'])
    align.main(primer_csv=primer_csv)
    blast.main(primer_csv=primer_csv, genome_name=params['genome'])
    select.main(result_dir=params['output_dir'], algorithm="Greedy")


cli.add_command(pipeline)


if __name__ == "__main__":
    cli()
