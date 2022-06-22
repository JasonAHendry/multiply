import os
import sys
import datetime


def print_header():
    """
    Print header for a script

    params
        None
    returns
        t0 : datetime.datetime
            Time when `print_header()` was run.

    """

    command = os.path.basename(sys.argv[0])
    args = " ".join(sys.argv[1:])

    t0 = datetime.datetime.now().replace(microsecond=0)
    print("=" * 80)
    print(f"Command: {command} {args}")
    print(f"Started at: {t0.strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 80)

    return t0


def print_footer(t0):
    """
    Print footer for a script

    params
        t0 : datetime object
            Time script was initialised, passed from
            `print_header`. Used to compute runtime.
    returns
        None

    """

    t1 = datetime.datetime.now().replace(microsecond=0)
    print("-" * 80)
    print("Finished at: %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Elapsed: %s" % (t1 - t0))
    print("=" * 80)


def print_parameters(design, params):
    """
    Print an overview of the parsed parameters

    """
    print("Design parameters")
    print(f"  Input file: {design}")
    print(f"  Genome: {params['genome']}")
    print(f"  Include region(s): {params['from_regions']}")
    if params["from_regions"]:
        print(f"    Region BED: {params['region_bed']}")
    print(f"  Include gene(s): {params['from_genes']}")
    if params["from_genes"]:
        print(f"    No. genes: {len(params['target_ids'])}")
        print(f"    Gene IDs: {', '.join(params['target_ids'][:3])}...")
    print(f"  Include primer tails: {params['include_tails']}")
    if params["include_tails"]:
        print(f"    F tail: {params['F_tail']}")
        print(f"    R tail: {params['R_tail']}")
    print(f"  Amplicon size range: {params['min_size_bp']}-{params['max_size_bp']}bp")
    print(f"  primer3 settings: {', '.join(params['primer3_settings'])}")
    print(f"  Output directory: {params['output_dir']}")
    print("Done.\n")
