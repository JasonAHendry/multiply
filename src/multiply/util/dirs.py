import os
import datetime


def produce_dir(*args):
    """
    Produce a new directory by concatenating `args`,
    if it does not already exist

    params
        *args: str1, str2, str3 ...
            Comma-separated strings which will
            be combined to produce the directory,
            e.g. str1/str2/str3

    returns
        dir_name: str
            Directory name created from *args.

    """

    # Define directory path
    dir_name = os.path.join(*args)

    # Create if doesn't exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    return dir_name


def create_output_directory(params):
    """Create output directroy for an experiment"""

    # Create a date stamped output directory
    today = datetime.datetime.today().strftime("%Y-%m-%d")
    output_dir = f"results/{today}_{params['output_name']}"
    produce_dir(output_dir)

    return output_dir
