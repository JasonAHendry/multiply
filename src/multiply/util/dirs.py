import os
import sys
import shutil


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

def check_output_dir_overwrite(output_dir):
    """
    Check if the `output_dir` already exists, 
    if so, ask before overwriting it or exit.
    
    """
    if os.path.exists(output_dir):
        
        print(f"Experiment directory with name {output_dir} already exists.")
        clean = input("Do you want to overwrite it? [Yes/No]: ")
       
        if clean == "Yes":
            shutil.rmtree(f"{output_dir}")
        else:
            print("Exiting. Change 'name' in [Output] section of design file before running again.")
            sys.exit()

