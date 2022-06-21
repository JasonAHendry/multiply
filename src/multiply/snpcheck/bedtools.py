import os
import subprocess


def bedtools_intersect(a, b, flags, output_path):
    """
    Run `bedtools intersect`

    """

    # Better to catch here, as subprocess hides traceback
    if not os.path.exists(a):
        raise FileNotFoundError(f"File -a {a} does not exist.")
    if not os.path.exists(b):
        raise FileNotFoundError(f"File -b {b} does not exist.")

    # Construct command
    cmd = "bedtools intersect"
    cmd += f" -a {a}"
    cmd += f" -b {b}"
    cmd += f" {' '.join(flags)}"
    cmd += f" > {output_path}"

    # Run
    subprocess.run(cmd, check=True, shell=True)
