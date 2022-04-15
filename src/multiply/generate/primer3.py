import json
import subprocess
import numpy as np


class Primer3Runner:

    settings_dir = "settings/primer3"

    def __init__(self):
        """
        Run primer3, including:

        - Loading primer3 settings
        - Preparing possible amplicon sizes
        - Choosing a target
        - Parsing outputs

        Want to include checks to make sure functions are run in a
        sensible order

        """

        # Ready state
        self.settings_loaded = False
        self.amplicon_sizes_set = False
        self.target_selected = False

    def load_primer3_settings(self, setting_name="default"):
        """
        Load primer3 settings stored as a .json file

        """

        # Load the settings
        self.setting_name = setting_name
        self.settings = json.load(open(f"{self.settings_dir}/{setting_name}.json", "r"))

        # Record
        self.settings_loaded = True

        return self

    def set_amplicon_size_ranges(self, min_size_bp, max_size_bp, n_intervals=3):
        """
        Set the size ranges of the amplicons, i.e. prepare the primer3
        string for 'PRIMER_PRODUCT_SIZE_RANGE'

        'PRIMER_PRODUCT_SIZE_RANGE' accepts windows of
        amplicon sizes, ordered by preference,
        e.g. "200-400 400-800 800-2000"

        By default, this function creates a string that
        prioritizes small amplicons. This can be reversed by
        setting `min_size_bp` to larger than `max_size_bp`.

        params
            min_size_bp : int
                Minimum amplicon size.
            max_size_bp : int
                Maximum amplicon size.
            n_intervals : int
                Number of size intervals.

        returns
            sizes : string
                Primer3 amplicon size string.

        """

        # Define cuts between interals
        cuts = np.linspace(min_size_bp, max_size_bp, n_intervals + 1, dtype="int")

        # Create sizes string for primer3
        sizes = ""
        for start, end in zip(cuts[:-1], cuts[1:]):
            sizes += "%d-%d " % (start, end)

        # Store
        self.sizes = sizes.rstrip()

        # Update settings
        self.settings["PRIMER_PRODUCT_SIZE_RANGE"] = self.sizes

        # Record
        self.amplicon_sizes_set = True

        return self

    def set_target(self, *, ID, seq, pad_start, start, length):
        """
        Set which target primers will be generated for, defining
        the revelant primer3 input fields.

        """
        # Store
        self.ID = ID
        self.seq = seq
        self.pad_start = pad_start
        self.start = start
        self.length = length

        # Compute start and end, relative to `seq`
        self.target_start = start - pad_start

        # Define target specific settings
        self.target_specific_settings = {
            "SEQUENCE_ID": self.ID,
            "SEQUENCE_TEMPLATE": self.seq,
            "SEQUENCE_TARGET": f"{self.target_start},{self.length}",
        }

        # Update settings
        self.settings.update(self.target_specific_settings)

        # Record
        self.target_selected = True

        return self

    def write_primer3_input(self, input_path):
        """
        Write a primer3 input file based on a dictionary
        of primer3 settings

        params
            primer3_settings : dict
                Dictionary of primer3 settings.
            input_path : str
                Write primer3 input file here.

        returns
            None

        """

        with open(input_path, "w+") as fn:
            fn.write("=\n")
            for k, v in self.settings.items():
                fn.write("%s=%s\n" % (k, str(v)))
            fn.write("=\n")

        return None

    def run(self, output_dir):
        """
        Run primer3 on defined settings

        """

        # Basic idea, but would want to improve
        if not (self.settings_loaded and self.amplicon_sizes_set and self.set_target):
            raise ValueError(
                "Ensure settings are loaded, amplicon sizes and target have been set."
            )

        # Prepare the input file
        self.input_path = f"{output_dir}/{self.ID}.primer3.input"
        self.output_path = self.input_path.replace("input", "output")
        self.write_primer3_input(self.input_path)

        # Run the primer3 as a subprocesss
        cmd = "primer3_core %s > %s" % (self.input_path, self.output_path)
        subprocess.run(cmd, shell=True, check=True)
