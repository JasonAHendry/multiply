import numpy as np
from dataclasses import dataclass, field
from multiply.generate.targets import Target


@dataclass
class Primer:
    seq: str
    direction: str
    start: int
    # end: int
    length: int
    tm: float
    gc: float
    # target: Target = None
    primer_name: str=""

    def give_name(self, name, primer_code, primer_ix):
        self.primer_name = f"{name}_{primer_code}{primer_ix:d}_{self.direction}"

    def add_tail(self, tail_seq):
        """
        Add a tail sequence to the 5' end of a Primer

        """
        
        assert isinstance(tail_seq, str), "Tail sequence must be a string."
        length_tail = len(tail_seq)
        self.seq = tail_seq + self.seq
        self.length += length_tail

@dataclass
class PrimerPair:
    """
    Define a pair of primers

    TODO:
    - A nice __repr__?
    - A joint name field?
    - Consider how to implement uniqueness; does target need to be defined?

    """

    F: Primer
    R: Primer
    product_bp: int
    pair_penalty: float
    pair_id: str = field(default="", repr=False)
    target: Target = None

    def __post_init__(self):
        """
        There might be a better way to do this; do I need target
        to ensure uniqueness?

        What about identical primers that anneal to different places?

        """

        F_info = f"{self.F.start}:{self.F.seq}"
        R_info = f"{self.R.start}:{self.R.seq}"
        self.pair_id = f"{F_info}+{R_info}"

    def get_primer_as_dict(self, direction, add_product_info=True, add_target_info=True, add_pair_id=False, add_pair_name=True):
        """
        Get either the forward or reverse primer, as a dictionary
        
        """

        if direction == "F":
            primer_info = self.F.__dict__.copy()
        elif direction == "R":
            primer_info = self.R.__dict__.copy()
        else:
            raise ValueError("Primer direction must be in ['F', 'R'].")

        if add_product_info:
            primer_info.update({
                "product_bp": self.product_bp,
                "pair_penalty": self.pair_penalty,
            })

        if add_target_info:
            primer_info.update({
                "target_id": self.target.ID,
                "target_name": self.target.name,
                "chrom": self.target.chrom,
            })
        
        if add_pair_id:
            primer_info["pair_id"] = self.pair_id

        if add_pair_name:
            primer_info["pair_name"] = self.pair_name

        return primer_info

    def give_primers_names(self, primer_code, primer_ix):
        """
        Give both primers names

        In theory, I think I violating some code coupling
        laws here
        
        """
        if self.target is None:
            raise ValueError("Must specify '.target' before giving primers names.")

        self.F.give_name(self.target.name, primer_code, primer_ix)
        self.R.give_name(self.target.name, primer_code, primer_ix)
        self.pair_name = f"{self.target.name}_{primer_code}{primer_ix:d}"

    # Allow set(), specifically on self.pair_id
    def __hash__(self):
        return hash(self.pair_id)

    def __eq__(self, other):
        if not isinstance(other, PrimerPair):
            return NotImplemented
        return self.pair_id == other.pair_id


def load_primer_pairs_from_primer3_output(primer3_output_path, add_target=None):
    """
    Given an output file from primer3, return a list of
    PrimerPair objects

    params
        primer3_output_path: str
            Path to an output file produced by primer3. This will
            contain information a series of primer pairs, in a format
            <key>=<value>.
        add_target: Target [optional]
            A Target object, containing information about on which target
            primer3 run.

    """

    # Parse primer3 output
    with open(primer3_output_path, "r") as f:
        # Iterate until determine number of primers returned
        for line in f:
            if line.startswith("PRIMER_PAIR_NUM_RETURNED"):
                n_returned = int(line.strip().split("=")[1])
                break

        # Return an empty list if no primers discovered
        if n_returned == 0:
            return []

        # Store remaining lines in a dictionary
        primer3_dt = {k: v.strip() for k, v in [l.split("=") for l in f]}

    # Define indexes and directions for primer pairs returned
    ixs = np.arange(n_returned)
    directions = ["LEFT", "RIGHT"]

    # Iterate over pairs
    primer_pairs = []
    for ix in ixs:

        # Get information about individual primers
        pair = {}
        for d in directions:

            primer_name = f"PRIMER_{d}_{ix}"
            s, l = primer3_dt[primer_name].split(",")
            s = int(s)
            l = int(l)

            # Note now that the start from primer3 is from
            # beginning of pad
            if add_target is not None:
                s += add_target.pad_start  # TODO: might be off by one

            pair[d] = Primer(
                seq=primer3_dt[f"{primer_name}_SEQUENCE"],
                direction="F" if d == "LEFT" else "R",
                start=s,
                length=l,
                tm=float(primer3_dt[f"{primer_name}_TM"]),
                gc=float(primer3_dt[f"{primer_name}_GC_PERCENT"]),
            )

        # Get information about pair
        pair_name = f"PRIMER_PAIR_{ix}"
        primer_pair = PrimerPair(
            F=pair["LEFT"],
            R=pair["RIGHT"],
            product_bp=int(primer3_dt[f"{pair_name}_PRODUCT_SIZE"]),
            pair_penalty=float(primer3_dt[f"{pair_name}_PENALTY"]),
            target=add_target
        )

        # Stores
        primer_pairs.append(primer_pair)

    return primer_pairs
