import pysam
import pandas as pd
from dataclasses import dataclass, field
from multiply.util.exceptions import NoTargetsFoundError, TargetSizeError, TargetPositionError


# ================================================================================
# Encapsulate a single multiplex PCR target
#
# ================================================================================


@dataclass(order=True)  # these get sorted
class Target:
    """
    Define a target for PCR

    """

    chrom: str
    start: int
    end: int
    ID: str = field(compare=False)
    name: str = field(default="", compare=False)
    strand: str = field(default=".", compare=False)
    length: int = field(default=0, compare=False, repr=False)
    pad_start: int = field(default=0, compare=False, repr=False)
    pad_end: int = field(default=0, compare=False, repr=False)
    seq: str = field(default="", compare=False, repr=False)

    @classmethod
    def from_series(cls, series):
        """
        Create Target from a pandas Series
        
        """

        return cls(
            ID=series["ID"],
            name=series["name"] if "name" in series else None,
            strand=series["strand"] if "strand" in series else ".",
            chrom=series["seqname"],
            start=series["start"],
            end=series["end"],
        )

    def __post_init__(self):
        """
        Compute the length of the target based on the start
        and end position

        Note that we require a name; these are used to
        define primer names later on.

        """

        self.length = self.end - self.start
        self.name = self.ID if not self.name else self.name

    def calc_pads(self, max_size_bp):
        """
        Compute the start and end position of the `pads`, these define
        the maximum extent of the region inside of which primers may be 
        found

        """

        pad = max_size_bp / 2
        self.pad_start = int(self.start - pad)
        self.pad_end = int(self.end + pad)

        return self

    def extract_seq(self, reference_fasta_path, include_pads=True):
        """
        Given a path to a reference genome .fasta file, `reference_fasta_path`,
        extract the sequence of the target

        """

        self.pads_included = include_pads
        with pysam.FastaFile(reference_fasta_path) as fasta:

            # Define start and end of sequence to extract
            if include_pads:
                if not self.pad_start or not self.pad_end:
                    raise ValueError(
                        "If `include_pads` is True, must run `.calc_pads()` first."
                    )
                start, end = self.pad_start, self.pad_end
            else:
                start, end = self.start, self.end

            self.seq = fasta.fetch(self.chrom, start, end)

        return self


# ================================================================================
# Encapsulate a set of Targets for multiplex PCR
#
# ================================================================================


class TargetSet:
    def __init__(self, targets):
        """
        Collect and co-ordinate information about a set of Target objects

        """

        if len(targets) == 0:
            raise NoTargetsFoundError(f"No targets found.")
        self.targets = targets.copy()
        self.targets.sort()

    def check_size_compatible(self, max_size_bp):
        """
        Check that all none of the targets are larger than
        the maximum amplicon size `max_size_bp`

        TODO:
        - Definitely can make this print out nicer, and I should
        because this will probably throw a lot
        - Custom exception, for example

        params
            max_size_bp: int
                Maximum size of multiplex PCR amplicons, in basepairs.

        returns
            None

        """
        # Store maximum amplicon size, used in `.calc_pads()`
        self.max_size_bp = max_size_bp

        # Find targets that are too large
        too_large = [
            target for target in self.targets if target.length > self.max_size_bp
        ]

        # Throw warning
        if too_large:
            msg_info = "; ".join([f"{target.ID} [{target.length}bp]" 
            for target in too_large])
            raise TargetSizeError(
                f"Maximum amplicon size is {max_size_bp}bp, "
                f"but {len(too_large)} targets are larger than this: "
                f"{msg_info}.\n"
            )

        return self

    def _adjust_overlapping_pads(self):
        """
        For any targets that have overlapping pads, adjust them
        such that the intervening difference is split equally

        NB: the algorithm below assumes the targets are sorted;
        but sorting occurs in initialisation method.

        """

        for (left, right) in zip(self.targets[:-1], self.targets[1:]):

            # No possibility of overlap if on different chromosomes
            if left.chrom != right.chrom:
                continue

            # Ensure targets themselves do not overlap
            bp_bw_targets = right.start - left.end
            if bp_bw_targets <= 0:
                raise TargetPositionError(
                    f"Targets {left.ID} and {right.ID} overlap. Cannot build multiplex."
                )

            # Check if pads overlap
            bp_bw_pads = right.pad_start - left.pad_end
            if bp_bw_pads <= 0:
                print(f"Pads overlap between {left.ID} and {right.ID}")
                print(
                    f"Automatically adjusting. Note this may compromise ability to find primers later on."
                )

                # If pads overlap, adjust to use split space b/w targets equally
                middle_point = int(left.end + bp_bw_targets / 2)
                right.pad_start = middle_point + 1
                left.pad_end = middle_point

    def calc_pads(self, max_size_bp=None):
        """Calculate pads consistently across all targets"""

        if max_size_bp is not None:
            self.max_size_bp = max_size_bp

        # First, calculate all pads independently
        for target in self.targets:
            target.calc_pads(max_size_bp=self.max_size_bp)

        # Adjust any overlapping pads, if necesssary
        self._adjust_overlapping_pads()

        return self

    def extract_seqs(self, reference_fasta_path, include_pads=True):
        """Extract sequences for every target in the set"""

        for target in self.targets:
            target.extract_seq(reference_fasta_path, include_pads)

        return self

    def to_csv(
        self, 
        csv_path,
        keep_columns=["ID", "name", "chrom", "start", "end", "length", "strand", "pad_start", "pad_end"]
    ):
        """
        Write all targets to an output .csv

        """

        self.targets_df = pd.DataFrame(self.targets)
        self.targets_df = self.targets_df[keep_columns]
        self.targets_df.to_csv(csv_path, index=False)

        return self

    def to_fasta(self, fasta_path):
        """ Write all targets to an ouptut .fasta """

        with open(fasta_path, "w") as fasta:
            for target in self.targets:
                header = f">ID={target.ID}|name={target.name}"
                header += f"|ORF={target.chrom}:{target.start}-{target.end}"
                header += f"|PRIMER_PAD_REGION={target.chrom}:{target.pad_start}-{target.pad_end}|"
                header += f"PADS_INCLUDED={target.pads_included}\n"
                seq = f"{target.seq}\n"
                fasta.write(header)
                fasta.write(seq)

        return self
