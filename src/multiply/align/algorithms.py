import json
from numba import njit
from numba.typed import Dict
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from .nn_model import create_nn_score_dt


# ================================================================================
# Define an alignment between two primers
#
# ================================================================================


@dataclass(order=True)
class PrimerAlignment:
    """
    Represent the alignment of two primers
    """

    primer1_name: str = field(compare=False)
    primer2_name: str = field(compare=False)
    primer1: str = field(compare=False)
    primer2: str = field(compare=False)
    score: float = field(compare=True)
    alignment: str = field(compare=False, repr=False)
    # index: int=field(compare=False, repr=False)


# ================================================================================
# Abstract base class for various primer alignment algorithms
#
# ================================================================================


class AlignmentAlgorithm(ABC):
    """
    Alignment algorithm for a pair of primers

    """

    def __init__(self):
        pass

    def set_primers(self, primer1, primer2, primer1_name, primer2_name):
        """
        Set a pair for primers to align

        """
        self.primer1 = primer1
        self.primer2 = primer2
        self.primer1_name = primer1_name
        self.primer2_name = primer2_name

        self.score = None  # reset

    @abstractmethod
    def load_parameters():
        pass

    @abstractmethod
    def align():
        """
        Align the primers

        """
        pass

    @abstractmethod
    def get_alignment_string():
        """
        Create an ASCII string representing the aligned primers

        """
        pass

    def print_alignment(self):
        """
        Print an ASCII view of the aligned primers

        """
        print(self.get_alignment_string())

    def get_primer_alignment(self):
        """
        Return an alignment object

        """
        return PrimerAlignment(
            primer1=self.primer1,
            primer2=self.primer2,
            primer1_name=self.primer1_name,
            primer2_name=self.primer2_name,
            score=self.score,
            alignment=self.get_alignment_string(),
        )


# ================================================================================
# Concrete primer alignment algorithms
#
# ================================================================================


class PrimerDimerAlgorithm(AlignmentAlgorithm):
    """
    Align two primers using ann algorithm similar to Primer Dimer
    from TK et al.

    RUNNING -- but needs a lot of testing / re-factoring
    Plus want a Numba version

    """

    param_path = "settings/alignment/primer_dimer/parameters.json"
    rc_map = {"A": "T", "T": "A", "C": "G", "G": "C"}  # Make this a global

    def load_parameters(self):
        """
        Load parameters necessary for Primer Dimer algorithm,
        and set as attributes
        
        """
        # Load parameter JSON
        params = json.load(open(self.param_path, "r"))

        # Load nearest neighbour model
        self.nn_scores = create_nn_score_dt(
            match_json=params["match_scores"],
            single_mismatch_json=params["single_mismatch_scores"]
        )

        # Load penalties
        self.end_length = params["end_length"]
        self.end_penalty = params["end_penalty"]
        self.end_bonus = params["end_bonus"]

    def align(self):
        """
        Calculate a score for a pair of primers that is
        related to their likelihood of forming a
        primer dimer, inspired by PrimerDimer

        A change in Gibb's free energy is calculated
        based on the extent of complementary base pairing

        End penalties / bonus are added to the delta G to
        return the final score.

        See: http://www.primer-dimer.com/

        """

        # Redfine sequences based on length
        n1 = len(self.primer1)
        n2 = len(self.primer2)

        if n1 > n2:
            l, nL = self.primer1, n1
            s, nS = self.primer2[::-1], n2
        else:
            l, nL = self.primer2, n2
            s, nS = self.primer1[::-1], n1

        # Iterate
        start_ix = 0
        dimer_score = 100
        for i in range(nL - 2 + 1):
            if i + nS > nL:
                w = nL - i - 2 + 1
            else:
                w = nS - 2 + 1

            current_score = 0
            penalties = 0
            for j in range(w):
                nn = l[(i + j) : (i + j + 2)]
                nn += "/"
                nn += s[j : (j + 2)]
                current_score += self.nn_scores[nn]

                # 3' end penalties
                if j < self.end_length or (nL - 1 - i - j) < self.end_length:
                    if self.rc_map[l[i + j]] != s[j]:
                        penalties += self.end_penalty
                if i + j == nL - 2:
                    if self.rc_map[l[i + j + 1]] != s[j + 1]:
                        penalties += self.end_penalty

            # Add penalties
            current_score += penalties

            # Add end bonus if no penalties
            if penalties == 0:
                current_score += self.end_bonus

            if current_score < dimer_score:
                dimer_score = current_score
                start_ix = i

        # SET
        self.score = dimer_score
        self.start_ix = start_ix

    def get_alignment_string(self):
        """
        Could definitely be cleaned up...
        """

        # Redfine sequences based on length
        n1 = len(self.primer1)
        n2 = len(self.primer2)

        if n1 > n2:
            l, nL, lname = self.primer1, n1, self.primer1_name
            s, nS, sname = self.primer2[::-1], n2, self.primer2_name
        else:
            l, nL, lname = self.primer2, n2, self.primer2_name
            s, nS, sname = self.primer1[::-1], n1, self.primer1_name

        # Define long and short primer strings
        lstr = "'5-" + l + "-3'"
        sstr = " " * self.start_ix + "'3-" + s + "-5'"

        # Determine string specifying match / mismatch
        match = " " * 3  # compensate for 3'-, 5'-
        match += " " * self.start_ix  # offset of alignment
        for s_ix, l_nt in enumerate(l[self.start_ix : (self.start_ix + nS)]):
            if self.rc_map[l_nt] == s[s_ix]:
                match += "|"
            else:
                match += " "

        # Sequence names
        if lname and sname:
            nc = max([len(lname), len(sname)])
        else:
            nc = 1
            lname = "L"
            sname = "S"

        # Template to allow for names of varying character length
        seq_str = "{:>%d}:    {}" % nc
        middle_str = "{:>%d}     {}" % nc

        # Final print
        self.alignment = f"Dimer Score: {self.score:.02f}\n"
        self.alignment += f"{seq_str.format(lname, lstr)}\n"
        self.alignment += f"{middle_str.format('', match)}\n"
        self.alignment += f"{seq_str.format(sname, sstr)}\n"

        return self.alignment


class PrimerDimerAlgorithmNumba(AlignmentAlgorithm):
    """
    A quick Numba implementation of the Primere Dimer
    Algorithm above

    A few points are that: (1) jitclass didn't work natively;
    (2) you can't use instances variables, hence make a static method;
    (3) Numba has it's own dictionary typing, which I convert to
    in load_parameters()

    I found a modest performance improvement of about 50%.

    """

    param_path = "settings/alignment/primer_dimer/parameters.json"
    rc_map = {"A": "T", "T": "A", "C": "G", "G": "C"}  # Make this a global

    def load_parameters(self):
        """
        Load parameters necessary for Primer Dimer algorithm,
        and set as attributes
        
        """
        # Load parameter JSON
        params = json.load(open(self.param_path, "r"))

        # Load nearest neighbour model
        self.nn_scores = create_nn_score_dt(
            match_json=params["match_scores"],
            single_mismatch_json=params["single_mismatch_scores"]
        )

        # Load penalties
        self.end_length = params["end_length"]
        self.end_penalty = params["end_penalty"]
        self.end_bonus = params["end_bonus"]
        
        # Create Numba versions of dictionaries
        self.nn_scores_numba = Dict()
        for k, v in self.nn_scores.items():
            self.nn_scores_numba[k] = v
            
        self.rc_map_numba = Dict()
        for k, v in self.rc_map.items():
            self.rc_map_numba[k] = v
        
        
    def align(self):
        """ Wrapper for jit implementation of _align() """
        
        
        self.score, self.start_ix = \
        self._align(
            primer1=self.primer1,
            primer2=self.primer2,
            rc_map=self.rc_map_numba,
            nn_scores=self.nn_scores_numba,
            end_length=self.end_length,
            end_penalty=self.end_penalty,
            end_bonus=self.end_bonus
        )
        
    @staticmethod
    @njit
    def _align(primer1, primer2, rc_map, nn_scores, end_length, end_penalty, end_bonus):
        """
        Calculate a score for a pair of primers that is
        related to their likelihood of forming a
        primer dimer, inspired by PrimerDimer

        A change in Gibb's free energy is calculated
        based on the extent of complementary base pairing

        End penalties / bonus are added to the delta G to
        return the final score.

        See: http://www.primer-dimer.com/

        """

        # Redfine sequences based on length
        n1 = len(primer1)
        n2 = len(primer2)

        if n1 > n2:
            l, nL = primer1, n1
            s, nS = primer2[::-1], n2
        else:
            l, nL = primer2, n2
            s, nS = primer1[::-1], n1

        # Iterate
        start_ix = 0
        dimer_score = 100
        for i in range(nL - 2 + 1):
            if i + nS > nL:
                w = nL - i - 2 + 1
            else:
                w = nS - 2 + 1

            current_score = 0
            penalties = 0
            for j in range(w):
                nn = l[(i + j) : (i + j + 2)]
                nn += "/"
                nn += s[j : (j + 2)]
                current_score += nn_scores[nn]

                # 3' end penalties
                if j < end_length or (nL - 1 - i - j) < end_length:
                    if rc_map[l[i + j]] != s[j]:
                        penalties += end_penalty
                if i + j == nL - 2:
                    if rc_map[l[i + j + 1]] != s[j + 1]:
                        penalties += end_penalty

            # Add penalties
            current_score += penalties

            # Add end bonus if no penalties
            if penalties == 0:
                current_score += end_bonus

            if current_score < dimer_score:
                dimer_score = current_score
                start_ix = i

        # SET
        score = dimer_score
        start_ix = start_ix
        
        return score, start_ix

    def get_alignment_string(self):
        """
        Could definitely be cleaned up...
        """

        # Redfine sequences based on length
        n1 = len(self.primer1)
        n2 = len(self.primer2)

        if n1 > n2:
            l, nL, lname = self.primer1, n1, self.primer1_name
            s, nS, sname = self.primer2[::-1], n2, self.primer2_name
        else:
            l, nL, lname = self.primer2, n2, self.primer2_name
            s, nS, sname = self.primer1[::-1], n1, self.primer1_name

        # Define long and short primer strings
        lstr = "'5-" + l + "-3'"
        sstr = " " * self.start_ix + "'3-" + s + "-5'"

        # Determine string specifying match / mismatch
        match = " " * 3  # compensate for 3'-, 5'-
        match += " " * self.start_ix  # offset of alignment
        for s_ix, l_nt in enumerate(l[self.start_ix : (self.start_ix + nS)]):
            if self.rc_map[l_nt] == s[s_ix]:
                match += "|"
            else:
                match += " "

        # Sequence names
        if lname and sname:
            nc = max([len(lname), len(sname)])
        else:
            nc = 1
            lname = "L"
            sname = "S"

        # Template to allow for names of varying character length
        seq_str = "{:>%d}:    {}" % nc
        middle_str = "{:>%d}     {}" % nc

        # Final print
        self.alignment = f"Dimer Score: {self.score:.02f}\n"
        self.alignment += f"{seq_str.format(lname, lstr)}\n"
        self.alignment += f"{middle_str.format('', match)}\n"
        self.alignment += f"{seq_str.format(sname, sstr)}\n"

        return self.alignment