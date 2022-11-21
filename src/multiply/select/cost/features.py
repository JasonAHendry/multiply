import pandas as pd
import numpy as np


class IndividualCosts:
    def __init__(self, cost_name, primer_values, weight):

        # Sanity check
        self._check_primer_values(primer_values)

        # Assigned at instantiation
        self.cost_name = cost_name
        self.primer_values = primer_values
        self.weight = weight

        # To be assigned
        self.primer_pair_values = None  # Raw values
        self.primer_pair_costs = None  # Post-normalisation

    @staticmethod
    def _check_primer_values(primer_values):
        assert isinstance(
            primer_values, pd.Series
        ), "`primer_values` must be a pandas series."
        assert all(
            [p.endswith("_F") or p.endswith("_R") for p in primer_values.index]
        ), "`primer_values` index must be primer names."

    def collapse_to_per_pair(self, collapse_func=sum):
        """
        Collapse per primer values to per primer pair values,
        by aggregating with a collapsing fuction,
        `collapse_func`

        """

        # Build into data frame
        pair_df = pd.DataFrame(
            {
                "primer_pair_name": [n[:-2] for n in self.primer_values.index],
                "primer_pair_values": self.primer_values,
            }
        )

        # Collapse
        self.primer_pair_values = pair_df.groupby(
            "primer_pair_name"
        ).primer_pair_values.apply(collapse_func)

        return self

    def normalise_costs(self):
        """
        Z-score standardise, and then multiply by weight

        NB:
        - It is critical this normalisation handles costs
        that have no variation across primer pairs

        """

        if self.primer_pair_values is None:
            raise ValueError("Run `.collapse_to_per_pair()` first.")

        # Compute mean and standard deviation
        mu = self.primer_pair_values.mean()
        std = self.primer_pair_values.std()

        if not std > 0:
            print(
                f"No variation in {self.cost_name} observed across primer pairs."
                "Will not contribute to scoring."
            )
            std = 1

        # Compute normalisation
        self.primer_pair_costs = self.weight * (self.primer_pair_values - mu) / std

        return self


class PairwiseCosts:
    def __init__(self, cost_name, primer_values, weight):

        # Sanity check
        self._check_primer_values(primer_values)

        # Assigned at instantiation
        self.cost_name = cost_name
        self.primer_values = primer_values
        self.weight = weight

        # To be assigned
        self.primer_pair_values = None
        self.primer_pair_costs = None

    @staticmethod
    def _check_primer_values(primer_values):
        """
        Input data structure `primer_values` is quite specific,
        namely it is a square pandas dataframe with an index
        and column giving the primer names

        Check that these conditions are met here.

        """
        assert isinstance(
            primer_values, pd.DataFrame
        ), "`primer_values` must be a pandas DataFrame."
        assert (
            primer_values.shape[0] == primer_values.shape[1]
        ), "`primer_values` is a pairwise matrix oof primers, it should be square."
        assert all(
            [p.endswith("_F") or p.endswith("_R") for p in primer_values.index]
        ), "`primer_values` index must be primer names."
        assert all(
            primer_values.index == primer_values.columns
        ), "`primer_values` index and column names must be the same."

    def collapse_to_per_pair(self, collapse_func=sum):
        """
        Collapse from a dataframe where each column and row
        corresponds to an individual primer, to a dataframe
        where they correspond to an individual primer pair

        Note that the way that I achieve this below (in two steps)
        imposes some restrictions on the `collapse_func`.

        """

        primer_pair_name = [n[:-2] for n in self.primer_values.columns]

        # Add pair name column
        collapse_df = self.primer_values.copy()
        collapse_df.insert(0, "primer_pair_name", primer_pair_name)

        # Collapse first axis
        self.primer_pair_values = (
            collapse_df.groupby("primer_pair_name").aggregate(collapse_func).transpose()
        )

        # Repeat with other axis
        self.primer_pair_values.insert(0, "primer_pair_name", primer_pair_name)
        self.primer_pair_values = self.primer_pair_values.groupby(
            "primer_pair_name"
        ).aggregate(collapse_func)

        return self

    def normalise_costs(self):
        """
        Normalise to Z-scores and then multiply by weight

        NB:
        - It is critical this normalisation handles costs
        that have no variation across primer pairs

        """
        if self.primer_pair_values is None:
            raise ValueError("Run `.collapse_to_per_pair()` first.")

        # Compute mean and standard devation
        arr = np.array(self.primer_pair_values)
        mu = arr.mean()
        std = arr.std()

        if not std > 0:
            print(
                f"No variation in {self.cost_name} observed across primer pairs."
                "Will not contribute to scoring."
            )
            std = 1

        # Normalise
        self.primer_pair_costs = self.weight * (self.primer_pair_values - mu) / std

        return self
