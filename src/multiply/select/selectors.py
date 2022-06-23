import sys
import random
from abc import ABC, abstractmethod
from .multiplex import Multiplex


# ================================================================================
# Abstract class for selection algorithm
#
# ================================================================================


class MultiplexSelector(ABC):
    def __init__(self, primer_df, cost_function):
        self.primer_df = primer_df
        self.cost_function = cost_function

    @abstractmethod
    def run(self):
        """
        Run the selection method

        """
        pass


# ================================================================================
# Concrete selection algorithms
#
# ================================================================================


class GreedySearch(MultiplexSelector):
    """
    Try to find the optimal multiplex using a greedy search algorithm

    """

    N = 50

    def run(self):
        """
        This runs but it seems to be insanely slow, even for an 8plex

        """

        # Get target pairs
        target_pairs = {
            target_id: target_df["pair_name"].tolist()
            for target_id, target_df in self.primer_df.groupby("target_id")
        }

        # IDs
        target_ids = list(target_pairs)

        # Iterate
        multiplexes = []
        sys.stdout.write(f"  Iterations complete: {0}/{self.N}")
        for ix in range(self.N):

            # Prepare empty new multiplex
            multiplex = []

            # Shuffle the target IDs in place
            random.shuffle(target_ids)

            # Compute scores of each possible pair
            for target_id in target_ids:
                costs = [
                    self.cost_function.calc_cost(multiplex + [primer_pair])
                    for primer_pair in target_pairs[target_id]
                ]

                # Add max scoring from this step
                idxmax = costs.index(min(costs))
                multiplex.append(target_pairs[target_id][idxmax])

            # Add to list of all multiplexes
            multiplexes.append(Multiplex(cost=min(costs), primer_pairs=multiplex))

            # Print
            sys.stdout.write("\r")
            sys.stdout.flush()
            sys.stdout.write(f"  Iterations complete: {ix+1}/{self.N}")
        print("\nDone.\n")

        return multiplexes


class BruteForce(MultiplexSelector):
    pass


class RandomSearch(MultiplexSelector):

    N = 50

    def run(self):
        """Run the random  selection algorithm"""
        # Get target pairs
        target_pairs = {
            target_id: target_df["pair_name"].tolist()
            for target_id, target_df in self.primer_df.groupby("target_id")
        }

        # Iterate
        multiplexes = []
        sys.stdout.write(f"  Iterations complete: {0}/{self.N}")
        for ix in range(self.N):

            # Randomly generate a multiplex
            multiplex = [random.choice(pairs) for _, pairs in target_pairs.items()]

            # Compute the cost
            cost = self.cost_function.calc_cost(multiplex)

            # Store
            multiplexes.append(Multiplex(cost=cost, primer_pairs=multiplex))

            # Print
            sys.stdout.write("\r")
            sys.stdout.flush()
            sys.stdout.write(f"  Iterations complete: {ix+1}/{self.N}")
        print("\nDone.\n")

        return multiplexes
