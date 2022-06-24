import sys
import random
from itertools import product
from functools import reduce
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

    Note, that it would be possible to compute exhaustively
    the number of possible permutations through the greedy algorithm;

    Alternatively; I could create a unique set of orders *instead*
    of shuffling; depending on the ratio of the number of iterations
    to the number of permutations, this would remove some redundant
    calculation

    """

    def run(self, N=10_000):
        """
        Run a greedy search algorithm for the lowest cost multiplex

        """

        # Get every UNIQUE primer pair, for each target
        # NB: from `primer_df` these are doubled, must use set()
        target_pairs = {
            target_id: list(set(target_df["pair_name"]))
            for target_id, target_df in self.primer_df.groupby("target_id")
        }

        # IDs
        target_ids = list(target_pairs)

        # Iterate
        multiplexes = []
        sys.stdout.write(f"  Iterations complete: {0}/{N}")
        for ix in range(N):

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
            sys.stdout.write(f"  Iterations complete: {ix+1}/{N}")
        print("\nDone.\n")

        return multiplexes


class BruteForce(MultiplexSelector):
    def run(self, store_maximum=200):
        """
        Run a brute force search for the highest scoring multiplex

        Note, we control the maximum number of multiplexes stored
        using the `store_maximum` argument; otherwise we can get
        exceptionally long lists of multiplexes

        """
        # Split target pairs into list of sets
        target_pairs = [
            set(target_df["pair_name"])
            for _, target_df in self.primer_df.groupby("target_id")
        ]

        # Compute number of iterations required
        total_N = reduce(lambda a, b: a * b, [len(t) for t in target_pairs])
        print(
            f"Found {int(self.primer_df.shape[0]/2)} primer pairs across {len(target_pairs)} targets."
        )
        print(f"A total of {total_N} possible multiplexes exist.")

        # Iterate over all possible multiplexes
        sys.stdout.write(f"  Iterations complete: {0}/{total_N}")
        stored_multiplexes = []
        stored_costs = []
        for ix, primer_pairs in enumerate(product(*target_pairs)):

            # Create the multiplex
            multiplex = Multiplex(
                primer_pairs=primer_pairs,
                cost=self.cost_function.calc_cost(primer_pairs),
            )

            # Store
            if len(stored_multiplexes) < store_maximum:
                stored_multiplexes.append(multiplex)
                stored_costs.append(multiplex.cost)
                highest_stored_cost = max(stored_costs)
            elif multiplex.cost < highest_stored_cost:
                stored_multiplexes.insert(
                    stored_costs.index(highest_stored_cost), multiplex
                )

            # Print
            sys.stdout.write("\r")
            sys.stdout.flush()
            sys.stdout.write(f"  Iterations complete: {ix+1}/{total_N}")

        print("\nDone.\n")

        return stored_multiplexes


class RandomSearch(MultiplexSelector):
    def run(self, N=10_000):
        """Run the random  selection algorithm"""
        # Get target pairs
        target_pairs = {
            target_id: list(set(target_df["pair_name"]))
            for target_id, target_df in self.primer_df.groupby("target_id")
        }

        # Iterate
        multiplexes = []
        sys.stdout.write(f"  Iterations complete: {0}/{N}")
        for ix in range(N):

            # Randomly generate a multiplex
            multiplex = [random.choice(pairs) for _, pairs in target_pairs.items()]

            # Compute the cost
            cost = self.cost_function.calc_cost(multiplex)

            # Store
            multiplexes.append(Multiplex(cost=cost, primer_pairs=multiplex))

            # Print
            sys.stdout.write("\r")
            sys.stdout.flush()
            sys.stdout.write(f"  Iterations complete: {ix+1}/{N}")
        print("\nDone.\n")

        return multiplexes


# ================================================================================
# Collection of selection algorithms
#
# ================================================================================


selector_collection = {
    "Greedy": GreedySearch,
    "Random": RandomSearch,
    "BruteForce": BruteForce,
}
