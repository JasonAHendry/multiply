import pandas as pd
from abc import ABC, abstractmethod


# ================================================================================
# Abstract implementation of a cost function
#
# ================================================================================


class CostFunction(ABC):
    def __init__(self, indv_costs, pairwise_costs):
        """
        Compute the cost of a multiplex, given a set of individual
        primer pair costs `indv_costs` and pairwise costs `pairwise_costs`

        Again, really need to be sure of ordering.

        """
        # Set an instantiation
        self.indv_costs = indv_costs
        self.pairwise_costs = pairwise_costs
        
        # Ordering and consistency
        self._check_cost_consistency()
        
        # Computed
        self.indv_combined = None
        self.pairwise_combined = None
        
        
    def _check_cost_consistency(self):
        """
        Check that the primer pairs for all of the individual
        and pairwise costs are consistent
        
        This ensures they are (1) the same and (2) in the
        same order.
        
        """
        
        # Benchmark against first example
        bench = self.indv_costs[0]
        self._primer_pairs = bench.primer_pair_costs.index.tolist()
        
        # Ensure the primer pairs are unique
        assert len(self._primer_pairs) == len(set(self._primer_pairs)), \
        "Primer pairs must be unique across all `indv_costs` and `pairwise_costs`."
        
        # Ensure primers the same for individual costs
        for indv_cost in self.indv_costs:
            assert indv_cost.primer_pair_costs.index.tolist() == self._primer_pairs, \
            f"Primer pairs must be the same across all `indv_costs`. Different primers in {indv_cost.cost_name}."
            
        # Ensure primers the same for pairwise costs
        for pair_cost in self.pairwise_costs:
            assert pair_cost.primer_pair_costs.index.tolist() == self._primer_pairs, \
            f"Index in {pair_cost.cost_name} has inconsistent primer pairs."
            assert pair_cost.primer_pair_costs.columns.tolist() == self._primer_pairs, \
            f"Columns in {pair_cost.cost_name} has inconsistent primer pairs."
            
    def combine_costs(self):
        """
        Combine individual and pariwise costs to facilitate evaluation

        NB:
        - Can sum across columns immediately

        """

        self._combine_individual_costs()
        self._combine_pairwise_costs()

    def _combine_individual_costs(self):
        """ 
        Combine individual costs 
        
        Needs to change, like this it defeats the purpose of flexiblity
        in the cost function. Also note that have the same order is assumed,
        although we check this with self._check_cost_consistency().
        
        """
        self.indv_combined = sum([i.primer_pair_costs for i in self.indv_costs])

    def _combine_pairwise_costs(self):
        """
        Combine pairwise costs
        
        Same comment as above

        """
        self.pairwise_combined = sum([p.primer_pair_costs for p in self.pairwise_costs])

    @abstractmethod
    def calc_cost(self, primer_pairs):
        """
        Calculate the cost of a set of primer pairs

        """
        pass


# ================================================================================
# Concrete cost functions
#
# ================================================================================


class LinearCost(CostFunction):
    def calc_cost(self, primer_pairs):
        indv = self.indv_combined.loc[primer_pairs].sum()
        pairwise = self.pairwise_combined.loc[primer_pairs][primer_pairs].sum().sum()
        return indv + pairwise

