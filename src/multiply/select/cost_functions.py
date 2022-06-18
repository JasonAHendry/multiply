from abc import ABC, abstractmethod


# ================================================================================
# Abstract base class for multiplex cost function
#
# ================================================================================

class CostFunction(ABC):
    def __init__(self, indv_costs, pairwise_costs):
        """
        Compute the cost of a multiplex, given a set of individual
        primer pair costs `indv_costs` and pairwise costs `pairwise_costs`
        
        Again, really need to be sure of ordering.
        
        """
        self.indv_costs = indv_costs
        self.pairwise_costs = pairwise_costs
        
        self.indv_combined = None
        self.pairwise_combined = None
        
        
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
        ORDER IS ASSUMED! 
        
        
        """
        self.indv_combined = sum([i.primer_pair_costs for i in self.indv_costs])
        
    def _combine_pairwise_costs(self):
        """
        ORDER IS ASSUMED!
        
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
    
            

