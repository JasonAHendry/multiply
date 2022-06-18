import numpy as np
import pandas as pd

# TODO:
# - Remove `primer_names`, but be more strict about `primer_values` data
#  - e.g. it should be a pandas series or data frame, where the index is the primer name



class IndividualCosts:
    def __init__(self, cost_name, primer_names, primer_values, weight):
        """
        Encapsulate primer cost information for *individual* primers
        
        The class handles collapsing per primer costs (e.g. CRT1_v1_F) to
        per primer pair costs (e.g. CRT1_v1), and then normalising and
        weighting the collapsed values.
        
        TODO:
        - Careful defensive programming -- this is a key class
        
        params
            cost_name: str
                Name of the primer costs.
            primer_names: list, str, shape(n_primers,)
                Names of the individual primers (e.g. CRT1_v1_F).
            primer_values: list, float, shape(n_primers,)
                Input values associated with each primer.
            weight: float
                Value by which to multiply resulting costs.
        
        """
        
        self.cost_name = cost_name
        self.primer_names = primer_names
        self.primer_values = primer_values
        self.weight = weight
        
        self.primer_pair_values = None
        self.primer_pair_costs = None
    
    def collapse_to_per_pair(self, collapse_func=sum):
        """
        Collapse per primer values to per primer pair values,
        by aggregating with a collapsing fuction,
        `collapse_func`
        
        """
        
        # Build into data frame
        pair_df = pd.DataFrame({
            "primer_pair_name": [n[:-2] for n in self.primer_names],
            "primer_pair_values": self.primer_values
        })
        
        # Collapse
        self.primer_pair_values = (
            pair_df
            .groupby("primer_pair_name")
            .primer_pair_values
            .apply(collapse_func)
        )
        
        return self
    
    def normalise_costs(self):
        """
        Z-score standardise, and then multiply by weight
        
        """
        
        if self.primer_pair_values is None:
            raise ValueError("Run `.collapse_to_per_pair()` first.")
            
        # Normalise
        mu = self.primer_pair_values.mean()
        std = self.primer_pair_values.std()
        self.primer_pair_costs = self.weight * (self.primer_pair_values - mu) / std
        
        return self
        
    def __repr__(self):
        return f"IndividualCosts(cost_name={self.cost_name}, "\
        f"weight={self.weight}, "\
        f"primer_names={','.join(self.primer_names.values[:3])}..., "\
        f"primer_values={','.join([str(v) for v in self.primer_values.values[:3]])}...)"



class PairwiseCosts:
    def __init__(self, cost_name, primer_names, primer_values, weight):
        """
        Create a pairwise costs
        
        """
                
        self.cost_name = cost_name
        self.primer_names = primer_names
        self.primer_values = primer_values
        self.weight = weight
        
        self.primer_pair_values = None
        self.primer_pair_costs = None
    
    def collapse_to_per_pair(self, collapse_func=sum):
        """
        Collapse to per pair
        
        Form of collapse func has some restrictions here
        
        """
        
        primer_pair_name = [n[:-2] for n in self.primer_names]

        
        # Add pair name column
        collapse_df = self.primer_values.copy()
        collapse_df.insert(0, "primer_pair_name", primer_pair_name)
        
        # Collapse first axis
        self.primer_pair_values = (
            collapse_df
            .groupby("primer_pair_name")
            .aggregate(collapse_func)
            .transpose()
        )
        
        # Repeat with other axis
        self.primer_pair_values.insert(0, "primer_pair_name", primer_pair_name)
        self.primer_pair_values = (
            self.primer_pair_values
            .groupby("primer_pair_name")
            .aggregate(collapse_func)
        )
        
        return self
    
    def normalise_costs(self):
        """
        Iclude re-weighting here, for simplicity
        
        """
        if self.primer_pair_values is None:
            raise ValueError("Run `.collapse_to_per_pair()` first.")
        
        # Normalise
        arr = np.array(self.primer_pair_values)
        mu = arr.mean()
        std = arr.std()
        self.primer_pair_costs = self.weight * (self.primer_pair_values - mu) / std
