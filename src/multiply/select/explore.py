import pandas as pd


class MultiplexExplorer:
    def __init__(self, primer_df, multiplexes):
        """
        Explore `multiplexes` derived from a given `primer_df`

        params
            primer_df: pandas DataFrame
                Each row is an individual primer.
            multiplexes: List[Multiplex]
                List of multiplex objects

        """

        # Set at instantiation
        self.primer_df = primer_df
        self.multiplexes = multiplexes

        # Compute
        self._get_unique_multiplexes_and_sort()

    def _get_unique_multiplexes_and_sort(self):
        """
        Reduce list of multiplexes to a unique set,
        and sort

        """

        self.uniq_multiplexes = sorted(set(self.multiplexes), key=lambda m: m.cost)

    def _extract_from_df(self, df):
        """
        Extract multiplex rows from a dataframe

        """

        edf = df.copy()
        edf.index = df["primer_name"]

        # Extract associated dataframes
        multiplex_dfs = [edf.loc[m.get_primer_names()] for m in self.top_multiplexes]

        return multiplex_dfs

    def set_top_multiplexes(self, top_N=3):
        """
        Create a smaller list of the `top_N` multiplexes

        """

        # Limit to `top_N` multiplexes
        self.top_N = top_N
        self.top_multiplexes = self.uniq_multiplexes[:top_N]

    def get_union_dataframe(self, output_path=None):
        """
        Write a union dataframe giving information about union of
        primers across all multiplexes

        """

        # Extract multiplex dataframes
        dfs = self._extract_from_df(df=self.primer_df)

        # Combine
        union_df = (
            pd.concat(dfs)
            .reset_index(drop=True)
            .fillna(False)
            .drop_duplicates("primer_name")
            .sort_values("primer_name")
        )

        # Add inclusion column
        for ix in range(self.top_N):

            in_multiplex = [
                p in self.top_multiplexes[ix].primer_pairs
                for p in union_df["pair_name"]
            ]

            union_df.insert(ix + 4, f"in_multiplex{ix:02d}", in_multiplex)

        # Store
        self.union_df = union_df

        # Optionally write
        if output_path is not None:
            self.union_df.to_csv(output_path, index=False)

    def get_order_dataframe(self, output_path=None):
        """Write dataframe for ordering"""
        assert self.union_df is not None, "Must run `.get_union_dataframe()` first."

        # Create smaller datafram with columns for IDT order of primers
        order_df = self.union_df[["primer_name", "seq"]]
        order_df.insert(2, "concentration", "25nm")
        order_df.insert(3, "purification", "STD")

        # Store
        self.order_df = order_df

        if output_path is not None:
            self.order_df.to_csv(output_path, index=False)
