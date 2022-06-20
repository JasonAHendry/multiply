import click
import pandas as pd
from multiply.util.dirs import produce_dir
from .costs import IndividualCosts, PairwiseCosts
from .cost_functions import LinearCost
from .selectors import GreedySearch


@click.command(short_help="Select optimal multiplex(es).")
@click.option(
    "-r",
    "--result_dir",
    type=click.Path(exists=True),
    required=True,
    help="Path to results directory for multiplex design (e.g. `results/2022-06-11_pf-default`).",
)
def select(result_dir):
    """
    Select optimal multiplex(es) by minimising a cost function incorporating
    information about primer quality (`multiply generate`), primer dimers (`multiply align`),
    off-target amplicons (`multiply blast`)

    """

    # PARSE CLI
    output_dir = produce_dir(result_dir, "select")

    # LOAD DATA
    # Primer data
    primer_df = pd.read_csv(f"{result_dir}/table.candidate_primers.csv")
    primer_df.index = primer_df["primer_name"]

    # Blast data
    blast_df = pd.read_csv(
        f"{result_dir}/blast/table.blast.candidate_primers.summary.csv"
    )

    # Primer-dimer data
    align_df = pd.read_csv(
        f"{result_dir}/align/matrix.pairwise_scores.csv", index_col=0
    )

    # CONSTRUCT INDIVIDUAL COSTS
    indv_costs = [
        IndividualCosts(
            "blast", blast_df["primer_name"], blast_df["predicted_bound"], weight=1.0
        ),
        IndividualCosts(
            "primer3", primer_df["primer_name"], primer_df["pair_penalty"], weight=1.0
        ),
    ]
    # Normalise
    for indv_cost in indv_costs:
        indv_cost.collapse_to_per_pair(collapse_func=sum)
        indv_cost.normalise_costs()

    # CONSTRUCT PAIRWISE COSTS
    # NB: Weight is negative, because we want to *maximise* the dimer score
    pairwise_costs = [
        PairwiseCosts("dimer", align_df.columns, align_df, weight=-1.0),
    ]
    # Normalise
    for pairwise_cost in pairwise_costs:
        pairwise_cost.collapse_to_per_pair(collapse_func=sum)
        pairwise_cost.normalise_costs()

    # SET COST FUNCTION
    cost_function = LinearCost(indv_costs=indv_costs, pairwise_costs=pairwise_costs)
    cost_function.combine_costs()

    # SET SELECTION ALGORITHM
    selector = GreedySearch(primer_df, cost_function)
    multiplexes = selector.run()

    # FORMAT OUTPUTSs
    # Reduce to unique and sort
    final_multiplexes = sorted(set(multiplexes))

    # Create a joint dataframe, and save
    dfs = []
    final_columns = primer_df.columns.tolist()
    for ix in range(3):
        column_name = f"in_multiplex{ix+1:02d}"
        df = primer_df.loc[final_multiplexes[ix].get_primer_names()]
        df.insert(df.shape[1], column_name, True)
        dfs.append(df)
        final_columns.insert(4 + ix, column_name)
    top_multiplex_df = (
        pd.concat(dfs)
        .reset_index(drop=True)
        .fillna(False)
        .drop_duplicates("primer_name")
        .sort_values("primer_name")
    )[final_columns]
    top_multiplex_df.to_csv(f"{output_dir}/table.multiplexes_overview.csv")

    # Sketch of OO for future
    # output_formatter = MultiplexFormatter(primer_df, multiplexes)
    # output_formatter.set_top_multiplexes(N=3)
    # output_formatter.write_joint_dataframe(f"{output_dir}/select/table.multiplexes_overview.csv")
    # output_formatter.write_order_dataframe(f"{output_dir}/select/table.multiplexes_order.csv")
