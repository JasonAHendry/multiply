import click
import pandas as pd
from multiply.util.dirs import produce_dir
from multiply.util.printing import print_header, print_footer
from .cost.factories import IndividualCostFactory, PairwiseCostFactory
from .cost.functions import LinearCost
from .selectors import GreedySearch


# PARAMETERS
INDV_INI_PATH = "settings/select/individual_costs.ini"
PAIR_INI_PATH = "settings/select/pairwise_costs.ini"
N_SELECT = 3


# ================================================================================
# Main function wrapped for Click CLI
#
# ================================================================================


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
    Select optimal multiplex(es) from a set of candidate primers

    It is assumed that the following commands have been run:
    `multiply generate`, `multiply align`, `multiply blast`.

    Information about primer quality, primer dimers, and off-target
    binding sites are fed into a cost function, which is then minimised.

    """
    main(result_dir)


# ================================================================================
# Main function, unwrapped
#
# ================================================================================


def main(result_dir):
    # PARSE CLI
    t0 = print_header()
    print("Parsing inputs...")
    output_dir = produce_dir(result_dir, "select")
    primer_df = pd.read_csv(f"{result_dir}/table.candidate_primers.csv")
    primer_df.index = primer_df["primer_name"]
    print(f"  Results directory: {result_dir}")
    print(f"  Primer CSV: {result_dir}/table.candidate_primers.csv")
    print(f"  Output directory: {output_dir}")
    print("Done.\n")

    # CREATE INDIVIDUAL COSTS
    print("Preparing inputs to cost function...")
    indv_factory = IndividualCostFactory(INDV_INI_PATH, result_dir)
    indv_costs = [
        indv_cost
        .collapse_to_per_pair()
        .normalise_costs()
        for indv_cost in indv_factory.get_individual_costs()
    ]
    print(f"  Individual costs: {' ,'.join([i.cost_name for i in indv_costs])}")

    # CREATE PAIRWISE COSTS
    pairwise_factory = PairwiseCostFactory(PAIR_INI_PATH, result_dir)
    pairwise_costs = [
        pair_cost
        .collapse_to_per_pair()
        .normalise_costs()
        for pair_cost in pairwise_factory.get_pairwise_costs()
    ]
    print(f"  Pairwise costs: {' ,'.join([i.cost_name for i in indv_costs])}")

    # SET COST FUNCTION
    print("Building cost function...")
    cost_function = LinearCost(indv_costs=indv_costs, pairwise_costs=pairwise_costs)
    cost_function.combine_costs()
    print("Done.\n")

    # SET SELECTION ALGORITHM
    print(f"Seaching for top {N_SELECT} optimal multiplexes...")
    print(f"  Search algorithm: Greedy")
    selector = GreedySearch(primer_df, cost_function)
    multiplexes = selector.run()

    # FORMAT OUTPUTS
    print("Writing output...")
    # Reduce to unique and sort
    final_multiplexes = sorted(set(multiplexes))

    # Create a summary data frame
    dfs = [
        primer_df.loc[final_multiplexes[ix].get_primer_names()]
        for ix in range(N_SELECT)
    ]
    top_multiplexes = (
        pd.concat(dfs)
        .reset_index(drop=True)
        .fillna(False)
        .drop_duplicates("primer_name")
        .sort_values("primer_name")
    )
    for ix in range(N_SELECT):
        top_multiplexes.insert(
            ix+4, 
            f"in_multiplex{ix:02d}", 
            [p in final_multiplexes[ix].primer_pairs 
            for p in top_multiplexes["pair_name"]])
    overview_csv = f"{output_dir}/table.multiplexes_overview.csv"
    top_multiplexes.to_csv(overview_csv)
    print(f"  to: {overview_csv}")
    print("Done.\n")

    print_footer(t0)

    # Sketch of OO for future
    # output_formatter = MultiplexFormatter(primer_df, multiplexes)
    # output_formatter.set_top_multiplexes(N=3)
    # output_formatter.write_joint_dataframe(f"{output_dir}/select/table.multiplexes_overview.csv")
    # output_formatter.write_order_dataframe(f"{output_dir}/select/table.multiplexes_order.csv")
