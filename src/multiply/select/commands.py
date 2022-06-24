import click
import pandas as pd
from multiply.util.dirs import produce_dir
from multiply.util.printing import print_header, print_footer
from .cost.factories import IndividualCostFactory, PairwiseCostFactory
from .cost.functions import LinearCost
from .selectors import GreedySearch
from .explore import MultiplexExplorer

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
    print(f"  Individual costs: {', '.join([i.cost_name for i in indv_costs])}")

    # CREATE PAIRWISE COSTS
    pairwise_factory = PairwiseCostFactory(PAIR_INI_PATH, result_dir)
    pairwise_costs = [
        pair_cost
        .collapse_to_per_pair()
        .normalise_costs()
        for pair_cost in pairwise_factory.get_pairwise_costs()
    ]
    print(f"  Pairwise costs: {', '.join([i.cost_name for i in pairwise_costs])}")

    # SET COST FUNCTION
    print("Building cost function...")
    cost_function = LinearCost(indv_costs=indv_costs, pairwise_costs=pairwise_costs)
    cost_function.combine_costs()
    print("Done.\n")

    # SET SELECTION ALGORITHM
    print(f"Seaching for optimal multiplexes...")
    print(f" Algorithm: Greedy")
    selector = GreedySearch(primer_df, cost_function)
    multiplexes = selector.run()

    # EXPLORE OUTPUTS
    print("Exploring resulting multiplexes...")
    explorer = MultiplexExplorer(primer_df, multiplexes)
    explorer.set_top_multiplexes(top_N=N_SELECT)
    print(f"  Total multiplexes generated: {len(explorer.multiplexes)}")
    print(f"  No. unique: {len(explorer.uniq_multiplexes)}")
    print(f"  Top {len(explorer.top_multiplexes)} retained.")

    # WRITE OUTPUTS
    print("Writing outputs...")
    info_path = f"{output_dir}/table.multiplexes.information.csv"
    order_path = f"{output_dir}/table.multiplexes.order.csv"
    explorer.get_union_dataframe(f"{output_dir}/table.multiplexes_overview.csv")
    explorer.get_order_dataframe(f"{output_dir}/table.multiplexes_order.csv")
    print(f"  Information about primers in multiplexes: {info_path}")
    print(f"  Primer ordering table: {order_path}")
    print("Done.\n")

    print_footer(t0)
