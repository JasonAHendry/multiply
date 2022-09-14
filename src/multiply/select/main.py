import pandas as pd
import numpy as np
from multiply.util.dirs import produce_dir
from multiply.util.printing import print_header, print_footer
from multiply.util.plot import visualise_pairwise_costs
from .cost.factories import IndividualCostFactory, PairwiseCostFactory
from .cost.functions import LinearCost
from .selectors import selector_collection
from .explore import MultiplexExplorer
from .plot import plot_explorer_costs


# PARAMETERS
INDV_INI_PATH = "settings/select/individual_costs.ini"
PAIR_INI_PATH = "settings/select/pairwise_costs.ini"
N_SELECT = 3


def select(result_dir, algorithm):
    """
    Select optimal multiplex(es) from a set of candidate primers

    It is assumed that the following commands have been run:
    `multiply generate`, `multiply align`, `multiply blast`.

    Information about primer quality, primer dimers, and off-target
    binding sites are fed into a cost function, which is then minimised.

    """
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
        indv_cost.collapse_to_per_pair().normalise_costs()
        for indv_cost in indv_factory.get_individual_costs()
    ]
    print(f"  Individual costs: {', '.join([i.cost_name for i in indv_costs])}")

    # CREATE PAIRWISE COSTS
    pairwise_factory = PairwiseCostFactory(PAIR_INI_PATH, result_dir)
    pairwise_costs = [
        pair_cost.collapse_to_per_pair().normalise_costs()
        for pair_cost in pairwise_factory.get_pairwise_costs()
    ]
    print(f"  Pairwise costs: {', '.join([i.cost_name for i in pairwise_costs])}")

    # SET COST FUNCTION
    print("Building cost function...")
    cost_function = LinearCost(indv_costs=indv_costs, pairwise_costs=pairwise_costs)
    cost_function.combine_costs()
    print("Done.\n")

    # RUN SELECTION ALGORITHM
    print("Searching for optimal multiplex(es)...")
    print(f"  Using algorithm: {algorithm}")
    # Run
    selector = selector_collection[algorithm](primer_df, cost_function)
    multiplexes = selector.run()

    # Prepare to explore results
    explorer = MultiplexExplorer(primer_df, multiplexes)
    explorer.set_top_multiplexes(top_N=N_SELECT)

    # Benchmark with random algorithm
    print("Benchmarking search performance with random search...")
    rnd_selector = selector_collection["Random"](primer_df, cost_function)
    rnd_multiplexes = rnd_selector.run(N=len(explorer.uniq_multiplexes))

    # Prepare to explore random results
    rnd_explorer = MultiplexExplorer(primer_df, rnd_multiplexes)
    rnd_explorer.set_top_multiplexes(top_N=N_SELECT)

    # EXPLORING RESULTS
    print("Exploring search results...")
    algo_costs = [
        m.cost for m in explorer.uniq_multiplexes
    ]  # NB: only looking at unique
    rnd_costs = [m.cost for m in rnd_explorer.uniq_multiplexes]
    print(f"  {'Algorithm':>10}  {'Mean Cost':>10}  {'Lowest Cost':>10}")
    print(
        f"  {algorithm:>10}  {np.mean(algo_costs):>10.3f}  {np.min(algo_costs):>10.3f}"
    )
    print(
        f"  {'Random':>10}  {np.mean(rnd_costs):>10.3f}  {np.min(rnd_costs):>10.3f}\n"
    )
    plot_explorer_costs(
        algorithm,
        algo_costs,
        rnd_costs,
        output_path=f"{output_dir}/search.diagnostics.png",
    )

    # WRITE OUTPUTS
    print("Writing outputs...")
    info_path = f"{output_dir}/table.multiplexes.information.csv"
    order_path = f"{output_dir}/table.multiplexes.order.csv"
    explorer.get_union_dataframe(f"{output_dir}/table.multiplexes_overview.csv")
    explorer.get_order_dataframe(f"{output_dir}/table.multiplexes_order.csv")
    print(f"  Information about primers in multiplexes: {info_path}")
    print(f"  Primer ordering table: {order_path}")
    print("Done.\n")

    # PLOTTING
    # - This is just a sketch, needs cleaning / encapsulation
    alignments_df = pd.read_csv(f"{result_dir}/align/table.alignment_scores.csv")
    for ix, multiplex in enumerate(explorer.top_multiplexes):

        # Make output directory
        multiplex_name = f"multiplex{ix:02d}"
        multiplex_output_dir = produce_dir(output_dir, multiplex_name)
        primer_names = multiplex.get_primer_names()

        # Write multiplex df
        multiplex_df = primer_df.loc[primer_names]
        multiplex_df.to_csv(
            f"{multiplex_output_dir}/table.{multiplex_name}_overview.csv", index=False
        )

        # Get pairwise dataframe
        pairwise_df = pairwise_costs[0].primer_values.loc[primer_names, primer_names]
        visualise_pairwise_costs(
            pairwise_df,
            cbar_title="Interaction Score",
            output_path=f"{multiplex_output_dir}/plot.{multiplex_name}.primer_interactions.pdf",
        )
        multiplex_align_df = alignments_df.query(
            "primer1_name in @primer_names and primer2_name in @primer_names"
        )

        with open(f"{multiplex_output_dir}/diagrams.{multiplex_name}.txt", "w") as fn:
            for _, row in multiplex_align_df.iterrows():
                fn.write(f"Overall rank: {row['rank']}\n")
                fn.write(row["alignment"])

    print_footer(t0)
