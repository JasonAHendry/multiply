import numpy as np
import matplotlib.pyplot as plt


def plot_explorer_costs(algorithm, algo_costs, rnd_costs, output_path=None):
    """
    Plot the cost distributions of the user selected algorithm,
    versus a random search algorithm
    
    """

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    
    # Distribution
    ax.hist([algo_costs, rnd_costs], color=["skyblue", "orange"],
            stacked=False, bins=40, label=[algorithm, "Random"])
    ax.legend()

    # Summary statistics
    ax.axvline(np.mean(algo_costs), color='skyblue', ls='dashed')
    ax.axvline(np.mean(rnd_costs), color='orange', ls='dashed')

    # Labels
    title = "Cost Distribution of Unique Multiplexes\n"
    title += f"{algorithm} lowest score: {algo_costs[0]:.02f}\n"
    title += f"Random lowest cost: {rnd_costs[0]:.02f}"
    ax.set_title(title, loc="left")
    ax.set_xlabel("Multiplex Cost\n[Lower = Better]")
    ax.set_ylabel("No. of Unique Candidate Multiplexes")
    
    if output_path is not None:
        fig.savefig(output_path, pad_inches=0.5, bbox_inches="tight")
        plt.close(fig)