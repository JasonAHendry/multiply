import numpy as np
import matplotlib.pyplot as plt


def visualise_pairwise_costs(pairwise_df, cmap="Spectral", cbar_title=None, output_path=None, **kwargs):
    """
    Plot a matrix of pairwise costs
    
    params
        pairwise_df: pandas DataFrame shape(n, m)
            Pandas dataframe of pairwise costs; columns and indices
            are used to label axes.
        cmap: str [optional]
            Color map to use for values.
        cbar_title: str [optional]
            Title to colorbar.
        output_path: str [optional]
            If set, figure is saved at `output_path`.
    
    returns
        None
            
    """
    
    # Prepare size
    RESCALE = 0.6
    c, r = pairwise_df.shape
    
    # Set canvas
    fig, ax = plt.subplots(1, 1, figsize=(c*RESCALE, r*RESCALE))
    
    # Plot
    cax = ax.matshow(pairwise_df, cmap=cmap)

    # Colorbar
    if cbar_title is not None:
        fig.colorbar(
            cax,
            shrink=0.85 if c > 10 else 0.8,
            pad=0.01 if c > 10 else 0.05,
            label=cbar_title,
        )

    # Axis, Ticks, Grid
    ax.xaxis.set_label_position("top")
    ax.set_yticks(range(r))
    ax.set_xticks(range(c))
    ax.set_yticks(np.arange(0.5, r + 0.5), minor=True)
    ax.set_xticks(np.arange(0.5, c + 0.5), minor=True)
    ax.set_yticklabels(pairwise_df.index)
    ax.set_xticklabels(pairwise_df.columns, rotation=90)
    ax.grid(which="minor", ls="dotted")
    
    # Optionally save figure
    if output_path is not None:
        fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5, dpi=300)
        plt.close(fig)