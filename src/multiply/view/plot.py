from multiply.util.statistics import (
    get_homopolymer_runs,
    calc_sliding_percentGC,
    get_array_encoding,
)
from collections import namedtuple
from functools import partial
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


class SequencePlotter:
    def __init__(self, seq):
        """
        Plot various summary statistics of a nucleotide sequennce

        """
        self.seq = seq
        self._calc_summaries()

    def _calc_summaries(self):
        """
        Calculate summary statistics of the sequence

        """

        self.hp_runs = get_homopolymer_runs(self.seq)
        self.per_gc = calc_sliding_percentGC(self.seq, window=20)
        self.seq_array = get_array_encoding(self.seq)

    def plot_sequence_array(self, ax, start, end):
        """
        Plot array giving nucleotide composition

        """
        # Plot
        ax.imshow(
            self.seq_array, cmap="Blues", aspect="auto", extent=(start, end, 3.5, -0.5)
        )

        # Ticks
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels(["A", "T", "C", "G"])

        # Labels
        ax.set_xlabel("Position (bp)")
        #ax.label_outer()

        return None

    def plot_sequence_complexity(self, ax, start, end):
        """
        Plot homopolymer run length and AT (%) on a single
        axis

        """

        HP_COL = "firebrick"
        GC_COL = "forestgreen"

        # Define x values
        xs = np.arange(start, end)

        # Homopolymers
        ax.plot(xs, self.hp_runs, lw=1, color=HP_COL, label="Homopolymer Length (bp)")
        ax.set_ylabel("Homopolymer \nLength (bp)", color=HP_COL)

        # Limits
        ax.set_ylim((0, ax.get_ylim()[1]))
        ax.set_xlim(start, end)

        # GC
        ax.patch.set_visible(False)
        axm = ax.twinx()
        axm.set_zorder(ax.get_zorder() - 1)
        axm.fill_between(
            x=xs,
            y1=0,
            y2=100 * (1 - self.per_gc),
            alpha=0.5,
            color=GC_COL,
            label="% AT",
        )
        axm.set_ylim(ax.get_ylim()[0], 100)
        axm.set_ylabel("AT (%)\n[20bp sliding average]", color=GC_COL)

        # Clean axis
        axm.xaxis.set_visible(False)
        plt.setp(ax.get_xticklabels(), visible=False)


class GffPlotter:

    gff_features = ["protein_coding_gene", "CDS"]

    def __init__(self, gff, chrom, start, end):
        """
        Plot information from a GFF over a defined
        region

        """
        self.gff = gff
        self.chrom = chrom
        self.start = start
        self.end = end

        self.plot_gff = self._get_relevant_gff()

    def _get_relevant_gff(self):
        """Get the relevant proportion of the .gff for plotting"""

        qry = f"(seqname == '{self.chrom}')"
        qry += f" and ({self.start} <= start <= {self.end}"
        qry += f" or {self.start} <= end <= {self.end}"
        qry += f" or (start <= {self.start} and {self.end} <= end))"
        self._qry = qry

        # Filter to rows for plotting
        plot_gff = self.gff.query(qry).query("feature in @self.gff_features")

        # Ensure there are some regions
        assert plot_gff.shape[0] > 0, "No features in this region."

        return plot_gff

    def plot_gff_features(self, ax):
        """Plot features of gff in this region"""
        
        PLUS_STRAND_Y = 3/4
        NEG_STRAND_Y = 1/4

        # Plot features
        for _, row in self.plot_gff.iterrows():

            # Define color and size from feature
            if row["feature"] == "protein_coding_gene":
                lw = 3
                color = "darkgrey"
            elif row["feature"] == "CDS":
                lw = 8
                color = "teal"

            # Define y position from strand
            if row["strand"] == "+":
                ypos = PLUS_STRAND_Y
            elif row["strand"] == "-":
                ypos = NEG_STRAND_Y

            # Plot the feature
            ax.plot([row["start"], row["end"]], [ypos, ypos], lw=lw, color=color)

            # Could add annotation text
            
        # Indicate strands themselves
        ax.plot([self.start, self.end], 
                [PLUS_STRAND_Y, PLUS_STRAND_Y],
                lw=1, color='lightgrey', zorder=-10)
        ax.annotate(xy=(self.end, PLUS_STRAND_Y), 
                    xycoords="data",
                    ha="left", va="center",
                    fontsize=8,
                    text=" ($+$) strand")
        
        # Indicate strands themselves
        ax.plot([self.start, self.end], 
                [NEG_STRAND_Y, NEG_STRAND_Y],
                lw=1, color='lightgrey', zorder=-10)
        ax.annotate(xy=(self.end, NEG_STRAND_Y), 
                    xycoords="data",
                    ha="left", va="center",
                    fontsize=8,
                    text=" ($-$) strand")

        # Clean axis ticks
        ax.axis("off")
        ax.set_xlim((self.start, self.end))
        ax.set_ylim((0, 1))
        ax.label_outer()

        return None


class PrimerPlotter:
    def __init__(self, primer_df):

        # Set primer data frame
        self.primer_df = primer_df

        # Group into pairs, get colors
        self._get_primer_pair_groups()
        self._set_colors()

    def _get_primer_pair_groups(self):
        self.grps = self.primer_df.groupby("pair_name")
        self.n_grps = len(self.grps)
        self.group_names = list(self.grps.groups.keys())

    def _set_colors(self, cmap="plasma"):
        self.col_dt = dict(zip(self.group_names, sns.color_palette(cmap, self.n_grps)))

    def plot(self, ax, start, end):
        """
        Plot, scaled

        """

        for ix, (pair_name, pair_df) in enumerate(self.grps):

            # Extract information
            F_row = pair_df.query("direction == 'F'")
            R_row = pair_df.query("direction == 'R'")

            F_start = F_row["start"]
            F_end = F_start + F_row["length"]

            R_start = R_row["start"]
            R_end = R_start - R_row["length"]

            # Forward
            ax.plot([F_start, F_end], [ix, ix], lw=4, color=self.col_dt[pair_name])
            ax.scatter(x=F_end, y=ix, marker=9, s=30, color=self.col_dt[pair_name])

            # Reverse
            ax.plot([R_start, R_end], [ix, ix], lw=4, color=self.col_dt[pair_name])
            ax.scatter(x=R_end, y=ix, marker=8, s=30, color=self.col_dt[pair_name])

            # Amplicon
            ax.plot([F_start, R_start], [ix, ix], color=self.col_dt[pair_name])

            # Annotate
            ax.annotate(
                xy=(R_start, ix),
                ha="left",
                va="center",
                fontsize=6,
                text=f"   {pair_name} {F_row['product_bp'].values[0]}bp",
            )
            # color=col_dt[pair_name])

            # Clean ticks
            ax.set_xlim(start, end)
            for s in ["top", "right", "bottom", "left"]:
                ax.spines[s].set_visible(False)
            ax.get_yaxis().set_ticks([])
            ax.label_outer()


class CombinedPlotter:
    """Runs, but super ugly"""

    def __init__(self, sequence_plotter, gff_plotter, primer_plotter):
        self.seq_plotter = sequence_plotter
        self.gff_plotter = gff_plotter
        self.primer_plotter = primer_plotter

    def plot(self, start, end, title=None, output_path=None):
        """
        Create a combined plot of sequence composition, complexity,
        gene locations, and candidate primers
        
        """

        # Constants
        SCALING = 0.08
        ROWS_PER_PRIMER = 2
        N_PRIMERS = self.primer_plotter.n_grps

        # Axes order and sizing
        AxesFrame = namedtuple("AxesFrame", ["name", "rows", "plot_func"])
        axes_order = [
            AxesFrame("genes", 6, self.gff_plotter.plot_gff_features),
            AxesFrame("complexity", 32, partial(self.seq_plotter.plot_sequence_complexity, start=start, end=end)),
            AxesFrame("primers", N_PRIMERS*ROWS_PER_PRIMER, partial(self.primer_plotter.plot, start=start, end=end)),
            AxesFrame("seq", 10, partial(self.seq_plotter.plot_sequence_array, start=start, end=end))
        ]

        # Figure size
        total_rows = sum([a.rows for a in axes_order]) + len(axes_order)
        height = total_rows * SCALING
        width = 10

        # Create figure
        fig = plt.figure(figsize=(width, height))
        fig.subplots_adjust(hspace=0.2)

        # Create grid
        gs = GridSpec(nrows=total_rows, ncols=1)

        # Iterate over axes and plot
        l = 0
        for i, axis_item in enumerate(axes_order):
            
            # Prepare axis
            ax = plt.subplot(gs[l:(l+axis_item.rows)])
            
            # Plot
            axis_item.plot_func(ax)
            
            # Optionally add title
            if i == 0 and title is not None:
                ax.set_title(title, loc="left")
            
            # Define boundary of next plot
            l = l + axis_item.rows + 1

        # Optionally write
        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5)