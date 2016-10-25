###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
from collections import defaultdict, namedtuple

from phylorank.rel_dist import RelativeDistance
from phylorank.common import read_taxa_file, filter_taxa_for_dist_inference

from biolib.taxonomy import Taxonomy
from biolib.plots.abstract_plot import AbstractPlot

from numpy import (mean as np_mean,
                   std as np_std,
                   median as np_median,
                   arange as np_arange,
                   linspace as np_linspace,
                   percentile as np_percentile)

from scipy.stats import norm

import mpld3

import dendropy


class DistributionPlot(AbstractPlot):
    """Plot distribution of taxa in each taxonomic rank."""

    def __init__(self):
        """Initialize."""
        Options = namedtuple('Options', 'width height font_size dpi')
        options = Options(6, 6, 12, 96)

        AbstractPlot.__init__(self, options)

    def _percent_correct_plot(self, rel_dists, taxa_for_dist_inference, output_prefix):
        """Create plots showing correctly classified taxa for different relative distance values.

        Parameters
        ----------
        rel_dists : d[rank_index][taxon] -> relative divergence
            Relative divergence of taxa at each rank.
        taxa_for_dist_inference : iterable
            Taxa to consider when inferring relative divergence thresholds.
        output_prefix : str
            Prefix for plots.
        """

        print ''
        print '  Relative divergence thresholds (rank, threshold, parent taxa, child taxa):'

        ranks = sorted(rel_dists.keys())
        rel_dist_thresholds = []
        for i in xrange(ranks[0], ranks[-1]):
            parent_rank = i
            child_rank = i + 1

            # determine classification results for relative divergence
            # values between the medians of adjacent taxonomic ranks
            parent_rds = []
            for taxa, rd in rel_dists[parent_rank].iteritems():
                if taxa in taxa_for_dist_inference:
                    parent_rds.append(rd)
            parent_p50 = np_percentile(parent_rds, 50)

            child_rds = []
            for taxa, rd in rel_dists[child_rank].iteritems():
                if taxa in taxa_for_dist_inference:
                    child_rds.append(rd)

            child_p50 = np_percentile(child_rds, 50)

            r = []
            y_parent = []
            y_child = []
            y_mean_corr = []
            for test_r in np_linspace(parent_p50, child_p50, 100):
                parent_cor = float(sum([1 for rd in parent_rds if rd <= test_r])) / len(parent_rds)
                child_cor = float(sum([1 for rd in  child_rds if rd > test_r])) / len(child_rds)

                r.append(test_r)
                y_parent.append(parent_cor)
                y_child.append(child_cor)
                y_mean_corr.append(0.5 * parent_cor + 0.5 * child_cor)

            # create plot of correctly classified taxa
            self.fig.clear()
            self.fig.set_size_inches(6, 6)
            ax = self.fig.add_subplot(111)

            ax.plot(r, y_parent, 'k--', label=Taxonomy.rank_labels[i])
            ax.plot(r, y_child, 'k:', label=Taxonomy.rank_labels[i + 1])
            ax.plot(r, y_mean_corr, 'r-', label='mean')

            legend = ax.legend(loc='upper left')
            legend.draw_frame(False)

            # find maximum of mean correct classification
            max_mean = max(y_mean_corr)
            r_max_values = [r[i] for i, rd in enumerate(y_mean_corr) if rd == max_mean]
            r_max_value = np_mean(r_max_values)  # Note: this will fail if there are multiple local maxima
            print '    %s\t%.3f\t%d\t%d' % (Taxonomy.rank_labels[parent_rank], r_max_value, len(parent_rds), len(child_rds))

            # check that there is a single local maximum
            rd_indices = [i for i, rd in enumerate(y_mean_corr) if rd == max_mean]
            for rd_index in xrange(0, len(rd_indices) - 1):
                if rd_indices[rd_index] != rd_indices[rd_index + 1] - 1:
                    print '[Warning] There are multiple local maxima, so estimated relative divergence threshold will be invalid.'

            rel_dist_thresholds.append(r_max_value)

            y_min, _y_max = ax.get_ylim()
            ax.axvline(x=r_max_value, ymin=0, ymax=1, color='r', ls='--')
            ax.text(r_max_value + 0.001, y_min + 0.01, '%.3f' % r_max_value, horizontalalignment='left')

            ax.set_xlabel('relative distance')
            ax.set_ylabel('% taxa correctly classified')

            self.prettify(ax)

            self.fig.tight_layout(pad=1)
            self.fig.savefig(output_prefix + '.%s_%s.png' % (Taxonomy.rank_labels[parent_rank], Taxonomy.rank_labels[child_rank]), dpi=96)

        print ''

        return rel_dist_thresholds

    def _distribution_plot(self, rel_dists, rel_dist_thresholds, taxa_for_dist_inference, distribution_table, plot_file):
        """Create plot showing the distribution of taxa at each taxonomic rank.

        Parameters
        ----------
        rel_dists: d[rank_index][taxon] -> relative divergence
            Relative divergence of taxa at each rank.
        rel_dist_thresholds: list
            Relative distances cutoffs for defining ranks.
        taxa_for_dist_inference : iterable
            Taxa to considered when inferring distributions.
        distribution_table : str
            Desired name of output table with distribution information.
        plot_file : str
            Desired name of output plot.
        """

        self.fig.clear()
        self.fig.set_size_inches(12, 6)
        ax = self.fig.add_subplot(111)

        # create normal distributions
        for i, rank in enumerate(sorted(rel_dists.keys())):
            v = [dist for taxa, dist in rel_dists[rank].iteritems() if taxa in taxa_for_dist_inference]
            u = np_mean(v)
            rv = norm(loc=u, scale=np_std(v))
            x = np_linspace(rv.ppf(0.001), rv.ppf(0.999), 1000)
            nd = rv.pdf(x)
            ax.plot(x, 0.75 * (nd / max(nd)) + i, 'b-', alpha=0.6, zorder=2)
            ax.plot((u, u), (i, i + 0.5), 'b-', zorder=2)

        # create percentile lines
        percentiles = {}
        for i, rank in enumerate(sorted(rel_dists.keys())):
            v = [dist for taxa, dist in rel_dists[rank].iteritems() if taxa in taxa_for_dist_inference]
            p10, p50, p90 = np_percentile(v, [10, 50, 90])
            ax.plot((p10, p10), (i, i + 0.5), 'r-', zorder=2)
            ax.plot((p50, p50), (i, i + 0.5), 'r-', zorder=2)
            ax.plot((p90, p90), (i, i + 0.5), 'r-', zorder=2)

            percentiles[i] = [p10, p50, p90]

        # create scatter plot and results table
        fout = open(distribution_table, 'w')
        fout.write('Taxa\tRelative Distance\tRank cutoff\tRank outlier\tP10\tMedian\tP90\tPercentile outlier\n')
        x = []
        y = []
        c = []
        labels = []
        rank_labels = []
        rel_dist_thresholds += [1.0]  # append boundry for species
        for i, rank in enumerate(sorted(rel_dists.keys())):
            rank_label = Taxonomy.rank_labels[rank]
            rank_labels.append(rank_label + ' (%d)' % len(rel_dists[rank]))

            for clade_label, dist in rel_dists[rank].iteritems():
                x.append(dist)
                y.append(i)
                labels.append(clade_label)

                if clade_label in taxa_for_dist_inference:
                    c.append((0.0, 0.0, 0.5))
                else:
                    c.append((0.5, 0.5, 0.5))

                p10, p50, p90 = percentiles[i]
                percentile_outlier = not (dist >= p10 and dist <= p90)

                if i == 0:
                    rank_cutoff = rel_dist_thresholds[i]
                    rank_outlier = dist > rank_cutoff
                else:
                    rank_cutoff = rel_dist_thresholds[i]
                    upper_rank_cutoff = rel_dist_thresholds[i - 1]
                    rank_outlier = not (dist >= upper_rank_cutoff and dist <= rank_cutoff)

                v = [clade_label, dist, rank_cutoff, str(rank_outlier)]
                v += percentiles[i] + [str(percentile_outlier)]
                fout.write('%s\t%.2f\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%s\n' % tuple(v))
        fout.close()

        scatter = ax.scatter(x, y, alpha=0.5, s=48, c=c, zorder=1)

        # set plot elements
        ax.grid(color=(0.8, 0.8, 0.8), linestyle='dashed')

        ax.set_xlabel('relative distance')
        ax.set_xticks(np_arange(0, 1.05, 0.1))
        ax.set_xlim([-0.05, 1.05])

        ax.set_ylabel('rank (no. taxa)')
        ax.set_yticks(xrange(0, len(rel_dists)))
        ax.set_ylim([-0.2, len(rel_dists) - 0.01])
        ax.set_yticklabels(rank_labels)

        self.prettify(ax)

        # plot relative divergence threshold lines
        y_min, y_max = ax.get_ylim()
        for threshold in rel_dist_thresholds[0:-1]:  # don't draw species boundary
            ax.plot((threshold, threshold), (y_min, y_max), color='r', ls='--')
            ax.text(threshold + 0.001, y_max, '%.3f' % threshold, horizontalalignment='center')

        # make plot interactive
        mpld3.plugins.connect(self.fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=10))
        mpld3.save_html(self.fig, plot_file[0:plot_file.rfind('.')] + '.html')

        self.fig.tight_layout(pad=1)
        self.fig.savefig(plot_file, dpi=96)

    def run(self, input_tree, output_prefix, plot_taxa_file, trusted_taxa_file, min_children, min_support):
        """Determine distribution of taxa at each taxonomic rank.

        Parameters
        ----------
        input_tree : str
            Name of input tree.
        output_prefix : str
            Desired prefix for generated files.
        plot_taxa_file : str
            File specifying taxa to plot. Set to None to consider all taxa.
        trusted_taxa_file : str
            File specifying trusted taxa to consider when inferring distribution. Set to None to consider all taxa.
        min_children : int
            Only consider taxa with at least the specified number of children taxa when inferring distribution.
        min_support : float
            Only consider taxa with at least this level of support when inferring distribution.
        """

        # read tree
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)

        # read taxa to plot
        taxa_to_plot = None
        if plot_taxa_file:
            taxa_to_plot = read_taxa_file(plot_taxa_file)

        # read trusted taxa
        trusted_taxa = None
        if trusted_taxa_file:
            trusted_taxa = read_taxa_file(trusted_taxa_file)

        # determine taxa to be used for inferring distribution
        taxa_for_dist_inference = filter_taxa_for_dist_inference(tree, trusted_taxa, min_children, min_support)

        # calculate relative distance to taxa
        rd = RelativeDistance()
        rel_dists = rd.rel_dist_to_named_clades(tree, taxa_to_plot)

        # report number of taxa at each rank
        print ''
        print '  Number of taxa plotted at each taxonomic rank:'
        for rank, taxa in rel_dists.iteritems():
            print '    %s\t%d' % (Taxonomy.rank_labels[rank], len(taxa))

        # create performance plots
        rel_dist_thresholds = self._percent_correct_plot(rel_dists, taxa_for_dist_inference, output_prefix)

        # create distribution plot
        distribution_table = output_prefix + '.tsv'
        plot_file = output_prefix + '.png'
        self._distribution_plot(rel_dists, rel_dist_thresholds, taxa_for_dist_inference, distribution_table, plot_file)
