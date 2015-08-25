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
from collections import defaultdict

from phylorank.infer_rank import InferRank
from phylorank.newick import parse_label

from skbio import TreeNode

from biolib.taxonomy import Taxonomy
from biolib.plots.abstract_plot import AbstractPlot

from numpy import (mean as np_mean,
                   std as np_std,
                   arange as np_arange,
                   linspace as np_linspace,
                   percentile as np_percentile)

from scipy.stats import norm

import mpld3


class DistributionPlot(AbstractPlot):
    """Plot distribution of groups in each taxonomic rank."""

    def __init__(self):
        """Initialize."""
        AbstractPlot.__init__(self, None)

    def rel_dist_to_named_clades(self, tree, start_label, taxa_to_consider, min_support):
        """Determine relative distance to named clades.

        Parameters
        ----------
        tree : str or TreeNode
            A newick string or a TreeNode
        start_label : str
            Label of node to start marking taxonomic ranks.
        taxa_to_consider : set
            Taxonomic groups to consider.
        min_support : float
            Required support to consider taxon.

        Returns
        -------
        dict : d[rank_index][taxon] -> relative distance to root
        """

        # make sure we have a TreeNode object
        root = tree
        if not isinstance(root, TreeNode):
            root = TreeNode.read(root, convert_underscores=False)

        # calculate relative distance for all nodes
        infer_rank = InferRank()
        infer_rank.decorate_rel_dist(root)

        # find specified starting node
        start_node = None
        for n in root.preorder():
            if n.name and start_label in n.name:
                start_node = n

        if not start_node:
            print 'Unable to locate node with label: ' + start_label
            return None

        # assign internal nodes with ranks from
        rel_dists = defaultdict(dict)
        for node in start_node.preorder(include_self=False):
            if not node.name or node.is_tip():
                continue

            # check for support value
            support, taxon_name = parse_label(node.name)

            if not taxon_name:
                continue

            if support and float(support) < min_support:
                continue
            elif not support and min_support > 0:
                # no support value, so inform user if they were trying to filter on this property
                print '[Error] Tree does not contain support values. As such, --min_support must be set to 0.'
                continue

            # get most-specific rank if a node represents multiple ranks
            if ';' in taxon_name:
                taxon_name = taxon_name.split(';')[-1].strip()

            if taxon_name not in taxa_to_consider:
                continue

            most_specific_rank = taxon_name[0:3]
            rel_dists[Taxonomy.rank_index[most_specific_rank]][taxon_name] = node.rel_dist

        return rel_dists

    def percent_correct_plot(self, rel_dists, output_prefix):
        """Create plots showing correctly classified taxa for different relative distance values."""

        print ''
        print '  Relative divergence thresholds:'

        ranks = sorted(rel_dists.keys())
        rel_dist_thresholds = []
        for i in xrange(ranks[0], ranks[-1]):
            parent_rank = i
            child_rank = i + 1

            # determine classification results for relative divergence
            # values between the medians of adjacent taxonomic ranks
            parent_rds = rel_dists[parent_rank].values()
            parent_p50 = np_percentile(parent_rds, 50)

            child_rds = rel_dists[child_rank].values()
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
            print '    %s\t%.3f' % (Taxonomy.rank_labels[parent_rank], r_max_value)

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
            self.fig.savefig(output_prefix + '.%s_%s.png' % (Taxonomy.rank_labels[parent_rank], Taxonomy.rank_labels[child_rank]), dpi=300)

        return rel_dist_thresholds

    def distribution_plot(self, rel_dists, rel_dist_thresholds, title, output_prefix):
        """Create plot showing the distribution of taxa at each taxonomic rank."""

        self.fig.clear()
        self.fig.set_size_inches(12, 4.5)
        ax = self.fig.add_subplot(111)

        # create normal distributions
        for i, rank in enumerate(sorted(rel_dists.keys())):
            v = rel_dists[rank].values()
            u = np_mean(v)
            rv = norm(loc=u, scale=np_std(v))
            x = np_linspace(rv.ppf(0.001), rv.ppf(0.999), 1000)
            nd = rv.pdf(x)
            ax.plot(x, 0.75 * (nd / max(nd)) + i, 'b-', alpha=0.6, zorder=2)
            ax.plot((u, u), (i, i + 0.5), 'b-', zorder=2)

        # create percentile lines
        for i, rank in enumerate(sorted(rel_dists.keys())):
            p10, p50, p90 = np_percentile(rel_dists[rank].values(), [5, 50, 95])
            ax.plot((p10, p10), (i, i + 0.5), 'r-', zorder=2)
            ax.plot((p50, p50), (i, i + 0.5), 'r-', zorder=2)
            ax.plot((p90, p90), (i, i + 0.5), 'r-', zorder=2)

            if rank == Taxonomy.rank_index['p__']:
                rel_dist_thresholds = [p10] + rel_dist_thresholds

        # create scatter plot and results table
        fout = open(output_prefix + '.tsv', 'w')
        x = []
        y = []
        labels = []
        rank_labels = []
        for i, rank in enumerate(sorted(rel_dists.keys())):
            rank_label = Taxonomy.rank_labels[rank]
            rank_labels.append(rank_label + ' (%d)' % len(rel_dists[rank]))

            fout.write(rank_label)
            for clade_label, dist in rel_dists[rank].iteritems():
                x.append(dist)
                y.append(i)
                labels.append(clade_label)

                fout.write('\t%s\t%.3f\n' % (clade_label, dist))
        fout.close()

        scatter = ax.scatter(x, y, alpha=0.5, s=48, c=(0.5, 0.5, 0.5), zorder=1)

        # set plot elements
        ax.grid(color=(0.8, 0.8, 0.8), linestyle='dashed')
        if title:
            ax.set_title(title, size=12)

        ax.set_xlabel('relative distance')
        ax.set_xticks(np_arange(0, 1.05, 0.1))
        ax.set_xlim([-0.05, 1.05])

        ax.set_ylabel('rank (no. taxa)')
        ax.set_yticks(xrange(0, len(rel_dists)))
        ax.set_ylim([-0.2, len(rel_dists) - 0.01])
        ax.set_yticklabels(rank_labels)

        self.prettify(ax)

        # plot relative divergence threshold lines
        _y_min, y_max = ax.get_ylim()
        for threshold in rel_dist_thresholds:
            ax.axvline(x=threshold, ymin=0, ymax=1, color='r', ls='--')
            ax.text(threshold + 0.001, y_max, '%.3f' % threshold, horizontalalignment='center')

        # make plot interactive
        mpld3.plugins.connect(self.fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=12))
        mpld3.save_html(self.fig, output_prefix + '.html')

        plot_file = output_prefix + '.png'
        self.fig.tight_layout(pad=1)
        self.fig.savefig(plot_file, dpi=300)

        return plot_file

    def run(self, input_tree, output_prefix, min_children, min_support, title):

        # read taxonomy and determine children taxa for each named group
        taxonomy = Taxonomy().read_from_tree(input_tree)
        taxon_children = Taxonomy().taxon_children(taxonomy)

        # read tree
        tree = TreeNode.read(input_tree, convert_underscores=False)

        # determine taxa with at least the specified number of children
        taxa_to_consider = set()
        for taxon, children_taxa in taxon_children.iteritems():
            if len(children_taxa) >= min_children:
                taxa_to_consider.add(taxon)

            if 's__' in taxon:
                taxa_to_consider.add(taxon)

        # calculate relative distance to taxa
        rel_dists = self.rel_dist_to_named_clades(tree, 'd__Bacteria', taxa_to_consider, min_support)

        # report number of taxa at each rank
        print ''
        print '  Number of taxa considered at each taxonomic rank:'
        for rank, taxa in rel_dists.iteritems():
            print '    %s\t%d' % (Taxonomy.rank_labels[rank], len(taxa))

        # create performance plots
        rel_dist_thresholds = self.percent_correct_plot(rel_dists, output_prefix)

        # create distribution plot
        plot_file = self.distribution_plot(rel_dists, rel_dist_thresholds, title, output_prefix)

        return plot_file
