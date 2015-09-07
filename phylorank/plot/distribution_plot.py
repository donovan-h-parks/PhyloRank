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
from phylorank.newick import parse_label, read_from_tree

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
    """Plot distribution of taxa in each taxonomic rank."""

    def __init__(self):
        """Initialize."""
        AbstractPlot.__init__(self, None)

    def _rel_dist_to_named_clades(self, root, taxa_to_consider):
        """Determine relative distance to specific taxa.

        Parameters
        ----------
        root : TreeNode
            Root of tree.
        taxa_to_consider : set
            Named taxonomic groups to consider.

        Returns
        -------
        dict : d[rank_index][taxon] -> relative divergence
        """

        # calculate relative distance for all nodes
        infer_rank = InferRank()
        infer_rank.decorate_rel_dist(root)

        # assign internal nodes with ranks from
        rel_dists = defaultdict(dict)
        for node in root.preorder(include_self=False):
            if not node.name or node.is_tip():
                continue

            # check for support value
            _support, taxon_name, _auxiliary_info = parse_label(node.name)

            if not taxon_name:
                continue

            # get most-specific rank if a node represents multiple ranks
            if ';' in taxon_name:
                taxon_name = taxon_name.split(';')[-1].strip()

            if taxon_name not in taxa_to_consider:
                continue

            most_specific_rank = taxon_name[0:3]
            rel_dists[Taxonomy.rank_index[most_specific_rank]][taxon_name] = node.rel_dist

        return rel_dists

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
            self.fig.savefig(output_prefix + '.%s_%s.png' % (Taxonomy.rank_labels[parent_rank], Taxonomy.rank_labels[child_rank]), dpi=300)

        return rel_dist_thresholds

    def _distribution_plot(self, rel_dists, rel_dist_thresholds, distribution_table, plot_file):
        """Create plot showing the distribution of taxa at each taxonomic rank.

        Parameters
        ----------
        rel_dists: d[rank_index][taxon] -> relative divergence
            Relative divergence of taxa at each rank.
        rel_dist_thresholds: ?
            ???
        distribution_table : str
            Desired name of output table with distribution information.
        plot_file : str
            Desired name of output plot.
        """

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
        fout = open(distribution_table, 'w')
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
        for threshold in rel_dist_thresholds:
            ax.plot((threshold, threshold), (y_min, y_max), color='r', ls='--')
            ax.text(threshold + 0.001, y_max, '%.3f' % threshold, horizontalalignment='center')

        # make plot interactive
        mpld3.plugins.connect(self.fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        mpld3.plugins.connect(self.fig, mpld3.plugins.MousePosition(fontsize=12))
        mpld3.save_html(self.fig, plot_file[0:plot_file.rfind('.')] + '.html')

        self.fig.tight_layout(pad=1)
        self.fig.savefig(plot_file, dpi=300)

    def _read_taxa_file(self, taxa_file):
        """Read taxa from file.

        Parameters
        ----------
        taxa_file : str
            File specifying taxa to consider. One per line.

        Returns
        -------
        set
            Taxa.
        """

        taxa = set()
        for line in open(taxa_file):
            taxa.add(line.strip())

        return taxa

    def _taxa_for_dist_inference(self, tree, trusted_taxa, min_children, min_support):
        """Determine taxa to use for inferring distribution of relative divergences.

        Parameters
        ----------
        tree : TreeNode
            Phylogenetic tree.
        trusted_taxa : iterable
            Trusted taxa to consider when inferring distribution.
        min_children : int
            Only consider taxa with at least the specified number of children taxa when inferring distribution.
        min_support : float
            Only consider taxa with at least this level of support when inferring distribution.
        """

        # read taxonomy and determine children taxa for each named group
        taxonomy = Taxonomy()
        t = read_from_tree(tree)
        taxon_children = taxonomy.taxon_children(t)

        # get all named groups
        taxa_for_dist_inference = set()
        for taxon_id, taxa in t.iteritems():
            for taxon in taxa:
                taxa_for_dist_inference.add(taxon)

        # sanity check species names as these are a common problem
        for taxon_id, taxa in t.iteritems():
            if len(taxa) > Taxonomy.rank_index['s__']:
                species_name = taxa[Taxonomy.rank_index['s__']]
                valid, error_msg = taxonomy.validate_species_name(species_name, require_full=True, require_prefix=True)
                if not valid:
                    print '[Error] Species name %s for %s is invalid: %s' % (species_name, taxon_id, error_msg)
                    sys.exit(-1)

        # restrict taxa to those with a sufficient number of named children
        # Note: a taxonomic group with no children will not end up in the
        # taxon_children data structure so care must be taken when applying
        # this filtering criteria.
        if min_children > 0:
            valid_taxa = set()
            for taxon, children_taxa in taxon_children.iteritems():
                if len(children_taxa) >= min_children:
                    valid_taxa.add(taxon)

            taxa_for_dist_inference.intersection_update(valid_taxa)

        # restrict taxa used for inferring distribution to those with sufficient support
        if min_support > 0:
            for node in tree.preorder(include_self=False):
                if not node.name or node.is_tip():
                    continue

                # check for support value
                support, taxon_name, _auxiliary_info = parse_label(node.name)

                if not taxon_name:
                    continue

                if support and float(support) < min_support:
                    taxa_for_dist_inference.difference_update([taxon_name])
                elif not support and min_support > 0:
                    # no support value, so inform user if they were trying to filter on this property
                    print '[Error] Tree does not contain support values. As such, --min_support should be set to 0.'
                    continue

        # restrict taxa used for inferring distribution to the trusted set
        if trusted_taxa:
            taxa_for_dist_inference = trusted_taxa.intersection(taxa_for_dist_inference)

        return taxa_for_dist_inference

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
        tree = TreeNode.read(input_tree, convert_underscores=False)

        # read taxa to plot
        taxa_to_plot = None
        if plot_taxa_file:
            taxa_to_plot = self._read_taxa_file(plot_taxa_file)

        # read trusted taxa
        trusted_taxa = None
        if trusted_taxa_file:
            trusted_taxa = self._read_taxa_file(trusted_taxa_file)

        # determine taxa to be used for inferring distribution
        taxa_for_dist_inference = self._taxa_for_dist_inference(tree, trusted_taxa, min_children, min_support)

        # calculate relative distance to taxa
        rel_dists = self._rel_dist_to_named_clades(tree, taxa_to_plot)

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
        self._distribution_plot(rel_dists, rel_dist_thresholds, distribution_table, plot_file)
