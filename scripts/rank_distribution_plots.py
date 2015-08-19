#!/usr/bin/env python

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

__prog_name__ = 'rank_distribution_plots'
__prog_desc__ = 'create plots showing the distribution of taxonomic ranks'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import argparse
from collections import defaultdict

from skbio import TreeNode

from numpy import (mean as np_mean,
                   std as np_std,
                   arange as np_arange,
                   linspace as np_linspace,
                   percentile as np_percentile)

from scipy.stats import norm

import matplotlib.pyplot as plt
import mpld3


class RankDistributionPlots(object):
    """Plot relative distance to taxonomic ranks.

    Relative distance is calculated as:
        dist_to_root / (dist_to_root + average_dist)

        where dist_to_root is the distance from a
        named clade to the specified starting node,
        and average_dist is the average distance from
        the named clade to all leaf nodes.
    """

    def __init__(self):
        """Initialize."""

        self.rank_prefixes = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')
        self.rank_labels = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
        self.rank_index = {'d__': 0, 'p__': 1, 'c__': 2, 'o__': 3, 'f__': 4, 'g__': 5, 's__': 6}

    def _mean_dist_to_leaves(self, node):
        """Calculate mean distance from an internal node to its leaves.

        Parameters
        ----------
        node : TreeNode
            Internal node of interest.

        Returns
        -------
        float
            Mean distance to leaves.
        """

        if node.is_tip():
            return 0

        dist = []
        for leaf in node.tips():
            dist.append(leaf.accumulate_to_ancestor(node))

        return np_mean(dist)

    def rel_dist_to_named_clades(self, tree, start_label, taxa_to_consider):
        """Determine relative distance to named clades.

        Relative distance is calculated as:
            dist_to_root / (dist_to_root + average_dist)

            where dist_to_root is the distance from a
            named clade to the specified starting node,
            and average_dist is the average distance from
            the named clade to all leaf nodes.

        Parameters
        ----------
        tree : str or TreeNode
            A newick string or a TreeNode
        start_label: str
            Label of node to start marking taxonomic ranks.
        taxa_to_consider: set
            Taxonomic groups to consider.

        Returns
        -------
        dict : d[rank_index][label] -> relative distance to root
        """

        # make sure we have a TreeNode object
        root = tree
        if not isinstance(root, TreeNode):
            root = TreeNode.read(root, convert_underscores=False)

        # find specified starting node
        start_node = None
        for n in root.preorder():
            if n.name and start_label in n.name:
                start_node = n

        if not start_node:
            print 'Unable to locate node with label: ' + start_label
            return None

        print 'Tips below start node %s: %d' % (start_label, len(list(start_node.tips())))

        # assign internal nodes with ranks from
        output = defaultdict(dict)
        for node in start_node.preorder(include_self=False):
            if not node.name or node.is_tip():
                continue

            # get most-specific rank if a node represents multiple ranks
            if ';' in node.name:
                node.name = node.name.split(';')[-1].strip()

            if node.name not in taxa_to_consider:
                continue

            most_specific_rank = node.name[0:3]

            # get relative distance from root to named child clade
            average_dist = self._mean_dist_to_leaves(node)
            dist_to_root = node.accumulate_to_ancestor(start_node)
            rel_dist_to_root = dist_to_root / (dist_to_root + average_dist)

            if ':' in node.name:
                node.name = node.name[node.name.find(':') + 1:]
            output[self.rank_index[most_specific_rank]][node.name] = rel_dist_to_root

        return output

    def prettify(self, axis):
        """Modify axis properties to make a cleaner plot."""

        axes_colour = (0.5, 0.5, 0.5)

        for a in axis.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in axis.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axis.yaxis.get_ticklines():
            line.set_color(axes_colour)

        for line in axis.xaxis.get_ticklines():
            line.set_color(axes_colour)

        for loc, spine in axis.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(axes_colour)

    def run(self, tree_file, start_label, output_prefix, min_children):
        """Plot relative distance to named groups."""

        # read tree
        tree = TreeNode.read(tree_file, convert_underscores=False)

        # read taxonomy from tree
        print 'Reading taxonomy from tree.'
        taxonomy = {}
        for tip in tree.tips():
            taxa = []

            node = tip.parent
            while node:
                if node.name:
                    taxa = [x.strip() for x in node.name.split(';')] + taxa
                node = node.parent

            taxonomy[tip.name] = taxa
            if len(taxonomy[tip.name]) == len(self.rank_prefixes):
                # all ranks are present so mark the genome id
                taxonomy[tip.name] += [tip.name]

        # determine named taxa with at least the specified number of children
        print 'Determining taxa with sufficient named children.'
        taxon_children = defaultdict(set)
        for taxa in taxonomy.values():
            for i, taxon in enumerate(taxa):
                if len(taxa) > i + 1:
                    taxon_children[taxon].add(taxa[i + 1])

        taxa_to_consider = set()
        for taxon, children_taxa in taxon_children.iteritems():
            if len(children_taxa) >= min_children:
                taxa_to_consider.add(taxon)

        # calculate relative distance to named taxa
        print 'Calculating relative distance to taxa.'
        rel_dist = self.rel_dist_to_named_clades(tree, start_label, taxa_to_consider)
        if not rel_dist:
            return

        fig, ax = plt.subplots(figsize=(12, 4.5), tight_layout=True, frameon=False)

        # create normal distributions
        for i, rank in enumerate(sorted(rel_dist.keys())):
            v = rel_dist[rank].values()
            u = np_mean(v)
            rv = norm(loc=u, scale=np_std(v))
            x = np_linspace(rv.ppf(0.001), rv.ppf(0.999), 1000)
            nd = rv.pdf(x)
            ax.plot(x, 0.75 * (nd / max(nd)) + i, 'b-', alpha=0.6)
            ax.plot((u, u), (i, i + 0.75), 'b--')

        # create percentile lines
        for i, rank in enumerate(sorted(rel_dist.keys())):
            p5, p50, p95 = np_percentile(rel_dist[rank].values(), [5, 50, 95])
            ax.plot((p5, p5), (i, i + 0.5), 'r--')
            ax.plot((p50, p50), (i, i + 0.75), 'r--')
            ax.plot((p95, p95), (i, i + 0.5), 'r--')

        # create scatter plot and results table
        fout = open(output_prefix + '.tsv', 'w')
        x = []
        y = []
        labels = []
        rank_labels = []
        for i, rank in enumerate(sorted(rel_dist.keys())):
            rank_label = self.rank_labels[rank]
            rank_labels.append(rank_label)

            fout.write(rank_label)
            for clade_label, dist in rel_dist[rank].iteritems():
                x.append(dist)
                y.append(i)
                labels.append(clade_label)

                fout.write('\t%s\t%.3f\n' % (clade_label, dist))
        fout.close()

        scatter = ax.scatter(x, y, alpha=0.5, s=48, c=(0.5, 0.5, 0.5))

        # set plot elements
        ax.grid(color=(0.8, 0.8, 0.8), linestyle='dashed')
        ax.set_title("Distribution for %s" % start_label, size=12)

        ax.set_xlabel('relative distance')
        ax.set_xticks(np_arange(0, 1.05, 0.1))
        ax.set_xlim([-0.05, 1.05])

        ax.set_ylabel('rank')
        ax.set_yticks(xrange(0, len(rel_dist)))
        ax.set_ylim([-0.2, len(rel_dist) - 0.01])
        ax.set_yticklabels(rank_labels)

        self.prettify(ax)

        # make plot interactive
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        mpld3.plugins.connect(fig, mpld3.plugins.MousePosition(fontsize=12))

        mpld3.save_html(fig, output_prefix + '.html')
        fig.savefig(output_prefix + '.png', dpi=300)

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('tree', help='input tree file with decorated internal nodes')
    parser.add_argument('start_label', help='label of node to use a root for distance calculations')
    parser.add_argument('output_prefix', help='prefix for output files')
    parser.add_argument('-m', '--min_children', help='minimum named child taxa to consider taxa', type=int, default=2)

    args = parser.parse_args()

    try:
        p = RankDistributionPlots()
        p.run(args.tree, args.start_label, args.output_prefix, args.min_children)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
