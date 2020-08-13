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

import logging
from collections import defaultdict

from biolib.taxonomy import Taxonomy

from phylorank.newick import parse_label


class RelativeDistance(object):
    """Determine relative rates of evolutionary divergence."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def _avg_descendant_rate(self, tree):
        """Calculate average rate of divergence for each nodes in a tree.

        The average rate is the arithmetic mean of the
        branch length to all descendant taxa.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
        """

        # calculate the mean branch length to extant taxa
        for node in tree.postorder_node_iter():
            avg_div = 0
            if node.is_leaf():
                node.mean_dist = 0.0
                node.num_taxa = 1
            else:
                node.num_taxa = sum([1 for _ in node.leaf_iter()])
                for c in node.child_node_iter():
                    num_tips = c.num_taxa
                    avg_div += (float(c.num_taxa) / node.num_taxa) * (c.mean_dist + c.edge_length)

            node.mean_dist = avg_div

    def decorate_rel_dist(self, tree, mblet=False):
        """Calculate relative distance to each internal node.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_taxa: number of terminal taxa
          rel_dists: relative distance of node between root and extant organisms
        """

        self._avg_descendant_rate(tree)
        if mblet:
            for node in tree.preorder_node_iter():
                node.rel_dist = node.mean_dist
        else:
            for node in tree.preorder_node_iter():
                if node == tree.seed_node:
                    node.rel_dist = 0.0
                elif node.is_leaf():
                    node.rel_dist = 1.0
                else:
                    a = node.edge_length
                    b = node.mean_dist
                    x = node.parent_node.rel_dist

                    if (a + b) != 0:
                        rel_dist = x + (a / (a + b)) * (1.0 - x)
                    else:
                        # internal node has zero length to parent,
                        # so should have the same relative distance
                        # as the parent node
                        rel_dist = x

                    node.rel_dist = rel_dist

    def rel_dist_to_named_clades(self, tree, mblet=False):
        """Determine relative distance to specific taxa.

        Parameters
        ----------
        tree : Dendropy Tree
            Phylogenetic tree.

        Returns
        -------
        dict : d[rank_index][taxon] -> relative divergence
        """

        # calculate relative distance for all nodes
        self.decorate_rel_dist(tree, mblet)

        # tabulate values for internal nodes with ranks
        rel_dists = defaultdict(dict)
        for node in tree.preorder_node_iter(lambda n: n != tree.seed_node):
            if not node.label or node.is_leaf():
                continue

            _support, taxon_name, _auxiliary_info = parse_label(node.label)
            if not taxon_name:
                continue

            # get most-specific rank if a node represents multiple ranks
            if ';' in taxon_name:
                taxon_name = taxon_name.split(';')[-1].strip()

            most_specific_rank = taxon_name[0:3]
            rel_dists[Taxonomy.rank_index[most_specific_rank]][taxon_name] = node.rel_dist

        return rel_dists
