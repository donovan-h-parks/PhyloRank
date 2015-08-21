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
import logging

from numpy import (mean as np_mean,
                   std as np_std,
                   arange as np_arange,
                   percentile as np_percentile)


class InferRank():
    """Infer taxonomic ranks from relative rates of evolutionary divergence."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        # branch length used to define an OTU
        self.otu_branch_length = 0.01

    def _avg_descendant_rate(self, tree):
        """Calculate average rate of divergence for each nodes in a tree.

        The average rate is weighted by the number of OTUs in each child
        lineage. A node is classified as an OTU if it is the most proximal
        node in a lineage where the sum of the distal branch length and
        average branch length to extant organisms is greater than a
        defined threshold. For OTUs and nodes more proximal than the
        defined OTU, the average rate is the arithmetic average of the
        branch length to all descendant taxa.

        Parameters
        ----------
        tree : TreeNode
            Root node of (sub)tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_otus: number of OTUs descendant from node
          weighted_dist: OTU-weighted average distance
        """

        # calculate the mean branch length to extant taxa
        for node in tree.postorder():
            avg_div = 0
            if not node.is_tip():
                total_tips = len(list(node.tips()))
                for c in node.children:
                    num_tips = 1
                    if not c.is_tip():
                        num_tips = len(list(c.tips()))
                    avg_div += (float(num_tips) / total_tips) * (c.mean_dist + c.length)

            node.mean_dist = avg_div

        # calculate the OTU-weighted average distance for all nodes
        for node in tree.postorder():
            node.num_otus = 0
            node.weighted_dist = 0

            if node.is_tip():
                continue

            sum_otus = sum([max(c.num_otus, 1) for c in node.children])
            weighted_dist = 0
            for c in node.children:
                otus = c.num_otus
                if c.mean_dist <= self.otu_branch_length:
                    # use mean distance as child is more distal than
                    # the mean branch length use to define an OTU
                    dist = c.mean_dist + c.length
                else:
                    # use OTU-weighted average distance
                    dist = c.weighted_dist + c.length

                weighted_dist += (float(max(otus, 1)) / sum_otus) * dist

                # accumulate the number of OTUs below this node
                if otus != 0:
                    node.num_otus += otus
                elif dist > self.otu_branch_length:
                    # child lineage meets the requirements of an OTU
                    node.num_otus += 1

            node.weighted_dist = weighted_dist

    def decorate_rel_dist(self, tree):
        """Calculate relative distance to each internal node.

        Parameters
        ----------
        tree : TreeNode
            Root node of (sub)tree.

        Returns
        -------
        The following attributes are added to each node:
          mean_dist: mean distance to tips
          num_otus: number of OTUs descendant from node
          weighted_dist: OTU-weighted average distance
          rel_dists: relative distance of node between root and extant organisms
        """

        self._avg_descendant_rate(tree)

        for node in tree.preorder():
            if node.is_root():
                node.rel_dist = 0.0
            elif node.is_tip():
                node.rel_dist = 1.0
            else:
                a = node.length
                b = node.weighted_dist
                x = node.parent.rel_dist
                rel_dist = x + (a / (a + b)) * (1.0 - x)

                node.rel_dist = rel_dist
