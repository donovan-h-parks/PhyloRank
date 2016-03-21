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

from collections import defaultdict

from phylorank.infer_rank import InferRank
from phylorank.newick import parse_label, read_from_tree

from biolib.taxonomy import Taxonomy

from skbio import TreeNode


def rel_dist_to_named_clades(root, taxa_to_consider):
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

        if taxa_to_consider and taxon_name not in taxa_to_consider:
            continue

        most_specific_rank = taxon_name[0:3]
        rel_dists[Taxonomy.rank_index[most_specific_rank]][taxon_name] = node.rel_dist

    return rel_dists
