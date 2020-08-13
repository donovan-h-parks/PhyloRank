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

import dendropy

from phylorank.newick import parse_label
from phylorank.rel_dist import RelativeDistance

'''
To do:
  - should produce a flat file indicating existing taxa labels,
    predicted taxa labels, relative divergence, and percentiles.
'''


class MarkTree(object):
    """Mark nodes with taxonomic ranks inferred from evolutionary divergence."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        self.rank_prefixes = ['D__', 'P__', 'C__', 'O__', 'F__', 'G__', 'S__', 'ST__']
        self.rank_designators = ['d', 'p', 'c', 'o', 'f', 'g', 's', 'st']
        self.highly_basal_designator = 'X__'

    def run(self, input_tree,
            output_tree,
            min_support,
            only_named_clades,
            min_length,
            show_percentiles,
            show_relative_divergence,
            show_prediction,
            thresholds):
        """Read distribution file.

        Parameters
        ----------
        input_tree : str
            Name of input tree.
        output_tree : str
            Name of output tree.
        min_support : int
            Only decorate nodes above specified support value.
        only_named_clades : boolean
            Only decorate nodes with existing labels.
        min_length : float
            Only decorate nodes above specified length.
        show_percentiles : bool
            Flag indicating if percentiles should be placed on nodes.
        show_relative_divergence : bool
            Flag indicating if relative divergences should be placed on nodes.
        show_prediction : bool
            Flag indicating if predicate ranks should be placed on nodes.
        thresholds : d[rank] -> threshold
            Relative divergence threshold for defining taxonomic ranks.
        """

        # make sure we have a TreeNode object
        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        # calculate relative distance for all nodes
        rd = RelativeDistance()
        rd.decorate_rel_dist(tree)

        # decorate nodes based on specified criteria
        self.logger.info('')
        self.logger.info('  %s\t%s' % ('Rank', 'Prediction results'))

        correct = defaultdict(int)
        incorrect = defaultdict(int)

        fout = open(output_tree + '.info', 'w')
        fout.write(
            'Taxon name\tPredicted rank\tRelative divergence\tCurrent rank percentile\tPredicted rank percentile\n')
        for n in tree.preorder_node_iter():
            if n.is_leaf():
                continue

            if n.edge_length and n.edge_length < min_length:
                continue

            # parse taxon name and support value from node label
            if n.label:
                support, taxon_name, _auxiliary_info = parse_label(n.label)
                n.label += '|'
            else:
                support = 100
                taxon_name = None
                n.label = ''

            if support and float(support) < min_support:
                continue

            if only_named_clades and not taxon_name:
                continue

            # Decorate node with predicted rank prefix. Nodes with
            # a relative divergence greater than the genus threshold
            # are a species. Nodes with a relative divergence less than
            # the domain threshold have no real prediction, so are marked
            # with an 'X__', All other nodes will be assigned an intermediate
            # rank based on the threshold values.
            if show_prediction:
                # calculate distance to each median threshold
                min_dist = 1e6
                predicted_rank = None
                for rank, threshold in thresholds.items():
                    d = abs(n.rel_dist - threshold)
                    if d < min_dist:
                        min_dist = d
                        rank_index = self.rank_designators.index(rank)
                        predicted_rank = self.rank_prefixes[rank_index]

                n.label += predicted_rank

            if show_relative_divergence:
                n.label += '[rd=%.2f]' % n.rel_dist

            if taxon_name and predicted_rank != self.highly_basal_designator:
                # tabulate number of correct and incorrect predictions
                named_rank = taxon_name.split(';')[-1][0:3]
                if named_rank == predicted_rank.lower():
                    correct[named_rank] += 1
                else:
                    incorrect[named_rank] += 1

            if taxon_name:
                fout.write('%s\t%s\t%.3f\n' % (taxon_name, predicted_rank, n.rel_dist))

        fout.close()
        tree.write_to_path(output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)

        for rank_prefix in self.rank_prefixes[1:7]:
            correct_taxa = correct[rank_prefix.lower()]
            incorrect_taxa = incorrect[rank_prefix.lower()]
            total_taxa = correct_taxa + incorrect_taxa
            if total_taxa > 0:
                self.logger.info('  {}\t{:,} of {:,} ({:.2f}%)'.format(rank_prefix.lower(), correct_taxa, total_taxa,
                                                                       correct_taxa * 100.0 / total_taxa))
            else:
                self.logger.info('  No taxa at rank {}.'.format(rank_prefix))
