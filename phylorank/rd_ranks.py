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
import os
from collections import defaultdict

import dendropy
from biolib.taxonomy import Taxonomy
from numpy import (mean as np_mean,
                   std as np_std)

from phylorank.common import get_phyla_lineages
from phylorank.newick import parse_label
from phylorank.rel_dist import RelativeDistance


class RdRanks(object):
    """Calculate number of taxa for specified relative divergence thresholds.

        <blah>
    """

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger()

    def write_rank_count(self, ranks_below_taxon, results_table):
        """Write table indicating number of ranks below each taxa.

        Parameters
        ----------
        ranks_below_taxon : d[taxon][rank prefix] -> count, or list of counts
            Number of ranks below named taxon.
        results_table : str
            Desired output file.
        """

        # determine if count is a scalar or vectors
        taxon = ranks_below_taxon.keys()[0]
        rank_prefix = ranks_below_taxon[taxon].keys()[0]
        count = ranks_below_taxon[taxon][rank_prefix]

        count_is_scalar = True
        if isinstance(count, (list, tuple)):
            count_is_scalar = False

        # write out results sorted by taxonomic rank        
        sorted_taxon = []
        for rank_prefix in (['root'] + list(Taxonomy.rank_prefixes) + ['RS_', 'GB_', 'U_']):
            taxa_at_rank = []
            for taxon in ranks_below_taxon:
                if taxon.startswith(rank_prefix):
                    taxa_at_rank.append(taxon)

            sorted_taxon += sorted(taxa_at_rank)

        fout = open(results_table, 'w')
        fout.write('Taxon')
        for rank_prefix in Taxonomy.rank_prefixes:
            if count_is_scalar:
                fout.write('\t%s' % rank_prefix.capitalize())
            else:
                fout.write('\t%s\t%s\t%s\t%s' % ('Mean: ' + rank_prefix.capitalize(),
                                                 'Std: ' + rank_prefix.capitalize(),
                                                 'Min: ' + rank_prefix.capitalize(),
                                                 'Max: ' + rank_prefix.capitalize()))
        fout.write('\n')

        for taxon in sorted_taxon:
            fout.write(taxon)

            for rank_prefix in Taxonomy.rank_prefixes:
                count = ranks_below_taxon[taxon][rank_prefix.capitalize()]
                if count_is_scalar:
                    fout.write('\t%d' % count)
                else:
                    if len(count) > 0:
                        fout.write('\t%.1f\t%.2f\t%d\t%d' % (np_mean(count), np_std(count), min(count), max(count)))
                    else:
                        fout.write('\t%d\t%d\t%d\t%d' % (0, 0, 0, 0))

            fout.write('\n')

        fout.close()

    def run(self, input_tree, rd_thresholds, output_dir):
        """Calculate number of taxa for specified relative divergence thresholds.

        Parameters
        ----------
        input_tree : str
            Name of input tree.
        rd_thresholds : d[rank] -> threshold
            Relative divergence threshold for defining taxonomic ranks.
        output_dir : str
            Desired output directory.
        """

        # get list of phyla level lineages
        tree = tree = dendropy.Tree.get_from_path(input_tree,
                                                  schema='newick',
                                                  rooting='force-rooted',
                                                  preserve_underscores=True)
        phyla = get_phyla_lineages(tree)
        self.logger.info('Identified %d phyla for rooting.' % len(phyla))

        self.logger.info('Reading taxonomy from tree.')
        taxonomy_file = os.path.join(output_dir, 'taxonomy.tsv')
        taxonomy = Taxonomy().read_from_tree(input_tree)
        Taxonomy().write(taxonomy, taxonomy_file)

        rd = RelativeDistance()
        overall_ranks_below_taxon = defaultdict(lambda: defaultdict(list))
        for p in phyla:
            phylum_children = Taxonomy().children(p, taxonomy)
            phylum = p.replace('p__', '')
            self.logger.info('Calculating information with rooting on %s.' % phylum)

            phylum_dir = os.path.join(output_dir, phylum)
            if not os.path.exists(phylum_dir):
                os.makedirs(phylum_dir)

            output_tree = os.path.join(phylum_dir, 'rerooted.tree')
            os.system('genometreetk outgroup %s %s %s %s' % (input_tree, taxonomy_file, p, output_tree))

            # calculate relative distance for all nodes
            cur_tree = dendropy.Tree.get_from_path(output_tree,
                                                   schema='newick',
                                                   rooting='force-rooted',
                                                   preserve_underscores=True)
            rd.decorate_rel_dist(cur_tree)

            # determine ranks
            for n in cur_tree.postorder_node_iter(lambda n: n != tree.seed_node):
                ranks = []
                for rank_prefix, threshold in rd_thresholds.items():
                    if n.rel_dist >= threshold and n.parent_node.rel_dist < threshold:
                        ranks.append(rank_prefix.capitalize() + '__')

                if ranks:
                    if not n.label:
                        n.label = '|%s [rd=%.2f]' % (';'.join(ranks), n.rel_dist)
                    else:
                        n.label += '|%s [rd=%.2f]' % (';'.join(ranks), n.rel_dist)

            cur_tree.write_to_path(os.path.join(phylum_dir, 'rd_ranks.tree'),
                                   schema='newick',
                                   suppress_rooting=True,
                                   unquoted_underscores=True)

            # determine number of ranks below root and all named nodes
            ranks_below_taxon = defaultdict(lambda: defaultdict(int))
            for cur_node in cur_tree.postorder_node_iter():
                if cur_node == cur_tree.seed_node:
                    cur_taxon = 'root'
                elif cur_node.label:
                    _support, cur_taxon, _auxiliary_info = parse_label(cur_node.label)
                    if not cur_taxon or cur_taxon.strip() == '':
                        continue
                else:
                    continue

                for n in cur_node.postorder_iter():
                    if not n.label:
                        continue

                    _support, _taxon, auxiliary_info = parse_label(n.label)
                    if auxiliary_info:
                        ranks = auxiliary_info[0:auxiliary_info.rfind('[')]
                        ranks = [r.strip() for r in ranks.split(';')]

                        for r in ranks:
                            ranks_below_taxon[cur_taxon][r] += 1

            for taxon in ranks_below_taxon:
                if taxon == p or taxon in phylum_children:
                    # do not record results for named groups in the lineage 
                    # used for rooting
                    continue

                for rank, count in ranks_below_taxon[taxon].items():
                    overall_ranks_below_taxon[taxon][rank].append(count)

            results_table = os.path.join(phylum_dir, 'rd_ranks.tsv')
            self.write_rank_count(ranks_below_taxon, results_table)

        results_table = os.path.join(output_dir, 'mean_rd_ranks.tsv')
        self.write_rank_count(overall_ranks_below_taxon, results_table)
