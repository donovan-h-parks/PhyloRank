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
from collections import defaultdict, namedtuple

from phylorank.newick import parse_label
from phylorank.common import get_phyla_lineages

from biolib.taxonomy import Taxonomy

from skbio import TreeNode

from numpy import (mean as np_mean,
                   std as np_std,
                   percentile as np_percentile)


class BranchLengthDistribution():
    """Calculate distribution of branch lengths at each taxonomic rank."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger()

    def run(self, input_tree, output_dir):
        """Calculate distribution of branch lengths at each taxonomic rank.

        Parameters
        ----------
        input_tree : str
            Name of input tree.
        output_dir : str
            Desired output directory.
        """

        # get list of phyla level lineages
        tree = TreeNode.read(input_tree, convert_underscores=False)
       
        # determine ranks
        rank_bl_dist = defaultdict(list)
        taxa_bl_dist = defaultdict(list)
        for node in tree.postorder():
            if node.is_tip() or not node.name:
                continue
                
            _support, taxon, _auxiliary_info = parse_label(node.name)
            if not taxon:
                continue
                
            # get most specific rank in multi-rank taxa string
            taxon = [t.strip() for t in taxon.split(';')][-1]
                
            for t in node.tips():
                bl = t.accumulate_to_ancestor(node)
                taxa_bl_dist[taxon].append(bl)
                
                if bl > 0:
                    # ignore zero length branches as these are indicative
                    # of a species with only identical (or extremely similar) genomes
                    rank = Taxonomy.rank_labels[Taxonomy.rank_index[taxon[0:3]]]
                    if rank != 'species' or Taxonomy().validate_species_name(taxon):
                        rank_bl_dist[rank].append(bl)
                    
        # report results sorted by rank
        sorted_taxon = []
        for rank_prefix in Taxonomy.rank_prefixes:
            taxa_at_rank = []
            for taxon in taxa_bl_dist:
                if taxon.startswith(rank_prefix):
                    taxa_at_rank.append(taxon)
                    
            sorted_taxon += sorted(taxa_at_rank)
                
        # report results for each named group
        taxa_file = os.path.join(output_dir, 'taxa_bl_dist.tsv')
        fout = open(taxa_file, 'w')
        fout.write('Taxa\tMean\tStd\t5th\t10th\t50th\t90th\t95\n')
        for taxon in sorted_taxon:
            dist = taxa_bl_dist[taxon]
            p = np_percentile(dist, [5, 10, 50, 90, 95])
            fout.write('%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (taxon,
                                                                np_mean(dist),
                                                                np_std(dist),
                                                                p[0], p[1], p[2], p[3], p[4]))
        fout.close()
        
        # report results for each taxonomic rank
        rank_file = os.path.join(output_dir, 'rank_bl_dist.tsv')
        fout = open(rank_file, 'w')
        fout.write('Rank\tMean\tStd\t5th\t10th\t50th\t90th\t95\n')
        for rank in Taxonomy.rank_labels:
            dist = rank_bl_dist[rank]
            p = np_percentile(dist, [5, 10, 50, 90, 95])
            fout.write('%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (rank,
                                                                np_mean(dist),
                                                                np_std(dist),
                                                                p[0], p[1], p[2], p[3], p[4]))
        fout.close()
        