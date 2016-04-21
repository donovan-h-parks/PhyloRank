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

from phylorank.rel_dist import RelativeDistance
from phylorank.newick import parse_label

from biolib.taxonomy import Taxonomy

from skbio import TreeNode

from numpy import (mean as np_mean,
                   std as np_std,
                   median as np_median,
                   abs as np_abs)


class RdRanks():
    """ Calculate number of taxa for specified relative divergence thresholds.

        <blah>
    """

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger()

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

        # midpoint root tree
        # midpoint_tree = os.path.join(output_dir, 'midpoint.tree')
        # os.system('genometreetk midpoint %s %s' % (input_tree, midpoint_tree))

        # read tree
        self.logger.info('Reading tree.')
        tree = TreeNode.read(input_tree, convert_underscores=False)

        # calculate relative distance for all nodes
        rd = RelativeDistance()
        rd.decorate_rel_dist(tree)

        # determine ranks
        for n in tree.postorder():
            if n.is_root():
                continue
                
            ranks = []
            for rank_prefix, threshold in rd_thresholds.iteritems():
                if n.rel_dist >= threshold and n.parent.rel_dist < threshold:
                    ranks.append(rank_prefix.capitalize() + '__')
                    
            if ranks:
                if not n.name:
                    n.name = '|%s [rd=%.2f]' % (';'.join(ranks), n.rel_dist)
                else:
                    n.name += '|%s [rd=%.2f]' % (';'.join(ranks), n.rel_dist)

        tree.write(os.path.join(output_dir, 'rd_ranks.tree'))
        
        # determine number of ranks below root and all named nodes
        ranks_below_taxon = defaultdict(lambda: defaultdict(int))
        for cur_node in tree.postorder():
            if cur_node.is_root():
                taxon = 'root'
            elif cur_node.name:
                _support, taxon, _auxiliary_info = parse_label(cur_node.name)
                if not taxon:
                    continue
                    
            for n in cur_node.postorder():
                if not n.name:
                    continue
                    
                _support, _taxon, auxiliary_info = parse_label(n.name)
                if auxiliary_info:
                    ranks = auxiliary_info[0:auxiliary_info.rfind('[')]
                    ranks = [r.strip() for r in ranks.split(';')]
                    
                    for r in ranks:
                        ranks_below_taxon[taxon][r] += 1
                        
        fout = open('rd_ranks.tsv', 'w')
        fout.write('Taxon')
        for rank_prefix in Taxonomy.rank_prefixes:
            fout.write('\t' + rank_prefix.capitalize())
        fout.write('\n')
            
        for taxon in ranks_below_taxon:
            fout.write(taxon)
            
            for rank_prefix in Taxonomy.rank_prefixes:
                fout.write('\t%d' % ranks_below_taxon[taxon].get(rank_prefix.capitalize(), 0))
                
        fout.close()
