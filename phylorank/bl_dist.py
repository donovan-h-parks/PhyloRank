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
from phylorank.common import (get_phyla_lineages, filter_taxa_for_dist_inference)

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

    def run(self, input_tree, trusted_taxa_file, min_children, taxonomy_file, output_dir):
        """Calculate distribution of branch lengths at each taxonomic rank.

        Parameters
        ----------
        input_tree : str
            Name of input tree.
        trusted_taxa_file : str
            File specifying trusted taxa to consider when inferring distribution. Set to None to consider all taxa.
        min_children : int
            Only consider taxa with at least the specified number of children taxa when inferring distribution.
        taxonomy_file : str
            File containing taxonomic information for leaf nodes (if NULL, read taxonomy from tree).
        output_dir : str
            Desired output directory.
        """

        # get list of phyla level lineages
        tree = TreeNode.read(input_tree, convert_underscores=False)
        
        # pull taxonomy from tree
        if not taxonomy_file:
            self.logger.info('Reading taxonomy from tree.')
            taxonomy_file = os.path.join(output_dir, 'taxonomy.tsv')
            taxonomy = Taxonomy().read_from_tree(input_tree)
            Taxonomy().write(taxonomy, taxonomy_file)
        else:
            self.logger.info('Reading taxonomy from file.')
            taxonomy = Taxonomy().read(taxonomy_file)
            
        # read trusted taxa
        trusted_taxa = None
        if trusted_taxa_file:
            trusted_taxa = read_taxa_file(trusted_taxa_file)
        
        # determine taxa to be used for inferring distribution
        taxa_for_dist_inference = filter_taxa_for_dist_inference(tree, taxonomy, set(), min_children, -1)
        
        # determine branch lengths to leaves for named lineages
        rank_bl_dist = defaultdict(list)
        taxa_bl_dist = defaultdict(list)
        taxa_at_rank = defaultdict(list)
        for node in tree.postorder():
            if node.is_tip() or not node.name:
                continue
                
            _support, taxon, _auxiliary_info = parse_label(node.name)
            if not taxon:
                continue
                
            # get most specific rank in multi-rank taxa string
            taxon = [t.strip() for t in taxon.split(';')][-1]
            
            most_specific_rank = taxon[0:3]
            taxa_at_rank[Taxonomy.rank_index[most_specific_rank]].append(taxon)
                
            for t in node.tips():
                bl = t.accumulate_to_ancestor(node)
                taxa_bl_dist[taxon].append(bl)
                
                if bl > 0:
                    # ignore zero length branches as these are indicative
                    # of a species with only identical (or extremely similar) genomes
                    rank = Taxonomy.rank_labels[Taxonomy.rank_index[taxon[0:3]]]
                    if rank != 'species' or Taxonomy().validate_species_name(taxon):
                        if taxon in taxa_for_dist_inference:
                            rank_bl_dist[rank].append(bl)
                            
        # report number of taxa at each rank
        print ''
        print 'Rank\tTaxa\tTaxa for Inference'
        for rank, taxa in taxa_at_rank.iteritems():
            taxa_for_inference = [x for x in taxa if x in taxa_for_dist_inference]
            print '%s\t%d\t%d' % (Taxonomy.rank_labels[rank], len(taxa), len(taxa_for_inference))
        print ''
                    
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
        fout.write('Taxa\tUsed for Inference\tMean\tStd\t5th\t10th\t50th\t90th\t95th\n')
        for taxon in sorted_taxon:
            dist = taxa_bl_dist[taxon]
            p = np_percentile(dist, [5, 10, 50, 90, 95])
            fout.write('%s\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (taxon,
                                                                str(taxon in taxa_for_dist_inference),
                                                                np_mean(dist),
                                                                np_std(dist),
                                                                p[0], p[1], p[2], p[3], p[4]))
        fout.close()
        
        # report results for each taxonomic rank
        rank_file = os.path.join(output_dir, 'rank_bl_dist.tsv')
        fout = open(rank_file, 'w')
        fout.write('Rank\tMean\tStd\t5th\t10th\t50th\t90th\t95th\n')
        for rank in Taxonomy.rank_labels:
            dist = rank_bl_dist[rank]
            p = np_percentile(dist, [5, 10, 50, 90, 95])
            fout.write('%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (rank,
                                                                np_mean(dist),
                                                                np_std(dist),
                                                                p[0], p[1], p[2], p[3], p[4]))
        fout.close()
        