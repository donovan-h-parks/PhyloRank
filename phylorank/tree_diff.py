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
from collections import defaultdict

from phylorank.newick import parse_label

from biolib.taxonomy import Taxonomy

import dendropy


class TreeDiff():
    """Quantify supported topological differences between trees.
    
    
    """

    def __init__(self, dpi=96):
        """Initialize."""
        self.logger = logging.getLogger()
  
  
    def run(self, tree1_file, tree2_file, output_file, min_support, min_taxa, named_only):
        """Calculate supported topological differences between trees.
        
        Parameters
        ----------
        tree1_file : str
            File with tree in Newick format.
        tree2_file : str
            File with tree in Newick format.
        output_file : str
            Output file.
        min_support : float
            Minimum value to consider a lineage well supported.
        min_taxa : int
            Only consider lineage with sufficient number of taxa.
        named_only : boolean
            Only consider named lineages.  
        """
        
        if not named_only:
            self.logger.error("This command currently assumes the 'named_only' flag will be thrown.")
            sys.exit()
        
        tree1 = dendropy.Tree.get_from_path(tree1_file, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
                                            
        tree2 = dendropy.Tree.get_from_path(tree2_file, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        
        # prune both trees to the set of common taxa
        taxa1 = set()
        for t in tree1.leaf_node_iter():
            taxa1.add(t.taxon.label)
            
        taxa2 = set()
        for t in tree2.leaf_node_iter():
            taxa2.add(t.taxon.label)
            
        taxa_in_common = taxa1.intersection(taxa2)
        self.logger.info('Tree 1 contains %d taxa.' % len(taxa1))
        self.logger.info('Tree 2 contains %d taxa.' % len(taxa2))
        self.logger.info('Pruning trees to the %d taxa in common.' % len(taxa_in_common))
        
        tree1.retain_taxa_with_labels(taxa_in_common)
        tree2.retain_taxa_with_labels(taxa_in_common)
        
        # identify nodes meeting specified criteria
        tree1_nodes = {}
        tree2_nodes = {}
        for tree, tree_nodes in ([tree1, tree1_nodes],[tree2, tree2_nodes]):
            for n in tree.preorder_internal_node_iter():
                support, taxon_name, _auxiliary_info = parse_label(n.label)
                if named_only and not taxon_name:
                    continue
                    
                num_taxa = sum([1 for _ in n.leaf_iter()])
                if support >= min_support and num_taxa >= min_taxa:
                    tree_nodes[taxon_name] = [support, num_taxa, n]
                    
        self.logger.info('Tree 1 has %d supported nodes.' % len(tree1_nodes))
        self.logger.info('Tree 2 has %d supported nodes.' % len(tree2_nodes))
        
        # identify supported nodes the differ in the two trees
        diffs = {}
        for taxon, data1 in tree1_nodes.iteritems():
            support1, num_taxa1, node1 = data1
            
            if taxon in tree2_nodes:
                support2, num_taxa2, node2 = tree2_nodes[taxon]
                
                taxa1 = set([t.taxon.label for t in node1.leaf_iter()])
                taxa2 = set([t.taxon.label for t in node2.leaf_iter()])
                
                diff_taxa = taxa1.symmetric_difference(taxa2)
                
                if len(diff_taxa) > 0:
                    diffs[taxon] = [len(diff_taxa), ','.join(taxa1 - taxa2), ','.join(taxa2- taxa1)]
                    
        fout = open(output_file, 'w')
        fout.write('Taxon\tNo. Incongruent Taxa\tTree1 - Tree2\tTree2 - Tree1\n')
        for taxon in Taxonomy().sort_taxa(diffs.keys()):
            num_diffs, t12_diff_str, t21_diff_str = diffs[taxon]
            fout.write('%s\t%d\t%s\t%s\n' % (taxon,
                                                num_diffs,
                                                t12_diff_str,
                                                t21_diff_str))
        
        fout.close()
        