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
  
  
    def run(self, tree1_file, tree2_file, output_dir, min_support, min_taxa, named_only):
        """Calculate supported topological differences between trees.
        
        Parameters
        ----------
        tree1_file : str
            File with tree in Newick format.
        tree2_file : str
            File with tree in Newick format.
        output_dir : str
            Output directory.
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
            
        tree1_name = os.path.splitext(os.path.basename(tree1_file))[0]
        tree2_name = os.path.splitext(os.path.basename(tree2_file))[0]
        
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
        node_support1 = {}
        node_support2 = {}
        for tree, tree_nodes, support_values in ([tree1, tree1_nodes, node_support1],[tree2, tree2_nodes, node_support2]):
            for n in tree.preorder_internal_node_iter():
                support, taxon_name, _auxiliary_info = parse_label(n.label)
                if named_only and not taxon_name:
                    continue
                    
                if not support:
                    continue
                    
                support = int(support)
                support_values[taxon_name] = support
                
                num_taxa = sum([1 for _ in n.leaf_iter()])
                if support >= min_support and num_taxa >= min_taxa:
                    tree_nodes[taxon_name] = [support, num_taxa, n]
                    
        self.logger.info('Tree 1 has %d supported nodes.' % len(tree1_nodes))
        self.logger.info('Tree 2 has %d supported nodes.' % len(tree2_nodes))
        
        # compare supported nodes between the two trees
        diffs = {}
        congruent_taxa = defaultdict(list)       # same node bootstrap supported in both trees
        incongruent_taxa = defaultdict(list)     # node supported in both trees, but have different extant taxa
        unresolved_taxa = defaultdict(list)      # supported node in one tree is not present and/or well support in the other tree

        for taxon, data1 in tree1_nodes.iteritems():
            most_specific_taxon = taxon.split(';')[-1].strip()
            rank_index = Taxonomy.rank_prefixes.index(most_specific_taxon[0:3])
            support1, num_taxa1, node1 = data1
            
            if taxon in tree2_nodes:
                support2, num_taxa2, node2 = tree2_nodes[taxon]
                
                taxa1 = set([t.taxon.label for t in node1.leaf_iter()])
                taxa2 = set([t.taxon.label for t in node2.leaf_iter()])
                
                diff_taxa = taxa1.symmetric_difference(taxa2)
                
                if len(diff_taxa) > 0:
                    diffs[taxon] = [len(diff_taxa), ','.join(taxa1 - taxa2), ','.join(taxa2- taxa1)]
                    incongruent_taxa[rank_index].append((taxon, len(diff_taxa)))
                else:
                    congruent_taxa[rank_index].append((taxon, support1, support2))
            else:
                unresolved_taxa[rank_index].append((taxon, tree1_name, support1, tree2_name, node_support2.get(taxon, -1)))
                
        # identify unresolved taxa in tree 2
        for taxon, data2 in tree2_nodes.iteritems():
            support2, num_taxa2, node2 = data1
            if taxon not in tree1_nodes:
                unresolved_taxa[rank_index].append((taxon, tree2_name, support2, tree1_name, node_support1.get(taxon, -1)))
        
        # write out difference in extant taxa for incongruent taxa
        tax_diff_file = os.path.join(output_dir, 'incongruent_taxa.tsv')
        fout = open(tax_diff_file, 'w')
        fout.write('Taxon\tNo. Incongruent Taxa\tTree1 - Tree2\tTree2 - Tree1\n')
        for taxon in Taxonomy().sort_taxa(diffs.keys()):
            num_diffs, t12_diff_str, t21_diff_str = diffs[taxon]
            fout.write('%s\t%d\t%s\t%s\n' % (taxon,
                                                num_diffs,
                                                t12_diff_str,
                                                t21_diff_str))
        
        fout.close()
        
        # write out classification of each node
        classification_file = os.path.join(output_dir, 'taxon_classification.tsv')
        fout_classification = open(classification_file, 'w')
        fout_classification.write('Rank\tTaxon\tClassification\tDescription\n')
        
        stats_file = os.path.join(output_dir, 'tree_diff_stats.tsv')
        fout_stats = open(stats_file, 'w')
        fout_stats.write('Rank\tCongruent\tIncongruent\tUnresolved for %s\tUnresolved for %s\n' % (tree1_name, tree2_name))
        for rank, rank_label in enumerate(Taxonomy.rank_labels):
            for info in congruent_taxa[rank]:
                taxon, support1, support2 = info
                
                desc = 'Taxon is congruent with %d and %d support.' % (support1, support2)
                fout_classification.write('%s\t%s\t%s\t%s\n' % (rank_label, taxon, 'congruent', desc))
                
            for info in incongruent_taxa[rank]:
                taxon, num_diff_taxa = info
                desc = 'Taxon has %d extant taxa in disagreement.' % num_diff_taxa
                fout_classification.write('%s\t%s\t%s\t%s\n' % (rank_label, taxon, 'incongruent', desc))
                
            unresolved1 = 0
            unresolved2 = 0
            for info in unresolved_taxa[rank]:
                taxon, supported_tree_name, support1, unsupported_tree_name, support2 = info
                desc = 'Taxon is supported in %s (%d), but not in %s (%d)' % (supported_tree_name, support1, unsupported_tree_name, support2)
                fout_classification.write('%s\t%s\t%s\t%s\n' % (rank_label, taxon, 'incongruent', desc))
                
                if supported_tree_name == tree1_name:
                    unresolved1 += 1
                else:
                    unresolved2 += 1
                
            fout_stats.write('%s\t%d\t%d\t%s\t%s\n' % (rank_label, 
                                                        len(congruent_taxa[rank]),
                                                        len(incongruent_taxa[rank]), 
                                                        unresolved1,
                                                        unresolved2))
                
        fout_classification.close()
        fout_stats.close()