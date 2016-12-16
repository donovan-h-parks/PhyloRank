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

import sys
import logging

import dendropy
from numpy import (median as np_median)

from biolib.taxonomy import Taxonomy
from biolib.newick import parse_label, create_label

from phylorank.common import (read_taxa_file,
                              filter_taxa_for_dist_inference)
from phylorank.outliers import Outliers


class Decorate():
    """Decorate internal nodes with taxa labels."""

    def __init__(self):
        """Initialize."""
        
        self.logger = logging.getLogger()
        
    def _fmeasure(self, tree, taxonomy):
        """Find node with highest F-measure for each taxon.
        
        Finds best placement for each taxon label
        by calculating the F-measure for every taxon
        at every node.
        
        Parameters
        ----------
        tree : Tree
          Dendropy Tree.
        taxonomy : d[extent_taxon_id] -> taxa list
          Taxon labels for extant taxa.
          
        Returns
        -------
        d[taxon] -> [(Node, F-measure, precision, recall_, ...]
            Node(s) with highest F-measure for each taxon.
        """
    
        # get named lineages/taxa at each taxonomic rank
        taxa_at_rank = Taxonomy().named_lineages_at_rank(taxonomy)
        
        # get extant taxa for each taxon label
        extent_taxa_with_label = {}
        for i, rank in enumerate(Taxonomy.rank_labels):
            extent_taxa_with_label[i] = Taxonomy().extant_taxa_for_rank(rank, taxonomy)
            
        # get parent taxon for each taxon:
        taxon_parents = Taxonomy().parents(taxonomy)
    
        # find node with best F-measure for each taxon
        fmeasure_for_taxa = {}
        for rank_index in xrange(0, len(Taxonomy.rank_labels)):
            self.logger.info('Processing %d taxa at %s rank.' % (len(taxa_at_rank[rank_index]),
                                                                    Taxonomy.rank_labels[rank_index].capitalize()))
            
            for taxon in taxa_at_rank[rank_index]:
                if rank_index == 0:
                    # processing taxa at the domain is a special case
                    taxon_parent_node = tree.seed_node
                else:
                    parent_taxon = taxon_parents[taxon][-1]
                    if parent_taxon in fmeasure_for_taxa:
                        # only need to process the lineage below the parent node
                        taxon_parent_node = fmeasure_for_taxa[parent_taxon][0][0]
                    else:
                        # the parent for this taxon was not placed so
                        # it can be ignored (e.g., bacterial phylum in archaeal tree)
                        continue
                    
                cur_taxon_fmeasure = -1
                for node in taxon_parent_node.preorder_iter():
                    total_taxa = len(extent_taxa_with_label[rank_index][taxon])
                    
                    taxa_in_lineage = 0
                    num_leaves_with_taxa = 0
                    for leaf in node.leaf_iter():
                        leaf_taxon = taxonomy[leaf.taxon.label][rank_index]
                        if leaf_taxon == taxon:
                            taxa_in_lineage += 1
                        if leaf_taxon != Taxonomy.rank_prefixes[rank_index]:
                            num_leaves_with_taxa += 1
                    
                    if taxa_in_lineage != 0 and num_leaves_with_taxa != 0:
                        precision = float(taxa_in_lineage) / num_leaves_with_taxa
                        recall = float(taxa_in_lineage) / total_taxa
                        fmeasure = (2*precision*recall) / (precision + recall)

                        if fmeasure > cur_taxon_fmeasure:
                            cur_taxon_fmeasure = fmeasure
                            fmeasure_for_taxa[taxon] = [(node, fmeasure, precision, recall)]
                        elif fmeasure == cur_taxon_fmeasure:
                            fmeasure_for_taxa[taxon].append((node, fmeasure, precision, recall))
                          
        return fmeasure_for_taxa
        
    def _strip_taxon_labels(self, tree):
        """Remove any previous taxon labels.
        
        Parameters
        ----------
        tree : Tree
          Dendropy Tree.
        """
        
        for node in tree.internal_nodes():
            support, _taxon, _aux_info = parse_label(node.label)
            if support:
                node.label = create_label(support, None, None)
                    
    def _assign_taxon_labels(self, fmeasure_for_taxa):
        """Assign taxon labels to nodes.
        
        Parameters
        ----------
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall), ...]
          Node with highest F-measure for each taxon.
          
        Returns
        -------
        set
            Taxon labels placed in tree.
        """

        placed_taxon = set()
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys()):
            if len(fmeasure_for_taxa[taxon]) == 1:
                placed_taxon.add(taxon)
                node, fmeasure, precision, recall = fmeasure_for_taxa[taxon][0]
                
                support, taxon_label, aux_info = parse_label(node.label)
                if taxon_label:
                    taxon_label += ';' + taxon
                else:
                    taxon_label = taxon
                node.label = create_label(support, taxon_label, aux_info)

        return placed_taxon
        
    def _write_statistics_table(self, fmeasure_for_taxa, out_table):
        """Write table containing statistics for each taxon.
        
        Parameters
        ----------
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall)]
          Node with highest F-measure for each taxon.
        out_table : str
          Output table to write statistics for assigned labels.  
        """
    
        fout_table = open(out_table, 'w')
        fout_table.write('Taxon\tF-measure\tPrecision\tRecall\n')
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys()):
            if len(fmeasure_for_taxa[taxon]) != 1:
                self.logger.error('Multiple positions specified for taxon label.')
                sys.exit()
                
            node, fmeasure, precision, recall = fmeasure_for_taxa[taxon][0]
            fout_table.write('%s\t%.2f\t%.2f\t%.2f\n' % (taxon, fmeasure, precision, recall))
                
        fout_table.close()
        
    def _write_taxonomy(self, tree, out_taxonomy):
        """Write taxonomy decorated on tree to file.
        
        Parameters
        ----------
        tree : Tree
          Dendropy Tree.
        out_taxonomy : str
          Output file.
        """
        
        fout = open(out_taxonomy, 'w')
        for leaf in tree.leaf_node_iter():
            leaf_taxa = []
            
            parent = leaf
            while parent:
                _support, taxon, _aux_info = parse_label(parent.label)
                
                if taxon:
                    for t in taxon.split(';')[::-1]:
                        leaf_taxa.append(t)
                        
                parent = parent.parent_node 
                    
            ordered_taxa = leaf_taxa[::-1]
            filled_ordered_taxa = Taxonomy().fill_missing_ranks(ordered_taxa)
                  
            fout.write('%s\t%s\n' % (leaf.taxon.label, ';'.join(filled_ordered_taxa)))
        
        fout.close()
        
    def _median_rank_rd(self, 
                            tree, 
                            placed_taxon, 
                            taxonomy,
                            trusted_taxa_file, 
                            min_children, 
                            min_support):
        """Calculate median relative divergence to each node and thresholds for each taxonomic rank.
        
        Parameters
        ----------
        tree : Tree
          Dendropy Tree.
        placed_taxon : set
          Taxon currently placed in tree which can be used for relative divergence inference.
        taxonomy: d[taxon_id] -> taxonomy info
          Taxonomic information for extant taxa.
        trusted_taxa_file : str
          File specifying trusted taxa to consider when inferring distribution. Set to None to consider all taxa.
        min_children : int
          Only consider taxa with at least the specified number of children taxa when inferring distribution.
        min_support : float
          Only consider taxa with at least this level of support when inferring distribution.
        
        Returns
        -------
        d[rank_index] -> float
          Median relative divergence for each taxonomic rank.
        """
                      
        # read trusted taxa
        trusted_taxa = None
        if trusted_taxa_file:
            trusted_taxa = read_taxa_file(trusted_taxa_file)
            
        # determine taxa to be used for inferring distribution
        taxa_for_dist_inference = filter_taxa_for_dist_inference(tree, 
                                                                    taxonomy, 
                                                                    trusted_taxa, 
                                                                    min_children, 
                                                                    min_support)
        taxa_for_dist_inference.intersection_update(placed_taxon)
 
        # infer distribution                                        
        outliers = Outliers()
        phylum_rel_dists, rel_node_dists = outliers.median_rd_over_phyla(tree, 
                                                                            taxa_for_dist_inference, 
                                                                            taxonomy)    
        median_for_rank = outliers.rank_median_rd(phylum_rel_dists, 
                                                    taxa_for_dist_inference)
                                                    
        # set edge lengths to median value over all rootings
        tree.seed_node.rel_dist = 0.0
        for n in tree.preorder_node_iter(lambda n: n != tree.seed_node):
            n.rel_dist = np_median(rel_node_dists[n.id])
            
        return median_for_rank
        
    def _resolve_ambiguous_placements(self, tree, fmeasure_for_taxa, median_rank_rd):
        """Resolve ambiguous taxon label placements using median relative divergences.
        
        Parameters
        ----------
        tree : Tree
          Dendropy tree.
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall)]
          Node with highest F-measure for each taxon.
        median_rank_rd : d[rank_index] -> float
          Median relative divergence for each taxonomic rank.
        """
        
        # For ambiguous nodes place them closest to median for rank 
        # and within accept relative divergence distance. Taxon labels
        # are placed in reverse taxonomic order (species to domain) and
        # this ordering used to ensure taxonomic consistency.
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys(), reverse=True):
            if len(fmeasure_for_taxa[taxon]) == 1:
                continue
                                
            rank_prefix = taxon[0:3]
            rank_index = Taxonomy.rank_prefixes.index(rank_prefix)
            rd = median_rank_rd[rank_index]

            # Find node closest to median distance, but making sure
            # taxon is not placed below a more specific taxon label.
            # The 'fmeasure_for_taxa' stores node information in preorder.
            closest_index = None
            closest_dist = 1e9
            closest_node = None
            for i, d in enumerate(fmeasure_for_taxa[taxon]):
                cur_node = d[0]    

                cur_rank_index = -1
                _support, cur_taxon, _aux_info = parse_label(cur_node.label)
                if cur_taxon:
                    cur_prefix = cur_taxon.split(';')[-1][0:3]
                    cur_rank_index = Taxonomy.rank_prefixes.index(cur_prefix)
                    
                if cur_rank_index > rank_index:
                    # reached a node with a more specific label so
                    # label should be appended to this node or
                    # placed above it
                    if closest_index is None:
                        closest_index = i
                        closest_node = cur_node
                    break
                    
                rd_diff = abs(rd - cur_node.rel_dist)
                if rd_diff < 0.1 and rd_diff < closest_dist:
                    closest_dist = rd_diff
                    closest_index = i
                    closest_node = cur_node
                    
            if closest_index is None:
                # no node is within an acceptable relative divergence distance 
                # for this label so it should be placed at the most extant node
                # which should be a leaf node
                closest_index = len(fmeasure_for_taxa[taxon]) - 1
                closest_node = fmeasure_for_taxa[taxon][closest_index][0]
                
                if not closest_node.is_leaf():
                    self.logger.error('Leaf node expected!')
                    sys.exit()
                    
            # add label to node
            support, cur_taxon, aux_info = parse_label(closest_node.label)
            if not cur_taxon:
                taxa_str = taxon
            else:
                taxa = cur_taxon.split(';') + [taxon]
                taxa_str = ';'.join(Taxonomy().sort_taxa(taxa))
                
            closest_node.label = create_label(support, taxa_str, aux_info)
                    
            # remove other potential node assignments
            fmeasure_for_taxa[taxon] = [fmeasure_for_taxa[taxon][closest_index]]
                    
    def run(self, 
                input_tree, 
                taxonomy_file, 
                trusted_taxa_file, 
                min_children, 
                min_support,
                output_tree):
        """Decorate internal nodes with taxa labels.

        Parameters
        ----------
        input_tree : str
          Tree to decorate
        taxonomy_file : str
          File indicating taxonomic information for extant taxa.
        trusted_taxa_file : str
          File specifying trusted taxa to consider when inferring distribution. Set to None to consider all taxa.
        min_children : int
          Only consider taxa with at least the specified number of children taxa when inferring distribution.
        min_support : float
          Only consider taxa with at least this level of support when inferring distribution.
        output_tree: str
          Name of output tree.
        """
        
        # read tree
        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
                                            
        # remove any previous taxon labels
        self.logger.info('Removing any previous internal node labels.')
        self._strip_taxon_labels(tree)
                                   
        # read taxonomy and trim to taxa in tree
        self.logger.info('Reading taxonomy.')
        full_taxonomy = Taxonomy().read(taxonomy_file)
        
        taxonomy = {}
        for leaf in tree.leaf_node_iter():
            taxonomy[leaf.taxon.label] = full_taxonomy.get(leaf.taxon.label, Taxonomy.rank_prefixes)

        # find best placement for each taxon based 
        # on the F-measure statistic
        self.logger.info('Calculating F-measure statistic for each taxa.')
        fmeasure_for_taxa = self._fmeasure(tree, taxonomy)
        
        # place labels with only one acceptable position and calculate
        # the relative divergence thresholds from these as a guide for
        # placing the remaining labels
        self.logger.info('Placing labels with unambiguous position in tree.')
        placed_taxon = self._assign_taxon_labels(fmeasure_for_taxa)

        # calculating relative 
        self.logger.info('Establishing median relative divergence for taxonomic ranks.')
        median_rank_rd = self._median_rank_rd(tree, 
                                                placed_taxon, 
                                                taxonomy,
                                                trusted_taxa_file, 
                                                min_children, 
                                                min_support)
                                                                                      
        # resolve ambiguous position in tree
        self.logger.info('Resolving ambiguous taxon label placements using median relative divergences.')
        self._resolve_ambiguous_placements(tree, fmeasure_for_taxa, median_rank_rd)
       
        # write statistics for placed taxon labels
        self.logger.info('Writing out statistics for taxa.')
        out_table = output_tree + '-table'
        self._write_statistics_table(fmeasure_for_taxa, out_table)
                                          
        # output taxonomy of extant taxa on tree
        self.logger.info('Writing out taxonomy for extant taxa.')
        out_taxonomy = output_tree + '-taxonomy'
        self._write_taxonomy(tree, out_taxonomy)
        
        # validate taxonomy
        self.logger.info('Validating taxonomy for extant taxa.')
        tree_taxonomy = Taxonomy().read(out_taxonomy)
        Taxonomy().validate(tree_taxonomy,
                          check_prefixes=True,
                          check_ranks=True,
                          check_hierarchy=True,
                          check_species=True,
                          check_group_names=True,
                          check_duplicate_names=True,
                          report_errors=True)
                                                                                  
        # output decorated tree
        self.logger.info('Writing out decorated tree.')
        tree.write_to_path(output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)
