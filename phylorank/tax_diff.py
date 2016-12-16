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


class TaxDiff():
    """Tabulate differences between two taxonomies.
    
    Every taxon label in the first taxonomy is classified
    to indicate how it has changed relative to the
    second taxonomy.
    
    This implies that every label in the first taxonomy is
    given precisely one label and that labels in the second 
    taxonomy may often be seen simply as deprecated. 
    
    Taxon are classified as follows in the order listed:
    
    retained: the exact same taxon label is in the GTDB and NCBI taxonomies
    corrected: only the suffix of the taxon has changed
    more specific: taxon has been moved down 1 or more ranks (e.g., phylum to class)
    more general: taxon has been moved up 1 or more ranks (e.g., class to phylum)
    new: the label has been introduced by the GTDB taxonomy
    deprecated: the label has been removed in the GTDB taxonomy
    
    E.g.,
    The first taxonomy has the label o__Archaeoglobales,
    but no label c__Archaeoglobi. The second taxonomy has
    the labels c__Archaeoglobi and o__Archaeoglobales. AssertionError
    such, o__Archaeoglobales will be classified as retained
    while c__Archaeoglobi will be classified as deprecated
    (and not 'more specific').
    """

    def __init__(self, dpi=96):
        """Initialize."""
        self.logger = logging.getLogger()
        
        # c__Halobacteria => o__Halobacteriales
        # c__Hadesarchaea => c__Hadesarchaeia
        # c__Archaeoglobi => o__Archaeoglobales
        
        # p__Pacearchaeota => f__Pacearchaeaceae
        
        # domain = 'a', phylum = 'aeota', class = 'ia', order = 'ales', family = 'aceae'
        self.rank_suffices = ['a', 'aeota', 'ia', 'ales', 'aceae', '']
        
    def _renamed(self, taxon, cur_rank, taxa_at_rank):
        """Determine if taxon has been renamed."""
        
        # remove paraphyletic suffix if one exists
        taxon = taxon[0:3] + taxon[3:].split('_')[0]
        
        # determine order of ranks to examine in order
        # to first identify 'more specific' and then
        # 'more general' label changes
        rank_order = [r for r in xrange(cur_rank-1, -1, -1)]
        rank_order += [r for r in xrange(cur_rank+1, len(Taxonomy.rank_prefixes))]
     
        # check if taxon has been moved up or down in rank
        # as determined by removing suffix characters and 
        # appending canonical rank prefixes
        for old_rank in rank_order:
            if old_rank >= len(self.rank_suffices):
                return None, None

            # eat suffix character by character and append canonical suffix
            old_rank_prefix = Taxonomy.rank_prefixes[old_rank]
            cur_suffix = self.rank_suffices[cur_rank]
            old_suffix = self.rank_suffices[old_rank]
            for i in xrange(1, max([len(s) for s in self.rank_suffices]) + 1): # eat taxon name
                if i == 0:
                    modified_taxon = taxon[3:]
                else:
                    modified_taxon = taxon[3:-i]
                    
                for j in xrange(0, len(old_suffix)): # eat suffix
                    old_taxon = old_rank_prefix + modified_taxon + old_suffix[j:]
                    if old_taxon in taxa_at_rank[old_rank]:
                        return old_taxon, old_rank
                        
                    old_taxon = old_rank_prefix + modified_taxon + old_suffix[0:j]
                    if old_taxon in taxa_at_rank[old_rank]:
                        return old_taxon, old_rank
        
        return None, None
        
    def _change_suffix(self, taxon, cur_rank, taxa_at_rank):
        """Change suffix of taxon."""
        
        # check if taxon name has been corrected
        for old_rank, old_suffix in enumerate(self.rank_suffices):
            # eat suffix character by character and append canonical suffix
            for i in xrange(0, max([len(s) for s in self.rank_suffices])): # eat taxon name
                if i == 0:
                    modified_taxon = taxon
                else:
                    modified_taxon = taxon[0:-i]
                    
                for j in xrange(0, len(old_suffix)): # eat suffix
                    old_taxon = modified_taxon + old_suffix[j:]
                    if old_taxon in taxa_at_rank[cur_rank]:
                        return old_taxon
  
                    old_taxon = modified_taxon + old_suffix[0:j]
                    if old_taxon in taxa_at_rank[cur_rank]:
                        return old_taxon
                        
        return None
        
    def _tax_diff_table(self, tax1, tax2, output_table):
        """Tabulate incongruency of taxonomy strings at each rank."""
        
        fout = open(output_table, 'w')
        fout.write('Lineage\tNo. Extent Taxa')
        for rank_label in Taxonomy.rank_labels:
            fout.write('\t%s (%%)' % rank_label.title())
        fout.write('\n')
        
        taxonomy = Taxonomy()
        named_lineages_at_rank = taxonomy.named_lineages_at_rank(tax1)
        for rank, taxa in named_lineages_at_rank.iteritems():
            rank_label = Taxonomy.rank_labels[rank]
            if rank_label == 'species':
                continue
                
            extant_taxa_for_rank = taxonomy.extant_taxa_for_rank(rank_label, tax1)
            
            for taxon in taxa:
                extent_taxa = extant_taxa_for_rank[taxon]
                fout.write('%s\t%d' % (taxon, len(extent_taxa)))
                
                row = defaultdict(list)
                for genome_id in extent_taxa:
                    taxa1 = tax1[genome_id]
                    taxa2 = tax2[genome_id]
                    
                    for cur_rank, (taxa1, taxa2) in enumerate(zip(taxa1, taxa2)):
                         row[cur_rank].append(taxa1 == taxa2)
                         
                for cur_rank, matches in row.iteritems():
                    if cur_rank <= rank:
                        fout.write('\t-')
                    else:
                        perc_match = sum(matches) * 100.0 / len(matches)
                        fout.write('\t%.1f' % (100.0 - perc_match))
                fout.write('\n')
        fout.close()
        
    def tree_tax_diff(self, tree1_file, tree2_file, output_dir):
        """Tabulate differences between two taxonomies on a tree.
        
        Parameters
        ----------
        tree1_file : str
            File with tree in Newick format.
        tree2_file : str
            File with tree in Newick format.
        output_dir : str
            Output directory.
        """
        
        tree1 = dendropy.Tree.get_from_path(tree1_file, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
                                            
        tree2 = dendropy.Tree.get_from_path(tree2_file, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
        
        # prune both trees to a set of common taxa
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
        
        # get named lineages at each taxonomic rank
        taxonomy = Taxonomy()
        tax1 = taxonomy.read_from_tree(tree1)
        tax2 = taxonomy.read_from_tree(tree2)
        
        taxa_at_rank1 = taxonomy.named_lineages_at_rank(tax1)
        taxa_at_rank2 = taxonomy.named_lineages_at_rank(tax2)

        # identify retained taxonomic names
        tax_file_name = os.path.splitext(os.path.basename(tree1_file))[0]
        output_file = os.path.join(output_dir, '%s.taxa_diff.tsv' % tax_file_name)
        fout = open(output_file, 'w')
        fout.write('Rank\tClassification\tTaxonomy 1\tTaxonomy 2\n')
        taxon2_accounted_for = defaultdict(set)
        for rank, rank_label in enumerate(Taxonomy.rank_labels[0:-1]):
            for taxon in taxa_at_rank1[rank]: 
                # check if taxon has been retained
                if taxon in taxa_at_rank2[rank]:
                    fout.write('%s\t%s\t%s\t%s\n' % (rank_label, 'retained', taxon, taxon))
                    taxon2_accounted_for[rank].add(taxon)
                    continue
                    
                # check if name was simply corrected by changing suffix
                old_taxon = self._change_suffix(taxon, rank, taxa_at_rank2)  
                if old_taxon:
                    fout.write('%s\t%s\t%s\t%s\n' % (rank_label, 'corrected', taxon, old_taxon))
                    taxon2_accounted_for[rank].add(old_taxon)
                    continue
                                         
                # check if taxon has been moved up or down in rank
                old_taxon, old_rank = self._renamed(taxon, rank, taxa_at_rank2)
                if old_taxon:
                    if rank < old_rank:
                        fout.write('%s\t%s\t%s\t%s\n' % (rank_label, 'more general', taxon, old_taxon))
                    elif rank == old_rank:
                        fout.write('%s\t%s\t%s\t%s\n' % (rank_label, 'corrected', taxon, old_taxon))
                    else:
                        fout.write('%s\t%s\t%s\t%s\n' % (rank_label, 'more specific', taxon, old_taxon))
                    
                    taxon2_accounted_for[old_rank].add(old_taxon)   
                    continue
                          
                # otherwise, the taxon appears to be new
                fout.write('%s\t%s\t%s\t%s\n' % (rank_label, 'new', taxon, 'NA'))
               
        # report deprecated taxa
        for rank, rank_label in enumerate(Taxonomy.rank_labels[0:-1]):
            for taxon in taxa_at_rank2[rank]:
                if taxon not in taxon2_accounted_for[rank]:
                    fout.write('%s\t%s\t%s\t%s\n' % (rank_label, 'deprecated', 'NA', taxon))

        fout.close()
        
        # tabulate congruence of taxonomy strings
        output_table = os.path.join(output_dir, '%s.perc_diff.tsv' % tax_file_name)
        self._tax_diff_table(tax1, tax2, output_table)
    
    def tax_diff(self, tax1_file, tax2_file, include_user_taxa, output_dir):
        """Tabulate differences between two taxonomies.
        
        Parameters
        ----------
        tax1_file : str
            First taxonomy file.
        tax2_file : str
            Second taxonomy file.
        include_user_taxa : boolean
            Flag indicating if User genomes should be considered.
        output_dir : str
            Output directory.
        """
        
        tax1 = Taxonomy().read(tax1_file)
        tax2 = Taxonomy().read(tax2_file)
        
        if not include_user_taxa:
            new_tax1 = {}
            for genome_id, taxonomy in tax1.iteritems():
                if not genome_id.startswith('U_'):
                    new_tax1[genome_id] = taxonomy
            tax1 = new_tax1
            
            new_tax2 = {}
            for genome_id, taxonomy in tax2.iteritems():
                if not genome_id.startswith('U_'):
                    new_tax2[genome_id] = taxonomy
            tax2 = new_tax2
        
        common_taxa = set(tax1.keys()).intersection(tax2.keys())
        
        self.logger.info('First taxonomy contains %d taxa.' % len(tax1))
        self.logger.info('Second taxonomy contains %d taxa.' % len(tax2))
        self.logger.info('Taxonomies have %d taxa in common.' % len(common_taxa))
        
        # identify differences between taxonomies
        tax_file_name1 = os.path.splitext(os.path.basename(tax1_file))[0]
        tax_file_name2 = os.path.splitext(os.path.basename(tax2_file))[0]
        output_table = os.path.join(output_dir, '%s.tax_diff.tsv' % tax_file_name1)
        
        fout = open(output_table, 'w')
        fout.write('Genome ID\tChange\tRank\t%s\t%s\n' % (tax_file_name1, tax_file_name2))
        
        unchanged = defaultdict(int)           # T2 = g__Bob -> T1 = g__Bob, or T2 = g__ -> T1 = g__
        active_change = defaultdict(int)       # T2 = g__Bob -> T1 = g__Jane, or T2 = g__Bob -> T1 = g__Bob_A
        passive_change = defaultdict(int)      # T2 = g__??? -> T1 = g__Jane
        unresolved_change = defaultdict(int)   # T2 = g__Box -> T1 = g__???
        for taxa in common_taxa:
            t1 = tax1[taxa]
            t2 = tax2[taxa]
            
            for rank, (taxon1, taxon2) in enumerate(zip(t1, t2)):
                if taxon1 == taxon2:
                    unchanged[rank] += 1
                elif taxon1 != Taxonomy.rank_prefixes[rank] and taxon2 != Taxonomy.rank_prefixes[rank]:
                    active_change[rank] += 1
                    fout.write('%s\t%s\t%s\t%s\t%s\n' % (taxa, 'active', Taxonomy.rank_labels[rank], ';'.join(t1), ';'.join(t2)))
                elif taxon2 == Taxonomy.rank_prefixes[rank]:
                    passive_change[rank] += 1
                    fout.write('%s\t%s\t%s\t%s\t%s\n' % (taxa, 'passive', Taxonomy.rank_labels[rank], ';'.join(t1), ';'.join(t2)))
                elif taxon1 == Taxonomy.rank_prefixes[rank]:
                    unresolved_change[rank] += 1
                    fout.write('%s\t%s\t%s\t%s\t%s\n' % (taxa, 'unresolved', Taxonomy.rank_labels[rank], ';'.join(t1), ';'.join(t2)))
                    
        fout.close()
  
        # report results
        output_table = os.path.join(output_dir, '%s.tax_diff_summary.tsv' % tax_file_name1)
        
        fout = open(output_table, 'w')
        fout.write('Rank\tUnchanged\tUnchanged (%)\tActive\t Active (%)\tPassive\tPassive (%)\tUnresolved\tUnresolved (%)\n')
        print 'Rank\tUnchanged\tActive\tPassive\tUnresolved\tTotal'
        for rank in xrange(0, len(Taxonomy.rank_prefixes)):
            total = unchanged[rank] + active_change[rank] + passive_change[rank] + unresolved_change[rank]
            if total != 0:
                fout.write('%s\t%d\t%.1f\t%d\t%.1f\t%d\t%.1f\t%d\t%.1f\n' %
                                    (Taxonomy.rank_labels[rank],
                                    unchanged[rank], unchanged[rank] * 100.0 / total,
                                    active_change[rank], active_change[rank] * 100.0 / total,
                                    passive_change[rank], passive_change[rank] * 100.0 / total,
                                    unresolved_change[rank], unresolved_change[rank] * 100.0 / total))
                print '%s\t%d\t%d\t%d\t%d\t%d' % (Taxonomy.rank_labels[rank],
                                                    unchanged[rank],
                                                    active_change[rank],
                                                    passive_change[rank],
                                                    unresolved_change[rank],
                                                    total)
                                
            
                    
                
                
            