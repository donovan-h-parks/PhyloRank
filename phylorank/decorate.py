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
import sys
import time
from collections import defaultdict, namedtuple

import dendropy
from biolib.newick import parse_label, create_label
from biolib.taxonomy import Taxonomy
from numpy import (median as np_median)

from phylorank.common import (read_taxa_file,
                              filter_taxa_for_dist_inference)
from phylorank.outliers import Outliers
from phylorank.viral_taxonomy import (translate_viral_taxonomy,
                                      rev_translate_output_file)


class Decorate(object):
    """Decorate internal nodes with taxa labels."""

    def __init__(self):
        """Initialize."""

        self.logger = logging.getLogger()

        self.StatsTable = namedtuple('StatsTable',
                                     'node fmeasure precision recall taxa_in_lineage total_taxa num_leaves_with_taxa rogue_out rogue_in')

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

        # get number of leaves and taxon in each lineage
        self.logger.info('Calculating taxa within each lineage.')
        for node in tree.preorder_node_iter():
            gids = set()
            num_leaves = 0
            taxa_count = defaultdict(lambda: defaultdict(int))
            for leaf in node.leaf_iter():
                gids.add(leaf.taxon.label)
                num_leaves += 1
                for rank_index, taxon in enumerate(taxonomy[leaf.taxon.label]):
                    if taxon != Taxonomy.rank_prefixes[rank_index]:
                        taxa_count[rank_index][taxon] += 1

            node.num_leaves = num_leaves
            node.taxa_count = taxa_count
            node.descendant_gids = gids

        taxa_in_tree = defaultdict(int)
        for leaf in tree.leaf_node_iter():
            for taxon in taxonomy[leaf.taxon.label]:
                taxa_in_tree[taxon] += 1

        # find node with best F-measure for each taxon
        fmeasure_for_taxa = {}
        for rank_index in range(0, len(Taxonomy.rank_labels)):
            start = time.time()
            self.logger.info('Processing {:,} taxa at {} rank.'.format(
                len(taxa_at_rank[rank_index]),
                Taxonomy.rank_labels[rank_index].capitalize()))

            for idx, taxon in enumerate(taxa_at_rank[rank_index]):
                if rank_index == 0:
                    # processing taxa at the domain is a special case
                    taxon_parent_node = tree.seed_node
                else:
                    # find first named parent 
                    # e.g., Cyanobacteria for Synechococcales in d__Bacteria;p__Cyanobacteria;c__;o__Synechococcales
                    parent_taxon = 'x__'
                    parent_index = rank_index - 1
                    while len(parent_taxon) == 3 and parent_index != -1:
                        parent_taxon = taxon_parents[taxon][parent_index]
                        parent_index = parent_index - 1

                    if parent_taxon in fmeasure_for_taxa:
                        # only need to process the lineage below the parent node,
                        # but must take the MRCA if the placement of the parent
                        # taxon is unresolved
                        parent_nodes = []
                        for stat_table in fmeasure_for_taxa[parent_taxon]:
                            parent_nodes.append(stat_table.node)

                        if len(parent_nodes) == 1:
                            taxon_parent_node = parent_nodes[0]
                        else:
                            taxa = []
                            for p in parent_nodes:
                                taxa += [leaf.taxon for leaf in p.leaf_iter()]
                            taxon_parent_node = tree.mrca(taxa=taxa)

                        if taxon_parent_node.taxa_count[rank_index][taxon] < 0.5 * taxa_in_tree[taxon]:
                            # substantial portion of genomes for taxon fall outside 
                            # the parent lineages so best search the entire tree
                            taxon_parent_node = tree.seed_node
                    else:
                        # the parent for this taxon was not placed so
                        # it can be ignored (e.g., bacterial phylum in archaeal tree)
                        continue

                cur_taxon_fmeasure = -1
                cur_gids = set(extent_taxa_with_label[rank_index][taxon])

                for node in taxon_parent_node.preorder_iter():
                    taxa_in_lineage = node.taxa_count[rank_index][taxon]
                    num_leaves_with_taxa = sum(node.taxa_count[rank_index].values())

                    if taxa_in_lineage != 0 and num_leaves_with_taxa != 0:
                        precision = float(taxa_in_lineage) / num_leaves_with_taxa
                        recall = float(taxa_in_lineage) / len(cur_gids)
                        fmeasure = (2 * precision * recall) / (precision + recall)

                        if fmeasure >= cur_taxon_fmeasure:
                            descendant_gids = node.descendant_gids
                            rogue_out = cur_gids - descendant_gids
                            rogue_in = [gid
                                        for gid in descendant_gids - cur_gids 
                                        if taxonomy[gid][rank_index] != Taxonomy.rank_prefixes[rank_index]]

                            if fmeasure >= cur_taxon_fmeasure:
                                stat_table = self.StatsTable(node=node,
                                                             fmeasure=fmeasure,
                                                             precision=precision,
                                                             recall=recall,
                                                             taxa_in_lineage=taxa_in_lineage,
                                                             total_taxa=len(cur_gids),
                                                             num_leaves_with_taxa=num_leaves_with_taxa,
                                                             rogue_out=rogue_out,
                                                             rogue_in=rogue_in)

                                if fmeasure > cur_taxon_fmeasure:
                                    cur_taxon_fmeasure = fmeasure
                                    fmeasure_for_taxa[taxon] = [stat_table]
                                elif fmeasure == cur_taxon_fmeasure:
                                    fmeasure_for_taxa[taxon].append(stat_table)

                statusStr = ' - processed {:,} of {:,} ({:.2f}%) taxa.'.format(
                                    idx+1, 
                                    len(taxa_at_rank[rank_index]), 
                                    float(idx*100)/len(taxa_at_rank[rank_index])).ljust(86)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            sys.stdout.write('\n')
            end = time.time()
            print(' - elapsed time: {:.2f} minutes'.format((end-start)/60.0))

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
            if support is not None:
                node.label = create_label(support, None, None)
            else:
                node.label = None

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

                stat_table = fmeasure_for_taxa[taxon][0]
                node = stat_table.node
                fmeasure = stat_table.fmeasure
                precision = stat_table.precision
                recall = stat_table.recall

                support, taxon_label, aux_info = parse_label(node.label)
                if taxon_label:
                    taxon_label += '; ' + taxon
                else:
                    taxon_label = taxon
                node.label = create_label(support, taxon_label, aux_info)

        return placed_taxon

    def _write_statistics_table(self, fmeasure_for_taxa, taxonomy, out_table):
        """Write table containing statistics for each taxon.
        
        Parameters
        ----------
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall)]
          Node with highest F-measure for each taxon.
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
          Taxonomic information for taxa in tree of interest.
        out_table : str
          Output table to write statistics for assigned labels.  
        """

        # get extent taxa
        extant_taxa = Taxonomy().extant_taxa(taxonomy)

        fout_table = open(out_table, 'w')
        fout_table.write('Taxon\tNo. Expected in Tree\tF-measure\tPrecision\tRecall')
        fout_table.write('\tNo. Genomes from Taxon\tNo. Genome In Lineage')
        fout_table.write('\tRogue out\tRogue in\n')
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys()):
            if len(fmeasure_for_taxa[taxon]) != 1:
                self.logger.error('Multiple positions specified for taxon label.')
                sys.exit()

            num_genomes = len(extant_taxa[taxon])

            stat_table = fmeasure_for_taxa[taxon][0]
            fout_table.write('%s\t%d\t%.4f\t%.4f\t%.4f\t%d\t%d\t%s\t%s\n' % (
                taxon,
                num_genomes,
                stat_table.fmeasure,
                stat_table.precision,
                stat_table.recall,
                stat_table.taxa_in_lineage,
                stat_table.num_leaves_with_taxa,
                ','.join(stat_table.rogue_out),
                ','.join(stat_table.rogue_in)))

        fout_table.close()
        
    def _write_summary_table(self, fmeasure_for_taxa, taxonomy, summary_table):
        """Write table containing statistics for each taxonomic rank.
        
        Parameters
        ----------
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall)]
          Node with highest F-measure for each taxon.
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
          Taxonomic information for taxa in tree of interest.
        summary_table : str
          Output table to write statistics for assigned labels.  
        """

        # get number of monophyletic, operationally monophyletic, and polyphyletic
        # taxa at each taxonomic rank
        taxon_count = defaultdict(int)
        mono = defaultdict(int)
        op_mono = defaultdict(int)
        poly = defaultdict(int)
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys()):
            if len(fmeasure_for_taxa[taxon]) != 1:
                self.logger.error('Multiple positions specified for taxon label.')
                sys.exit()
                
            rank_prefix = taxon[0:3]
            taxon_count[rank_prefix] += 1

            stat_table = fmeasure_for_taxa[taxon][0]
            if stat_table.fmeasure == 1.0:
                mono[rank_prefix] += 1
            elif stat_table.fmeasure >= 0.95:
                op_mono[rank_prefix] += 1
            else:
                poly[rank_prefix] += 1

        fout = open(summary_table, 'w')
        fout.write('Rank\tNo. taxon')
        fout.write('\tNo. monophyletic\tNo. operationally monophyletic\tNo. polyphyletic')
        fout.write('\tMonophyletic (%)\tOperationally monophyletic (%)\tPolyphyletic (%)\n')
        for idx, rank_prefix in enumerate(Taxonomy.rank_prefixes):
            fout.write('{}\t{}'.format(Taxonomy.rank_labels[idx],
                                        taxon_count[rank_prefix]))
            fout.write('\t{}\t{}\t{}'.format(mono[rank_prefix],
                                                op_mono[rank_prefix],
                                                poly[rank_prefix]))
            fout.write('\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(mono[rank_prefix]*100.0/taxon_count[rank_prefix],
                                                            op_mono[rank_prefix]*100.0/taxon_count[rank_prefix],
                                                            poly[rank_prefix]*100.0/taxon_count[rank_prefix]))
        fout.close()

    def _leaf_taxa(self, leaf):
        """Get taxonomic information for leaf node.
        
        Parameters
        ----------
        leaf : Node
          Node in tree.
          
        Returns
        -------
        list
          Taxa for leaf in rank order.
        """

        leaf_taxa = []

        parent = leaf
        while parent:
            _support, taxon, _aux_info = parse_label(parent.label)

            if taxon:
                for t in taxon.split(';')[::-1]:
                    leaf_taxa.append(t.strip())

            parent = parent.parent_node

        ordered_taxa = leaf_taxa[::-1]

        # fill in missing ranks
        last_rank = ordered_taxa[-1][0:3]
        for i in range(Taxonomy.rank_prefixes.index(last_rank) + 1, len(Taxonomy.rank_prefixes)):
            ordered_taxa.append(Taxonomy.rank_prefixes[i])

        return ordered_taxa

    def _resolve_missing_taxa(self, taxa_list):
        """Resolve missing taxa."""

        for rank_index, taxon in enumerate(taxa_list[0:Taxonomy.rank_labels.index('species')]):
            if taxon == Taxonomy.rank_prefixes[rank_index]:
                # taxon is missing, check if any child taxon are specified
                specified = False
                for r in range(rank_index + 1, Taxonomy.rank_labels.index('species') + 1):
                    if taxa_list[r] != Taxonomy.rank_prefixes[r]:
                        specified = True
                        break

                if specified:
                    # fill in missing ranks
                    for i in range(rank_index, r):
                        taxon_label = '%s{unclassified %s}' % (Taxonomy.rank_prefixes[i],
                                                               taxa_list[rank_index - 1][3:])
                        taxa_list[i] = taxon_label

                        print(taxa_list)

        return taxa_list

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
            taxa = self._leaf_taxa(leaf)
            # taxa = self._resolve_missing_taxa(taxa)
            fout.write('%s\t%s\n' % (leaf.taxon.label, '; '.join(taxa)))

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

    def _resolve_ambiguous_placements(self,
                                      fmeasure_for_taxa,
                                      median_rank_rd,
                                      max_rd_diff=0.1):
        """Resolve ambiguous taxon label placements using median relative divergences.
        
        Parameters
        ----------
        fmeasure_for_taxa : d[taxon] -> [(Node, F-measure, precision, recall)]
          Node with highest F-measure for each taxon.
        median_rank_rd : d[rank_index] -> float
          Median relative divergence for each taxonomic rank.
        max_rd_diff : float
          Maximum difference in relative divergence for assigning a taxonomic label.
        """

        # For ambiguous nodes place them closest to median for rank 
        # and within accepted relative divergence distance. Taxon labels
        # are placed in reverse taxonomic order (species to domain) and
        # this ordering used to ensure taxonomic consistency.
        for taxon in Taxonomy().sort_taxa(fmeasure_for_taxa.keys(), reverse=True):
            if len(fmeasure_for_taxa[taxon]) == 1:
                continue

            rank_prefix = taxon[0:3]
            rank_index = Taxonomy.rank_prefixes.index(rank_prefix)
            if rank_index not in median_rank_rd:
                del fmeasure_for_taxa[taxon]
                continue  # handles trees without any defined taxa at a given rank (e.g. species)
            rd = median_rank_rd[rank_index]

            # Find node closest to median distance, but making sure
            # taxon is not placed below a more specific taxon label.
            # The 'fmeasure_for_taxa' stores node information in preorder.
            closest_index = None
            closest_dist = 1e9
            closest_node = None
            for i, stat_table in enumerate(fmeasure_for_taxa[taxon]):
                cur_node = stat_table.node

                cur_rank_index = -1
                _support, cur_taxon, _aux_info = parse_label(cur_node.label)
                if cur_taxon:
                    cur_prefix = cur_taxon.split(';')[-1].strip()[0:3]
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
                if rd_diff > max_rd_diff:
                    continue

                if rd_diff < closest_dist:
                    closest_dist = rd_diff
                    closest_index = i
                    closest_node = cur_node

            if closest_index is None:
                # no node is within an acceptable relative divergence distance 
                # for this label so it should be placed at the most extant node
                # in order to be conservative
                closest_index = len(fmeasure_for_taxa[taxon]) - 1
                closest_node = fmeasure_for_taxa[taxon][closest_index].node

            # add label to node
            support, cur_taxon, aux_info = parse_label(closest_node.label)
            if not cur_taxon:
                taxa_str = taxon
            else:
                taxa = [t.strip() for t in cur_taxon.split(';')] + [taxon]
                taxa_str = '; '.join(Taxonomy().sort_taxa(taxa))

            closest_node.label = create_label(support, taxa_str, aux_info)

            # remove other potential node assignments
            fmeasure_for_taxa[taxon] = [fmeasure_for_taxa[taxon][closest_index]]

    def run(self,
            input_tree,
            taxonomy_file,
            viral,
            trusted_taxa_file,
            min_children,
            min_support,
            skip_rd_refine,
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
        skip_rd_refine : boolean
          Skip refinement of taxonomy based on relative divergence information.
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

        if viral:
            self.logger.info('Translating viral prefixes.')
            full_taxonomy = translate_viral_taxonomy(full_taxonomy)

        taxonomy = {}
        for leaf in tree.leaf_node_iter():
            taxonomy[leaf.taxon.label] = full_taxonomy.get(leaf.taxon.label, Taxonomy.rank_prefixes)

        # find best placement for each taxon based 
        # on the F-measure statistic
        self.logger.info('Calculating F-measure statistic for each taxa.')
        fmeasure_for_taxa = self._fmeasure(tree, taxonomy)

        # calculating relative
        if not skip_rd_refine:
            # place labels with only one acceptable position and calculate
            # the relative divergence thresholds from these as a guide for
            # placing the remaining labels
            self.logger.info('Placing labels with unambiguous position in tree.')
            placed_taxon = self._assign_taxon_labels(fmeasure_for_taxa)

            self.logger.info('Establishing median relative divergence for taxonomic ranks.')
            median_rank_rd = self._median_rank_rd(tree,
                                                  placed_taxon,
                                                  taxonomy,
                                                  trusted_taxa_file,
                                                  min_children,
                                                  min_support)

            # resolve ambiguous position in tree
            self.logger.info('Resolving ambiguous taxon label placements using median relative divergences.')
            self._resolve_ambiguous_placements(fmeasure_for_taxa, median_rank_rd)
        else:
            # simply select most terminal placement in order to be conservative
            ambiguous_placements = set()
            for taxon, fmeasures in fmeasure_for_taxa.items():
                if len(fmeasures) != 1:
                    ambiguous_placements.add(taxon)
                    fmeasure_for_taxa[taxon] = [fmeasures[-1]]

            if len(ambiguous_placements) > 0:
                self.logger.warning('There are %d taxon with multiple placements of equal quality.' % len(ambiguous_placements))
                self.logger.warning('These were resolved by placing the label at a terminal position.')

            # place all labels on tree
            self.logger.info('Placing labels on tree.')
            placed_taxon = self._assign_taxon_labels(fmeasure_for_taxa)

        # write statistics for placed taxon labels
        self.logger.info('Writing out statistics for taxa.')
        out_table = output_tree + '-table'
        self._write_statistics_table(fmeasure_for_taxa, taxonomy, out_table)
        
        summary_table = output_tree + '-summary'
        self._write_summary_table(fmeasure_for_taxa, taxonomy, summary_table)

        # output taxonomy of extant taxa on tree
        self.logger.info('Writing out taxonomy for extant taxa.')
        out_taxonomy = output_tree + '-taxonomy'
        self._write_taxonomy(tree, out_taxonomy)

        # output decorated tree
        self.logger.info('Writing out decorated tree.')
        tree.write_to_path(output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)

        if viral:
            self.logger.info('Translating output files to viral prefixes.')
            rev_translate_output_file(out_table)
            rev_translate_output_file(out_taxonomy)
            rev_translate_output_file(output_tree)

