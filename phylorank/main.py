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
import sys
from collections import defaultdict

import dendropy
from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           check_file_exists)
from biolib.external.execute import check_dependencies
from biolib.misc.time_keeper import TimeKeeper
from biolib.taxonomy import Taxonomy

from phylorank.bl_dist import BranchLengthDistribution
from phylorank.decorate import Decorate
from phylorank.mark_tree import MarkTree
from phylorank.newick import parse_label
from phylorank.outliers import Outliers
from phylorank.plot.robustness_plot import RobustnessPlot
from phylorank.rd_ranks import RdRanks
from phylorank.rel_dist import RelativeDistance
from phylorank.rogue_test import RogueTest
from phylorank.viral_taxonomy import (VIRAL_RANK_LABELS,
                                      VIRAL_RANK_PREFIXES,
                                      sort_viral_taxa)


class OptionsParser(object):
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def outliers(self, options):
        """Create information for identifying taxnomic outliers"""

        check_file_exists(options.input_tree)
        check_file_exists(options.taxonomy_file)

        if options.plot_taxa_file:
            check_file_exists(options.plot_taxa_file)

        if options.trusted_taxa_file:
            check_file_exists(options.trusted_taxa_file)

        if not os.path.exists(options.output_dir):
            os.makedirs(options.output_dir)

        if options.highlight_polyphyly and not options.fmeasure_table:
            self.logger.error("The '--highlight_polyphyly' flag must be used with the '--fmeasure_table' flag.")
            return

        o = Outliers(options.skip_mpld3, options.dpi, options.output_dir)
        o.run(options.input_tree,
              options.taxonomy_file,
              options.viral,
              options.plot_taxa_file,
              options.plot_dist_taxa_only,
              options.plot_domain,
              options.highlight_polyphyly,
              options.highlight_taxa_file,
              options.trusted_taxa_file,
              options.fixed_root,
              options.min_children,
              options.min_support,
              options.mblet,
              options.fmeasure_table,
              options.min_fmeasure,
              options.fmeasure_mono,
              options.verbose_table)

        self.logger.info('Done.')

    def scale_tree(self, options):
        """Scale a rooted tree based on RED."""

        check_file_exists(options.input_tree)

        self.logger.info('Reading tree.')
        tree = dendropy.Tree.get_from_path(options.input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        self.logger.info('Scaling tree based on RED.')
        rd = RelativeDistance()
        rd.decorate_rel_dist(tree)
        for n in tree.preorder_node_iter(lambda n: n != tree.seed_node):
            rd_to_parent = n.rel_dist - n.parent_node.rel_dist
            n.edge_length = rd_to_parent

        tree.write_to_path(options.output_tree,
                           schema='newick',
                           suppress_rooting=True,
                           unquoted_underscores=True)

        self.logger.info('Done.')

    def compare_red(self, options):
        """Compare RED values of taxa calculated over different trees."""

        check_file_exists(options.red_table1)
        check_file_exists(options.red_table2)
        check_file_exists(options.red_dict2)

        median_reds = eval(open(options.red_dict2).readline())

        red1 = {}
        red2 = {}
        lineage = {}
        for d, red_file in [(red1, options.red_table1), (red2, options.red_table2)]:
            with open(red_file) as f:
                f.readline()

                for line in f:
                    line_split = line.strip().split('\t')
                    taxon = line_split[0]
                    median_red = float(line_split[2])
                    d[taxon] = median_red

                    if d == red1:
                        lineage[taxon] = line_split[1]

        red1_label = os.path.splitext(os.path.basename(options.red_table1))[0]
        red2_label = os.path.splitext(os.path.basename(options.red_table2))[0]

        fout = open(options.output_table, 'w')
        fout.write('Taxon\tLineage\t%s\t%s\tDifference\tAbs. Difference\tChanged rank\n' % (red1_label, red2_label))
        if options.viral:
            sorted_taxa = sort_viral_taxa(set(red1.keys()).union(red2.keys()))
        else:
            sorted_taxa = Taxonomy().sort_taxa(set(red1.keys()).union(red2.keys()))

        for taxon in sorted_taxa:
            r1 = red1.get(taxon, 'NA')
            r2 = red2.get(taxon, 'NA')
            if r1 == 'NA':
                fout.write('%s\t%s\t%s\t%.3f\t%s\t%s' % (taxon, 'NA', 'NA', r2, 'NA', 'NA'))
            elif r2 == 'NA':
                fout.write('%s\t%s\t%.3f\t%s\t%s\t%s\t%s\n' % (taxon, lineage[taxon], r1, 'NA', 'NA', 'NA', 'NA'))
            else:
                fout.write('%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f' % (taxon, lineage[taxon], r1, r2, r1 - r2, abs(r1 - r2)))

            if r2 != 'NA':
                rank_prefix = taxon[0:3]
                if rank_prefix == 'd__':
                    continue

                if options.viral:
                    rank_label = VIRAL_RANK_LABELS[VIRAL_RANK_PREFIXES.index(rank_prefix)]
                else:
                    rank_label = Taxonomy.rank_labels[Taxonomy.rank_prefixes.index(rank_prefix)]
                rank_median = median_reds[rank_label]

                closest_rank = rank_label
                closest_dist = 1e6
                if r2 < rank_median - 0.1 or r2 > rank_median + 0.1:
                    for rank, median_red in median_reds.items():
                        d = abs(r2 - median_red)
                        if d < closest_dist:
                            closest_dist = d
                            closest_rank = rank

                if rank_label != closest_rank:
                    fout.write('\tTrue (%s: %.3f)' % (closest_rank, closest_dist))
                else:
                    fout.write('\tFalse')
                fout.write('\n')

        fout.close()

    def mark_tree(self, options):
        """Mark tree command."""

        check_file_exists(options.input_tree)

        mt = MarkTree()
        mt.run(options.input_tree,
               options.output_tree,
               options.min_support,
               options.only_named_clades,
               options.min_length,
               not options.no_percentile,
               not options.no_relative_divergence,
               not options.no_prediction,
               options.thresholds)

        self.logger.info('Marked tree written to: %s' % options.output_tree)

    def rogue_test(self, options):
        """Rogue taxa command."""

        check_dir_exists(options.input_tree_dir)
        check_file_exists(options.taxonomy_file)
        make_sure_path_exists(options.output_dir)

        if options.decorate:
            check_dependencies(['genometreetk'])

        rt = RogueTest()
        rt.run(options.input_tree_dir,
               options.taxonomy_file,
               options.outgroup_taxon,
               options.decorate,
               options.output_dir)

        self.logger.info('Finished rogue taxa test.')

    def decorate(self, options):
        """Place internal taxonomic labels on tree."""

        check_file_exists(options.input_tree)
        check_file_exists(options.taxonomy_file)

        decorate = Decorate()
        decorate.run(options.input_tree,
                     options.taxonomy_file,
                     options.viral,
                     options.trusted_taxa_file,
                     options.min_children,
                     options.min_support,
                     options.skip_rd_refine,
                     options.output_tree)

        self.logger.info('Finished decorating tree.')

    def taxon_stats(self, options):
        """Taxon stats command"""

        check_file_exists(options.taxonomy_file)

        taxonomy = Taxonomy().read(options.taxonomy_file)
        taxon_children = Taxonomy().taxon_children(taxonomy)

        fout = open(options.output_file, 'w')
        fout.write('Taxa')
        for rank in Taxonomy.rank_labels[1:]:
            fout.write('\t# named %s' % rank)
        fout.write('\t# extant taxon with complete taxonomy')
        fout.write('\n')

        for rank_prefix in Taxonomy.rank_prefixes:
            # find taxon at the specified rank
            cur_taxa = []
            for taxon in taxon_children:
                if taxon.startswith(rank_prefix):
                    cur_taxa.append(taxon)

            cur_taxa.sort()

            for taxon in cur_taxa:
                fout.write(taxon)
                fout.write('\t-' * Taxonomy.rank_index[rank_prefix])

                next_taxa = [taxon]
                for _ in range(Taxonomy.rank_index[rank_prefix], Taxonomy.rank_index['s__'] + 1):
                    children_taxa = set()
                    for t in next_taxa:
                        children_taxa.update(taxon_children[t])

                    fout.write('\t%d' % len(children_taxa))
                    next_taxa = children_taxa
                fout.write('\n')

        fout.close()

        self.logger.info('Summary statistics written to: %s' % options.output_file)

    def robustness_plot(self, options):
        """Robustness plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - robustness_plot] Plotting distances across a set of tree.')
        self.logger.info('*******************************************************************************')

        robustness_plot = RobustnessPlot()
        robustness_plot.run(options.rank,
                            options.input_tree_dir,
                            options.full_tree_file,
                            options.derep_tree_file,
                            options.taxonomy_file,
                            options.output_prefix,
                            options.min_children,
                            options.title)

        self.time_keeper.print_time_stamp()

    def rd_ranks(self, options):
        """Calculate number of taxa for specified rd thresholds."""

        check_file_exists(options.input_tree)
        make_sure_path_exists(options.output_dir)

        r = RdRanks()
        r.run(options.input_tree,
              options.thresholds,
              options.output_dir)

        self.logger.info('Done.')

    def bl_dist(self, options):
        """Calculate distribution of branch lengths at each taxonomic rank."""

        check_file_exists(options.input_tree)
        make_sure_path_exists(options.output_dir)

        b = BranchLengthDistribution()
        b.run(options.input_tree,
              options.trusted_taxa_file,
              options.min_children,
              options.taxonomy_file,
              options.output_dir)

        self.logger.info('Done.')

    def bl_optimal(self, options):
        """Determine branch length for best congruency with existing taxonomy."""

        b = BranchLengthDistribution()
        optimal_bl, correct_taxa, incorrect_taxa = b.optimal(options.input_tree,
                                                             options.rank,
                                                             options.min_dist,
                                                             options.max_dist,
                                                             options.step_size,
                                                             options.output_table)

        prec = float(correct_taxa) / (correct_taxa + incorrect_taxa)

        self.logger.info('Optimal branch length is %f.' % optimal_bl)
        self.logger.info('This results in %d correct and %d incorrect taxa (precision = %.2f).' % (
        correct_taxa, incorrect_taxa, prec))

    def bl_decorate(self, options):
        """Decorate tree based using a mean branch length criterion."""

        check_file_exists(options.input_tree)

        b = BranchLengthDistribution()
        b.decorate(options.input_tree,
                   options.taxonomy_file,
                   options.threshold,
                   options.rank,
                   options.retain_named_lineages,
                   options.keep_labels,
                   options.prune,
                   options.output_tree)

        self.logger.info('Done.')

    def bl_table(self, options):
        """Produce table with number of lineage for increasing mean branch lengths."""

        check_file_exists(options.input_tree)
        check_file_exists(options.taxon_category)

        b = BranchLengthDistribution()
        b.table(options.input_tree,
                options.taxon_category,
                options.step_size,
                options.output_table)

        self.logger.info('Done.')

    def rank_res(self, options):
        """Calculate taxonomic resolution at each rank."""

        check_file_exists(options.input_tree)
        check_file_exists(options.taxonomy_file)

        if options.taxa_file:
            taxa_out = open(options.taxa_file, 'w')
            taxa_out.write('Rank\tLowest Rank\tTaxon\n')

        # determine taxonomic resolution of named groups
        tree = dendropy.Tree.get_from_path(options.input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        rank_res = defaultdict(lambda: defaultdict(int))
        for node in tree.preorder_node_iter(lambda n: n != tree.seed_node):
            if not node.label or node.is_leaf():
                continue

            _support, taxon_name, _auxiliary_info = parse_label(node.label)

            if taxon_name:
                lowest_rank = [x.strip() for x in taxon_name.split(';')][-1][0:3]
                for rank_prefix in Taxonomy.rank_prefixes:
                    if rank_prefix in taxon_name:
                        rank_res[rank_prefix][lowest_rank] += 1
                        if options.taxa_file:
                            rank_prefix_name = Taxonomy.rank_labels[Taxonomy.rank_index[rank_prefix]]
                            lowest_rank_name = Taxonomy.rank_labels[Taxonomy.rank_index[lowest_rank]]
                            taxa_out.write('%s\t%s\t%s\n' % (rank_prefix_name, lowest_rank_name, taxon_name))

        # identify any singleton taxa which are treated as having species level resolution
        for line in open(options.taxonomy_file):
            line_split = line.split('\t')
            genome_id = line_split[0]
            taxonomy = line_split[1].split(';')

            for i, rank_prefix in enumerate(Taxonomy.rank_prefixes):
                if taxonomy[i] == rank_prefix:
                    # this taxa is undefined at the specified rank so
                    # must be the sole representative; e.g., a p__
                    # indicates a taxon that represents a novel phyla
                    rank_res[rank_prefix]['s__'] += 1
                    if options.taxa_file:
                        rank_prefix_name = Taxonomy.rank_labels[Taxonomy.rank_index[rank_prefix]]
                        taxa_out.write('%s\t%s\t%s (%s)\n' % (rank_prefix_name, 'species', taxonomy[i], genome_id))
        if options.taxa_file:
            taxa_out.close()

        # write out results
        fout = open(options.output_file, 'w')
        fout.write('Category')
        for rank in Taxonomy.rank_labels[1:]:
            fout.write('\t' + rank)
        fout.write('\n')

        for i, rank_prefix in enumerate(Taxonomy.rank_prefixes[1:]):
            fout.write(Taxonomy.rank_labels[i + 1])

            for j, r in enumerate(Taxonomy.rank_prefixes[1:]):
                if i >= j:
                    fout.write('\t' + str(rank_res[r].get(rank_prefix, 0)))
                else:
                    fout.write('\t-')
            fout.write('\n')
        fout.close()

        self.logger.info('Done.')

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        # check_dependencies(('diamond', 'ktImportText'))

        if options.subparser_name == 'outliers':
            self.outliers(options)
        elif options.subparser_name == 'scale_tree':
            self.scale_tree(options)
        elif options.subparser_name == 'compare_red':
            self.compare_red(options)
        elif options.subparser_name == 'mark_tree':
            self.mark_tree(options)
        elif options.subparser_name == 'rogue_test':
            self.rogue_test(options)
        elif options.subparser_name == 'decorate':
            self.decorate(options)
        elif options.subparser_name == 'taxon_stats':
            self.taxon_stats(options)
        elif options.subparser_name == 'robustness_plot':
            self.robustness_plot(options)
        elif options.subparser_name == 'rd_ranks':
            self.rd_ranks(options)
        elif options.subparser_name == 'bl_dist':
            self.bl_dist(options)
        elif options.subparser_name == 'bl_optimal':
            self.bl_optimal(options)
        elif options.subparser_name == 'bl_decorate':
            self.bl_decorate(options)
        elif options.subparser_name == 'bl_table':
            self.bl_table(options)
        elif options.subparser_name == 'rank_res':
            self.rank_res(options)
        else:
            self.logger.error('  [Error] Unknown PhyloRank command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
