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

from phylorank.decorate import Decorate
from phylorank.newick import read_from_tree
from phylorank.plot.robustness_plot import RobustnessPlot
from phylorank.plot.distribution_plot import DistributionPlot

from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           check_file_exists)
from biolib.taxonomy import Taxonomy
from biolib.misc.time_keeper import TimeKeeper
from biolib.external.execute import check_dependencies

from skbio import TreeNode


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def dist_plot(self, options):
        """Distribution plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - dist_plot] Plotting distribution of each taxonomic rank.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)

        if options.trusted_taxa_file:
            check_file_exists(options.trusted_taxa_file)

        dist_plot = DistributionPlot()
        dist_plot.run(options.input_tree,
                            options.output_prefix,
                            options.trusted_taxa_file,
                            options.min_children,
                            options.min_support)

        self.time_keeper.print_time_stamp()

    def decorate(self, options):
        """Decorate command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - decorate] Decorating nodes with distribution and rank info.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)

        decorate = Decorate()
        decorate.run(options.input_tree,
                        options.output_tree,
                        options.min_support,
                        options.only_named_clades,
                        options.min_length,
                        not options.no_percentile,
                        not options.no_relative_divergence,
                        not options.no_prediction,
                        options.thresholds)

        self.logger.info('')
        self.logger.info('  Decorated tree written to: %s' % options.output_tree)

        self.time_keeper.print_time_stamp()

    def pull(self, options):
        """Pull command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - pull] Extracting taxonomy from tree.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)

        t = read_from_tree(options.input_tree)
        if not options.no_rank_fill:
            for taxon_id, taxa in t.iteritems():
                t[taxon_id] = Taxonomy().fill_missing_ranks(taxa)

        Taxonomy().write(t, options.output_file)

        self.logger.info('')
        self.logger.info('  Taxonomy strings written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

    def validate(self, options):
        """Validate command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - validate] Validating consistency of taxonomy.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.taxonomy_file)

        taxonomy = Taxonomy()
        t = taxonomy.read(options.taxonomy_file)

        errors = taxonomy.validate(t,
                                     not options.no_prefix,
                                     not options.no_all_ranks,
                                     not options.no_hierarhcy,
                                     not options.no_species,
                                     True)

        invalid_ranks, invalid_prefixes, invalid_species_name, invalid_hierarchies = errors

        self.logger.info('')
        if sum([len(e) for e in errors]) == 0:
            self.logger.info('  No errors identified in taxonomy file.')
        else:
            self.logger.info('  Identified %d incomplete taxonomy strings.' % len(invalid_ranks))
            self.logger.info('  Identified %d rank prefix errors.' % len(invalid_prefixes))
            self.logger.info('  Identified %d invalid species names.' % len(invalid_species_name))
            self.logger.info('  Identified %d taxa with multiple parents.' % len(invalid_hierarchies))

        self.time_keeper.print_time_stamp()

    def append(self, options):
        """Append command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - append] Appending taxonomy to extant tree labels.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)
        check_file_exists(options.taxonomy_file)

        taxonomy = Taxonomy().read(options.taxonomy_file)

        root = TreeNode.read(options.input_tree, convert_underscores=False)
        for n in root.tips():
            taxa_str = taxonomy.get(n.name, None)
            if taxa_str == None:
                self.logger.error('  [Error] Taxonomy file does not contain an entry for %s.' % n.name)
                sys.exit(-1)
            n.name = n.name + '|' + ';'.join(taxonomy[n.name])

        root.write(options.output_tree)

        self.logger.info('')
        self.logger.info('  Decorated tree written to: %s' % options.output_tree)

        self.time_keeper.print_time_stamp()

    def taxon_stats(self, options):
        """Taxon stats command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - taxon_stats] Summarizing statistics of taxonomic groups.')
        self.logger.info('*******************************************************************************')

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
                for _ in xrange(Taxonomy.rank_index[rank_prefix], Taxonomy.rank_index['s__'] + 1):
                    children_taxa = set()
                    for t in next_taxa:
                        children_taxa.update(taxon_children[t])

                    fout.write('\t%d' % len(children_taxa))
                    next_taxa = children_taxa
                fout.write('\n')

        fout.close()

        self.logger.info('')
        self.logger.info('  Summary statistics written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

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

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        # check_dependencies(('diamond', 'ktImportText'))

        if(options.subparser_name == 'decorate'):
            self.decorate(options)
        elif(options.subparser_name == 'pull'):
            self.pull(options)
        elif(options.subparser_name == 'validate'):
            self.validate(options)
        elif(options.subparser_name == 'append'):
            self.append(options)
        elif(options.subparser_name == 'taxon_stats'):
            self.taxon_stats(options)
        elif(options.subparser_name == 'robustness_plot'):
            self.robustness_plot(options)
        elif(options.subparser_name == 'dist_plot'):
            self.dist_plot(options)
        else:
            self.logger.error('  [Error] Unknown PhyloRank command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
