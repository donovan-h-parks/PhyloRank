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
from phylorank.plot.robustness_plot import RobustnessPlot

from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           check_file_exists)
from biolib.misc.time_keeper import TimeKeeper
from biolib.external.execute import check_dependencies


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def decorate(self, options):
        """Decorate command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [PhyloRank - decorate] Decorating nodes with inferred taxonomic ranks.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.input_tree)

        decorate = Decorate()
        decorate.run(options.input_tree,
                        options.output_tree,
                        options.min_support,
                        options.named_clades,
                        options.only_named_clades,
                        options.min_length)

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
        elif(options.subparser_name == 'robustness_plot'):
            self.robustness_plot(options)
        else:
            self.logger.error('  [Error] Unknown AutoRank command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
