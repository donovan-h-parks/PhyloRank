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

from autorank.decorate import Decorate

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
        self.logger.info(' [Autorank - decorate] Decorating nodes with inferred taxonomic ranks.')
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

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        # check_dependencies(('diamond', 'ktImportText'))

        if(options.subparser_name == 'decorate'):
            self.decorate(options)
        else:
            self.logger.error('  [Error] Unknown AutoRank command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
