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

import argparse
import json
import logging
import ntpath
import os
import sys

from biolib.common import make_sure_path_exists
from biolib.misc.custom_help_formatter import CustomHelpFormatter

from phylorank import __version__
from phylorank.main import OptionsParser


def print_help():
    """Help menu."""

    # Deprecated:
    # rd_ranks    -> Calculate number of taxa for specified RED thresholds
    # dist_plot   -> Plot distribution of taxa in each taxonomic rank
    # robustness_plot -> Plot relative distance of groups across a set of trees.

    print('')
    print('                ...::: PhyloRank v' + __version__ + ' :::...''')
    print('''\

  Curation methods:
    outliers    -> Create RED table, scaled tree, and plot useful for identifying taxonomic outliers
    scale_tree  -> Scale a rooted tree based on RED
    compare_red -> Compare RED values of taxa calculated over different trees
    mark_tree   -> Mark nodes with distribution information and predicted taxonomic ranks
    rogue_test  -> Index indicating the incongruence of genomes over a set of tree

  Taxonomy decoration and inspection methods:
    decorate    -> Place internal taxonomic labels on tree
    taxon_stats -> Summary statistics of taxonomic groups
    rank_res    -> Calculate taxonomic resolution at each rank

  Mean branch length to extant taxa methods:
    bl_dist     -> Calculate distribution of branch lengths at each taxonomic rank
    bl_optimal  -> Determine branch length for best congruency with existing taxonomy
    bl_decorate -> Decorate tree using a mean branch length criterion
    bl_table    -> Produce table with number of lineage for increasing mean branch lengths

  Use: phylorank <command> -h for command specific help.

  Feature requests or bug reports can be sent to Donovan Parks (donovan.parks@gmail.com)
    or posted on GitHub (https://github.com/dparks1134/phylorank.)
    ''')


def logger_setup(log_file, silent):
    """Set logging for application.

    Parameters
    ----------
    log_file : str
        Name of log file.
    silent : boolean
        Flag indicating if output to stdout should be suppressed.
    """

    # setup general properties of logger
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")

    # setup logging to console
    if not silent:
        stream_logger = logging.StreamHandler(sys.stdout)
        stream_logger.setFormatter(log_format)
        stream_logger.setLevel(logging.INFO)
        logger.addHandler(stream_logger)

    if log_file:
        file_logger = logging.FileHandler(log_file, 'a')
        file_logger.setFormatter(log_format)
        logger.addHandler(file_logger)

    logger.info('PhyloRank v%s' % __version__)
    logger.info(ntpath.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))


def main(args=None):
    # initialize the options parser
    parser = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    # create table and plot useful for identifying taxonomic outliers.
    outliers_parser = subparsers.add_parser('outliers',
                                            formatter_class=CustomHelpFormatter,
                                            description='Create information for identifying taxonomic outliers')

    outliers_parser.add_argument('input_tree', help="decorated tree for inferring RED outliers")
    outliers_parser.add_argument('taxonomy_file', help='taxonomy file for inferring RED outliers', default=None)
    outliers_parser.add_argument('output_dir', help="desired output directory for generated files")
    outliers_parser.add_argument('--viral', action="store_true", help='indicates a viral input tree and taxonomy')
    outliers_parser.add_argument('--fixed_root', action="store_true", help='use single fixed root to infer outliers')
    outliers_parser.add_argument('-t', '--trusted_taxa_file',
                                 help="file indicating trusted taxonomic groups to use for inferring distribution (default: all taxa)",
                                 default=None)
    outliers_parser.add_argument('-m', '--min_children',
                                 help='minimum required child taxa to consider taxa when inferring distribution',
                                 type=int, default=2)
    outliers_parser.add_argument('-s', '--min_support',
                                 help="minimum support value to consider taxa when inferring distribution (default: 0)",
                                 type=float, default=0.0)
    outliers_parser.add_argument('--fmeasure_table', help="table indicating F-measure score for each taxa")
    outliers_parser.add_argument('--min_fmeasure',
                                 help="minimum F-measure to consider taxa when inferring distribution", type=float,
                                 default=0.95)
    outliers_parser.add_argument('--fmeasure_mono', help="minimum F-measure to consider taxa monophyletic", type=float,
                                 default=0.95)
    outliers_parser.add_argument('--highlight_polyphyly',
                                 help='highlight taxa with an F-measure less than --fmeasure_mono', action="store_true")
    outliers_parser.add_argument('--mblet', action="store_true",
                                 help="calculate 'mean branch length to extent taxa' instead of 'relative evolutionary distances'")
    outliers_parser.add_argument('-p', '--plot_taxa_file',
                                 help="file indicating taxonomic groups to plot (default: all taxa)", default=None)
    outliers_parser.add_argument('--plot_domain', action="store_true", help='show domain rank in plot')
    outliers_parser.add_argument('--plot_dist_taxa_only', help='only plot taxa used to infer distribution',
                                 action="store_true")
    outliers_parser.add_argument('--highlight_taxa_file', help='file indicating taxa to highlight')
    outliers_parser.add_argument('--dpi', help='DPI of plots', type=int, default=96)
    outliers_parser.add_argument('--verbose_table', action="store_true", help='add additional columns to output table')
    outliers_parser.add_argument('--skip_mpld3', action="store_true", help='skip plots requiring mpld3')

    # create table and plot useful for identifying taxonomic outliers.
    scale_tree_parser = subparsers.add_parser('scale_tree',
                                              formatter_class=CustomHelpFormatter,
                                              description='Scale a rooted tree based on RED')

    scale_tree_parser.add_argument('input_tree', help="rooted tree to scale")
    scale_tree_parser.add_argument('output_tree', help="tree scaled by RED")

    # Compare RED values of taxa calculated over different trees
    compare_red_parser = subparsers.add_parser('compare_red',
                                               formatter_class=CustomHelpFormatter,
                                               description='Compare RED values of taxa calculated over different trees')
    compare_red_parser.add_argument('red_table1', help="RED table calculated by 'outlier' command.")
    compare_red_parser.add_argument('red_table2', help="RED table calculated by 'outlier' command.")
    compare_red_parser.add_argument('red_dict2', help="Median RED dictionary calculated by 'outlier' command.")
    compare_red_parser.add_argument('output_table', help='output table')
    compare_red_parser.add_argument('--viral', action="store_true", help='indicates a viral input tree and taxonomy')

    # plot distribution of groups in each taxonomic rank
    dist_plot_parser = subparsers.add_parser('dist_plot',
                                             formatter_class=CustomHelpFormatter,
                                             description='Plot distribution of taxa in each taxonomic rank')

    dist_plot_parser.add_argument('input_tree',
                                  help="decorated tree for establishing relative divergence distributions")
    dist_plot_parser.add_argument('output_prefix', help="output prefix for generated files")
    dist_plot_parser.add_argument('-p', '--plot_taxa_file',
                                  help="file indicating taxonomic groups to plot (default: all taxa)", default=None)
    dist_plot_parser.add_argument('-t', '--trusted_taxa_file',
                                  help="file indicating trusted taxonomic groups to use for inferring distribution (default: all taxa)",
                                  default=None)
    dist_plot_parser.add_argument('-m', '--min_children',
                                  help='minimum required child taxa to consider taxa when inferring distribution  (default: 0)',
                                  type=int, default=0)
    dist_plot_parser.add_argument('-s', '--min_support',
                                  help="minimum support value to consider taxa when inferring distribution (default: 0)",
                                  type=float, default=0.0)

    # decorate nodes with inferred taxonomic ranks

    # ******** MAYBE THIS SHOULD JUST TAKE A 'DISTRIBUTIONS FILE' produce by 'dist_plot'
    # ************************************************************************************
    mark_tree_parser = subparsers.add_parser('mark_tree',
                                             formatter_class=CustomHelpFormatter,
                                             description='Mark nodes with distribution information and predicted taxonomic ranks.')

    mark_tree_parser.add_argument('input_tree', help="input tree to mark")
    mark_tree_parser.add_argument('output_tree', help="output tree with assigned taxonomic ranks")
    mark_tree_parser.add_argument('-t', '--thresholds', help="relative divergence thresholds for taxonomic ranks",
                                  type=json.loads,
                                  default='{"d": 0.33, "p": 0.56, "c": 0.65, "o": 0.78, "f": 0.92, "g": 0.99}')
    mark_tree_parser.add_argument('-s', '--min_support',
                                  help="only mark nodes above the specified support value (default=0)", type=float,
                                  default=0)
    mark_tree_parser.add_argument('-n', '--only_named_clades', help="only mark nodes with an existing label",
                                  action='store_true')
    mark_tree_parser.add_argument('-l', '--min_length',
                                  help="only mark nodes with a parent branch above the specified length (default=0)",
                                  type=float, default=0.0)
    mark_tree_parser.add_argument('--no_percentile', action="store_true",
                                  help="do not mark with percentile information")
    mark_tree_parser.add_argument('--no_relative_divergence', action="store_true",
                                  help="do not mark with relative divergence information")
    mark_tree_parser.add_argument('--no_prediction', action="store_true",
                                  help="do not mark with predicted rank information")

    # rogue test
    rogue_test_parser = subparsers.add_parser('rogue_test',
                                              formatter_class=CustomHelpFormatter,
                                              description='Index indicating the incongruence of genomes over a set of tree.')

    rogue_test_parser.add_argument('input_tree_dir', help="directory containing trees to assess incongruence over")
    rogue_test_parser.add_argument('taxonomy_file', help='file indicating taxonomy of extant taxa')
    rogue_test_parser.add_argument('output_dir', help="desired output directory for generated files")
    rogue_test_parser.add_argument('--outgroup_taxon',
                                   help='taxon to use as outgroup (e.g., d__Archaea); imples tree should be rooted')
    rogue_test_parser.add_argument('--decorate', help='indicates trees should be decorated', action='store_true')

    # decorate ree
    decorate_parser = subparsers.add_parser('decorate',
                                            formatter_class=CustomHelpFormatter,
                                            description='Place internal taxonomic labels on tree.')
    decorate_parser.add_argument('input_tree', help='tree to decorate')
    decorate_parser.add_argument('taxonomy_file', help='file indicating taxonomy of extant taxa')
    decorate_parser.add_argument('output_tree', help='decorated tree')
    decorate_parser.add_argument('--viral', action="store_true", help='indicates a viral input tree and taxonomy')
    decorate_parser.add_argument('--skip_rd_refine',
                                 help="skip refinement of taxonomy based on relative divergence information",
                                 action='store_true')
    decorate_parser.add_argument('-t', '--trusted_taxa_file',
                                 help="file indicating trusted taxonomic groups to use for inferring distribution (default: all taxa)",
                                 default=None)
    decorate_parser.add_argument('-m', '--min_children',
                                 help='minimum required child taxa to consider taxa when inferring distribution',
                                 type=int, default=2)
    decorate_parser.add_argument('-s', '--min_support',
                                 help="minimum support value to consider taxa when inferring distribution (default: 0)",
                                 type=float, default=0.0)

    # pull taxonomy strings from tree
    pull_parser = subparsers.add_parser('pull',
                                        formatter_class=CustomHelpFormatter,
                                        description='Pull taxonomy information from tree.')

    pull_parser.add_argument('input_tree', help="input tree to extract taxonomy from")
    pull_parser.add_argument('output_file', help="file to contain taxonomy strings for each extant taxon")
    pull_parser.add_argument('--no_rank_fill', action="store_true", help="do not automatically fill in missing ranks")

    # validate consistency of taxonomy
    validate_parser = subparsers.add_parser('validate',
                                            formatter_class=CustomHelpFormatter,
                                            description='Validate consistency of taxonomy.')

    validate_parser.add_argument('taxonomy_file', help="file with taxonomy for extant taxa")
    validate_parser.add_argument('--no_prefix', action="store_true", help="do not check taxon prefixes")
    validate_parser.add_argument('--no_all_ranks', action="store_true",
                                 help="do not check for the presence of all ranks")
    validate_parser.add_argument('--no_hierarhcy', action="store_true",
                                 help="do not check for inconsistencies in the taxonomic hierarchy")
    validate_parser.add_argument('--no_species', action="store_true",
                                 help="do not check for hierarchical inconsistencies with named species")

    # summary statistics of taxonomic groups
    taxon_stats_parser = subparsers.add_parser('taxon_stats',
                                               formatter_class=CustomHelpFormatter,
                                               description='Summary statistics of taxonomic groups.')

    taxon_stats_parser.add_argument('taxonomy_file', help="file with taxonomy for extant taxa")
    taxon_stats_parser.add_argument('output_file', help="output file with summary statistics")

    # plot relative distance of groups across a set of trees.
    robustness_plot_parser = subparsers.add_parser('robustness_plot',
                                                   formatter_class=CustomHelpFormatter,
                                                   description='Plot relative divergence of groups across a set of trees')

    robustness_plot_parser.add_argument('rank', help="taxonomic rank of named groups to plot", type=int,
                                        choices=[1, 2, 3, 4, 5, 6])
    robustness_plot_parser.add_argument('input_tree_dir',
                                        help="directory containing trees to inferred relative divergence across")
    robustness_plot_parser.add_argument('full_tree_file',
                                        help="unmodified tree to include in plot; must be decorate with taxonomy")
    robustness_plot_parser.add_argument('derep_tree_file', help="dereplicated tree to include in plot")
    robustness_plot_parser.add_argument('taxonomy_file', help="file indicating taxonomy string for each genome")
    robustness_plot_parser.add_argument('output_prefix', help="output prefix for generated files")
    robustness_plot_parser.add_argument('-m', '--min_children', help='minimum named child taxa to consider taxa',
                                        type=int, default=2)
    robustness_plot_parser.add_argument('-t', '--title', help='title of plot', default=None)

    rd_ranks_parser = subparsers.add_parser('rd_ranks',
                                            formatter_class=CustomHelpFormatter,
                                            description='Calculate number of taxa for specified rd thresholds.')

    rd_ranks_parser.add_argument('input_tree', help="input tree to calculate ranks over")
    rd_ranks_parser.add_argument('output_dir', help="desired output directory for generated files")
    rd_ranks_parser.add_argument('-t', '--thresholds', help="relative divergence thresholds for taxonomic ranks",
                                 type=json.loads,
                                 default='{"p": 0.35, "c": 0.52, "o": 0.67, "f": 0.79, "g": 0.94, "s":0.996}')

    bl_dist_parser = subparsers.add_parser('bl_dist',
                                           formatter_class=CustomHelpFormatter,
                                           description='Calculate distribution of branch lengths at each taxonomic rank.')

    bl_dist_parser.add_argument('input_tree', help="input tree to calculate branch length distributions")
    bl_dist_parser.add_argument('output_dir', help="desired output directory for generated files")
    bl_dist_parser.add_argument('-t', '--trusted_taxa_file',
                                help="file indicating trusted taxonomic groups to use for inferring distribution (default: all taxa)",
                                default=None)
    bl_dist_parser.add_argument('-m', '--min_children',
                                help='minimum required child taxa to consider taxa when inferring distribution',
                                type=int, default=2)
    bl_dist_parser.add_argument('--taxonomy_file', help='read taxonomy from this file instead of directly from tree',
                                default=None)

    bl_optimal_parser = subparsers.add_parser('bl_optimal',
                                              formatter_class=CustomHelpFormatter,
                                              description='Determine branch length for best congruency with existing taxonomy.')

    bl_optimal_parser.add_argument('input_tree', help="input tree to calculate branch length distributions")
    bl_optimal_parser.add_argument('rank', help="rank of labels", type=int, choices=[1, 2, 3, 4, 5, 6])
    bl_optimal_parser.add_argument('output_table', help="desired named of output table")
    bl_optimal_parser.add_argument('--min_dist', help='minimum mean branch length value to evaluate', type=float,
                                   default=0.5)
    bl_optimal_parser.add_argument('--max_dist', help='maximum mean branch length value to evaluate', type=float,
                                   default=1.2)
    bl_optimal_parser.add_argument('--step_size', help='step size of mean branch length values', type=float,
                                   default=0.025)

    bl_decorate_parser = subparsers.add_parser('bl_decorate',
                                               formatter_class=CustomHelpFormatter,
                                               description='Decorate tree using a mean branch length criterion.')

    bl_decorate_parser.add_argument('input_tree', help="input tree to decorate")
    bl_decorate_parser.add_argument('taxonomy_file', help="file with taxonomic information for each taxon")
    bl_decorate_parser.add_argument('threshold', help="mean branch length threshold", type=float)
    bl_decorate_parser.add_argument('rank', help="rank of labels", type=int, choices=[1, 2, 3, 4, 5, 6])
    bl_decorate_parser.add_argument('output_tree', help="decorate tree")
    bl_decorate_parser.add_argument('--retain_named_lineages', action="store_true",
                                    help='retain existing named lineages at the specified rank')
    bl_decorate_parser.add_argument('--keep_labels', action="store_true", help='keep all existing internal labels')
    bl_decorate_parser.add_argument('--prune', action="store_true",
                                    help='prune tree to preserve only the shallowest and deepest taxa in each child lineage from newly decorated nodes')

    bl_table_parser = subparsers.add_parser('bl_table',
                                            formatter_class=CustomHelpFormatter,
                                            description='Produce table with number of lineage for increasing mean branch lengths.')

    bl_table_parser.add_argument('input_tree', help="input tree to calculate branch length distributions")
    bl_table_parser.add_argument('taxon_category', help="file indicating category for each taxon in the tree")
    bl_table_parser.add_argument('output_table', help="desired named of output table")
    bl_table_parser.add_argument('--step_size', help="step size for mean branch length criterion", type=float,
                                 default=0.01)

    rank_res_parser = subparsers.add_parser('rank_res',
                                            formatter_class=CustomHelpFormatter,
                                            description='Calculate taxonomic resolution at each rank.')

    rank_res_parser.add_argument('input_tree', help="decorated tree")
    rank_res_parser.add_argument('taxonomy_file', help="file with taxonomy for extant taxa")
    rank_res_parser.add_argument('output_file', help="output file with resolution of taxa at each rank")
    rank_res_parser.add_argument('--taxa_file', help="output file indicating taxa within each resolution category",
                                 default=None)

    # get and check options
    if len(sys.argv) == 1 or sys.argv[1] in {'-h', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()

    if hasattr(args, 'output_dir'):
        make_sure_path_exists(args.output_dir)
        logger_setup(os.path.join(args.output_dir, 'phylorank.log'), False)
    elif hasattr(args, 'output_prefix'):
        output_dir, output_prefix = os.path.split(args.output_prefix)
        if output_dir:
            make_sure_path_exists(output_dir)
        logger_setup(os.path.join(output_dir, 'phylorank.log'), False)
    else:
        logger_setup('phylorank.log', False)

    # do what we came here to do
    try:
        parser = OptionsParser()
        if (False):
            # import pstats
            # p = pstats.Stats('prof')
            # p.sort_stats('cumulative').print_stats(10)
            # p.sort_stats('time').print_stats(10)
            import cProfile
            cProfile.run('parser.parse_options(args)', 'prof')
        elif False:
            import pdb
            pdb.run(parser.parse_options(args))
        else:
            parser.parse_options(args)
    except SystemExit:
        print("\n  Controlled exit resulting from an unrecoverable error or warning.")
        sys.exit(1)
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise


if __name__ == "__main__":
    main()
