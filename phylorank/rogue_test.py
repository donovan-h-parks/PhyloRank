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
import operator
import os
from collections import defaultdict

import dendropy
from biolib.common import make_sure_path_exists
from biolib.taxonomy import Taxonomy
from numpy import mean as np_mean

from phylorank.decorate import Decorate


class RogueTest(object):
    """Index indicating the incongruence of genomes over a set of tree."""

    def __init__(self):
        """Initialize."""

        self.logger = logging.getLogger()

    def _root(self, input_tree_dir, outgroup_taxon, taxonomy_file, output_dir):
        """Root trees in input specified directory."""

        make_sure_path_exists(output_dir)

        for f in os.listdir(input_tree_dir):
            if f.endswith('.tre') or f.endswith('.tree') or f.endswith('.nwk'):
                input_tree = os.path.join(input_tree_dir, f)
                tree_ext = os.path.splitext(f)[1]
                out_tree = os.path.join(output_dir, f.replace(tree_ext, '.rooted' + tree_ext))
                cmd = 'genometreetk outgroup %s %s %s %s' % (input_tree,
                                                             taxonomy_file,
                                                             outgroup_taxon,
                                                             out_tree)
                os.system(cmd)

    def _decorate(self, input_tree_dir, taxonomy_file, output_dir):
        """Decorate trees in specified directory."""

        make_sure_path_exists(output_dir)

        for f in os.listdir(input_tree_dir):
            if f.endswith('.tre') or f.endswith('.tree') or f.endswith('.nwk'):
                decorate = Decorate()

                input_tree = os.path.join(input_tree_dir, f)
                tree_ext = os.path.splitext(f)[1]
                output_tree = os.path.join(output_dir, f.replace(tree_ext, '.decorated' + tree_ext))
                decorate.run(input_tree,
                             taxonomy_file,
                             None,
                             2,
                             0,
                             True,
                             output_tree)

    def run(self, input_tree_dir, taxonomy_file, outgroup_taxon, decorate, output_dir):
        """Calculate rogue taxa index for each genome."""

        taxonomy = Taxonomy().read(taxonomy_file)

        # check if trees need to be rooted
        if outgroup_taxon:
            self.logger.info('Rooting trees.')
            root_dir = os.path.join(output_dir, 'rooted')
            self._root(input_tree_dir, outgroup_taxon, taxonomy_file, root_dir)

        # check if trees need to be decorated
        if decorate:
            self.logger.info('Decorating trees.')
            tree_dir = input_tree_dir
            if outgroup_taxon:
                tree_dir = root_dir

            decorate_dir = os.path.join(output_dir, 'decorated')
            self._decorate(tree_dir, taxonomy_file, decorate_dir)

        # count occurrence of genomes in tree
        # (not all genome need be in all trees!)
        self.logger.info('Counting occurrence of each genome across the trees.')
        genome_count = defaultdict(int)
        tree_count = 0
        for f in os.listdir(input_tree_dir):
            if f.endswith('.tre') or f.endswith('.tree') or f.endswith('.nwk'):
                input_tree = os.path.join(input_tree_dir, f)
                tree = dendropy.Tree.get_from_path(input_tree,
                                                   schema='newick',
                                                   rooting='force-rooted',
                                                   preserve_underscores=True)
                tree_count += 1
                for n in tree.leaf_node_iter():
                    genome_count[n.taxon.label] += 1
        self.logger.info('Genomes were identified in %.1f%% of the %d input trees on average.' % (
            np_mean(list(genome_count.values())) * 100.0 / tree_count,
            tree_count))

        # determine rogue out and rogue in genomes
        self.logger.info('Determining rogue genomes.')
        table_dir = input_tree_dir
        if decorate:
            table_dir = decorate_dir

        rogue_out = defaultdict(int)
        rogue_out_rank = defaultdict(lambda: defaultdict(int))
        rogue_in = defaultdict(lambda: defaultdict(int))
        for table_file in os.listdir(table_dir):
            if not table_file.endswith('-table'):
                continue

            rogue_out_in_tree = set()
            with open(os.path.join(table_dir, table_file)) as f:
                header = f.readline().strip().split('\t')

                rogue_out_index = header.index('Rogue out')
                rogue_in_index = header.index('Rogue in')

                for line in f:
                    line_split = line.rstrip('\r\n').split('\t')

                    taxon = line_split[0]

                    for gid in line_split[rogue_out_index].split(','):
                        rogue_out_in_tree.add(gid)
                        rank_label = Taxonomy.rank_labels[Taxonomy.rank_index[taxon[0:3]]]
                        rogue_out_rank[gid][rank_label] += 1

                    for gid in line_split[rogue_in_index].split(','):
                        rogue_in[gid][taxon] += 1

            for gid in rogue_out_in_tree:
                rogue_out[gid] += 1

        # calculate rogue index
        self.logger.info('Calculating rogue index.')
        rogue_index = defaultdict(float)
        for gid in genome_count:
            r_index = float(rogue_out.get(gid, 0)) / genome_count[gid]
            rogue_index[gid] = r_index

        # write out results
        fout = open(os.path.join(output_dir, 'rogue_index.tsv'), 'w')
        fout.write('Accession\tGTDB taxonomy')
        fout.write('\tRogue out count\tRogue out index')
        for rank in Taxonomy.rank_labels:
            fout.write('\tRogue out (%s)' % rank)
        fout.write('\tRogue in count\tRogue in taxonomic summary\n')
        rogue_index_sorted = sorted(rogue_index.items(), key=operator.itemgetter(1), reverse=True)
        for gid, r_index in rogue_index_sorted:
            rogue_in_summary = []
            rogue_in_count = 0
            if gid in rogue_in:
                rogue_in_sorted = sorted(rogue_in[gid].items(), key=operator.itemgetter(1), reverse=True)
                for taxon, count in rogue_in_sorted:
                    rogue_in_summary.append('%s: %d' % (taxon, count))
                    rogue_in_count += count

            fout.write('%s\t%s\t%d\t%.3f' % (gid,
                                             '; '.join(taxonomy[gid]),
                                             rogue_out.get(gid, 0),
                                             r_index))

            for rank_label in Taxonomy.rank_labels:
                if gid in rogue_out_rank:
                    fout.write('\t%.3f' % (float(rogue_out_rank[gid].get(rank_label, 0)) / genome_count[gid]))
                else:
                    fout.write('\t0')

            fout.write('\t%d\t%s\n' % (rogue_in_count,
                                       ', '.join(rogue_in_summary)))
        fout.close()
