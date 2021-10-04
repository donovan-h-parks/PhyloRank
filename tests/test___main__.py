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
import shutil
import subprocess
import tempfile
import unittest


class TestPhyloRank(unittest.TestCase):

    def setUp(self):
        self.dir_tmp = tempfile.mkdtemp(prefix='phylorank_tmp_')
        self.path_ar122_r95_tree = os.path.join(os.getcwd(), 'data', 'ar122_r95.tree')
        self.path_ar122_tax_r95 = os.path.join(os.getcwd(), 'data', 'ar122_taxonomy_r95.tsv')

    def tearDown(self):
        shutil.rmtree(self.dir_tmp)

    def test_outliers(self):
        args = ['python', '-m', 'phylorank', 'outliers',
                self.path_ar122_r95_tree, self.path_ar122_tax_r95, self.dir_tmp]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_scale_tree(self):
        tree_out = os.path.join(self.dir_tmp, 'output.tree')
        args = ['python', '-m', 'phylorank', 'scale_tree',
                self.path_ar122_r95_tree, tree_out]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_compare_red(self):
        outlier_dir = os.path.join(self.dir_tmp, 'outliers')
        args = ['python', '-m', 'phylorank', 'outliers',
                self.path_ar122_r95_tree, self.path_ar122_tax_r95, outlier_dir]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

        red_table1 = os.path.join(outlier_dir, 'ar122_r95.tsv')
        red_table2 = red_table1
        red_dict2 = os.path.join(outlier_dir, 'ar122_r95.dict')
        output_table = os.path.join(self.dir_tmp, 'output_table.tsv')

        args = ['python', '-m', 'phylorank', 'compare_red',
                red_table1, red_table2, red_dict2, output_table]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_mark_tree(self):
        tree_out = os.path.join(self.dir_tmp, 'output.tree')
        args = ['python', '-m', 'phylorank', 'mark_tree',
                self.path_ar122_r95_tree, tree_out]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_rogue_test(self):
        input_tree_dir = os.path.join(self.dir_tmp, 'input_trees')
        os.makedirs(input_tree_dir)
        shutil.copy(self.path_ar122_r95_tree, input_tree_dir)
        dir_out = os.path.join(self.dir_tmp, 'output')
        args = ['python', '-m', 'phylorank', 'rogue_test',
                input_tree_dir, self.path_ar122_tax_r95, dir_out]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_decorate(self):
        output_tree = os.path.join(self.dir_tmp, 'output.tree')
        args = ['python', '-m', 'phylorank', 'decorate',
                self.path_ar122_r95_tree, self.path_ar122_tax_r95, output_tree]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_taxon_stats(self):
        output_file = os.path.join(self.dir_tmp, 'output.tsv')
        args = ['python', '-m', 'phylorank', 'taxon_stats',
                self.path_ar122_tax_r95, output_file]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_rank_res(self):
        output_file = os.path.join(self.dir_tmp, 'output.tsv')
        args = ['python', '-m', 'phylorank', 'rank_res',
                self.path_ar122_r95_tree, self.path_ar122_tax_r95, output_file]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_bl_dist(self):
        args = ['python', '-m', 'phylorank', 'bl_dist',
                self.path_ar122_r95_tree, self.dir_tmp]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_bl_optimal(self):
        output_file = os.path.join(self.dir_tmp, 'output.tsv')
        args = ['python', '-m', 'phylorank', 'bl_optimal',
                self.path_ar122_r95_tree, '3', output_file]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_bl_decorate(self):
        output_tree = os.path.join(self.dir_tmp, 'output.tree')
        args = ['python', '-m', 'phylorank', 'bl_decorate',
                self.path_ar122_r95_tree, self.path_ar122_tax_r95, '0.5',
                '3', output_tree]
        proc = subprocess.Popen(args, encoding='utf-8')
        proc.communicate()
        self.assertEqual(proc.returncode, 0)

    def test_bl_table(self):
        pass
