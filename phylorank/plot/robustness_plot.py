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
from collections import defaultdict

import dendropy
from biolib.plots.abstract_plot import AbstractPlot
from numpy import (mean as np_mean,
                   arange as np_arange,
                   percentile as np_percentile)

from phylorank.rel_dist import RelativeDistance


class RobustnessPlot(AbstractPlot):
    """Plot relative distance of named groups across a set of trees."""

    def __init__(self):
        """Initialize."""
        AbstractPlot.__init__(self, None)

    def rel_dist_to_specified_groups(self, tree_file, groups_to_consider, groups):
        """Determine relative distance to specified named clades.

        Parameters
        ----------
        tree_file : str
          File containing a tree in Newick format.
        groups_to_consider: set
          Taxonomic groups to consider.
        groups : d[taxon] -> list of children
          Children within named taxonomic groups.

        Returns
        -------
        dict : d[taxon] -> relative distance to root
        """

        tree = dendropy.Tree.get_from_path(tree_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        # calculate relative distance for all nodes
        rd = RelativeDistance()
        rd.decorate_rel_dist(tree)

        # gather information for nodes of interest
        rel_dists_to_taxon = {}
        dist_components_taxon = {}
        polyphyletic = set()
        for taxon, taxa_ids in groups.items():
            if taxon not in groups_to_consider:
                continue

            tips = []
            for t in taxa_ids:
                try:
                    tip = tree.find(t)
                    tips.append(tip)
                except:
                    continue

            if len(tips) == 0:
                # group is within the phylum removed from the tree
                continue

            lca_node = tree.lca(tips)

            if len(list(lca_node.tips())) != len(tips):
                print('  [Warning] Group is not monophyletic %s' % taxon)
                polyphyletic.add(taxon)
                continue

            # get relative distance from root to named child clade
            rel_dists_to_taxon[taxon] = lca_node.rel_dist
            dist_components_taxon[taxon] = [lca_node.parent.rel_dist, lca_node.length, lca_node.weighted_dist]

        return rel_dists_to_taxon, dist_components_taxon, polyphyletic

    def run(self, rank, input_tree_dir, full_tree_file, derep_tree_file, taxonomy_file,
            output_prefix, min_children, title):

        # determine named clades in full tree
        named_clades = set()
        tree = dendropy.Tree.get_from_path(full_tree_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        for node in tree.preorder_node_iter():
            if node.label:
                taxonomy = node.label.split(';')
                named_clades.add(taxonomy[-1].strip().split(':')[-1])

        print('Identified %d named clades in full tree.' % len(named_clades))

        # determine named groups with at least the specified number of children
        print('Determining taxa with sufficient named children lineages.')
        taxon_children = defaultdict(set)
        groups = defaultdict(list)
        for line in open(taxonomy_file):
            line_split = line.replace('; ', ';').split()
            genome_id = line_split[0]
            taxonomy = [x.strip() for x in line_split[1].split(';')]

            if len(taxonomy) > rank + 1:
                taxon_children[taxonomy[rank]].add(taxonomy[rank + 1])

            if len(taxonomy) > rank:
                groups[taxonomy[rank]].append(genome_id)

        groups_to_consider = set()
        for taxon, children_taxa in taxon_children.items():
            if len(children_taxa) >= min_children and taxon in named_clades:
                groups_to_consider.add(taxon)

        print('Assessing distribution over %d groups.' % len(groups_to_consider))

        # calculate RED for full tree
        print('')
        print('Calculating RED over full tree.')
        tree = dendropy.Tree.get_from_path(full_tree_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        full_rel_dist, _full_dist_components, polyphyletic = self.rel_dist_to_specified_groups(tree,
                                                                                               groups_to_consider,
                                                                                               groups)
        if len(polyphyletic) > 0:
            print('')
            print('[Warning] Full tree contains polyphyletic groups.')

        # calculate RED for dereplicated tree
        print('')
        print('Calculating RED over dereplicated tree.')
        tree = dendropy.Tree.get_from_path(derep_tree_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        derep_rel_dist, derep_dist_components, polyphyletic = self.rel_dist_to_specified_groups(tree,
                                                                                                groups_to_consider,
                                                                                                groups)

        groups_to_consider = groups_to_consider - polyphyletic
        print('Assessing distribution over %d groups after removing polyphyletic groups in original trees.' % len(
            groups_to_consider))

        # calculate RED to each group in each tree
        print('')
        rel_dists = defaultdict(list)
        dist_components = defaultdict(list)
        for f in os.listdir(input_tree_dir):
            if not f.endswith('.rooted.tree'):
                continue

            print(f)

            tree_file = os.path.join(input_tree_dir, f)
            tree = dendropy.Tree.get_from_path(tree_file,
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)

            # calculate relative distance to named taxa
            rel_dist, components, _polyphyletic = self.rel_dist_to_specified_groups(tree, groups_to_consider, groups)

            for taxon, dist in rel_dist.items():
                rel_dists[taxon].append(dist)
                dist_components[taxon].append(components[taxon])

        # create scatter plot
        x = []
        y = []
        xDerep = []
        yDerep = []
        xFull = []
        yFull = []
        perc10 = []
        perc90 = []
        labels = []
        fout = open(output_prefix + '.tsv', 'w')
        fout.write('Taxon\tP10\tP90\tP90-P10\tMean RED\tMean dist to parent\tMean dist to leaves\t'
                   'Original RED\tOrigial dist to parent\tOriginal dist to leaves\n')
        for i, taxon in enumerate(sorted(rel_dists.keys(), reverse=True)):
            labels.append(taxon + ' (%d)' % (len(rel_dists[taxon])))

            rd = rel_dists[taxon]
            for d in rd:
                x.append(d)
                y.append(i + 0.2)

            p10, p90 = np_percentile(rd, [10, 90])
            perc10.append(p10)
            perc90.append(p90)

            mean_x, mean_a, mean_b = np_mean(dist_components[taxon], axis=0)
            derep_x, derep_a, derep_b = derep_dist_components[taxon]
            fout.write('%s\t%.2f\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' % (
            taxon, p10, p90, p90 - p10, mean_x, mean_a, mean_b, derep_x, derep_a, derep_b))

            xDerep.append(derep_rel_dist[taxon])
            yDerep.append(i)

            xFull.append(full_rel_dist[taxon])
            yFull.append(i)
        fout.close()

        self.fig.clear()
        self.fig.set_size_inches(8, len(rel_dists) * 0.4)
        ax = self.fig.add_subplot(111)

        ax.scatter(x, y, alpha=0.5, s=24, c=(0.5, 0.5, 0.5), marker='s')
        ax.scatter(xDerep, yDerep, alpha=1.0, s=24, c=(1.0, 0.0, 0.0), marker='s')
        ax.scatter(xFull, yFull, alpha=1.0, s=24, c=(0.0, 0.0, 1.0), marker='*')

        for i in range(len(labels)):
            ax.plot((perc10[i], perc10[i]), (i, i + 0.4), 'r-')
            ax.plot((perc90[i], perc90[i]), (i, i + 0.4), 'r-')

        # set plot elements
        ax.grid(color=(0.8, 0.8, 0.8), linestyle='dashed')
        if title:
            ax.set_title(title, size=12)

        ax.set_xlabel('relative distance')
        ax.set_xticks(np_arange(0, 1.05, 0.1))
        ax.set_xlim([-0.05, 1.05])

        ax.set_ylabel('taxa')
        ax.set_yticks(range(0, len(rel_dists)))
        ax.set_ylim([-0.2, len(rel_dists) - 0.01])
        ax.set_yticklabels(labels)

        self.prettify(ax)

        # make plot interactive
        # mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(scatter, labels=labels))
        # mpld3.plugins.connect(fig, mpld3.plugins.MousePosition(fontsize=12))

        # mpld3.save_html(fig, output_prefix + '.html')
        self.fig.tight_layout(pad=1)
        self.fig.savefig(output_prefix + '.png', dpi=300)
