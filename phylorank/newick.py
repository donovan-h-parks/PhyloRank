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


from biolib.common import is_float

from skbio import TreeNode


'''Helper functions for parsing Newick information.'''


def parse_label(label):
    """Parse a Newick label which may contain a support value, taxon, and/or auxiliary information.

    Parameters
    ----------
    label : str
        Internal label in a Newick tree.

    Returns
    -------
    float
        Support value specified by label, or None
    str
        Taxon specified by label, or None
    str
        Auxiliary information, on None
    """

    support = None
    taxon = None
    auxiliary_info = None

    if '|' in label:
        label, auxiliary_info = label.split('|')

    if ':' in label:
        support, taxon = label.split(':')
        support = float(support)
    else:
        if is_float(label):
            support = float(label)
        else:
            taxon = label

    return support, taxon, auxiliary_info


def read_from_tree(tree, modify_tree=True):
    """Obtain the taxonomy for extant taxa as specified by internal tree labels.

    Special care is taken for species names to ensure they contain both a
    'specific epithet' (genus) and a 'specific name' (species). The
    modify_tree flag indicates if node names should be updated to ensure
    binomial species names.

    Parameters
    ----------
    tree : str or TreeNode
        Filename of newick tree or TreeNode object.
    modify_tree : bool
        Flag indicating if tree should be modified.

    Returns
    -------
    dict : d[unique_id] -> [d__<taxon>, ..., s__<taxon>]
        Taxa indexed by unique ids.
    """

    if not isinstance(tree, TreeNode):
        tree = TreeNode.read(tree, convert_underscores=False)

    taxonomy = {}
    for leaf in tree.tips():
        taxa = []

        node = leaf.parent
        species_node = None
        while node:
            if node.name:
                _support, taxon, _auxiliary_info = parse_label(node.name)

                if taxon:
                    taxa = [x.strip() for x in taxon.split(';')] + taxa

                    if 's__' in taxon:
                        species_node = node
            node = node.parent

        # check if genus name should be appended to species label
        if len(taxa) == 7:
            genus = taxa[5][3:]
            species = taxa[6][3:]
            if genus not in species:
                binomial_name = 's__' + genus + ' ' + species
                taxa[6] = binomial_name

                if modify_tree:
                    species_node.name = binomial_name

        taxonomy[leaf.name] = taxa

    return taxonomy
