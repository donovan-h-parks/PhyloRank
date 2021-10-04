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

import dendropy
from biolib.common import is_float

from phylorank.newick import parse_label, create_label

VIRAL_RANK_LABELS = ['phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species']
VIRAL_RANK_PREFIXES = ['P__', 'C__', 'O__', 'F__', 'f__', 'G__', 'S__']
VIRAL_PREFIX_TRANSLATION = {'P__': 'd__',
                            'C__': 'p__',
                            'O__': 'c__',
                            'F__': 'o__',
                            'f__': 'f__',
                            'G__': 'g__',
                            'S__': 's__'}
VIRAL_PREFIX_TRANSLATION_REV = {'d__': 'P__',
                                'p__': 'C__',
                                'c__': 'O__',
                                'o__': 'F__',
                                'f__': 'f__',
                                'g__': 'G__',
                                's__': 'S__'}


def translate_viral_taxonomy(taxonomy):
    """Translate prefixes of viral taxonomy to prokaryotic prefixes."""

    translated = {}
    for gid, taxa in taxonomy.items():
        translated_taxa = []
        for taxon in taxa:
            prefix = taxon[0:3]
            if prefix not in VIRAL_PREFIX_TRANSLATION:
                print('Unrecognized viral prefix for {}: {}'.format(gid, prefix))
                sys.exit(1)

            translated_taxa.append(taxon.replace(prefix, VIRAL_PREFIX_TRANSLATION[prefix]))

        translated[gid] = translated_taxa

    return translated


def translate_viral_tree(tree):
    """Translate prefixes of viral taxonomy in tree to prokaryotic prefixes."""

    if isinstance(tree, str):
        tree = dendropy.Tree.get_from_path(tree,
                                           schema='newick',
                                           rooting="force-rooted",
                                           preserve_underscores=True)

    for node in tree.preorder_node_iter():
        if not node.label or node.is_leaf():
            continue

        support, taxa, auxiliary_info = parse_label(node.label)
        if not taxa:
            continue

        translated_taxa = []
        for taxon in [t.strip() for t in taxa.split(';')]:
            prefix = taxon[0:3]
            if prefix not in VIRAL_PREFIX_TRANSLATION:
                print('Unrecognized viral prefix for {}: {}'.format(taxon, prefix))
                sys.exit(1)

            translated_taxa.append(taxon.replace(prefix, VIRAL_PREFIX_TRANSLATION[prefix]))

        taxa_str = ';'.join(translated_taxa)
        node.label = create_label(support, taxa_str, auxiliary_info)


def rev_translate_output_file(output_file):
    """Translate output file to viral prefixes."""

    if not os.path.exists(output_file):
        print('File does not exist: {}'.format(output_file))
        sys.exit(1)

    with open(output_file) as f:
        content = f.read()

        for prefix, viral_prefix in VIRAL_PREFIX_TRANSLATION_REV.items():
            content = content.replace(prefix, viral_prefix)

    fout = open(output_file, 'w')
    fout.write(content)
    fout.close()


def read_viral_taxonomy_from_tree(tree):
    """Obtain the taxonomy for each extant taxa as specified by internal tree labels.

    Parameters
    ----------
    tree : str or dendropy.Tree
        Filename of newick tree or dendropy tree object.

    Returns
    -------
    dict : d[unique_id] -> [d__<taxon>, ..., s__<taxon>]
        Taxa indexed by unique ids.
    """

    if isinstance(tree, str):
        tree = dendropy.Tree.get_from_path(tree,
                                           schema='newick',
                                           rooting="force-rooted",
                                           preserve_underscores=True)

    taxonomy = {}
    for leaf in tree.leaf_node_iter():
        taxa = []

        node = leaf.parent_node
        while node:
            if node.label:
                taxa_str = node.label
                if ':' in taxa_str:
                    taxa_str = taxa_str.split(':')[1]

                if not is_float(taxa_str):
                    if taxa_str[-1] == ';':
                        taxa_str = taxa_str[:-1]

                    # appears to be an internal label and not simply a support value
                    taxa = [x.strip() for x in taxa_str.split(';')] + taxa
            node = node.parent_node

        if len(taxa) > 7:
            logger = logging.getLogger()
            logger.error('Invalid taxonomy string read from tree for taxon %s: %s' % (leaf.taxon.label, ';'.join(taxa)))
            sys.exit(1)

        for taxon in taxa:
            prefix = taxon[0:3]
            if prefix not in VIRAL_PREFIX_TRANSLATION:
                print('Unrecognized viral prefix for {}: {}'.format(taxon, prefix))
                sys.exit(1)

        # fill missing ranks
        try:
            last_rank = VIRAL_RANK_PREFIXES.index(taxa[-1][0:3])
        except:
            logger = logging.getLogger()
            logger.error('Taxon {} is missing rank prefix: {}'.format(leaf.taxon.label, ';'.join(taxa)))
            sys.exit(1)

        for i in range(last_rank + 1, len(VIRAL_RANK_PREFIXES)):
            taxa.append(VIRAL_RANK_PREFIXES[i])

        taxonomy[leaf.taxon.label] = taxa

    return taxonomy


def sort_viral_taxa(taxa, reverse=False):
    """Sort taxa by rank and then alphabetically.

    Parameters
    ----------
    taxa_list : iterable
        Taxa with rank prefixes.
    
    Returns
    -------
    list
        Taxa sorted by rank and alphabetically within each rank.
    """

    ordered_taxa = []
    for rank_prefix in VIRAL_RANK_PREFIXES:
        rank_taxa = []
        for taxon in taxa:
            if taxon.startswith(rank_prefix):
                rank_taxa.append(taxon)

        ordered_taxa.extend(sorted(rank_taxa))

    if reverse:
        ordered_taxa = ordered_taxa[::-1]

    return ordered_taxa
