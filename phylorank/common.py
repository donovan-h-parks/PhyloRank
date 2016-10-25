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

import sys

from biolib.taxonomy import Taxonomy

from phylorank.newick import parse_label


def is_integer(s):
    """Test if a string represents an integer."""
    try:
        int(s)
        return True
    except ValueError:
        return False


def read_taxa_file(taxa_file):
    """Read taxa from file.

    Parameters
    ----------
    taxa_file : str
        File specifying taxa to consider. One per line.

    Returns
    -------
    set
        Taxa.
    """

    taxa = set()
    for line in open(taxa_file):
        taxa.add(line.rstrip().split('\t')[0])

    return taxa

def filter_taxa_for_dist_inference(tree, taxonomy, trusted_taxa, min_children, min_support):
    """Determine taxa to use for inferring distribution of relative divergences.

    Parameters
    ----------
    tree : Dendropy Tree
        Phylogenetic tree.
    taxonomy : d[taxon ID] -> [d__x; p__y; ...]
        Taxonomy for each taxon.
    trusted_taxa : iterable
        Trusted taxa to consider when inferring distribution.
    min_children : int
        Only consider taxa with at least the specified number of children taxa when inferring distribution.
    min_support : float
        Only consider taxa with at least this level of support when inferring distribution.
    """

    # determine children taxa for each named group
    taxon_children = Taxonomy().taxon_children(taxonomy)

    # get all named groups
    taxa_for_dist_inference = set()
    for taxon_id, taxa in taxonomy.iteritems():
        for taxon in taxa:
            taxa_for_dist_inference.add(taxon)

    # sanity check species names as these are a common problem
    species = set()
    for taxon_id, taxa in taxonomy.iteritems():
        if len(taxa) > Taxonomy.rank_index['s__']:
            species_name = taxa[Taxonomy.rank_index['s__']]
            valid, error_msg = True, None
            if species_name != 's__':
                valid, error_msg = Taxonomy().validate_species_name(species_name, require_full=True, require_prefix=True)
            if not valid:
                print '[Warning] Species name %s for %s is invalid: %s' % (species_name, taxon_id, error_msg)
                continue
                
            species.add(species_name)

    # restrict taxa to those with a sufficient number of named children
    # Note: a taxonomic group with no children will not end up in the
    # taxon_children data structure so care must be taken when applying
    # this filtering criteria.
    if min_children > 0:
        valid_taxa = set()
        for taxon, children_taxa in taxon_children.iteritems():
            if len(children_taxa) >= min_children:
                valid_taxa.add(taxon)

        taxa_for_dist_inference.intersection_update(valid_taxa)

        # explicitly add in the species since they have no
        # children and thus be absent from the taxon_child dictionary
        taxa_for_dist_inference.update(species)

    # restrict taxa used for inferring distribution to those with sufficient support
    if min_support > 0:
        for node in tree.preorder_node_iter():
            if not node.label or node.is_leaf():
                continue

            # check for support value
            support, taxon_name, _auxiliary_info = parse_label(node.label)

            if not taxon_name:
                continue

            if support and float(support) < min_support:
                taxa_for_dist_inference.difference_update([taxon_name])
            elif not support and min_support > 0:
                # no support value, so inform user if they were trying to filter on this property
                print '[Error] Tree does not contain support values. As such, --min_support should be set to 0.'
                continue

    # restrict taxa used for inferring distribution to the trusted set
    if trusted_taxa:
        taxa_for_dist_inference = trusted_taxa.intersection(taxa_for_dist_inference)

    return taxa_for_dist_inference
    
    
def get_phyla_lineages(tree):
    """Get list of phyla level lineages.

    Parameters
    ----------
    tree : Dendropy Tree
        Phylogenetic tree.

    Returns
    -------
    list
        List of phyla level lineages.
    """
    phyla = []
    for node in tree.preorder_node_iter():
        if not node.label or node.is_leaf():
            continue

        _support, taxon_name, _auxiliary_info = parse_label(node.label)
        if taxon_name:
            taxa = [x.strip() for x in taxon_name.split(';')]
            if taxa[-1].startswith('p__'):
                phyla.append(taxa[-1])
                
    return phyla
