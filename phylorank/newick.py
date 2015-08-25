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

'''Helper functions for parsing Newick information.'''


def parse_label(label):
    """Parse a Newick label which may contain a support value, taxon information, or both.

    Parameters
    ----------
    label : str
        Internal label in a Newick tree.

    Returns
    -------
    support : float
        Support value specified by label, or None
    taxon : str
        Taxon specified by label, or Non
    """

    support = None
    taxon = None
    if ':' in label:
        support, taxon = label.split(':')
        support = float(support)
    else:
        if is_float(label):
            support = float(label)
        else:
            taxon = label

    return support, taxon
