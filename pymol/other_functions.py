#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>

from pymol import cmd, stored

def coords(sel1="all", state1=1):
        """
DESCRIPTION

    Prints the coordinates of the selected elements

USAGE

    coords selection

EXAMPLES

    coords sele
    coords active_site_residues, 1
    coords protein
        """
        cmd.iterate_state(state1, sel1, 'print("{:<4}{:<7}[{:8.3f}, {:8.3f}, {:8.3f}]".format(name, ID, x, y, z))')

# let pymol know about the function
cmd.extend("coords", coords)
