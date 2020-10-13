#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>

from pymol import cmd, stored

def coords(sel1="all", state1=1):
    """
DESCRIPTION

    Prints the coordinates of the selected elements. The frame can be specified with the frame number after a comma

USAGE

    coords selection
    coords selection, frame

EXAMPLES

    coords sele
    coords active_site_residues, 1
    coords protein
    """
    cmd.iterate_state(state1, sel1, 'print("{:<4}{:<7}[{:8.3f}, {:8.3f}, {:8.3f}]".format(name, ID, x, y, z))')


def pyrama(sel1):
    from pyrama import calc_ramachandran, plot_ramachandran
    """
DESCRIPTION

    Generates a ramachandran plot using pyrama script.
    Written by Gábor Erdős, 2017
    Contact info: gerdos[at]caesar.elte.hu

    The preferences were calculated from the following artice:
    Lovell et al. Structure validation by Calpha geometry: phi,psi and Cbeta deviation. 2003
    DOI: 10.1002/prot.10286

USAGE

    pyrama selection

EXAMPLES

    pyrana protName
    """
    cmd.save("/tmp/{}.pdb".format(sel1), selection=sel1, state=-1, format="pdb")
    normals, outliers = calc_ramachandran(["/tmp/{}.pdb".format(sel1)])
    plot_ramachandran(normals, outliers)

# let pymol know about the function
cmd.extend("coords", coords)
cmd.extend("pyrama", pyrama)
