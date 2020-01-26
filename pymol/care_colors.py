#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>

from pymol import cmd

def care_color_HH(selection):
    """
DESCRIPTION

    Changes the color scheme of the amino acids depending on the 
    physicochemic properties: Acid, Basic, Polar and Non-Polar

USAGE

    care_color_HH selection

EXAMPLES

    care_color_HH 5DO7
    """
    acids = "{} and (resn ASP or resn ASH or resn GLU or resn GLH)".format(selection)
    bases = "{} and (resn HIS or resn HIE or resn HID or resn HIP or resn ARG \
             or resn LYS or resn LYN)".format(selection)
    polars = "{} and (resn CYS or resn CYX or resn GLN or resn ASN or resn SER \
             or resn TYR or resn THR)".format(selection)
    nonpolars = "{} and (resn GLY or resn ALA or resn LEU or resn ILE or resn PHE \
                 or resn TRP or resn MET or resn PRO or resn VAL)".format(selection)
    cmd.color("firebrick", acids)
    cmd.color("deepteal", bases)
    cmd.color("tv_orange", polars)
    cmd.color("smudge", nonpolars)
    util.cnc(selection)
# let pymol know about the function
cmd.extend("care_color_HH", care_color_HH)

def care_color_membrane(selection='all'):
    """
DESCRIPTION

    Changes the color scheme of a membrane, by default it selects all.
    Currently only supports POPC, POPE and POPA lipids.
    Adding more lipids should be done easily.

USAGE

    care_color_membrane selection

EXAMPLES

    care_color_membrane
    care_color_membrane system_MEMBR
    """
    membrane = "{} and (resn CHL or resn LA or resn MY or resn OL or resn PA \
                or resn PC or resn PE or resn POPC)".format(selection)
    cmd.color("gray50", membrane)
    util.cnc(selection)
# let pymol know about the function
cmd.extend("care_color_membrane", care_color_membrane)
