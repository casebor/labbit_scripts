from pymol import cmd, stored

def sirah_memb(sel1):
    """
DESCRIPTION

    Makes a visualization of the sirah ff for membranes

USAGE

    sirah_memb selection

EXAMPLES

    sirah_memb memProtName
    """
    # get the names of the proteins in the selection

    protein = "{} and (resn sA or resn sR or resn sN or resn sD or resn sDh or \
               resn sC or resn sX or resn sQ or resn sE or resn sEh or resn sG \
               or resn sHd or resn sHe or resn sI or resn sL or resn sK or resn \
               sM or resn sF or resn sP or resn sS or resn sT or resn sW or resn \
               sY or resn sV)".format(sel1)
    nucleic_acid = "{} and (resn DAX or resn AX3 or resn AX5 or resn DTX or resn \
                    TX3 or resn TX5 or resn DGX or resn GX3 or resn GX5 or resn \
                    DCX or resn CX3 or resn CX5)".format(sel1)
    membrane = "{} and (resn xPA or resn xPC or resn xOL)".format(sel1)
    ions = "{} and (resn KW or resn NaW or resn ClW or resn CaX)".format(sel1)
    water = "{} and (resn WT4 or resn WLS)".format(sel1)
    heads = "{} and (name BFO or name BCO or name BPE or name BPSO or name BPSN) \
             and not {}_prot".format(sel1, sel1)
    cmd.select("{}_prot".format(sel1), protein)
    cmd.select("{}_na".format(sel1), nucleic_acid)
    cmd.select("{}_memb".format(sel1), membrane)
    cmd.select("{}_ions".format(sel1), ions)
    cmd.select("{}_water".format(sel1), water)
    cmd.select("{}_heads".format(sel1), heads)
    cmd.hide("everything", sel1)
    cmd.show("sticks", protein)
    cmd.show("sticks", nucleic_acid)
    cmd.show("sticks", membrane)
    cmd.show("spheres", heads)
    cmd.color("gray50", membrane)
    cmd.color("tv_red", water)
    cmd.color("tv_blue", "{} and name BCO".format(heads))
    cmd.color("orange", "{} and name BFO".format(heads))
    cmd.color("firebrick", "{} and name BPSO".format(heads))
    cmd.color("deepblue", "{} and name BPSN".format(heads))
    cmd.color("magenta", "{} and name NaW".format(ions))
    cmd.color("limon", "{} and name ClW".format(ions))
    cmd.zoom(sel1)

# let pymol know about the function
cmd.extend("sirah_memb", sirah_memb)


