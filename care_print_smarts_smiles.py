#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>

import sys
from rdkit import Chem

def main(sdffile):
    """ Main function """
    sdf_mol = Chem.SDMolSupplier(sdffile, removeHs=False)
    # print("SMART", sdf_mol[0].GetProp("SMARTS"))
    # print("Molecule:", sdf_mol.GetItemText(0))
    m_smarts = Chem.rdmolfiles.MolToSmarts(sdf_mol[0])
    m_smiles = Chem.rdmolfiles.MolToSmiles(sdf_mol[0], canonical=False)
    print("Molecule smarts:", m_smarts)
    print("Molecule smiles:", m_smiles)
    q_smarts = "[#17]-[#6&!a]"
    patt = Chem.MolFromSmarts(q_smarts)
    print("Has SMARTS:", q_smarts, "?:", sdf_mol[0].HasSubstructMatch(patt))

if __name__ == '__main__':
    # python3 sdffile
    # python3 bromocyclohexane.sdf
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print("Usage:\n\t python3 bromocyclohexane.sdf")
        quit(0)
    main(sys.argv[1])
