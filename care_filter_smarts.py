#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>

import os
import sys
from rdkit import Chem
from rdkit.Chem import rdFMCS

def main(smarts):
    """ Main function """
    patt = Chem.MolFromSmarts(smarts)
    for filename in os.listdir('.'):
        if filename.endswith('.sdf'):
            sdf_mol = Chem.SDMolSupplier(filename, removeHs=False)
            if not sdf_mol[0].HasSubstructMatch(patt):
                # print("Molecule smarts:", Chem.rdmolfiles.MolToSmarts(sdf_mol[0]))
                print("Molecule:", filename, "is invalid!")
            else:
                print("Molecule:", filename, "is valid")

if __name__ == '__main__':
    # python3 smarts
    # python3 '[#17]-[#6&!a]'
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print("Usage:\n\t care_filter_smarts.py '[#17]-[#6&!a]'")
        quit(0)
    main(sys.argv[1])
