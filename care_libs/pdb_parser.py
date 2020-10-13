#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import numpy as np
from care_libs import Atom

def parse_pdb(pdb_file):
    with open(pdb_file, "r") as pdbstream:
        pdb_list = list(pdbstream)
    i = 0
    atoms = []
    geom_center = np.zeros((3))
    min_x, max_x = 9999.9, -9999.9
    min_y, max_y = 9999.9, -9999.9
    min_z, max_z = 9999.9, -9999.9
    while i < len(pdb_list):
        if pdb_list[i].startswith("ATOM"):
            index = int(pdb_list[i][6:11])
            name = pdb_list[i][12:16]
            res_name = pdb_list[i][17:20]
            chain = pdb_list[i][21]
            res_index = int(pdb_list[i][22:26])
            x = float(pdb_list[i][30:38])
            if x < min_x:
                min_x = x
            if x > max_x:
                max_x = x
            y = float(pdb_list[i][38:46])
            if y < min_y:
                min_y = y
            if y > max_y:
                max_y = y
            z = float(pdb_list[i][46:54])
            if z < min_z:
                min_z = z
            if z > max_z:
                max_z = z
            occupancy = float(pdb_list[i][54:60])
            bfactor = float(pdb_list[i][60:66])
            atom = Atom.Atom(name=name, index=index, pos=(x,y,z),
                             res_index=res_index, res_name=res_name,
                             chain=chain, occupancy=occupancy, bfactor=bfactor)
            atoms.append(atom)
            geom_center += atom.pos
        i += 1
    geom_center /= len(atoms)
    print(geom_center)
    print(min_x, max_x)
    print(min_y, max_y)
    print(min_z, max_z)
