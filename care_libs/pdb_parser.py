#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import numpy as np
from care_libs import Atom

def parse_pdb(pdb_file, with_hex=False):
    with open(pdb_file, "r") as pdbstream:
        pdb_list = list(pdbstream)
    i = 0
    atoms = []
    geom_center = np.zeros((3))
    min_x, max_x = 9999.9, -9999.9
    min_y, max_y = 9999.9, -9999.9
    min_z, max_z = 9999.9, -9999.9
    hexflags = [False, False]
    num_at_ids = 0
    num_res_ids = 0
    while i < len(pdb_list):
        if pdb_list[i].startswith("ATOM"):
            if num_at_ids == 99999:
                hexflags[0] = True
            if hexflags[0]:
                index = int(pdb_list[i][6:11], 16)
            else:
                index = int(pdb_list[i][6:11])
            num_at_ids += 1
            if num_res_ids == 9999:
                hexflags[1] = True
            if hexflags[1]:
                res_index = int(pdb_list[i][22:26], 16)
            else:
                res_index = int(pdb_list[i][22:26])
            num_res_ids += 1
            name = pdb_list[i][12:16].strip()
            res_name = pdb_list[i][17:20].strip()
            chain = pdb_list[i][21].strip()
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
    return atoms, geom_center, min_x, max_x, min_y, max_y, min_z, max_z
