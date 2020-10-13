#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import numpy as np

RADII = {"C": (1.0, 2.0, 3.0), "H": (1.0, 2.0, 3.0), "O": (1.0, 2.0, 3.0),
         "N": (1.0, 2.0, 3.0), "P": (1.0, 2.0, 3.0), "S": (1.0, 2.0, 3.0)}

def _get_radius(name):
    pass

def _get_vdw_rad(name):
    pass

def _get_cov_rad(name):
    pass

class Atom:
    """ Class doc """
    def __init__ (self, name='X', index=None, pos=None, res_index=None,
                  res_name=None, chain='X', atomic_number=0, occupancy=0.0,
                  bfactor=0.0, charge=0.0, bonds_indices=None):
        if pos is None:
            pos = (0.0, 0.0, 0.0)
        self.name = name
        self.index = index
        self.pos = np.array((pos[0], pos[1], pos[2]))
        self.res_index = res_index
        self.res_name = res_name
        self.chain = chain
        self.atomic_number = atomic_number
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.charge = charge
        self.radius = _get_radius(name)
        self.vdw_rad = _get_vdw_rad(name)
        self.cov_rad = _get_cov_rad(name)
    
    def coords(self):
        """ Function doc """
        return self.pos
    

