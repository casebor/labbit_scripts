#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import numpy as np

# For each atom, we have (atoms radius, covalent radius and Van der Waals radius)
RADII = {"C": (0.67, 0.76, 1.7), "H": (0.53, 0.31, 1.0), "O": (0.48, 0.66, 1.52),
         "N": (0.56, 0.71, 1.55), "P": (0.98, 1.07, 1.8), "S": (0.87, 1.05, 1.8),
         "Na": (1.9, 1.6, 2.27), "Mg": (1.45, 1.41, 1.73), "K": (2.43, 2.03, 2.75),
         "Ca": (1.94, 1.76, 2.0), "F": (0.42, 0.57, 1.47), "Cl": (0.79, 1.02, 1.75),
         "Br": (0.94, 1.2, 1.85), "I": (1.15, 1.39, 1.98)}

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
            pos = [0.0, 0.0, 0.0]
        self.name = name
        self.index = index
        self.pos = np.array(pos)
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
    
    def x(self):
        return self.pos[0]
    
    def y(self):
        return self.pos[1]
    
    def z(self):
        return self.pos[2]
    
    def radius(self):
        return self.radius
    
