#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import numpy as np

class Octree():
    """ Class doc """
    # def __init__ (self, atoms, origin=None, top_left_front=None,
    #               bottom_right_back=None):
    def __init__ (self, atoms=None, origin=None, border=None):
        if atoms is None:
            atoms = []
        self.to_insert = atoms
        if origin is None:
            origin = [0., 0., 0.]
        self.origin = np.array(origin)
        # if top_left_front is None:
        #     top_left_front = [999., 999., 999]
        # self.top_left_front = top_left_front
        # if bottom_right_back is None:
        #     bottom_right_back = [-999., -999., -999.]
        if border is None:
            border = 512.
        self.border = int(border) + 1
        # self.bottom_right_back = bottom_right_back
        self.elements = []
        self.children = [None] * 8
        self.is_leaf = False
    
    def initialize(self):
        if abs(self.border - self.origin[0]) <= 4.0 or \
           abs(self.border - self.origin[1]) <= 4.0 or \
           abs(self.border - self.origin[2]) <= 4.0:
            self.is_leaf = True
        if self.is_leaf:
            self.elements = [e for e in self.to_insert]
            # self.calculate_collisions()
        else:
            while self.to_insert:
                atom = self.to_insert.pop()
                if abs(atom.x() - self.origin[0]) < atom.radius():
                    self.elements.append(atom)
                elif abs(atom.y() - self.origin[1]) < atom.radius():
                    self.elements.append(atom)
                elif abs(atom.z() - self.origin[2]) < atom.radius():
                    self.elements.append(atom)
                else:
                    self.add_atom_to_children(atom)
    
    def add_atom_to_children(self, atom):
        child_pos = np.abs(atom.coords() - self.origin)
        child_pos = np.array((atom.coords() - self.origin)/child_pos, dtype=np.int)
        child_pos = np.array(child_pos+1, dtype=np.bool)
        if child_pos[0] and child_pos[1] and child_pos[2]:
            if self.children[0] is None:
                _origin = self.origin + self.border
                _border = self.border / 2.
                self.children[0] = Octree(origin=_origin, border=_border)
            self.children[0].add_temp(atom)
        elif child_pos[0] and child_pos[1] and not child_pos[2]:
            if self.children[1] is None:
                _origin = self.origin + [self.border, self.border, -self.border]
                _border = self.border / 2.
                self.children[1] = Octree(origin=_origin, border=_border)
            self.children[1].add_temp(atom)
        elif not child_pos[0] and child_pos[1] and not child_pos[2]:
            if self.children[2] is None:
                _origin = self.origin + [-self.border, self.border, -self.border]
                _border = self.border / 2.
                self.children[2] = Octree(origin=_origin, border=_border)
            self.children[2].add_temp(atom)
        elif not child_pos[0] and child_pos[1] and child_pos[2]:
            if self.children[3] is None:
                _origin = self.origin + [-self.border, self.border, self.border]
                _border = self.border / 2.
                self.children[3] = Octree(origin=_origin, border=_border)
            self.children[3].add_temp(atom)
        elif child_pos[0] and not child_pos[1] and child_pos[2]:
            if self.children[4] is None:
                _origin = self.origin + [self.border, -self.border, self.border]
                _border = self.border / 2.
                self.children[4] = Octree(origin=_origin, border=_border)
            self.children[4].add_temp(atom)
        elif child_pos[0] and not child_pos[1] and not child_pos[2]:
            if self.children[5] is None:
                _origin = self.origin + [self.border, -self.border, -self.border]
                _border = self.border / 2.
                self.children[5] = Octree(origin=_origin, border=_border)
            self.children[5].add_temp(atom)
        elif not child_pos[0] and not child_pos[1] and not child_pos[2]:
            if self.children[6] is None:
                _origin = self.origin + [-self.border, -self.border, -self.border]
                _border = self.border / 2.
                self.children[6] = Octree(origin=_origin, border=_border)
            self.children[6].add_temp(atom)
        elif not child_pos[0] and not child_pos[1] and child_pos[2]:
            if self.children[7] is None:
                _origin = self.origin + [-self.border, -self.border, self.border]
                _border = self.border / 2.
                self.children[7] = Octree(origin=_origin, border=_border)
            self.children[7].add_temp(atom)
    
    def initialize_children(self):
        # assert len(self.to_insert) == 0
        self.children = [child for child in self.children if child is not None]
        for child in self.children:
            child.initialize()
            child.initialize_children()
    
    def add_temp(self, atom):
        self.to_insert.append(atom)
    
    def calculate_collisions(self):
        pass
    
    

