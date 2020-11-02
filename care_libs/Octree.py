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
    def __init__ (self, elements=None, origin=None, border=None):
        if elements is None:
            elements = []
        self.to_insert = elements
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
        self.atoms = []
        self.atoms_pos = None
        self.children = [None] * 8
        self.is_leaf = False
    
    def initialize(self):
        if abs(self.border - self.origin[0]) <= 8.0 or \
           abs(self.border - self.origin[1]) <= 8.0 or \
           abs(self.border - self.origin[2]) <= 8.0:
            self.is_leaf = True
        if self.is_leaf:
            self.atoms = [e for e in self.to_insert]
            # self.calculate_collisions()
        else:
            while self.to_insert:
                atom = self.to_insert.pop()
                if abs(atom.x() - self.origin[0]) < atom.radius():
                    self.atoms.append(atom)
                elif abs(atom.y() - self.origin[1]) < atom.radius():
                    self.atoms.append(atom)
                elif abs(atom.z() - self.origin[2]) < atom.radius():
                    self.atoms.append(atom)
                else:
                    self.add_atom_to_children(atom)
        self.atoms_pos = np.zeros((len(self.atoms), 3))
        for i, atom in enumerate(self.atoms):
            self.atoms_pos[i,:] = atom.coords()
    
    def add_atom_to_children(self, atom):
        """ The position of the children will depend on their positions wrt
            the 'origin' of the tree. Using the right hand rule, with the
            middle finger as X axis, the thumb as Y axis and the index as
            Z axis, we have eight children distributed as:
            children   0   1   2   3   4   5   6   7
                X      +   +   -   -   +   +   -   -
                y      +   +   +   +   -   -   -   -
                Z      +   -   -   +   +   -   -   +
        """
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
        # self.children = [child for child in self.children if child is not None]
        for child in self.children:
            if child is not None:
                child.initialize()
                child.initialize_children()
    
    def add_temp(self, atom):
        self.to_insert.append(atom)
    
    def calculate_collisions(self):
        for i in range(self.atoms_pos.shape[0] - 1):
            dists = self.atoms_pos[i] - self.atoms_pos[i+1:]
            dists *= dists
            dists = dists.sum(axis=1)
            for j, d in enumerate(dists):
                if d <= 4.41:
                    pos = i + j + 1
                    self.atoms[i].add_connection(self.atoms[pos])
        box_atoms = self._get_atoms_for_contact(self.origin)
        box_coords = np.zeros((len(box_atoms), 3))
        for i, atom in enumerate(box_atoms):
            box_coords[i,:] = atom.coords()
        for atom in self.atoms:
            dists = atom.coords() - box_coords
            dists *= dists
            dists = dists.sum(axis=1)
            for i, d in enumerate(dists):
                if d <= 4.41:
                    atom.add_connection(box_atoms[i])
    
    def _get_atoms_for_contact(self, reference_point):
        if self.is_leaf:
            return self.atoms
        else:
            to_send = [e for e in self.atoms]
            indices_to_search = self._get_child_indexes_for_contact(reference_point)
            for i in indices_to_search:
                if self.children[i] is not None:
                    to_send.extend(self.children[i]._get_atoms_for_contact(reference_point))
            return to_send
    
    def _get_child_indexes_for_contact(self, reference_point):
        point = reference_point - self.origin
        to_send = set()
        if point[0] > 0:
            to_send.update([0,1,4,5])
        else:
            to_send.update([2,3,6,7])
        if point[1] > 0:
            to_send.update([0,1,2,3])
        else:
            to_send.update([4,5,6,7])
        if point[2] > 0:
            to_send.update([0,3,4,7])
        else:
            to_send.update([1,2,5,6])
        return to_send
    

