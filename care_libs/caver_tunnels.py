#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import numpy as np

class CaverTunnel(object):
    """docstring for CaverTunnel"""
    def __init__(self, id, **kwargs):
        super(CaverTunnel, self).__init__()
        self.id = id
        if "color" in kwargs:
            self.color = kwargs["color"]
        else:
            self.color = [0., 0., 0.]
        if "nodes" in kwargs:
            self.nodes = kwargs["nodes"]
        else:
            self.nodes = []
        if "connections" in kwargs:
            self.connections = kwargs["connections"]
        else:
            self.connections = []
    
    def __str__(self):
        # output = "MODEL {:>6d}\n".format(self.id)
        output = ""
        for node in self.nodes:
            output += node.__str__() + "\n"
        if len(self.connections) == 0:
            self.build_connections()
        for connection in self.connections:
            output += connection + "\n"
        # output += "ENDMDL\n"
        output += "\n"
        return output
    
    def add_node(self, node):
        self.nodes.append(node)
    
    def get_nodes(self):
        return self.nodes
    
    def get_node(self, pos):
        assert pos >= 0 and pos < len(self.nodes)
        return self.nodes[pos]
    
    def build_connections(self):
        self.connections = []
        # for i in range(len(self.nodes)-1):
        #     self.connections.append("CONECT{:>5d}{:>5d}".format(i, i+1))
        for node in self.nodes:
            self.connections.append("CONECT{:>5d}{:>5d}".format(node.id, node.id+1))
        self.connections.pop()
    
    def reverse_node_order(self):
        self.nodes.reverse()

class CaverNode(object):
    """docstring for CaverNode"""
    def __init__(self, id, x, y, z, radius, **kwargs):
        super(CaverNode, self).__init__()
        self.id = id
        self.x = x
        self.y = y
        self.z = z
        self.pos = (x, y, z)
        self.radius = radius
        if "color" in kwargs:
            self.color = kwargs["color"]
        else:
            self.color = [0., 0., 0.]
        if "res_name" in kwargs:
            self.res_name = kwargs["res_name"]
        else:
            self.res_name = "FIL"
        if "atom_name" in kwargs:
            self.atom_name = kwargs["atom_name"]
        else:
            self.atom_name = "H"
        if "res_id" in kwargs:
            self.res_id = kwargs["res_id"]
        else:
            self.res_id = 1
        if "chain" in kwargs:
            self.chain = kwargs["chain"]
        else:
            self.chain = "T"
    
    def __str__(self):
        return "ATOM  {:>5d}{:>3}{:>6}{:>2}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}     {:>7.2f}              ".format(
                self.id, self.atom_name, self.res_name, self.chain, self.res_id, self.x, self.y, self.z, self.radius)

class GridBox(object):
    """docstring for Grid"""
    def __init__(self, grid_id, size, x_ini, y_ini, z_ini):
        super(GridBox, self).__init__()
        self.size = size
        self.grid_id = grid_id
        self.x_ini = x_ini
        self.x_end = x_ini + size
        self.y_ini = y_ini
        self.y_end = y_ini + size
        self.z_ini = z_ini
        self.z_end = z_ini + size
        self.atoms = set()
    
    def add_atom(self, atom):
        self.atoms.add(atom)
    
    def __str__(self):
        out = "GridBox: {:>6d} Size: {:>8.3f}\n"
        out += "             X limits: {:>8.3f} {:>8.3f} Y limits: {:>8.3f} {:>8.3f} Z limits: {:>8.3f} {:>8.3f}\n"
        out = out.format(self.grid_id, self.size, self.x_ini, self.x_end, self.y_ini, self.y_end, self.z_ini, self.z_end)
        out += "             Atoms:\n"
        for atom in self.atoms:
            out += "                   {}\n".format(atom)
        return out

class GridList(object):
    """docstring for GridList"""
    def __init__(self, size, level, num_elements, ini_coords):
        super(GridList, self).__init__()
        self.size = size
        self.level = level
        self.elements = [None] * num_elements
        self.ini_coords = ini_coords
        self._initialize()
    
    def _initialize(self):
        def _get_GL(level, factor):
            new_pos = np.copy(self.ini_coords)
            new_pos[self.level] -= 0.1
            new_pos[self.level] += self.size*factor
            return GridList(self.size, level, self.num_elements, new_pos)
        def _get_GB(factor):
            new_pos = np.copy(self.ini_coords)
            new_pos[self.level] -= 0.1
            new_pos[self.level] += self.size*factor
            return GridBox(self.size, *new_pos)
        if self.level < 1:
            self.elements = [_get_GL(self.level+1, i) for i in range(len(self.elements))]
        elif self.level == 1:
            self.elements = [_get_GB(i) for i in range(len(self.elements))]
        else:
            raise NotImplementedError
    
    def add_atom(self, atom):
        pass

class ProteinGrid(object):
    """docstring for ProteinGrid"""
    def __init__(self, atoms, geom_center, min_x, max_x, min_y, max_y, min_z, max_z, grid_size=10.):
        super(ProteinGrid, self).__init__()
        self.atoms = atoms
        self.cog = geom_center
        self.grid_size = grid_size
        self.ini_coords = np.array([min_x, min_y, min_z])
        x_grids = int((max_x-min_x)/grid_size) + 1
        y_grids = int((max_y-min_y)/grid_size) + 1
        z_grids = int((max_z-min_z)/grid_size) + 1
        self.grid = [None] * x_grids
        for i in range(x_grids):
            self.grid[i] = [None] * y_grids
            for j in range(y_grids):
                self.grid[i][j] = [None] * z_grids
        grid_id = 0
        # for i in range(len(self.x_grids)):
        #     self.x_grids[i] = Grid(grid_id, self.size, (min_x-.1) + (i*self.size),
        #                            (min_y-.1) + (i*self.size), (min_z-.1) + (i*self.size))
        #     grid_id += 1
        pos = [None] * 3
        for i in range(x_grids):
            for j in range(y_grids):
                for k in range(z_grids):
                    pos[0] = min_x-0.1 + i*grid_size
                    pos[1] = min_y-0.1 + j*grid_size
                    pos[2] = min_z-0.1 + k*grid_size
                    self.grid[i][j][k] = GridBox(grid_id, self.grid_size, *pos)
                    grid_id += 1
        self.distribute_atoms()
    
    def distribute_atoms(self):
        for atom in self.atoms:
            pos = (atom.coords() - self.ini_coords) / self.grid_size
            pos = pos.astype(int)
            self.grid[pos[0]][pos[1]][pos[2]].add_atom(atom)
