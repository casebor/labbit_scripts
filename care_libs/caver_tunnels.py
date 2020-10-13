#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

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
