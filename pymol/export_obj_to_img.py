#!/home/carlos/Programs/pymol/bin/python
# -*- coding: utf-8 -*-
#
#  export_obj_to_img.py
#  
#  Copyright 2018 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import pymol
from pymol import cmd

def export_objects_to_png(selection='all'):
    objects = cmd.get_object_list(selection)
    cmd.disable(selection)
    for pobj in objects:
        cmd.enable(pobj)
        cmd.png('{}.png'.format(pobj), width=2340, height=1239, ray=1)
        cmd.disable(pobj)

cmd.extend("export_objects_to_png", export_objects_to_png)