#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  care_fit_by_residue.py
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

import sys, os

# autocompletion
import readline
import rlcompleter
readline.parse_and_bind('tab: complete')

# pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
import pymol
from pymol import cmd
args = sys.argv[1:]
pymol.pymol_argv = ['pymol','-qc'] + args
pymol.finish_launching()

if '-h' in args:
    print("help")
else:
    if ('-r' in args) and ('-i' in args) and ('-o' in args) and ('-n' in args) and (len(args) == 9):
        prot_a = args[args.index('-r')+1]
        prot_b = args[args.index('-i')+1]
        outfile = args[args.index('-o')+1]
        pos = int(args[args.index('-n')+1])
        
        cmd.load(protein)
        cmd.load(loop)
        a_name = os.path.basename(prot_a).split('.')[0]
        b_name = os.path.basename(prot_b).split('.')[0]
        cmd.select("a_s1", "/%s///%d/N"  %(a_name, pos))
        cmd.select("a_s2", "/%s///%d/CA" %(a_name, pos))
        cmd.select("a_s3", "/%s///%d/C"  %(a_name, pos))
        cmd.select("a_s1", "/%s///%d/N"  %(b_name, pos))
        cmd.select("a_s2", "/%s///%d/CA" %(b_name, pos))
        cmd.select("a_s3", "/%s///%d/C"  %(b_name, pos))
        out_loop = cmd.select("out_prot", "/%s////" %(b_name))
        cmd.pair_fit("b_s1","a_s1", "b_s2","a_s2", "b_s3","a_s3")
        cmd.save(outfile, "out_prot")
    else:
        print("help")

