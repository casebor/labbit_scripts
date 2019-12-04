#!/home/carlos/Programs/pymol/bin/python
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
    if ('-i' in args) and ('-o' in args) and (len(args) == 4):
        prot_a = args[args.index('-i')+1]
        outfile = args[args.index('-o')+1]
        
        cmd.load(prot_a)
        name_a = os.path.basename(prot_a).split('.')[0]
        cmd.select("residue_15", "/%s///%d/" %(name_a, 15))
        cmd.select("out_prot", "residue_15 expand %f" %(12.3))
        cmd.select("out_prot", "byresidue out_prot")
        cmd.save(outfile, "out_prot")
    else:
        print(len(args))
        print("Error in commands")

