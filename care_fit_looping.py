#!/home/carlos/Programs/pymol/bin/python
# -*- coding: utf-8 -*-
#
#  care_fit_looping.py
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
    if ('-r' in args) and ('-i' in args) and ('-o' in args) and ('-c' in args) and ('-b' in args) and (len(args) == 11):
        protein = args[args.index('-r')+1]
        loop = args[args.index('-i')+1]
        outfile = args[args.index('-o')+1]
        chain = args[args.index('-c')+1]
        ini = int(args[args.index('-b')+1])
        end = int(args[args.index('-b')+2])
        
        cmd.load(protein)
        cmd.load(loop)
        p_name = os.path.basename(protein).split('.')[0]
        l_name = os.path.basename(loop).split('.')[0]
        #l_size = len(cmd.get_model(l_name).get_residues())
        cmd.select("p_s1", "/%s//%s/%d/N"  %(p_name, chain, ini-2))
        cmd.select("p_s2", "/%s//%s/%d/CA" %(p_name, chain, ini-2))
        cmd.select("p_s3", "/%s//%s/%d/C"  %(p_name, chain, ini-2))
        cmd.select("p_s4", "/%s//%s/%d/N"  %(p_name, chain, ini-1))
        cmd.select("p_s5", "/%s//%s/%d/CA" %(p_name, chain, ini-1))
        cmd.select("p_s6", "/%s//%s/%d/C"  %(p_name, chain, ini-1))
        cmd.select("p_e1", "/%s//%s/%d/N"  %(p_name, chain, end+1))
        cmd.select("p_e2", "/%s//%s/%d/CA" %(p_name, chain, end+1))
        cmd.select("p_e3", "/%s//%s/%d/C"  %(p_name, chain, end+1))
        cmd.select("p_e4", "/%s//%s/%d/N"  %(p_name, chain, end+2))
        cmd.select("p_e5", "/%s//%s/%d/CA" %(p_name, chain, end+2))
        cmd.select("p_e6", "/%s//%s/%d/C"  %(p_name, chain, end+2))
        cmd.select("l_s1", "/%s///%d/N"  %(l_name, ini-2))
        cmd.select("l_s2", "/%s///%d/CA" %(l_name, ini-2))
        cmd.select("l_s3", "/%s///%d/C"  %(l_name, ini-2))
        cmd.select("l_s4", "/%s///%d/N"  %(l_name, ini-1))
        cmd.select("l_s5", "/%s///%d/CA" %(l_name, ini-1))
        cmd.select("l_s6", "/%s///%d/C"  %(l_name, ini-1))
        cmd.select("l_e1", "/%s///%d/N"  %(l_name, end+1))
        cmd.select("l_e2", "/%s///%d/CA" %(l_name, end+1))
        cmd.select("l_e3", "/%s///%d/C"  %(l_name, end+1))
        cmd.select("l_e4", "/%s///%d/N"  %(l_name, end+2))
        cmd.select("l_e5", "/%s///%d/CA" %(l_name, end+2))
        cmd.select("l_e6", "/%s///%d/C"  %(l_name, end+2))
        out_loop = cmd.select("out_loop", "/%s////" %(l_name))
        cmd.pair_fit("l_s1","p_s1", "l_s2","p_s2", "l_s3","p_s3", "l_s4","p_s4", "l_s5","p_s5", "l_s6","p_s6",
                     "l_e1","p_e1", "l_e2","p_e2", "l_e3","p_e3", "l_e4","p_e4", "l_e5","p_e5", "l_e6","p_e6")
        cmd.save(outfile, "out_loop")
    else:
        print("help")

