#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pymol_fit_looping.py
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

import argparse
from pymol import cmd

def main():
    """ Main function """
    protein = cmd.load(args.ref)
    loop = cmd.load(args.infile)
    
    prt_s_1 = 
    prot_s = cmd.select("prot_start", "/%s//%s/%d+%d/N+CA+C+O" %(args.ref.split(".")[0], args.chain, args.bounds[0]-2, args.bounds[0]-1))
    prot_e = cmd.select("prot_end",   "/%s//%s/%d+%d/N+CA+C+O" %(args.ref.split(".")[0], args.chain, args.bounds[1]+1, args.bounds[1]+2))
    loop_s = cmd.select("loop_start", "/%s///%d+%d/N+CA+C+O" %(args.ref.split(".")[0], 1, 2))
    loop_e = cmd.select("loop_end",   "/%s///%d+%d/N+CA+C+O" %(args.ref.split(".")[0], args.bounds[1]-args.bounds[0]+3, , args.bounds[1]-args.bounds[0]+4))
    cmd.pair_fit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Your script description')
    parser.add_argument('-r', '--ref', required=True, type=argparse.FileType('r'), dest='ref', help='Input PDB protein file to fit the loop pdb (-i option).')
    parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), dest='infile', help='Input PDB loop file to fit onto -r option.')
    parser.add_argument('-c', '--chain', default="A", dest='chain', help='Chain where the loop is missing in the protein (-r, --ref option).')
    parser.add_argument('-o', '--output', required=True, type=argparse.FileType('w'), dest='outfile', help='Defines the name of the modified molecule file.')
    parser.add_argument('-b', '--boundaries', required=True, type=int, nargs=2, dest='bounds', help='Here you define the residue number of the initial and last residue of the loop, i.e. the residues that are missing in the structure.')
    args = parser.parse_args()
    main()
