#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  care_loop_paster.py
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

def main():
    """ Main function """
    start_pat = args.chain + "{:>4}".format(str(args.bounds[0]))
    end_pat = args.chain + "{:>4}".format(str(args.bounds[1]))
    limit = args.bounds[1] - args.bounds[0]
    p_list = list(args.protein)
    limit = args.bounds[1] - args.bounds[0] + 1
    loop = []
    i = 0
    while i < len(p_list):
        if start_pat in p_list[i]:
            while (start_pat in p_list[i]):
                loop.append(p_list[i])
                i += 1
            for j in range(1, limit):
                t_pat = args.chain + "{:>4}".format(args.bounds[0]+j)
                while (t_pat in p_list[i]):
                    loop.append(p_list[i])
                    i += 1
        i += 1
    args.outfile.write(''.join(loop))
    print("Finished")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Your script description')
    parser.add_argument('-p', '--protein', required=True, type=argparse.FileType('r'), dest='protein', help='Input protein file to perform the loop adition.')
    parser.add_argument('-o', '--outfile', required=True, type=argparse.FileType('w'), dest='outfile', help='Output file name.')
    parser.add_argument('-b', '--bounds', type=int, dest='bounds', nargs=2, help='The indexes where the loop begins and ends. Here you have to put the missing residue numbers, not the anchors!!!')
    parser.add_argument('-c', '--chain', dest='chain', default='A', help='The chain where the loop is located. Default is A.')
    args = parser.parse_args()
    main()
