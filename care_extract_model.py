#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  care_extract_model.py
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
    traj = list(args.infile)
    i = 0
    result = []
    while i < len(traj):
        if "MODEL" in traj[i]:
            if int(traj[i][5:]) == args.model:
                while "ENDMDL" not in traj[i]:
                    result.append(traj[i])
                    i += 1
                i = len(traj)
        i += 1
    args.outfile.write(''.join(result))
    print("Finished")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Your script description')
    parser.add_argument('-i', '--infile', required=True, type=argparse.FileType('r'), dest='infile', help='Input file to perform the operations.')
    parser.add_argument('-o', '--outfile', required=True, type=argparse.FileType('w'), dest='outfile', help='Output file name.')
    parser.add_argument('-m', '--model', type=int, dest='model', help='Number of the model to extract.')
    args = parser.parse_args()
    main()
