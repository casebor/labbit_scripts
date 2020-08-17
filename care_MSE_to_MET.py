#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  care_MSE_to_MET.py
#  
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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
    in_list = list(args.infile)
    modified = []
    for line in in_list:
        if len(line) > 20:
            if line[17:20] == "MSE":
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if line[12:14] == "SE":
                        outline = "ATOM  {} SD  MET{} S{}".format(line[6:12], line[20:76], line[78:])
                    else:
                        outline = "ATOM  {}MET{}".format(line[6:17], line[20:])
                    modified.append(outline)
                elif line.startswith("ANISOU"):
                    if line[12:14] == "SE":
                        outline = "ANISOU{} SD  MET{} S{}".format(line[6:12], line[20:76], line[78:])
                    else:
                        outline = "ANISOU{}MET{}".format(line[6:17], line[20:])
                    modified.append(outline)
            else:
                modified.append(line)
        else:
            modified.append(line)
    args.outfile.write(''.join(modified))
    print("Finished")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Changes the residue names and atoms from MSE to MET.')
    parser.add_argument('-i', '--input', required=True, type=argparse.FileType('r'), dest='infile', help='Input protein file to perform the change of modified MSE to MET.')
    parser.add_argument('-o', '--output', required=True, type=argparse.FileType('w'), dest='outfile', help='Output file name.')
    args = parser.parse_args()
    main()
