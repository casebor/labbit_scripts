#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2020 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  

import argparse
import numpy as np
from Bio import SeqIO
from biostructmap import biostructmap as bmap


def main():
    """ Main function """
    ref_seq = SeqIO.read(args.pdb_file, "pdb-atom")
    original_pdb = bmap.Structure(pdbfile=args.pdb_file, pdbname="original")
    rmsf_info = np.loadtxt(args.rmsf_file, usecols=(1,))
    rmsf_info = np.insert(rmsf_info, 0, 0)
    rmsf_pdb = original_pdb.map(data={('A',): rmsf_info}, radius=0,
        method_params={"method": max})
    rmsf_pdb.write_data_to_pdb_b_factor(fileobj=args.outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the reuslts of RMSF to a pdb file. "
        "Make sure that the sequence is the same length as the PDB file.")
    parser.add_argument('-r', '--rmsf', required=True, dest='rmsf_file',
        help='The Amber RMSF file to be used.')
    parser.add_argument('-p', '--pdb', required=True, dest='pdb_file',
        help='The PDB file to be used.')
    parser.add_argument('-o', '--outfile', required=True, dest='outfile',
        help='Name for the PDB file with the RMSF data in the B-factor column.')
    args = parser.parse_args()
    main()
