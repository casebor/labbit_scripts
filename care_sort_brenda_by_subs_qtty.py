#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  care_sort_brenda_by_subs_qtty.py
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

import sys, warnings

class BrendaEntry:
    
    def __init__ (self):
        self.ec = None
        self.name = None
        self.crystal_qtty = 0
        self.substrates_qtty = 0
        self.nat_substrates_qtty = 0
        self.crystals = []
        self.references = []
        self.proteins = {}
        self.substrates = {}
        self.nat_substrates = {}
        self.inhibitors = {}
        self.pdbs = set([])
    
    def _endline(self, line):
        return line.startswith('\n') or line.startswith('\\\\\\')
    
    def parse_id(self, line):
        elements = line.split()
        self.ec = elements[1]
    
    def add_pdb(self, code):
        self.pdbs.add(code)
    
    def parse_name(self, instream):
        line = instream.readline()
        while not self._endline(line):
            if line.startswith('RN\t'):
                self.name = line.split('\t')[1].strip()
            line = instream.readline()
    
    def parse_proteins(self, instream):
        line = instream.readline()
        buffer = ''
        while not self._endline(line):
            if line.startswith('\t'):
                buffer += ' '+line.strip()
            elif line.startswith('PR\t'):
                if buffer != '':
                    prots = self._get_proteins(buffer)
                    protein_name = self._get_protein_name(buffer)
                    self.proteins[prots[0]] = [protein_name, 0]
                buffer = line.rstrip()
            line = instream.readline()
        prots = self._get_proteins(buffer)
        protein_name = self._get_protein_name(buffer)
        self.proteins[prots[0]] = [protein_name, 0]
    
    def parse_subtrates(self, instream):
        line = instream.readline()
        buffer = ''
        while not self._endline(line):
            if line.startswith('\t'):
                buffer += ' '+line.strip().rstrip()
            elif line.startswith('SP\t'):
                self.substrates_qtty += 1
                if buffer != '':
                    proteins = self._get_proteins(buffer)
                    substrates = self._get_substrates(buffer)
                    self.substrates[substrates] = proteins
                buffer = line.strip()
            line = instream.readline()
        if buffer != '':
            proteins = self._get_proteins(buffer)
            substrates = self._get_substrates(buffer)
            self.substrates[substrates] = proteins
    
    def parse_nat_substrates(self, instream):
        line = instream.readline()
        buffer = ''
        while not self._endline(line):
            if line.startswith('\t'):
                buffer += ' '+line.strip().rstrip()
            elif line.startswith('NSP\t'):
                self.nat_substrates_qtty += 1
                if buffer != '':
                    proteins = self._get_proteins(buffer)
                    substrates = self._get_substrates(buffer)
                    self.nat_substrates[substrates] = proteins
                buffer = line.strip()
            line = instream.readline()
        if buffer != '':
            proteins = self._get_proteins(buffer)
            substrates = self._get_substrates(buffer)
            self.nat_substrates[substrates] = proteins
    
    def parse_crystals(self, instream):
        line = instream.readline()
        buffer = ''
        while not self._endline(line):
            if line.startswith('\t'):
                buffer += ' '+line.strip().rstrip()
            elif line.startswith('CR\t'):
                self.crystal_qtty += 1
                if buffer != '':
                    proteins = self._get_proteins(buffer)
                    self.crystals.append(proteins[0])
                buffer = line.strip()
            line = instream.readline()
        if buffer != '':
            proteins = self._get_proteins(buffer)
            self.crystals.append(proteins[0])
        assert(len(self.crystals)==self.crystal_qtty)
    
    def parse_references(self, instream):
        line = instream.readline()
        buffer = ''
        while not self._endline(line):
            if line.startswith('\t'):
                buffer += ' ' + line.strip().rstrip()
            elif line.startswith('RF\t'):
                if buffer != '':
                    self.references.append(buffer)
                buffer = line.strip()
            line = instream.readline()
        if buffer != '':
            self.references.append(buffer)
    
    def parse_inhibitors(self, instream):
        line = instream.readline()
        buffer = ''
        while not self._endline(line):
            if line.startswith('\t'):
                buffer += ' '+line.strip().rstrip()
            elif line.startswith('IN\t'):
                if buffer != '':
                    proteins = self._get_proteins(buffer)
                    inhibitor = self._get_protein_name(buffer)
                    self.inhibitors[inhibitor] = proteins
                buffer = line.strip()
            line = instream.readline()
        if buffer != '':
            proteins = self._get_proteins(buffer)
            inhibitor = self._get_protein_name(buffer)
            self.inhibitors[inhibitor] = proteins
    
    def _get_proteins(self, line):
        proteins = []
        if line != '':
            stop_counter = 2
            parse_number = ''
            for char in line:
                if char in '<>()':
                    return proteins
                elif char == '#':
                    stop_counter -= 1
                else:
                    if stop_counter == 0:
                        proteins.append(int(parse_number))
                        break
                    elif stop_counter == 1:
                        if char == ',' or char == ' ':
                            proteins.append(int(parse_number))
                            parse_number = ''
                        else:
                            parse_number += char
        return tuple(proteins)
    
    def _get_protein_name(self, line):
        if line != '':
            stop_counter = 2
            parse_name = ''
            for char in line:
                if char in '<>()':
                    return parse_name.strip()
                elif char == '#':
                    stop_counter -= 1
                else:
                    if stop_counter == 0:
                        parse_name += char
    
    def _get_substrates(self, line):
        start = line.index('#') + 1
        start = line[start:].index('#') + start + 1
        end = line.index('=')
        return tuple(line[start:end].strip().split(' + '))
    
    def update_substrates_per_protein(self):
        output = ''
        if self.substrates_qtty:
            for group in self.substrates.values():
                for protein in group:
                    try:
                        self.proteins[protein][1] += 1
                    except KeyError:
                        output += "EC enzyme: %s has protein #%d not specified.\n"\
                              %(self.ec, protein)
        else:
            output += "EC enzyme: %s has no substrate information.\n" %(self.ec)
        if self.nat_substrates_qtty:
            for group in self.nat_substrates.values():
                for protein in group:
                    try:
                        self.proteins[protein][1] += 1
                    except KeyError:
                        output += "EC enzyme: %s has protein #%d not specified.\n"\
                              %(self.ec, protein)
        else:
            output += "EC enzyme: %s has no natural substrate information.\n" %(self.ec)
        return output
    
    def get_refs(self):
        return self.references
    
    def __str__(self):
        out_str = """Recommended Name: %s
EC Number: %s
Substrates in Brenda: %d
Natural substrates in Brenda: %d
""" %(self.name, self.ec, self.substrates_qtty, self.nat_substrates_qtty)
        out_str += "-"*80 + "\n"
        return out_str
    
    def get_info(self, extended=False):
        out_str = """Recommended Name: %s
EC Number: %s
Substrates in Brenda: %d
Natural substrates in Brenda: %d
Quantity of PDB structures: %d
""" %(self.name, self.ec, self.substrates_qtty, self.nat_substrates_qtty, len(self.pdbs))
        if extended:
            out_str += "List of Proteins and substrates used by them:\n"
            for protein in self.proteins.values():
                out_str += "\t%s: %d\n" %(protein[0], protein[1])
        out_str += "-"*80 + "\n"
        return out_str
    
    def print_proteins(self):
        for id,name in self.proteins.items():
            print("Protein %d: %s" %(id, name[0]))
    
    def __lt__(self, other):
        return self.substrates_qtty < other.substrates_qtty
    def __le__(self, other):
        return self.substrates_qtty <= other.substrates_qtty
    def __eq__(self, other):
        return self.substrates_qtty == other.substrates_qtty
    def __ne__(self, other):
        return self.substrates_qtty != other.substrates_qtty
    def __gt__(self, other):
        return self.substrates_qtty > other.substrates_qtty
    def __ge__(self, other):
        return self.substrates_qtty >= other.substrates_qtty
    

def add_pdbs_to_entries(instream, db):
    for line in instream:
        sections = line.strip().split('\t')
        if len(sections[0].split('.')) == 4:
            try:
                db[sections[0].strip()].add_pdb(sections[1].strip())
            except KeyError:
                pass

def main(args):
    brendafile = open(sys.argv[1], 'r')
    entries = []
    brenda_db = {}
    line = brendafile.readline()
    buffer_line = ""
    while line:
        if line.startswith('ID'):
            entry = BrendaEntry()
            entry.parse_id(line)
            entries.append(entry)
        elif line.startswith('CRYSTALLIZATION'):
            entry.parse_crystals(brendafile)
        elif line.startswith('INHIBITORS'):
            entry.parse_inhibitors(brendafile)
        elif line.startswith('PROTEIN'):
            entry.parse_proteins(brendafile)
        elif line.startswith('RECOMMENDED_NAME'):
            entry.parse_name(brendafile)
        elif line.startswith('SUBSTRATE_PRODUCT'):
            entry.parse_subtrates(brendafile)
        elif line.startswith('NATURAL_SUBSTRATE_PRODUCT'):
            entry.parse_nat_substrates(brendafile)
        elif line.startswith('REFERENCE'):
            entry.parse_references(brendafile)
        line = brendafile.readline()
    brendafile.close()
    entries = sorted(entries, reverse=True)
    for entry in entries:
        entry.update_substrates_per_protein()
        brenda_db[entry.ec] = entry
    for file in sys.argv[3:]:
        with open(file, 'r') as instream:
            add_pdbs_to_entries(instream, brenda_db)
    with open(sys.argv[2], 'w') as outstream:
        for entry in entries:
            outstream.write(entry.get_info(True))
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
