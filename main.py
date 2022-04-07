'''
Python code to calculate SASA from given PDB and then calculate surface properties
using HPS-Urry or user provided parameters 
Author: Roshan M Regy
Email ID: roshanm.regy@tamu.edu 
''' 
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import argparse 

# from Get Area http://curie.utmb.edu/area_man.html
random_coil_values = {
        "ALA":64.9,
        "ARG":195.5,
        "ASN":114.3,
        "ASP":113.0,
        "CYS":102.3,
        "GLN":143.7,
        "GLU":141.2,
        "HIS":154.6,
        "ILE":147.3,
        "GLY":87.2,
        "LEU":146.2,
        "LYS":164.5,
        "MET":158.3,
        "PHE":180.1,
        "PRO":105.2,
        "SER":77.4,
        "THR":106.2,
        "TRP":224.6,
        "TYR":193.1,
        "VAL":122.3
        }


p = PDBParser(QUIET=1)
struct = p.get_structure('SAM','only_sam.pdb')
sr = ShrakeRupley()
sr.compute(struct, level="R")
for m,model in enumerate(sr):
    for c,chain in enumerate(model):
        for r,residue in enumerate(chain):
            name[r] = residue.get_name()
            sasa[r] = residue.sasa()
            normsasa[r] = sasa/random_coil_values[name]


