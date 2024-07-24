#!/usr/bin/env python3

from pymol import cmd
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("native")
parser.add_argument("peptide")
parser.add_argument("ids", nargs="+")
args = parser.parse_args()

native = args.native
peptide = args.peptide
ids = args.ids

backbones = {i: glob(f"backbones/*_{i}.pdb") for i in ids}
structures = {i: glob(f"structures/*/*_{i}_sample*_model_0.cif") for i in ids}

for k,v in backbones.items():
   assert len(v) == 1, k
for k,v in structures.items():
   assert len(v) == 1, k

cmd.load(native, "native")
cmd.findseq(peptide, "native", "pep_native")

for i in ids:
   # load objects
   cmd.load(backbones[i], "backbone")
   cmd.load(structures[i], "structure")
   # get sequences
   cmd.findseq(peptide, "backbone", "pep_backbone")
   cmd.findseq(peptide, "structure", "pep_structure")

   # calculate values
   print(i)

   try:
      rmsd_backbone_ca = cmd.fit("bycalpha structure", "bycalpha backbone", matchmaker=-1)
   except Exception as e:
      print(f"{i}: rmsd_backbone_ca: {e}")
   try:
      rmsd_peptide_ca = cmd.fit("bycalpha pep_structure", "bycalpha pep_native", matchmaker=-1)
   except Exception as e:
      print(f"{i}: rmsd_peptide_ca : {e}")
   try:
      rmsd_peptide = cmd.fit("pep_structure", "pep_native", matchmaker=-1)
   except Exception as e:
      print(f"{i}: rmsd_peptide    : {e}")


   print("RMSD backbone (CA)", rmsd_backbone_ca)
   print("RMSD peptide (full atom)", rmsd_peptide)
   print("RMSD peptide (CA)", rmsd_peptide_ca)

   # unload objects
   cmd.delete("structure")
   cmd.delete("backbone")
