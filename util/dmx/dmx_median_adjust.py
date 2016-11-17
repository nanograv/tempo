#!/usr/bin/env python

import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Adjust DM and DMX_xxxx so median DMX value is zero")
parser.add_argument("infile",   help="input .par file")
parser.add_argument("outfile",  help="output .par file")

args = parser.parse_args()
infile   = args.infile
outfile  = args.outfile

with open(infile,"r") as fin:
  p=fin.readlines()

dm = False
dmxlist = []

for s in p:
  if s.startswith("DMX_"):
    dmxlist.append(float(s.split()[1].replace("D","E")))

dmxmed = np.round(np.median(dmxlist),decimals=6)

with open(outfile,"w") as fout:
  for s in p:
    if s.startswith("DM "):
      dm = float(s.split()[1].replace("D","E"))+dmxmed
      dms = ("%23.6f" % dm).replace("e","D")
      s = "DM "+dms+s[26:]
    if s.startswith("DMX_"):
      dmx = float(s.split()[1].replace("D","E"))-dmxmed
      dmxs = ("%15.8e" % dmx).replace("e","D")
      s = s[0:11]+dmxs+s[26:]
    fout.write(s)


