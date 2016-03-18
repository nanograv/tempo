#!/usr/bin/env python

import sys
import argparse, struct

parser = argparse.ArgumentParser(description="Sort DMX ranges in MJD order")
parser.add_argument("infile",   help="input .par file")
parser.add_argument("outfile",  help="output .par file")
 
args = parser.parse_args()
infile   = args.infile
outfile  = args.outfile

# read in the DMX ranges

fin  = open(infile,"r")
fout = open(outfile,"w")

dmx = {}
for s in fin:
  ss = s.split()
  if ss[0].startswith("DMX") and len(ss[0])>3:
    sss = ss[0].split("_")
    idx = int(sss[1])
    key = sss[0][3:]
    if not idx in dmx.keys():
      dmx[idx] = {}
    dmx[idx][key] = ss[1:]
  else:
    fout.write(s)

for idx in dmx.keys():
  if not "R1" in dmx[idx].keys():
    print "Error: no DMXR1 key for DMX range ",idx
    sys.exit() 
  if float(dmx[idx]["R1"][0])<1. and float(dmx[idx]["R2"][0])<1.:
    del dmx[idx]
k = dmx.keys()
k = sorted(k,key=lambda x:dmx[x]["R1"])

for i,idx in zip(range(len(k)),k):
  for key in dmx[idx].keys():
    fout.write("DMX%s_%4.4d %s\n" % (key,i+1," ".join(dmx[idx][key])))

fin.close()
fout.close()
