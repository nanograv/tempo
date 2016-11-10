#!/usr/bin/env python

import sys
import argparse, struct

try:
  import tempo_utils
except:
  print "tmepo_utils is required to run this script"
  print "make sure it is in your python path"
  print "or install it from https://github.com/demorest/tempo_utils"
  sys.exit()

parser = argparse.ArgumentParser(description="List infomration about TOAs in each DMX range")
parser.add_argument("parfile",  help="name of .par file")
parser.add_argument("timfile",  help="name of .tim file")
parser.add_argument("-u","--unsorted", action="store_true",
                     help="print list DMX number order (don't sort by MJD)")
 
args = parser.parse_args()
parfile  = args.parfile
timfile  = args.timfile
unsorted = args.unsorted

# read in the DMX ranges

fpar = open(parfile,"r")
dmxr = []
for s in fpar:
  if s.startswith("DMXR"):
    ss = s.split()
    parname = ss[0]
    if s.startswith("DMXR"):
      irange = int(ss[0].split("_")[-1])-1  # DMX ranges start with 1, arrays start with 0
      parval  = float(ss[1].replace("D","E"))
      while irange>len(dmxr)-1:
        dmxr.append([0.,0.])
      if s.startswith("DMXR1"):
        dmxr[irange][0] = parval
      elif s.startswith("DMXR2"):
        dmxr[irange][1] = parval
fpar.close()


# read in the MJDs and INFO flags

toas = tempo_utils.read_toa_file(timfile)

# loop through the MJDs and INFO flags,
# develop a list of INFO strings,
# and accumulate the number of TOAs in
# each DMX range segretaged by INFO string.

infolist = []
orphans = []
dmxn = [[] for i in range(len(dmxr))]
norange = []

for t in toas:
  if not t.is_toa():
    continue
  t.parse_line()
  mjd = t.mjd
  info = t.flags['f']

  if not info in infolist:
    infolist.append(info)
    for i in range(len(dmxn)):
      dmxn[i].append(0)
    norange.append(0)
  iinfo = infolist.index(info) 

  idx = -1
  for i,r in zip(range(len(dmxr)),dmxr):
    if mjd>=r[0] and mjd<=r[1]:
      dmxn[i][iinfo] += 1
      break
  else:
    norange[iinfo] += 1
    orphans.append(mjd)  

print "# information about DMX ranges"
print "#"
print "# parfile:  ",parfile
print "# timfile:  ",timfile 
print "#"
print "# Shorthand used in table below:"
for i,info in zip(range(len(infolist)),infolist):
  print "# I%2.2d ==  %s" % (i,info)
print "#"

if len(orphans)>0:
  print "# TOAs not in DMX ranges:"
  print ("# %4s %8s  %8s  %6s  "+("   I%2.2d ")*len(infolist)) % \
      tuple( [ "" , "","","" ] + range(len(infolist)) )
  print ("  %3s  %8s  %8s  %6s  "+(" %5d ")*len(norange)) % tuple( ["","","",""] + norange)
  print "# A full list of these TOAs is below the table"

print "#"


print ("# %4s %8s  %8s  %6s  "+("   I%2.2d ")*len(infolist)) % \
    tuple( [ "DMXn" , "  start","   end"," days" ] + range(len(infolist)) )

k = range(len(dmxr))
if not unsorted:
  k = sorted(k, key=lambda x: dmxr[x][0])

# for i,n,r in zip(range(len(dmxn)),dmxn,dmxr):
#   print ("  %3d  %8.2f  %8.2f  %6.2f  "+(" %5d ")*len(n)) % tuple( [i+1,r[0],r[1],r[1]-r[0]] + n)

for i in k:
  print ("  %3d  %8.2f  %8.2f  %6.2f  "+(" %5d ")*len(dmxn[i])) % \
            tuple( [i+1,dmxr[i][0],dmxr[i][1],dmxr[i][1]-dmxr[i][0]] + dmxn[i] )


if len(orphans)>0:
  print "#"
  print "# Dates of TOAs not in any DMX range:"
  for mjd in orphans:
    print "#    %f" %(mjd)


