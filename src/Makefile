#      $Id$

#  Makefile for tempo and related programs

#  Always run make in the "src" subdirectory.  Binaries will be
#  created in the tempo root directory.

#  To make programs:
#     make tempo
#     make matrix
#     make tdbtrans

#  To make an archive file:
#     make archive



# compilers and compiler flags

FC = f77
FFLAGS= -O3 

CC = cc
CCFLAGS = 



# suppress sccs retrieval
.SCCS_GET:





# source code for tempo itself.  
# (note: variables SRC3 and SRC4 are used for making archives)

SRC1 = a1ut1f.f a1utcf.f ang.f arrtim.f atimfake.f bigendian.f blkdbat.f \
        bnrybt.f bnrydd.f bnryddgr.f bnryddt.f bnryeh.f bnryell1.f bnrymss.f \
        bootmc.f bt2ell1.f chebpc.f citem.f clockcor.f damoyr.f dot.f \
        dskrd.f earth.f  ephinit.f ephread.f equ2ecl.f fit.f fk4tofk5.f \
        galco.f gasdev.f inpar.f interp.f lmst.f matinv.f mxprt.f newbin.f \
        newsrc.f newval.f nutation.f obsite.f outpar.f pcard.f pcshft.f \
        pospm_conv.f prcnut.f precession.f propmo.f radian.f ran1.f resid.f \
        setup.f sort.f space_motion.f tasc2t0.f tdbinit.f tdbread.f tempo.f \
        timcalc.f tpohdr.f tzfit.f tzinit.f ut1red.f vmemw.f  ztim.f 
SRC2 = cutil.c
SRC3 = acom.h dim.h bcom.h clocks.h dp.h eph.h ephcom.h glitch.h \
	orbit.h tdbcom.h trnsfr.h tz.h vcom.h 
SRC4 = Makefile

OBJS = $(SRC1:.f=.o) $(SRC2:.c=.o)



# source code for "matrix"

MATRIX_SRC = matrix.f
MATRIX_OBJS = $(MATRIX_SRC:.f=.o)



# source code  for "tdbgen"

TDBGEN_SRC_F = tdbgen.f tdb1ns.f bigendian.f
TDBGEN_SRC_C = cutil.c
TDBGEN_SRC = $(TDBGEN_SRC_F) $(TDBGEN_SRC_C)	
TDBGEN_OBJS = $(TDBGEN_SRC_F:.f=.o) $(TDBGEN_SRC_C:.c=.o)



# file lists for creating an archive of tempo

# files in the 'src' directory
ARCHIVE_S_0 = $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(MATRIX_SRC) $(TDBGEN_SRC) 
ARCHIVE_S = $(ARCHIVE_S_0:%=src/%)

# files in the 'test' directory
ARCHIVE_T_0 = 0437.par 0437.tim tempo.lis.original 0437-4715.par.original
ARCHIVE_T = $(ARCHIVE_T_0:%=test/%)

# files in the 'clock' directory
ARCHIVE_C_0 = bipmnist.96 leap.sec ut1.dat utccorr.tot time.dat
ARCHIVE_C = $(ARCHIVE_C_0:%=clock/%)

# include the entire 'tzpar' directory (primarily so it gets created)
ARCHIVE_Z = tzpar

# files in the root directory (except for ephemerdies)
ARCHIVE_1 = obsys.dat tempo.cfg tempo.hlp 

tempo:	$(OBJS)
	$(FC) $(FFLAGS) -o ../tempo $(OBJS)

archive:
	cd ..;  tar -chvf tempo.tar $(ARCHIVE_S) $(ARCHIVE_T) $(ARCHIVE_1) \
	    $(ARCHIVE_C) $(ARCHIVE_Z);  gzip tempo.tar

all:	tempo matrix

matrix:	$(MATRIX_OBJS) dim.h
	$(FC) $(FFLAGS) -o ../matrix $(MATRIX_OBJS)

tdbgen: $(TDBGEN_OBJS)
	$(FC) $(FFLAGS) -o ../tdbgen $(TDBGEN_OBJS)

clean:
	rm -f *.o *~



arrtim.o: arrtim.f acom.h orbit.h dp.h bcom.h trnsfr.h clocks.h dim.h glitch.h tz.h
	$(FC) -c $(FFLAGS) arrtim.f
atimfake.o: atimfake.f dim.h tz.h
	$(FC) -c $(FFLAGS) atimfake.f
bnrybt.o: bnrybt.f dp.h orbit.h dim.h
	$(FC) -c $(FFLAGS) bnrybt.f
bnrydd.o: bnrydd.f dp.h orbit.h dim.h
	$(FC) -c $(FFLAGS) bnrydd.f
bnryddgr.o: bnryddgr.f orbit.h dim.h
	$(FC) -c $(FFLAGS) bnryddgr.f
bnryddt.o: bnryddt.f orbit.h dim.h
	$(FC) -c $(FFLAGS) bnryddt.f
bnryeh.o: bnryeh.f dp.h orbit.h dim.h
	$(FC) -c $(FFLAGS) bnryeh.f
bnryell1.o: bnryell1.f dp.h orbit.h dim.h
	$(FC) -c $(FFLAGS) bnryell1.f
bnrymss.o: bnrymss.f dp.h orbit.h dim.h
	$(FC) -c $(FFLAGS) bnrymss.f
bootmc.o: bootmc.f dim.h
	$(FC) -c $(FFLAGS) bootmc.f
bt2ell1.o: bt2ell1.f orbit.h
	$(FC) -c $(FFLAGS) bt2ell1.f
clockcor.o: clockcor.f dim.h clocks.h
	$(FC) -c $(FFLAGS) clockcor.f
cutil.o: cutil.c
	$(CC) -c $(CCFLAGS) cutil.c
fit.o: fit.f acom.h bcom.h orbit.h dim.h tz.h
	$(FC) -c $(FFLAGS) fit.f
ephinit.o: ephinit.f ephcom.h
	$(FC) -c $(FFLAGS) ephinit.f
ephread.o: ephread.f ephcom.h
	$(FC) -c $(FFLAGS) ephread.f
mxprt.o: mxprt.f dim.h
	$(FC) -c $(FFLAGS) mxprt.f
inpar.o: inpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h
	$(FC) -c $(FFLAGS) inpar.f
newbin.o: newbin.f acom.h orbit.h dim.h
	$(FC) -c $(FFLAGS) newbin.f
newsrc.o: newsrc.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h
	$(FC) -c $(FFLAGS) newsrc.f
newval.o: newval.f acom.h bcom.h dp.h trnsfr.h orbit.h dim.h glitch.h
	$(FC) -c $(FFLAGS) newval.f
outpar.o: outpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h eph.h
	$(FC) -c $(FFLAGS) outpar.f
pcard.o: pcard.f acom.h dim.h
	$(FC) -c $(FFLAGS) pcard.f
resid.o: resid.f acom.h bcom.h orbit.h dp.h dim.h glitch.h
	$(FC) -c $(FFLAGS) resid.f
setup.o: setup.f acom.h dim.h
	$(FC) -c $(FFLAGS) setup.f
tasc2t0.o: tasc2t0.f orbit.h
	$(FC) -c $(FFLAGS) tasc2t0.f
tdbinit.o: tdbinit.f tdbcom.h
	$(FC) -c $(FFLAGS) tdbinit.f	
tdbread.o: tdbread.f tdbcom.h
	$(FC) -c $(FFLAGS) tdbread.f	
tempo.o: tempo.f acom.h bcom.h dp.h orbit.h dim.h tz.h
	$(FC) -c $(FFLAGS) tempo.f
timcalc.o: timcalc.f acom.h dim.h trnsfr.h
	$(FC) -c $(FFLAGS) timcalc.f
tzinit.o: tzinit.f dim.h tz.h
	$(FC) -c $(FFLAGS) tzinit.f
tpohdr.o: tpohdr.f dim.h tz.h
	$(FC) -c $(FFLAGS) tpohdr.f
tzfit.o: tzfit.f dim.h acom.h tz.h
	$(FC) -c $(FFLAGS) tzfit.f
ut1red.o: ut1red.f dim.h 
	$(FC) -c $(FFLAGS) ut1red.f
vmemw.o: vmemw.f dim.h vcom.h
	$(FC) -c $(FFLAGS) vmemw.f
ztim.o: ztim.f trnsfr.h
	$(FC) -c $(FFLAGS) ztim.f

