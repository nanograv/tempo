#      $Id$

#  Makefile for tempo and related programs

#  Always run make in the "src" subdirectory.  Binaries will be
#  created in the tempo root directory.

#  To make programs:
#     make tempo
#     make tdbtrans

#  To make an archive file (in Princeton):
#     make archive



# compilers and compiler flags

FC = f77
FFLAGS= -O3 

CC = cc
CCFLAGS = 


# source code for tempo itself.  
# (note: variables SRC3 and SRC4 are used for making archives)

SRC1 = a1ut1f.f a1utcf.f ang.f arrtim.f atimfake.f bigendian.f blkdbat.f \
        bnrybt.f bnrydd.f bnryddgr.f bnryddt.f bnryeh.f bnryell1.f bnrymss.f \
        bootmc.f bt2ell1.f chebpc.f citem.f clockcor.f covar.f damoyr.f dot.f \
        dskrd.f earth.f  ephinit.f ephread.f equ2ecl.f fit.f fk4tofk5.f \
        galco.f gasdev.f inpar.f interp.f lmst.f mass2dd.f matinv.f mxprt.f \
        newbin.f \
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

MATRIX_DIR = ../contrib/matrix/
MATRIX_SRC = $(MATRIX_DIR)matrix.f

# source code  for "tdbgen"

TDBGEN_SRC_F = tdbgen.f tdb1ns.f bigendian.f
TDBGEN_SRC_C = cutil.c
TDBGEN_SRC = $(TDBGEN_SRC_F) $(TDBGEN_SRC_C)	
TDBGEN_OBJS = $(TDBGEN_SRC_F:.f=.o) $(TDBGEN_SRC_C:.c=.o)




all:	tempo 

tempo:	$(OBJS)
	$(FC) $(FFLAGS) -o tempo $(OBJS)

tdbgen: $(TDBGEN_OBJS)
	$(FC) $(FFLAGS) -o tdbgen $(TDBGEN_OBJS)

clean:
	rm -f *.o *~

# dependencies of Tempo subroutines

arrtim.o: arrtim.f acom.h orbit.h dp.h bcom.h trnsfr.h clocks.h dim.h glitch.h tz.h
atimfake.o: atimfake.f dim.h tz.h
bnrybt.o: bnrybt.f dp.h orbit.h dim.h
bnrydd.o: bnrydd.f dp.h orbit.h dim.h
bnryddgr.o: bnryddgr.f orbit.h dim.h
bnryddt.o: bnryddt.f orbit.h dim.h
bnryeh.o: bnryeh.f dp.h orbit.h dim.h
bnryell1.o: bnryell1.f dp.h orbit.h dim.h
bnrymss.o: bnrymss.f dp.h orbit.h dim.h
bootmc.o: bootmc.f dim.h
bt2ell1.o: bt2ell1.f orbit.h
clockcor.o: clockcor.f dim.h clocks.h
cutil.o: cutil.c
fit.o: fit.f acom.h bcom.h orbit.h dim.h tz.h
ephinit.o: ephinit.f ephcom.h
ephread.o: ephread.f ephcom.h
mxprt.o: mxprt.f dim.h
inpar.o: inpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h
newbin.o: newbin.f acom.h orbit.h dim.h
newsrc.o: newsrc.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h
newval.o: newval.f acom.h bcom.h dp.h trnsfr.h orbit.h dim.h glitch.h
outpar.o: outpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h eph.h
pcard.o: pcard.f acom.h dim.h
resid.o: resid.f acom.h bcom.h orbit.h dp.h dim.h glitch.h
setup.o: setup.f acom.h dim.h
tasc2t0.o: tasc2t0.f orbit.h
tdbinit.o: tdbinit.f tdbcom.h
tdbread.o: tdbread.f tdbcom.h
tempo.o: tempo.f acom.h bcom.h dp.h orbit.h dim.h tz.h
timcalc.o: timcalc.f acom.h dim.h trnsfr.h
tzinit.o: tzinit.f dim.h tz.h
tpohdr.o: tpohdr.f dim.h tz.h
tzfit.o: tzfit.f dim.h acom.h tz.h
ut1red.o: ut1red.f dim.h 
vmemw.o: vmemw.f dim.h vcom.h
ztim.o: ztim.f trnsfr.h



# file lists for creating an archive of tempo

# files in the 'src' directory
ARCHIVE_S_0 = $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(MATRIX_SRC) $(TDBGEN_SRC) 
ARCHIVE_S = $(ARCHIVE_S_0:%=src/%)

# files in the 'test' directory
ARCHIVE_T_0 = 0437.par 0437.tim tempo.lis.original 0437-4715.par.original
ARCHIVE_T = $(ARCHIVE_T_0:%=test/%)

# files in the 'clock' directory
ARCHIVE_C_0 = bipmnist.96 leap.sec ut1.dat utccorr.tot 
# files in the 'clock_dist' directory
ARCHIVE_C_1 = time.dat
ARCHIVE_C = $(ARCHIVE_C_0:%=clock/%) $(ARCHIVE_C_1:%=clock/%)
ARCHIVE_C_L = $(ARCHIVE_C_0:%=../../clock/%) $(ARCHIVE_C_1:%=../../clock_dist/%)

# include one file from the 'tzpar' directory (primarily to create directory)
ARCHIVE_Z = tzpar/1937+21.par

# files in the root directory (except for ephemerdies)
ARCHIVE_1 = obsys.dat tempo.cfg tempo.hlp README
ARCHIVE_1_L = $(ARCHIVE_1:%=../%)

archive:
	cd ..; mkdir tmp; mkdir tmp/clock; \
	ln -s ../src ../test $(ARCHIVE_1_L) tmp; \
	ln -s ../tzpar_dist tmp/tzpar;   \
	ln -s $(ARCHIVE_C_L) tmp/clock;    \
	cd tmp; \
	chmod -R u+w *; \
	tar -chvf tempo.tar $(ARCHIVE_S) $(ARCHIVE_T) $(ARCHIVE_1) \
	    $(ARCHIVE_C) $(ARCHIVE_Z);  \
	gzip tempo.tar; \
	mv tempo.tar.gz ..; \
	cd ..; rm -r tmp

