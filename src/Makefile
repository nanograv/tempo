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
FFLAGS= -O

CC = cc
CFLAGS = -O


# source code for tempo itself.  
# (note: variables SRC3 and SRC4 are used for making archives)

SRC1 = a1ut1f.f a1utcf.f ang.f arrtim.f atimfake.f bigendian.f \
        blkdbat.f bnrybt.f bnrybtx.f bnrydd.f bnryddgr.f bnryddt.f  \
        bnryeh.f bnryell1.f bnrymss.f bootmc.f bt2ell1.f \
        chebpc.f citem.f clockcor.f covar.f damoyr.f dot.f \
        earth.f  ephinit.f ephread.f equ2ecl.f fit.f fk4tofk5.f galco.f \
        gasdev.f inpar.f interp.f lmst.f mass2dd.f matinv.f mxprt.f newbin.f \
        newsrc.f newval.f nutation.f obsite.f outpar.f pcard.f pcshft.f \
        pospm_conv.f prcnut.f precession.f propmo.f radian.f ran1.f resid.f \
        setup.f sitea2n.f siten2a.f \
	sort.f space_motion.f ell12bt.f tdbinit.f tdbinit2.f tdbread.f tempo.f \
        timcalc.f tmalloc.f tparin.f tpohdr.f tzfit.f tzinit.f ut1red.f \
        vmemw.f  ztim.f 
SRC2 = cutil.c
SRC3 = acom.h array.h dim.h bcom.h clocks.h dp.h eph.h ephcom.h glitch.h \
	orbit.h tdbcom.h trnsfr.h tz.h vcom.h toa.h 
SRC4 = Makefile

OBJS = $(SRC1:.f=.o) $(SRC2:.c=.o)


# source code for "matrix"

MATRIX_SRC = matrix.f
MATRIX_OBJS = $(MATRIX_SRC:.f=.o)

# source code  for "tdbgen"  
#   SRC_F_1 are unique to tdgben;  SRC_F_2 is shared with Tempo

TDBGEN_SRC_F_1 = tdbgen.f tdb1ns.f
TDBGEN_SRC_F_2 = bigendian.f
TDBGEN_SRC_C = cutil.c
TDBGEN_SRC_1 = $(TDBGEN_SRC_F_1) 
TDBGEN_OBJS = $(TDBGEN_SRC_F_1:.f=.o) $(TDBGEN_SRC_F_2:.f=.o) $(TDBGEN_SRC_C:.c=.o)




all:	tempo 

tempo:	$(OBJS)
	$(FC) $(FFLAGS) -o tempo $(OBJS)

tdbgen: $(TDBGEN_OBJS)
	$(FC) $(FFLAGS) -o tdbgen $(TDBGEN_OBJS)

matrix: $(MATRIX_OBJS)
	$(FC) $(FFLAGS) -o matrix $(MATRIX_OBJS)

clean:
	rm -f *.o *~

# dependencies of Tempo subroutines

arrtim.o: arrtim.f acom.h orbit.h dp.h bcom.h trnsfr.h clocks.h dim.h glitch.h tz.h
atimfake.o: atimfake.f dim.h tz.h toa.h
bnrybt.o: bnrybt.f dp.h orbit.h dim.h
bnrybtx.o: bnrybt.f dp.h orbit.h dim.h
bnrydd.o: bnrydd.f dp.h orbit.h dim.h
bnryddgr.o: bnryddgr.f orbit.h dim.h
bnryddt.o: bnryddt.f orbit.h dim.h
bnryeh.o: bnryeh.f dp.h orbit.h dim.h
bnryell1.o: bnryell1.f dp.h orbit.h dim.h
bnrymss.o: bnrymss.f dp.h orbit.h dim.h
bootmc.o: bootmc.f dim.h
bt2ell1.o: bt2ell1.f orbit.h dim.h
clockcor.o: clockcor.f dim.h clocks.h
covar.o: acom.h dim.h
cutil.o: cutil.c
ell12bt.o: ell12bt.f dim.h orbit.h
ephinit.o: ephinit.f ephcom.h
ephread.o: ephread.f ephcom.h
fit.o: fit.f acom.h bcom.h orbit.h dim.h tz.h
inpar.o: inpar.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h
matinv.o: dim.h
mxprt.o: mxprt.f dim.h
newbin.o: newbin.f acom.h orbit.h dim.h
newsrc.o: newsrc.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h
newval.o: newval.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h
outpar.o: outpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h eph.h
pcard.o: pcard.f acom.h dim.h dp.h
resid.o: resid.f acom.h bcom.h orbit.h dp.h dim.h glitch.h
setup.o: setup.f acom.h dim.h
tdbinit.o: tdbinit.f tdbcom.h
tdbinit2.o: tdbinit2.f tdbcom.h
tdbread.o: tdbread.f tdbcom.h
tempo.o: tempo.f acom.h bcom.h dp.h array.h clocks.h eph.h vcom.h orbit.h dim.h trnsfr.h tz.h
timcalc.o: timcalc.f acom.h dim.h dp.h trnsfr.h
tmalloc.o: array.h
tparin.o: tparin.f acom.h dim.h vcom.h
tpohdr.o: tpohdr.f dim.h tz.h
tzfit.o: tzfit.f dim.h acom.h trnsfr.h tz.h
tzinit.o: tzinit.f dim.h tz.h
ut1red.o: ut1red.f dim.h 
vmemw.o: vmemw.f dim.h vcom.h
ztim.o: ztim.f trnsfr.h

# dependency of matrix:
matrix.o: matrix.f dim.h

# file lists for creating an archive of tempo

# files in the 'src' directory
ARCHIVE_S_0 = $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(MATRIX_SRC) $(TDBGEN_SRC_1) 
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

archive-src:
	cd ..; tar -cvf temposrc.tar $(ARCHIVE_S); gzip temposrc.tar
