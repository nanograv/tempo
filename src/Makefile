#      $Id$

#  Makefile for tempo


FC = f77
FFLAGS= -O3 

CC = cc
CCFLAGS = 

# suppress sccs retrieval
.SCCS_GET:


OBJS =	a1ut1f.o a1utcf.o ang.o arrtim.o atimfake.o bigendian.o blkdbat.o \
        bnrybt.o bnrydd.o bnryddgr.o bnryddt.o bnryeh.o bnryell1.o bnrymss.o \
        bootmc.o bt2ell1.o chebpc.o citem.o clockcor.o cutil.o damoyr.o dot.o \
        dskrd.o earth.o  ephinit.o ephread.o equ2ecl.o fit.o fk4tofk5.o \
        galco.o gasdev.o inpar.o interp.o lmst.o matinv.o mxprt.o newbin.o \
        newsrc.o newval.o nutation.o obsite.o outpar.o pcard.o pcshft.o \
        pospm_conv.o prcnut.o precession.o propmo.o radian.o ran1.o resid.o \
        setup.o sort.o space_motion.o tdbinit.o tdbread.o tempo.o  timcalc.o \
        tpohdr.o tzfit.o tzinit.o ut1red.o vmemw.o  ztim.o 

all:	tempo matrix

tempo:	$(OBJS)
	$(FC) $(FFLAGS) -o tempo $(OBJS)

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
	$(FC) -c -O1 tzfit.f
vmemw.o: vmemw.f dim.h vcom.h
	$(FC) -c $(FFLAGS) vmemw.f
ztim.o: ztim.f trnsfr.h
	$(FC) -c $(FFLAGS) ztim.f

matrix:	matrix.f dim.h
	$(FC) -o matrix matrix.f
clean:
	rm -f *.o *~

TDBGEN_OBJS = tdbgen.o tdb1ns.o

tdbgen: $(TDBGEN_OBJS)
	$(FC) $(FFLAGS) -o tdbgen $(TDBGEN_OBJS)

