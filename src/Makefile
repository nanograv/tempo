#  Makefile for tempo
#      $Id$
FFLAGS= -O3 
# FFLAGS= -C
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
        tpohdr.o tzfit.o tzinit.o ut1red.o vmemw.o  ztim.o \

all:	tpo matrix

tpo:	$(OBJS)
	f77 $(FFLAGS) -o ../tpo $(OBJS)

arrtim.o: arrtim.f acom.h orbit.h dp.h bcom.h trnsfr.h clocks.h dim.h glitch.h tz.h
	f77 -c $(FFLAGS) arrtim.f
atimfake.o: atimfake.f dim.h tz.h
	f77 -c $(FFLAGS) atimfake.f
bnrybt.o: bnrybt.f dp.h orbit.h dim.h
	f77 -c $(FFLAGS) bnrybt.f
bnrydd.o: bnrydd.f dp.h orbit.h dim.h
	f77 -c $(FFLAGS) bnrydd.f
bnryddgr.o: bnryddgr.f orbit.h dim.h
	f77 -c $(FFLAGS) bnryddgr.f
bnryddt.o: bnryddt.f orbit.h dim.h
	f77 -c $(FFLAGS) bnryddt.f
bnryeh.o: bnryeh.f dp.h orbit.h dim.h
	f77 -c $(FFLAGS) bnryeh.f
bnryell1.o: bnryell1.f dp.h orbit.h dim.h
	f77 -c $(FFLAGS) bnryell1.f
bnrymss.o: bnrymss.f dp.h orbit.h dim.h
	f77 -c $(FFLAGS) bnrymss.f
bootmc.o: bootmc.f dim.h
	f77 -c $(FFLAGS) bootmc.f
bt2ell1.o: bt2ell1.f orbit.h
	f77 -c $(FFLAGS) bt2ell1.f
clockcor.o: clockcor.f dim.h clocks.h
	f77 -c $(FFLAGS) clockcor.f
cutil.o: cutil.c
	cc -c $(CCFLAGS) cutil.c
fit.o: fit.f acom.h bcom.h orbit.h dim.h tz.h
	f77 -c $(FFLAGS) fit.f
ephinit.o: ephinit.f ephcom.h
	f77 -c $(FFLAGS) ephinit.f
ephread.o: ephread.f ephcom.h
	f77 -c $(FFLAGS) ephread.f
mxprt.o: mxprt.f dim.h
	f77 -c $(FFLAGS) mxprt.f
inpar.o: inpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h
	f77 -c $(FFLAGS) inpar.f
newbin.o: newbin.f acom.h orbit.h dim.h
	f77 -c $(FFLAGS) newbin.f
newsrc.o: newsrc.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h
	f77 -c $(FFLAGS) newsrc.f
newval.o: newval.f acom.h bcom.h dp.h trnsfr.h orbit.h dim.h glitch.h
	f77 -c $(FFLAGS) newval.f
outpar.o: outpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h eph.h
	f77 -c $(FFLAGS) outpar.f
pcard.o: pcard.f acom.h dim.h
	f77 -c $(FFLAGS) pcard.f
resid.o: resid.f acom.h bcom.h orbit.h dp.h dim.h glitch.h
	f77 -c $(FFLAGS) resid.f
setup.o: setup.f acom.h dim.h
	f77 -c $(FFLAGS) setup.f
tdbinit.o: tdbinit.f tdbcom.h
	f77 -c $(FFLAGS) tdbinit.f	
tdbread.o: tdbread.f tdbcom.h
	f77 -c $(FFLAGS) tdbread.f	
tempo.o: tempo.f acom.h bcom.h dp.h orbit.h dim.h tz.h
	f77 -c $(FFLAGS) tempo.f
timcalc.o: timcalc.f acom.h dim.h trnsfr.h
	f77 -c $(FFLAGS) timcalc.f
tzinit.o: tzinit.f dim.h tz.h
	f77 -c $(FFLAGS) tzinit.f
tpohdr.o: tpohdr.f dim.h tz.h
	f77 -c $(FFLAGS) tpohdr.f
tzfit.o: tzfit.f dim.h acom.h tz.h
	f77 -c -O1 tzfit.f
vmemw.o: vmemw.f dim.h vcom.h
	f77 -c $(FFLAGS) vmemw.f
ztim.o: ztim.f trnsfr.h
	f77 -c $(FFLAGS) ztim.f

matrix:	matrix.f dim.h
	f77 -o matrix matrix.f
clean:
	rm -f *.o *~

TDBGEN_OBJS = tdbgen.o tdb1ns.o

tdbgen: $(TDBGEN_OBJS)
	f77 $(FFLAGS) -o tdbgen $(TDBGEN_OBJS)

