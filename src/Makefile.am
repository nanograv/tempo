bin_PROGRAMS = tempo

EXTRA_PROGRAMS = matrix tdbgen
CLEANFILES = $(EXTRA_PROGRAMS)

tempo_SOURCES = \
        a1ut1f.f a1utcf.f ang.f arrtim.f atimfake.f bigendian.f \
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
        vmemw.f  ztim.f \
        cutil.c \
        acom.h array.h dim.h bcom.h clocks.h dp.h eph.h ephcom.h glitch.h \
	orbit.h tdbcom.h trnsfr.h tz.h vcom.h toa.h 

matrix_SOURCES = matrix.f

tdbgen_SOURCES = tdbgen.f tdb1ns.f bigendian.f cutil.c

# dependencies of Tempo subroutines

arrtim.$(OBJEXT): arrtim.f acom.h orbit.h dp.h bcom.h trnsfr.h clocks.h dim.h glitch.h tz.h
atimfake.$(OBJEXT): atimfake.f dim.h tz.h toa.h
bnrybt.$(OBJEXT): bnrybt.f dp.h orbit.h dim.h
bnrybtx.$(OBJEXT): bnrybt.f dp.h orbit.h dim.h
bnrydd.$(OBJEXT): bnrydd.f dp.h orbit.h dim.h
bnryddgr.$(OBJEXT): bnryddgr.f orbit.h dim.h
bnryddt.$(OBJEXT): bnryddt.f orbit.h dim.h
bnryeh.$(OBJEXT): bnryeh.f dp.h orbit.h dim.h
bnryell1.$(OBJEXT): bnryell1.f dp.h orbit.h dim.h
bnrymss.$(OBJEXT): bnrymss.f dp.h orbit.h dim.h
bootmc.$(OBJEXT): bootmc.f dim.h
bt2ell1.$(OBJEXT): bt2ell1.f orbit.h dim.h
clockcor.$(OBJEXT): clockcor.f dim.h clocks.h
covar.$(OBJEXT): acom.h dim.h
cutil.$(OBJEXT): cutil.c
ell12bt.$(OBJEXT): ell12bt.f dim.h orbit.h
ephinit.$(OBJEXT): ephinit.f ephcom.h
ephread.$(OBJEXT): ephread.f ephcom.h
fit.$(OBJEXT): fit.f acom.h bcom.h orbit.h dim.h tz.h
inpar.$(OBJEXT): inpar.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h
matinv.$(OBJEXT): dim.h
mxprt.$(OBJEXT): mxprt.f dim.h
newbin.$(OBJEXT): newbin.f acom.h orbit.h dim.h
newsrc.$(OBJEXT): newsrc.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h
newval.$(OBJEXT): newval.f acom.h bcom.h dp.h eph.h trnsfr.h orbit.h clocks.h dim.h glitch.h
outpar.$(OBJEXT): outpar.f acom.h bcom.h dp.h trnsfr.h orbit.h clocks.h dim.h glitch.h tz.h eph.h
pcard.$(OBJEXT): pcard.f acom.h dim.h dp.h
resid.$(OBJEXT): resid.f acom.h bcom.h orbit.h dp.h dim.h glitch.h
setup.$(OBJEXT): setup.f acom.h dim.h
tdbinit.$(OBJEXT): tdbinit.f tdbcom.h
tdbinit2.$(OBJEXT): tdbinit2.f tdbcom.h
tdbread.$(OBJEXT): tdbread.f tdbcom.h
tempo.$(OBJEXT): tempo.f acom.h bcom.h dp.h array.h clocks.h eph.h vcom.h orbit.h dim.h trnsfr.h tz.h
timcalc.$(OBJEXT): timcalc.f acom.h dim.h dp.h trnsfr.h
tmalloc.$(OBJEXT): array.h
tparin.$(OBJEXT): tparin.f acom.h dim.h vcom.h
tpohdr.$(OBJEXT): tpohdr.f dim.h tz.h
tzfit.$(OBJEXT): tzfit.f dim.h acom.h trnsfr.h tz.h
tzinit.$(OBJEXT): tzinit.f dim.h tz.h
ut1red.$(OBJEXT): ut1red.f dim.h
vmemw.$(OBJEXT): vmemw.f dim.h vcom.h
ztim.$(OBJEXT): ztim.f trnsfr.h

# dependency of matrix:
matrix.$(OBJEXT): matrix.f dim.h


