c      $Id$
        parameter(NMODELS=11)
        character bmodel(0:NMODELS)*8

	common/orbitp/ a1(4),e(4),t0(4),pb(4),omz(4),
     +     eps1,eps2,eps1dot,eps2dot,t0asc,
     +	   omdot,gamma,pbdot,si,am,am2,dth,xomdot,xpbdot,dr,a0,b0,xk,
     +     bp,bpp,xdot,edot,a0aligned,afac,om2dot,x2dot
        common/orbitn/ nbin,nplanets,nell1
	common/bmodel/ bmodel
