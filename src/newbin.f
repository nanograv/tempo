c      $Id$
	subroutine newbin(nits,jits)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'acom.h'
	include 'orbit.h'
        logical err
cvk     Parameters for finding mass function
	parameter (TWOPI=6.28318530717958648d0)
      	parameter (gm=1.3271243999e26)
	parameter (cvel=2.99792458e10)
	data e12/1.d12/,e6/1.d6/
	nerr(xx,n)=nint(xx*10.d0**n)

c --- output of Keplerian parameters --- 

	if(nbin.ne.9)then
	   write(31,'(//5x,''A1 sin(i)'',8x,''E'',10x,''T0(MJD)'',14x,
     +          ''PB'',9x,''OMEGA''/)')
	else
	   write(31,'(//4x,''A1 sin(i)'',4x,''EPS1'',10x,''EPS2'',12x,
     +         ''PB(days)'',11x,''TASC(MJD)''/)')
	endif 

	do i=1,1+nplanets
	   if(nplanets.gt.0) write(31,1005) i
 1005	   format('Orbit # ',i2)
	   if(nbin.ne.9)then
	      write(31,1006) a1(i),e(i),t0(i),pb(i)/8.64d4,omz(i)
 1006	      format(f14.9,f14.10,f16.9,f18.12,f11.6)
	   else
	      write(31,2006) a1(i),eps1,eps2,pb(i)/8.64d4,t0asc
 2006         format(f14.9,f14.10,f14.10,f18.12,f18.11)
           endif
	      
	   j=0
	   if(i.ge.2) j=17
	   if(i.eq.3) j=17+5
	   dt0=freq(j+11)
	   et0=ferr(j+11)
	   dpb=freq(j+12)/8.64d4
	   epb=ferr(j+12)/8.64d4
	   dom=freq(j+13)*57.2957795131d0
	   eom=ferr(j+13)*57.2957795131d0

	   if(nbin.ne.9)then
	      write(31,1006) freq(j+9),freq(j+10),dt0,dpb,dom
	      write(31,1006) ferr(j+9),ferr(j+10),et0,epb,eom
	   else
              write(31,2006) freq(9),freq(10),freq(13),dpb,dt0
              write(31,2006) ferr(9),ferr(10),ferr(13),epb,et0
	   endif

C  Update parameters

	    a1(i) = a1(i)  + freq(j+9)
	    pb(i) = pb(i)  + freq(j+12)
	    if(nbin.ne.9)then
	       e(i)   = e(i)   + freq(j+10)
	       omz(i) = omz(i) + dom
	       t0(i)  = t0(i)  + dt0
	    else
	       eps1  = eps1   + freq(10)
	       eps2  = eps2   + freq(13)
	       t0asc = t0asc  + dt0
	    endif

	    if(nbin.ne.9)then
	       write(31,1006) a1(i),e(i),t0(i),pb(i)/8.64d4,omz(i)
	    else
	       write(31,2006) a1(i),eps1,eps2,pb(i)/8.64d4,t0asc

C  Print (calculated) eccentricity and omega (for ELL1 model) and 
C  maximum error of ELL1 model. Since ELL1 neglects effects of order
C  e^2 errors are .le. a1*e^2

	       ecal=DSQRT(eps1**2+eps2**2)
	       if(ecal.gt.0)then                                       
  	          eeps1=ferr(10)
	          eeps2=ferr(13)
		  eecal=((eps1*eeps1)**2+(eps2*eeps2)**2)**.5/ecal 
		  write(31,2007) ecal,eecal                        
 2007		  format(/'Eccentricity, E:              ',
     +                 f14.12,' +/- ',f14.12)       
		  omcal=360.d0/TWOPI*DATAN2(eps1,eps2)
		  if(omcal.lt.0) omcal=omcal+360                      
		  eomcal=((eps2*eeps1)**2+(eps1*eeps2)**2)**.5/ecal**2
		  eomcal=360.d0/TWOPI*eomcal                      
		  write(31,2008) omcal,eomcal                       
 2008		  format('Longitude of periastron, OM:  ',
     +                 f14.10,' +/- ',f14.10,' deg')
		  t0peri=t0asc+omcal/360.0*pb(i)/8.64d4
		  et0peri=eomcal/360.0*pb(i)/8.64d4
		  write(31,2009) t0peri,et0peri
 2009		  format('Time of periastron, T0:       ',
     +                 f14.8,' +/- ',f14.8,' MJD')
		  write(31,2010) ecal**2*a1(1)*1.d9  		  
 2010		  format(/'[ Difference ELL1-DD < ',f7.1,' ns ]')    
 	       endif
	    endif                                                        

C  Write out mass function and error

c	    fm=4.d0*((twopi/2.0d0)**2)*((a1(i)*cvel)**3)/((pb(i)**2)*gm)
c	    efm=(fm/(a1(i)*cvel))*sqrt(9.0d0*((ferr(j+9)*cvel)**2)+
c     +		((4.0d0*(a1(i)*cvel)**2/(pb(i))**2)*(ferr(j+12)**2)))

	    fm = twopi**2 * (a1(i)*cvel)**3 / (pb(i)**2 * gm)
            call covar (j+9,j+12,c,err)
            if (.not. err) then
              efm = fm * sqrt(    9.d0*(ferr(j+9)/a1(i))**2 
     +             +  4.d0*(ferr(j+12)/pb(i))**2  
     +             - 12.d0*c/(a1(i)*pb(i))        )
              write(31,1009) fm,efm
 1009         format(/'Mass function: ',f13.10,' +/- ',f13.10,
     +             ' solar masses')
            else
              write(31,1010) fm
 1010         format(/'Mass function: ',f13.10,
     +             ' solar masses')
            endif

	    pb(i)=pb(i)/8.64d4	! Convert back to days for iteration
	    ferr(j+12)=ferr(j+12)/8.64d4 ! and for output
	 enddo

c --- output of non-Keplerian parameters, I ---

	if(     omdot.eq.0.0   .and.gamma.eq.0.0   .and.pbdot.eq.0.0
     +     .and.si.eq.0.0      .and.am.eq.0.0      .and.am2.eq.0.0
     +     .and.dth.eq.0.0     .and.xdot.eq.0.0    .and.edot.eq.0.0
     +     .and.xomdot.eq.0.0  .and.xpbdot.eq.0.0
     +     .and.om2dot.eq.0.0  .and.x2dot.eq.0.0
     +     .and.ferr(14).eq.0.0.and.ferr(15).eq.0.0.and.ferr(18).eq.0.0
     +     .and.ferr(20).eq.0.0.and.ferr(21).eq.0.0.and.ferr(22).eq.0.0
     +     .and.ferr(23).eq.0.0.and.ferr(24).eq.0.0.and.ferr(25).eq.0.0
     +     .and.ferr(37).eq.0.0.and.ferr(38).eq.0.0
     +     .and.ferr(39).eq.0.0.and.ferr(40).eq.0.0)
     +     goto 100

        if(nbin.eq.3.or.nbin.eq.8)then
	   write(31,10503)
10503	   format(//'     OMDOT     GAMMA    PBDOT(-12)   sin(i)',
     +              '       M          m2      DTH(-6)'/)
           write(31,10523) omdot,gamma,pbdot*e12,si,am,am2,e6*dth
           write(31,10523) freq(14),freq(15),e6*freq(18),freq(20),
     +        freq(21),freq(22),freq(23)*e6
           write(31,10523) ferr(14),ferr(15),e6*ferr(18),ferr(20),
     +        ferr(21),ferr(22),ferr(23)*e6
10523      format(f11.7,f10.7,f14.6,f11.6,f11.7,2f11.4)
        else if(nbin.eq.4)then
           write(31,10504)
10504      format(//'    XOMDOT     XPBDOT(-12)     M         m2'/)
           write(31,10524) xomdot,xpbdot*e12,am,am2
           write(31,10524) freq(37),e6*freq(38),freq(21),freq(22)
           write(31,10524) ferr(37),e6*ferr(38),ferr(21),ferr(22)
10524	   format(f11.7,f14.6,f11.6,f11.6)
        else if(nbin.eq.7)then
           write(31,10507)
10507      format(//'     OMDOT     GAMMA      PBDOT(-12)   sin(i)',
     +              '       m1         m2     bp   bpp'/)
           write(31,10527) omdot,gamma,pbdot*e12,si,am,am2,bp,bpp
           write(31,10527) freq(14),freq(15),e6*freq(18),freq(20),
     +        freq(21),freq(22)
           write(31,10527) ferr(14),ferr(15),e6*ferr(18),ferr(20),
     +        ferr(21),ferr(22)
10527      format(f11.7,f10.7,f14.6,f11.6,f11.7,f11.4,f6.2,f5.1)
        else if(nbin.eq.9)then
	   if(nell1.eq.0)then
	      write(31,10509)
10509	      format(//'   XDOT(-12)   EPS1DOT(-12)  EPS2DOT(-12)',
     +           '    PBDOT(-12)   sin(i)    m2'/)
	      write(31,10529) xdot*e12,eps1dot*e12,eps2dot*e12,
     +	         pbdot*e12,si,am2
	      write(31,10529) freq(24)*e12,freq(39)*e12,freq(40)*e12,
     +           e6*freq(18),freq(20),freq(22)
	      write(31,10529) ferr(24)*e12,ferr(39)*e12,ferr(40)*e12,
     +           e6*ferr(18),ferr(20),ferr(22)
10529	      format(f13.9,f14.9,f14.9,f14.9,f10.6,f10.6)
	   else
	      write(31,10559)
10559	      format(//'    XDOT(-12)  OMDOT(deg/yr)   EDOT(-12) ',
     +           '    PBDOT(-12)    sin(i)        m2'/)
	      write(31,10579) xdot*e12,omdot,edot*e12,pbdot*e12,si,am2
	      write(31,10579) freq(24)*e12,freq(14),freq(25)*e12,
     +           e6*freq(18),freq(20),freq(22)
	      write(31,10579) ferr(24)*e12,ferr(14),ferr(25)*e12,
     +           e6*ferr(18),ferr(20),ferr(22)
10579	      format(f13.9,f14.9,f14.9,f14.9,f10.6,f10.6)
	   endif
	else
           write(31,10501)
10501      format(//'     OMDOT     GAMMA      PBDOT(-12)   sin(i)',
     +              '       m1         m2'/)
           write(31,10521) omdot,gamma,pbdot*e12,si,am,am2
           write(31,10527) freq(14),freq(15),e6*freq(18),freq(20),
     +        freq(21),freq(22)
           write(31,10527) ferr(14),ferr(15),e6*ferr(18),ferr(20),
     +        ferr(21),ferr(22)
10521      format(f11.7,f10.7,f14.6,f11.6,f11.7,f11.4)
	endif       
        
        if(sim) then
           write(70,1110) nint(pb(1)*8.64d4),a1,e,nint(omz(1)),omdot,
     +        nint(e6*gamma),pbdot*e12,
     +        si,e6*dth,nerr(ferr(14),8),nerr(ferr(15),8),
     +        nerr(e6*ferr(18),6),nerr(ferr(20),5),nerr(ferr(22),4),
     +        nerr(ferr(23),8),nerr(ferr(25),18)
 1110      format(i6,2f6.3,i4,f7.3,i5,f8.3,f6.3,f5.1,7i5)
	endif

c  Update parameters

        omdot  = omdot  + freq(14)
	gamma  = gamma  + freq(15)
        pbdot  = pbdot  + freq(18)/e6
	si     = si     + freq(20)
	am     = am     + freq(21)
	am2    = am2    + freq(22)
	dth    = dth    + freq(23)
	xpbdot = xpbdot + freq(38)/e6
	xomdot = xomdot + freq(37)
	if(nbin.eq.9)then
	   xdot    = xdot    + freq(24)
	   edot    = edot    + freq(25)
	   eps1dot = eps1dot + freq(39)
	   eps2dot = eps2dot + freq(40)
	endif

c  Print updated parameters

	if(nbin.eq.3.or.nbin.eq.8)then
	   write(31,10523) omdot,gamma,pbdot*e12,si,am,am2,e6*dth
	else if (nbin.eq.4)then
	   write(31,10524) xomdot,xpbdot*e12,am,am2
	else if (nbin.eq.7)then
	   write(31,10527) omdot,gamma,pbdot*e12,si,am,am2,bp,bpp
	else if (nbin.eq.9)then
	   if(nell1.eq.1)then
	      write(31,10579) xdot*e12,omdot,edot*e12,pbdot*e12,si,am2
	   else
	      write(31,10529) xdot*e12,eps1dot*e12,eps2dot*e12,
     +           pbdot*e12,si,am2
	   endif
	else
	   write(31,10521) omdot,gamma,pbdot*e12,si,am,am2
	endif

	if(nbin.eq.3 .and. si.ne.0. .and. am2.ne.0.) then
	  amtot = (am2*si/(a1(1)*cvel))**1.5 * gm**0.5 * 
     +		pb(1)*86400.d0/twopi
	  write (31,10510) amtot
10510	  format (/'MTOT derived from sin i, M2: ',f10.7)
	endif

	if(nbin.eq.4) then  !  write out calculated values for DDGR
	   omd=360.d0*365.25d0*xk/pb(1)
	   write(31,10511)
	   write(31,10512) omd
	   write(31,10513) gamma
	   write(31,10514) pbdot*e12
	   write(31,10515) si
	   write(31,10516) a0*e6
	   write(31,10517) a0aligned*e6
10511	   format (/' Calculated values (assuming GR):')
10512	   format (/' Omegadot:           ',f11.7)
10513	   format (' Gamma:               ',f10.7)
10514	   format (' Pbdot(-12):     ',f15.7)
10515	   format (' Sin i:              ',f11.7)
10516	   format (' A0 used (-6):       ',f11.7)
10517	   format (' A0 if aligned (-6): ',f11.7)
	endif

C  Print (calculated) edot and omegadot for ELL1 model, if nell1=0 
C  (i.e. fit for eps1dot and/or eps2dot)

	if(nbin.eq.9.and.nell1.eq.0)then
	   ee1dot=ferr(39)
	   ee2dot=ferr(40)
           if(ee1dot.ne.0.or.ee2dot.ne.0)then
              edotcal=(eps1*eps1dot+eps2*eps2dot)/ecal
              eedotcal=((eps1*ee1dot)**2+(eps2*ee2dot)**2)**.5/ecal
              write(31,20511) edotcal*e12,eedotcal*e12
20511	      format(/'edot (-12 s-1): ',f12.8,' +/- ',f12.8)
              omdotcal=(eps2*eps1dot-eps1*eps2dot)/ecal**2
              eomdotcal=((eps1*ee2dot)**2+(eps2*ee1dot)**2)**.5/ecal**2
              omdotcal=omdotcal*360.0/TWOPI*365.25*86400
              eomdotcal=eomdotcal*360.0/TWOPI*365.25*86400
              write(31,20512) omdotcal,eomdotcal
20512	      format('omdot (deg/yr): ',f12.8,' +/- ',f12.8)
           endif
        endif

	if(nplanets.gt.0) write(31,3001) freq(18),ferr(18)
 3001	format(/(f12.1))

c --- Output of non-Keplerian parameters, II ---

	if(nbin.eq.9) goto 100
	   
	if(     xdot.eq.0.0  .and.edot.eq.0.0
     +     .and.om2dot.eq.0.0.and.x2dot.eq.0.0
     +     .and.ferr(24).eq.0.and.ferr(25).eq.0
     +     .and.ferr(39).eq.0.and.ferr(40).eq.0)
     +     goto 100

	if(nbin.eq.8)then
	   write(31,10608)
10608	   format(//'    XDOT(-12)   EDOT(-12)  OM2DOT(s-2)  ',
     +          'X2DOT(s-1)'/)
	   write(31,10628) xdot*e12,edot*e12,om2dot,x2dot
	   write(31,10628) freq(24)*e12,freq(25)*e12,freq(39),freq(40)
	   write(31,10628) ferr(24)*e12,ferr(25)*e12,ferr(39),ferr(40)	   
10628	   format(f12.6,f12.6,3x,e10.4,3x,e10.4)
	else
	   write(31,10601)
10601	   format(//'    XDOT(-12)   EDOT(-12)'/)
	   write(31,10621) xdot*e12,edot*e12
	   write(31,10621) freq(24)*e12,freq(25)*e12
	   write(31,10621) ferr(24)*e12,ferr(25)*e12
10621	   format(f12.6,f12.6)
	endif

c  update parameters

	xdot   = xdot   + freq(24)
	edot   = edot   + freq(25)
	om2dot = om2dot + freq(39)
	x2dot  = x2dot  + freq(40)

c  print updated parameters

	if(nbin.eq.8)then
	   write(31,10628) xdot*e12,edot*e12,om2dot,x2dot 
	else
	   write(31,10621) xdot*e12,edot*e12  
	endif

c ------

 100	continue

	if(gro.and.(nits.eq.0.or.jits.eq.nits)) then
	   open(34,file='gro.2',status='unknown')
	   write(34,1100) psrname(1:8),pb(1)*8.64d4,a1(1),e(1),t0(1),
     +        omz(1),omdot,gamma,pbdot,obsflag
1100	   format(a8,f17.6,f12.7,f11.8,f15.8,f11.6,f8.5,f9.6,
     +        1p,d11.2,0p,1x,a1)
	   close(34)
	endif

	call outbinpar

	return
	end

