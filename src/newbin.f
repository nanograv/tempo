c      $Id$
	subroutine newbin(nits,jits)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'acom.h'
	include 'orbit.h'
cvk     Parameters for finding mass function
	parameter (TWOPI=6.28318530717958648d0)
      	parameter (gm=1.3271243999e26)
	parameter (cvel=2.99792458e10)
	data e12/1.d12/,e6/1.d6/
	nerr(xx,n)=nint(xx*10.d0**n)

	write(31,1005)
1005    format(//'     A1 sin(i)        E          T0(MJD)      ',
     +    '        PB         OMEGA'/)

	jbin=1
	if(nbin.ge.9) jbin=nbin-7
	do 10 i=1,jbin
	   if(jbin.gt.1) write(31,1007) i
 1007	   format('Orbit # ',i2)
	   write(31,1006) a1(i),e(i),t0(i),pb(i)/8.64d4,omz(i)
 1006	   format(f14.9,f14.10,f16.9,f18.12,f11.6)
	   j=0
	   if(i.ge.2) j=17
	   if(i.eq.3) j=17+5
	   dt0=freq(j+11)
	   et0=ferr(j+11)
	   dom=freq(j+13)*57.2957795131d0
	   eom=ferr(j+13)*57.2957795131d0
	   
	   write(31,1006) freq(j+9),freq(j+10),dt0,freq(j+12)/8.64d4,dom
	   write(31,1006) ferr(j+9),ferr(j+10),et0,ferr(j+12)/8.64d4,eom
	   a1(i)=a1(i)+freq(j+9)
	   e(i)=e(i)+freq(j+10)
	   t0(i)=t0(i)+dt0
	   pb(i)=pb(i)+freq(j+12)
	   omz(i)=omz(i)+dom
	   write(31,1006) a1(i),e(i),t0(i),pb(i)/8.64d4,omz(i)

C 	   Write out mass function and error
	   fm=4.0d0*((twopi/2.0d0)**2)*((a1(i)*cvel)**3)/((pb(i)**2)*gm)
	   efm=(fm/(a1(i)*cvel))*sqrt(9.0d0*((ferr(j+9)*cvel)**2)+
     +		((4.0d0*(a1(i)*cvel)**2/(pb(i))**2)*(ferr(j+12)**2)))
	   write(31,*)
	   write(31,1009) fm,efm
 1009	   format(' Mass function: ',f12.10,' +/- ',f12.10,
     +		' solar masses')
	   pb(i)=pb(i)/8.64d4	! Convert back to days for iteration
	   ferr(j+12)=ferr(j+12)/8.64d4  ! and for output
 10	continue

	if(ferr(14).eq.0.0.and.ferr(15).eq.0.0.and.ferr(18).eq.0.0
     +  .and.ferr(20).eq.0.0.and.ferr(21).eq.0.0.and.ferr(22).eq.0.0
     +  .and.ferr(23).eq.0.0.and.ferr(24).eq.0.0.and.ferr(25).eq.0.0
     +  .and.ferr(37).eq.0.0.and.ferr(38).eq.0.0)
     +  go to 100


	if(nbin.eq.3) write(31,1050)
1050  format(//'     OMDOT     GAMMA    PBDOT(-12)   sin(i)',
     +  '       M          m2      DTH(-6)'/)
	if(nbin.eq.4) write(31,10501)
10501 format(//'    XOMDOT     XPBDOT(-12)',
     +  '     M         m2'/)
	if(nbin.ne.3.and.nbin.ne.4) write(31,10502)
10502 format(//'     OMDOT     GAMMA      PBDOT(-12)   sin(i)',
     +  '       m1         m2     bp   bpp'/)

       if(nbin.eq.3) then
         write(31,1052) omdot,gamma,pbdot*e12,si,am,am2,e6*dth
 1052    format(f11.7,f10.7,f14.6,f11.6,f11.7,2f11.4)
         dth=dth+freq(23)
         freq(23)=e6*freq(23)
         ferr(23)=e6*ferr(23)
       else if(nbin.eq.4)then
         write(31,10521) xomdot,xpbdot*e12,am,am2
10521    format(f11.7,f14.6,f11.6,f11.6)
       else
         write(31,10522) omdot,gamma,pbdot*e12,si,am,am2,bp,bpp
10522    format(f11.7,f10.7,f14.6,f11.6,f11.7,f11.4,f6.2,f5.1)
       endif
       
       if(nbin.eq.3) then
         write(31,1052) freq(14),freq(15),e6*freq(18),freq(20),
     +        freq(21),freq(22),freq(23)
         write(31,1052) ferr(14),ferr(15),e6*ferr(18),ferr(20),
     +        ferr(21),ferr(22),ferr(23)
       else if (nbin.eq.4) then
         write(31,10521) freq(37),e6*freq(38),
     +        freq(21),freq(22)
         write(31,10521) ferr(37),e6*ferr(38),
     +        ferr(21),ferr(22)
       else
         write(31,1052) freq(14),freq(15),e6*freq(18),freq(20),
     +        freq(21),freq(22)
         write(31,1052) ferr(14),ferr(15),e6*ferr(18),ferr(20),
     +        ferr(21),ferr(22)
       endif
       
       if(sim) then
         write(70,1110) nint(pb(1)*8.64d4),a1,e,nint(omz(1)),omdot,
     +        nint(e6*gamma),pbdot*e12,
     +        si,e6*dth,nerr(ferr(14),8),nerr(ferr(15),8),
     +        nerr(e6*ferr(18),6),nerr(ferr(20),5),nerr(ferr(22),4),
     +        nerr(1.d-6*ferr(23),8),nerr(ferr(25),18)
 1110    format(i6,2f6.3,i4,f7.3,i5,f8.3,f6.3,f5.1,7i5)
       endif

c  update parameters
	gamma=gamma+freq(15)
	si=si+freq(20)
	am=am+freq(21)
	am2=am2+freq(22)
        pbdot = pbdot + freq(18)/e6
        omdot = omdot + freq(14)
        xpbdot = xpbdot + freq(38)/e6
        xomdot = xomdot + freq(37)

c   print updated parameters

	if (nbin.eq.3) then
	  write(31,1052) omdot,gamma,pbdot*e12,si,am,am2,e6*dth
	else if (nbin.eq.4) then
	  write(31,10521) xomdot,xpbdot*e12,am,am2
	else
	  write(31,10522) omdot,gamma,pbdot*e12,si,am,am2,bp,bpp
	endif


       if(nbin.eq.4) then       !  write out calculated values
         omd=360.d0*365.25d0*xk/pb(1)
	 write(31,10511)
	 write(31,10512) omd
	 write(31,10513) gamma
	 write(31,10514) pbdot*e12
	 write(31,10515) si
	 write(31,10516) a0*e6
	 write(31,10517) a0aligned*e6
10511	 format (/" Calculated values (assuming GR):")
10512	 format (/" Omegadot:           ",f11.7)
10513	 format (" Gamma:               ",f10.7)
10514	 format (" Pbdot(-12):      ",f14.6)
10515	 format (" Sin i:              ",f11.7)
10516	 format (" A0 used (-6):       ",f11.7)
10517	 format (" A0 if aligned (-6): ",f11.7)
       endif


	if(nbin.ge.9) write(31,3001) freq(18),ferr(18)
3001	format(/(f12.1))

	xdot0=xdot
	xdot=xdot+freq(24)
	edot0=edot
	edot=edot+freq(25)
	if (ferr(24).ne.0.or.ferr(25).ne.0)then
          write(31,1053) e12*xdot0,e12*edot0,e12*freq(24),e12*freq(25),
     +     e12*ferr(24),e12*ferr(25),e12*xdot,e12*edot
	endif
1053	format(/'    XDOT(-12)   EDOT(-12)'//(f12.6,f12.6))

 100	continue

	if(gro.and.(nits.eq.0.or.jits.eq.nits)) then
	  open(34,file='gro.2',status='unknown')
	  t0mjd=t0(1)+40000.d0-0.5d0
	  write(34,1100) psrname(1:8),pb(1)*8.64d4,a1(1),e(1),t0mjd,
     +      omz(1),omdot,gamma,pbdot,obsflag
1100	  format(a8,f17.6,f12.7,f11.8,f15.8,f11.6,f8.5,f9.6,
     +      1p,d11.2,0p,1x,a1)
	  close(34)
	endif

	call outbinpar

	return
	end

