c      $Id$
	subroutine newval(chisqr,nfree,rms0,rms1,nits,jits,wmax,nboot)

	implicit real*8 (a-h,o-z)
	parameter (TWOPI=6.28318530717958648d0)
	character*1 decsgn,binflag,label*6
	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'trnsfr.h'
	include 'dp.h'
	include 'orbit.h'
	include 'clocks.h'
	include 'eph.h'
	include 'glitch.h'

	data label/'Offset'/

	if(jits.lt.nits)jits=jits+1

	write(31,1039) psrname,ephfile(nephem)(1:5),clklbl(nclk),
     :     pepoch
1039	format(/'PSR ',a12,'  Ephem.: ',a,'  Clock: ',a12,
     +       '  Ref. MJD: ',f12.4)


	if (eclcoord) then
	  write (31,1040)
 1040     format (/5x,'LAMBDA',15x,'BETA',7x,'PM LAMBDA',2x,'PM BETA',
     +         4x,'PM RV',3x,'PARALLAX'/)
          c1 = 360.d0/TWOPI     ! radians to degrees
          c = 360.*60.*60./TWOPI !radians to arcsec
          cc =c*365.25*86400000. !radians/sec to mas/yr
          write (31,1041) c1*pra, c1*pdec, pmra, pmdec, pmrv, px
 1041     format (2f18.13,4f10.5)
          freq5 = freq(5)*c1
          freq6 = freq(6)*c1
          write (31,1041) freq6, freq5, 1.d-9*freq(8)*cc, 
     +         1.d-9*freq(7)*cc, 1.d1*freq(36)*c, freq(17)
          ferr(5) = ferr(5)*c1
          ferr(6) = ferr(6)*c1
          pmrerr=1.d-9*ferr(8)*cc
          pmderr=1.d-9*ferr(7)*cc
          write(31,1041) ferr(6),ferr(5),pmrerr,pmderr,
     +       1.d1*ferr(36)*c,ferr(17)
          pra=pra+freq(6)
          pdec=pdec+freq(5)
          pmra=pmra+1.d-9*freq(8)*cc
          pmdec=pmdec+1.d-9*freq(7)*cc
          pmrv=pmrv+1.d1*freq(36)*c
          px=px+freq(17)
          write (31,1041) c1*pra, c1*pdec, pmra, pmdec, pmrv, px
	else
	  write (31,1048)
 1048     format (/6x,'RA',16x,'DEC',12x,'PM RA',4x,'PM DEC',5x,'PM RV',
     +         3x,'PARALLAX'/)

          c=360.*60./TWOPI*60.
          cc=c*365.25*86400000.

          call radian(pra,irh,irm,rsec,123,1)
          call radian(pdec,idd,idm,dsec,1,1)
          decsgn=' '
          if(pdec.lt.0.) decsgn='-'
          idd=iabs(idd)
          
          write(31,1050) irh,irm,rsec,decsgn,idd,idm,dsec,pmra,pmdec,
     +         pmrv,px
 1050     format (i2.2,i3.2,f12.8,1x,a1,i2.2,i3.2,f11.7,4f10.4)

          freq5=freq(5)*c
          freq6=freq(6)*c/15.
          write(31,1051)freq6,freq5,1.d-9*freq(8)*cc,1.d-9*freq(7)*cc,
     +         1.d1*freq(36)*c,freq(17)
 1051     format (f17.8,f18.7,4f10.4)
          
          ferr(5)=ferr(5)*c
          ferr(6)=ferr(6)*c/15.
          pmrerr=1.d-9*ferr(8)*cc
          pmderr=1.d-9*ferr(7)*cc
          write(31,1051) ferr(6),ferr(5),pmrerr,pmderr,
     +         1.d1*ferr(36)*c,ferr(17)

          pra=pra+freq(6)
          pdec=pdec+freq(5)
          call radian(pra,irh,irm,rsec,123,1)
          call radian(pdec,idd,idm,dsec,1,1)
          decsgn=' '
          if(pdec.lt.0.) decsgn='-'
          idd=iabs(idd)
          pmra=pmra+1.d-9*freq(8)*cc
          pmdec=pmdec+1.d-9*freq(7)*cc
          pmrv=pmrv+1.d1*freq(36)*c ! Convert from rad/century to mas/yr
          px=px+freq(17)
          write(31,1050) irh,irm,rsec,decsgn,idd,idm,dsec,pmra,pmdec,
     +         pmrv,px
        endif

	p0z=1.d0/f0
	p1z=-1.d15*f1/f0**2
	f0z=f0
	f1z=f1
	f2z=f2
	f3z=f3
	df0=-freq(2)*1.d-9
	df1=-freq(3)*1.d-18
	df2=-freq(4)*1.d-27
	df3=-freq(51)*1.d-36
	f0=f0z+df0
	f1=f1z+df1
	f2=f2z+df2
	f3=f3z+df3
	kf1=0
	kf2=0
	kf3=0
	if(abs(f1).gt.0.)kf1=int(-log10(abs(f1)))+1
	if(abs(f2).gt.0.)kf2=int(-log10(abs(f2)))+1
	if(abs(f3).gt.0.)kf3=int(-log10(abs(f3)))+1
	sf1=10**(dfloat(kf1))
	sf2=10**(dfloat(kf2))
	sf3=10**(dfloat(kf3))

	write(31,1052)kf1,kf2,kf3
 1052	format(/10x,'F0',18x,'F1(D-',i2.2,')',7x,'F2(D-',i2.2,')',5x,
     +       'F3(D-',i2.2,')'/)

	write(31,1053)f0z,sf1*f1z,sf2*f2z,sf3*f3z
1053	format(f22.17,f18.12,f14.9,f12.6)

	write(31,1053)df0,sf1*df1,sf2*df2,sf3*df3

	ef0=ferr(2)*1.d-9
	ef1=ferr(3)*1.d-18
	ef2=ferr(4)*1.d-27
	ef3=ferr(51)*1.d-36
	write(31,1053)ef0,sf1*ef1,sf2*ef2,sf3*ef3

	write(31,1053)f0,sf1*f1,sf2*f2,sf3*f3

	write(31,1054)
 1054	format(/10x,'P0',18x,'P1(D-15)',9x,'DM',9x,'DM1',7x,'PPN GAM'/)

	ppng=1.d0
	write(31,1055)p0z,p1z,dm,dmcof(1),ppng
 1055	format(f22.19,f18.12,3f12.6)

	p0=1.d0/f0
	p1=-1.d15*f1/f0**2
	write(31,1055)p0-p0z,p1-p1z,freq(16),freq(41),freq(19)

	p0e=1.d-9*ferr(2)/f0**2
	p1e=1.d-3*dsqrt((2.d9*ferr(2)*f1/f0**3)**2+(ferr(3)/f0**2)**2)
	write(31,1055)p0e,p1e,ferr(16),ferr(41),ferr(19)

	dm=dm+freq(16)
	dmcof(1)=dmcof(1)+freq(41)
	ppng=ppng+freq(19)
	write(31,1055)p0,p1,dm,dmcof(1),ppng

	if((nfit(7).ne.0.or.nfit(8).ne.0) .and. (.not.eclcoord)) 
     +       call propmo(pmra,pmdec,pmrerr,pmderr,pra,pdec)

C  Compute braking index
	if(nfit(4).ne.0) then
	  brkind=f0*f2/(f1*f1)
	  brkerr=f0*1.d-27*ferr(4)/(f1*f1)
	  write(31,1080) brkind,brkerr
1080      format(/'Braking index:',f13.4,' +/-',f13.4)

	  if(nfcalc.ge.4) then
	     do j=4,nfcalc,3
		ia=j-3
		ib=min(ia+2,nfcalc-3)
		write(31,1081) j,j+1,j+2
 1081		format(/25x,'f',z1,15x,'f',z1,16x,'f',z1)
		write(31,1082) (f4(i),i=ia,ib)
		write(31,1082) (-freq(51+i)*(1.d-9)**(i+4),i=ia,ib)
		write(31,1082) (ferr(51+i)*(1.d-9)**(i+4),i=ia,ib)
		write(31,1082) (f4(i)-freq(51+i)*(1.d-9)**(i+4),i=ia,ib)
 1082		format(15x,1p,3d18.8,0p)
	     enddo

	     do i=1,nfcalc-3
		f4(i)=f4(i)-freq(51+i)*(1.d-9)**(i+4)
	     enddo
	  endif
	endif

	if(ndmcalc.ge.3) then
	   write(31,1083) (i,i=2,9)
 1083	   format(/8(5x,'DM',z1,2x))
	   write(31,1084) (dmcof(i),i=2,ndmcalc-1)
	   write(31,1084) (freq(40+i),i=2,ndmcalc-1)
	   write(31,1084) (ferr(40+i),i=2,ndmcalc-1)
	   write(31,1084) (dmcof(i)+freq(40+i),i=2,ndmcalc-1)
 1084	   format(8f10.6)
	   do j=2,ndmcalc-1
	      dmcof(j)=dmcof(j)+freq(40+j)
	   enddo
	endif

	if(ngl.gt.0)then
	  do i=1,ngl
	    glphz=glph(i)
	    glf0pz=glf0p(i)
	    glf1pz=glf1p(i)
	    glf0d1z=glf0d1(i)
	    gltd1z=gltd1(i)
	    glph(i)=glphz-freq(60+(i-1)*NGLP+1)
	    glf0p(i)=glf0pz-freq(60+(i-1)*NGLP+2)*1.d-9
	    glf1p(i)=glf1pz-freq(60+(i-1)*NGLP+3)*1.d-18
	    glf0d1(i)=glf0d1z-freq(60+(i-1)*NGLP+4)*1.d-9
	    gltd1(i)=gltd1z-freq(60+(i-1)*NGLP+5)/86400.d0
	    write(31,1065)i,glepoch(i)
	    write(31,1066)
	    write(31,1067)glphz,glf0pz,glf1pz,glf0d1z,gltd1z
	    write(31,1067)glph(i)-glphz,glf0p(i)-glf0pz,glf1p(i)-glf1pz,
     :        glf0d1(i)-glf0d1z,gltd1(i)-gltd1z
	    write(31,1067)ferr(60+(i-1)*NGLP+1),
     :        ferr(60+(i-1)*NGLP+2)*1.d-9,ferr(60+(i-1)*NGLP+3)*1.d-18,
     :        ferr(60+(i-1)*NGLP+4)*1.d-9,ferr(60+(i-1)*NGLP+5)/86400.d0
	    write(31,1067)glph(i),glf0p(i),glf1p(i),glf0d1(i),gltd1(i)
	    if(glf0p(i)+glf0d1(i).ne.0.d0)then
	      fac=1.d0/(glf0p(i)+glf0d1(i))/86400.d0
	      iph=nint(glph(i))
	      fph=glph(i)-iph
	      glep1z=glepoch(i)-fph*fac
	      glep2z=glepoch(i)-(fph-sign(1.d0,fph))*fac
	      glepe=ferr(60+(i-1)*NGLP+1)*abs(fac)
	      write(31,1068)glep1z,glep2z,glepe
	    endif
	  enddo
	endif
1065	format(/' Glitch',i2,'  MJD:',f14.6)
1066    format('    glph',7x,'glf0 (Hz)',5x,'glf1 (Hz/s)',4x,
     :       'glf0d (Hz)',7x,'gltd (d)')
1067	format(f10.6,1p,3d15.6,0p,f15.6)
1068	format(' MJD for zero phase:',f14.6,' or',f14.6,'  Error:',f10.6)

C Output new parameters
	rewind 71
	call outpar(nits,irh,irm,rsec,ferr(6),decsgn,idd,idm,
     +       dsec,ferr(5))

C Output binary parameters
	if(a1(1).ne.0.0) call newbin(nits,jits)

	if(nxoff.gt.0) then
	  koff=60+NGLT*NGLP
	  do 70 n=1,(nxoff+3)/4
	  ib=n*4
	  ia=ib-3
	  ib=min(ib,nxoff)
	  write(31,1059) (label,i,i=ia,ib)
1059	  format(/1x,5(7x,a6,i3))

	  write(31,1060) ((xjdoff(j,i),j=1,2),i=ia,ib)
1060	  format(/'MJD:',5(f8.1,'-',f7.1))
	  do 63 i = ia, ib
	  do 63 j = 1, 2
63	  if(xjdoff(j,i).lt.0.001) xjdoff(j,i)=0.

	  write(31,1061) (dct(i),i=ia,ib)
1061	  format(f17.8,4f16.8)
	  write(31,1061) (freq(k)*p0,k=koff+ia,koff+ib)
	  write(31,1061) (ferr(k)*p0,k=koff+ia,koff+ib)
	  do 65 i=ia,ib
65	  dct(i)=dct(i)+freq(koff+i)*p0
	  write(31,1061) (dct(i),i=ia,ib)
70	  continue
	endif



	asig=asig*p0*1000.
	if(nboot.gt.0) write(31,1085) nboot
1085	format('Uncertainties by bootstrap Monte Carlo:',
     +    i6,' iterations.')
	rms0=1000.d0*sigma1
	rms1=1000.d0*asig
	write(31,1100) rms0,rms1
1100	format(/'Weighted RMS residual: pre-fit',f10.3,
     +  ' us. Predicted post-fit',f10.3,' us.')
	write(*,1101)  rms0,rms1
1101	format(/' Weighted RMS residual: pre-fit',f10.3,
     +  ' us. Predicted post-fit',f10.3,' us.')
	if(chisqr.ne.0.d0) then
	  write(31,1110) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
	  write(*,1110) chisqr*nfree,nfree,chisqr,rms0/rms1,wmax
1110	  format(' Chisqr/nfree:',f9.2,'/',i5,' =',f9.4,
     +    '   pre/post:',f7.2,'   Wmax:',f7.1)
	endif

	if(gro.and.(nits.eq.0.or.jits.eq.nits)) then
	  open(33,file='gro.1',status='unknown')
	  if(dt2sec.gt.0.d0) dt2sec=dt2sec-p0
	  t0geo=t0geo-dt2sec/86400.d0
	  if(abs(t0geo-pepoch).gt.1.d-3) write(*,1112) t0geo,pepoch
1112	  format(/5x,'### Warning: t0geo=',f10.3,' and pepoch=',f10.3,
     +      ' do not match! ###')
	  dt=(t0geo-pepoch)*86400.d0
	  ff0=f0 + f1*dt + 0.5d0*f2*dt**2
	  ff1=f1 + f2*dt
	  pp0=1.d0/ff0
	  pp1=-ff1*ff0**(-2)
	  pp2=-f2*ff0**(-2) + 2.d0 * ff1**2 * ff0**(-3)
	  open(99,file='gro.99',status='unknown')
	  if(oldpar)then
	     write(99,1113) pp0,1.d15*pp1,t0geo,1.d30*pp2
 1113	     format('P',f18.16,1x,f12.5,8x,f11.5,9x,f12.1)
	  else
	     write(99,'(''F0'',f22.16)')ff0
	     write(99,'(''F1'',1p,d22.12)')ff1
	     write(99,'(''F2'',1p,d22.12)')f2
	     write(99,'(''PEPOCH'',f18.6)')t0geo
	  endif
	  close(99)
	  f0=ff0
	  f1=ff1
	  rmsmp=asig/p0
	  binflag=' '
	  if(a1(1).ne.0.0) binflag='*'
	  write(33,1120) psrname(1:8),irh,irm,rsec,decsgn,idd,idm,dsec,
     +      mjd1,mjd2,t0geo,f0,f1,f2,rmsmp,obsflag,binflag,
     +      ephfile(nephem)(1:5),psrname
1120	  format(a8,2i3.2,f7.3,1x,a1,i2.2,i3.2,f6.2,2i6,f16.9,
     +      f18.13,1p,d13.5,d11.2,0p,f5.1,2(1x,a1),1x,a,1x,a)
	  close(33)
	endif

	return
	end
