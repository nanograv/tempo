c      $Id$
      subroutine fit(npts,mode,chisqr,varfit,xmean,ymean,sum,nz,
     +     wmax,lw,ddmch,
     +     buf,npmsav,ksav)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'orbit.h'
	include 'tz.h'

c       DJN changes Sep '98
c       --reformat code (indent loops, end all loops with 'continue')
c       --get rid of unneeded 'rmul' calculation
c       --get rid of separate sigmax() calculation since sigmax(i)=array(i,i)
c       --don't calculate gcor()=array(i,i) just before matinv() since 
c         array(i,i)=1.0 always here
c       --add 120/121 loops to un-normalize array() rather than doing so
c         in sigma0 calculation (this way un-normalized array is eventually
c         saved in common block, as desired for covar().)

c       moved declaration of real*8 array(NPA,NPA) to acom.h, djn, 8 Sep 98
	real*8 xmean(NPA),sigmax(NPA),fctn(NPAP1)
	real*8 r(NPA),a(NPA),sigmaa(NPA),gcor(NPA)
	logical lw
	real*8 ddmch(*)
        real*8 buf(*)
        integer npmsav(*), ksav(*)
	character mark*1,adn*14,date*9,damoyr*9

	rewind 32
	mprt=10**nprnt
	sigma=0.
	chisq=0.

	nterms=nparam-1
	do 29 j=1,nterms
          sigmax(j)=0.
          r(j)=0.
          a(j)=0.
          sigmaa(j)=0.
          do 28 k=1,nterms
            array(j,k)=0.
 28       continue
 29     continue

        ymean=ymean/sum
	do 53 j=1,nterms
          xmean(j)=xmean(j)/sum
 53     continue
	fnpts=npts-nz
	wmean=sum/fnpts

	do 67 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          weight=weight/wmean
          sigma=sigma+weight*(y-ymean)**2
          do 66 j=1,nterms
            fx=fctn(j)-xmean(j)
            wfx=weight*fx
            sigmax(j)=sigmax(j)+wfx*fx
            r(j)=r(j)+wfx*(y-ymean)
            do 65 k=1,j
              array(j,k)=array(j,k) + wfx * (fctn(k)-xmean(k))
 65         continue
 66       continue
 67     continue

	do 68 i = 1, nterms
	  sigmax(i) = array(i,i)
 68	continue

	free1=fnpts-1
	sigma=dsqrt(sigma/free1)
	do 78 j=1,nterms
          sigmax(j)=dsqrt(sigmax(j)/free1)
          fs=free1*sigmax(j)
          r(j)=r(j)/(fs*sigma)
          do 77 k=1,j
            array(j,k)=array(j,k)/(fs*sigmax(k))
            array(k,j)=array(j,k)
 77       continue
 78     continue
	call matinv(array,nterms,det)
	if(det.eq.0.d0) stop 'fit 78'
	do 80 i=1,nterms
          gcor(i)=sqrt(abs(1.d0 - 1.d0/array(i,i)))
 80     continue
	aa0=ymean

	if(lw) then
	  write(31,1170)
1170	  format(//'Normalized covariance matrix in "number of 9''s"',
     +      ' format:'/)
	  call mxprt(array,gcor,nterms,mfit,nbin,eclcoord)
	  unit=1.d-6
	  if(p0firs.gt.0.1) unit=1.d-3
	  nunit=dlog10(unit)-0.0001
	  write(31,1050) nunit
1050	  format(//'   N     TOBS         Date',
     +     '      FREQ    WGT      NPULSE     R1(P)  R2(',i2,') PH/DT'/)
	endif

	do 106 j=1,nterms
          do 104 k=1,nterms
            a(j)=a(j)+r(k)*array(j,k)
 104      continue
          a(j)=a(j)*sigma/sigmax(j)
          aa0=aa0-a(j)*xmean(j)
 106    continue

	do 108 i=1,npts
          call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms,
     +     buf,npmsav,ksav)
          dt2=y-aa0
          weight=weight/wmean
          do 107 j=1,nterms
            dt2=dt2-a(j)*fctn(j)
 107      continue
          nct=ct
          mark=char(mod(nct,26)+65)
          date=damoyr(int(fmjd))
          dt2sec=dt2*p0firs
          phase=0.
          if(a1(1).ne.0.)
     +         phase=dmod((ct-t0(1))*86400.d0/pb(1)+1000000.d0,1.d0)

          if(nprnt.eq.0.or.mod(i,mprt).eq.1) then
            write(adn,1106) dn
 1106       format(f14.0)

            if(lw) write(31,1107) i,fmjd,mark,date,
     +           rfrq,weight,adn(1:13),y,dt2sec/unit,phase
 1107       format(i4,f12.5,1x,a1,1x,a9,f9.3,f6.3,a13,
     +           f8.3,f8.2,f7.3)
          endif

C Correct tz ref TOA
          if(i.eq.ntzref) tzrmjd=tzrmjd-dt2sec/8.64d4

          if((i.lt.npts.or.(.not.gro)).and.lw) write(32) ct,dt2,
     +         dt2sec,phase,frq,weight,terr,y,ddmch(i)
          wmax=max(wmax,weight)
          chisq=chisq+weight*dt2**2
 108    continue

	freen=fnpts-nterms-1
	chisqr=chisq*wmean/freen
	if(mode.eq.0) varnce=chisqr
	if(mode.ge.1) varnce=1./wmean
	varfit=chisq/fnpts
	if(mode.eq.0) chisqr=0.

	do 120 j = 1, nterms
	  do 121 k = 1, nterms
	    array(j,k) = array(j,k)*varnce/(free1*sigmax(j)*sigmax(k))
 121	  continue
 120	continue

	do 133 j=1,nterms
          sigmaa(j)=dsqrt(array(j,j))
 133    continue
	sigma0=varnce/fnpts
	do 145 j=1,nterms
          do 144 k=1,nterms
            sigma0=sigma0+xmean(j)*xmean(k)*array(j,k) 
 144      continue
 145    continue

	sigma0=dsqrt(sigma0)

	do 146 j=1,NPAP1
          freq(j)=0.
          ferr(j)=0.
 146      continue

	freq(1)=aa0
	ferr(1)=sigma0

	do 150 j=1,nterms
          freq(mfit(j+1))=a(j)
          ferr(mfit(j+1))=sigmaa(j)
 150    continue


	return
	end
