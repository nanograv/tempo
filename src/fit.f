c      $Id$
	subroutine fit(npts,mode,chisqr,varfit,xmean,ymean,sum,nz,wmax,lw)

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'acom.h'
	include 'bcom.h'
	include 'orbit.h'
	include 'tz.h'

	real*8 xmean(NPA),sigmax(NPA),array(NPA,NPA),fctn(NPAP1)
	real*8 r(NPA),a(NPA),sigmaa(NPA),gcor(NPA)
	logical lw
	common/dmch/ddmch(NPTSMAX)
	character mark*1,adn*14,date*9,damoyr*9

	rewind 32
	mprt=10**nprnt
	sigma=0.
	chisq=0.
	rmul=0.
	nterms=nparam-1
	do 28 j=1,nterms
	sigmax(j)=0.
	r(j)=0.
	a(j)=0.
	sigmaa(j)=0.
	do 28 k=1,nterms
28      array(j,k)=0.

	ymean=ymean/sum
	do 53 j=1,nterms
53	xmean(j)=xmean(j)/sum
	fnpts=npts-nz
	wmean=sum/fnpts

	do 67 i=1,npts
	call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms)
	weight=weight/wmean
	sigma=sigma+weight*(y-ymean)**2
	do 67 j=1,nterms
	fx=fctn(j)-xmean(j)
	wfx=weight*fx
	sigmax(j)=sigmax(j)+wfx*fx
	r(j)=r(j)+wfx*(y-ymean)
	do 67 k=1,j
67	array(j,k)=array(j,k) + wfx * (fctn(k)-xmean(k))

	free1=fnpts-1
	sigma=dsqrt(sigma/free1)
	do 78 j=1,nterms
	sigmax(j)=dsqrt(sigmax(j)/free1)
	fs=free1*sigmax(j)
	r(j)=r(j)/(fs*sigma)
	do 78 k=1,j
	array(j,k)=array(j,k)/(fs*sigmax(k))
78	array(k,j)=array(j,k)
	do 79 i=1,nterms
79	gcor(i)=array(i,i)
	call matinv(array,NPA,nterms,det)
	if(det.eq.0.d0) stop 'fit 78'
	do 80 i=1,nterms
80	gcor(i)=sqrt(abs(1.d0 - 1.d0/(gcor(i)*array(i,i))))
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
104	a(j)=a(j)+r(k)*array(j,k)
	a(j)=a(j)*sigma/sigmax(j)
106	aa0=aa0-a(j)*xmean(j)

	do 108 i=1,npts
	call vmemr(i,fctn,ct,y,weight,dn,terr,frq,fmjd,rfrq,nterms)
	dt2=y-aa0
	weight=weight/wmean
	do 107 j=1,nterms
107	dt2=dt2-a(j)*fctn(j)
	nct=ct
	mark=char(mod(nct,26)+65)
	date=damoyr(int(fmjd))
	dt2sec=dt2*p0firs
	phase=0.
	if(a1(1).ne.0.)
     +    phase=dmod((ct-t0(1))*86400.d0/pb(1)+1000000.d0,1.d0)

	if(nprnt.eq.0.or.mod(i,mprt).eq.1) then
	write(adn,1106) dn
1106	format(f14.0)

	if(lw) write(31,1107) i,fmjd,mark,date,
     +    rfrq,weight,adn(1:13),y,dt2sec/unit,phase
1107    format(i4,f12.5,1x,a1,1x,a9,f9.3,f6.3,a13,
     +    f8.3,f8.2,f7.3)
	endif

C Correct tz ref TOA
	if(i.eq.ntzref) tzrmjd=tzrmjd-dt2sec/8.64d4

	if((i.lt.npts.or.(.not.gro)).and.lw) write(32) ct,dt2,
     +    dt2sec,phase,frq,weight,terr,y,ddmch(i)
	wmax=max(wmax,weight)
108	chisq=chisq+weight*dt2**2

	freen=fnpts-nterms-1
	chisqr=chisq*wmean/freen
	if(mode.eq.0) varnce=chisqr
	if(mode.ge.1) varnce=1./wmean
	varfit=chisq/fnpts
	if(mode.eq.1) varfit=chisq/fnpts
	if(mode.eq.0) chisqr=0.
	do 133 j=1,nterms
	sigmaa(j)=dsqrt(array(j,j)*varnce/free1)/sigmax(j)
133	rmul=rmul+a(j)*r(j)*sigmax(j)/sigma
C	freej=nterms
C	ftest=(rmul/freej)/((1.-rmul)/freen)
	rmul=dsqrt(rmul)
	sigma0=varnce/fnpts
	do 145 j=1,nterms
	do 145 k=1,nterms
145	sigma0=sigma0+varnce*xmean(j)*xmean(k)*array(j,k) /
     +  (free1*sigmax(j)*sigmax(k))
	sigma0=dsqrt(sigma0)

	do 146 j=1,NPAP1
	freq(j)=0.
146	ferr(j)=0.

	freq(1)=aa0
	ferr(1)=sigma0

	do 150 j=1,nterms
	freq(mfit(j+1))=a(j)
150	ferr(mfit(j+1))=sigmaa(j)

	return
	end
