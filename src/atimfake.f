c      $Id$
C	@(#)atimfake.f	9.19 2/3/94

	subroutine atimfake(oldpar,afmjd,nbin,nt,ncoord,alng,ipsr)
	implicit real*8 (a-h,o-z)
	parameter (n=31,pi=3.141592653589793d0)
	real*8 x(n)
	logical oldpar

	include 'dim.h'
	include 'tz.h'

	if(oldpar)then
	   write(50,1012) nbin,ncoord
 1012	   format('2',25x,i1,18x,'1',i1)
	   if(nbin.ge.1) then
	      np=6
	   else
	      np=4
	   end if
	   write(50,1000) (params(i),i=1,np)
 1000	   format(a80)
	   read(params(np),1013) tzrfrq
 1013	   format(15x,f9.0)
	else
	   i=0                                         ! Write ref toa line
	   write(50,1040)tzrsite,i,pname(1:8),tzrfrq,tzrmjd
	endif

	nspan=nsp(ipsr)
	maxha=mxha(ipsr)
	tzfreq=tzof(ipsr)
	if(tzfreq.lt.0.)tzfreq=tzrfrq
 
	hlst=24.d0*dmod(1.002737909d0*afmjd+0.154374d0 -
     +    alng/6.2831853071795864d0,1.d0)	
	read(pname,1030) irah,iram		
1030	format(2i2)
	rax=irah+(iram+2)/60.			
	wait=(rax-hlst)*0.99727				!Solar to sidereal
	if(wait.lt.-maxha) wait=wait + 23.9345
	fmjd1=afmjd+(wait-maxha)/24.+nspan/2880.		!Compute start time
	fmjd1=nint(48*fmjd1)/48.d0	
	nsets=(120*maxha+nspan-1)/nspan
 
	b=nspan/2 + 5
	a=-b
	bma=0.5*(b-a)
	bpa=0.5*(b+a)
	do 30 k=1,n
30	x(k)=cos(pi*(k-0.5)/n)*bma+bpa
 
	i=0
	do 50 j=1,nsets				
	fmjd2=fmjd1+(j-1)*nspan/1440.d0
	do 40 k=1,n
	i=i+1
	if(i.gt.800) stop ' Nspan too small'
	tfmjd(i)=fmjd2+x(k)/1440.d0
	tmin(i)=1440.d0*(tfmjd(i)-tfmjd(1))
	write(50,1040) tzsite,i,pname(1:8),tzfreq,tfmjd(i)
 40	continue
 50	continue
1040	format(a1,i5,1x,a,f9.3,f20.13)
 
	nt=i
	return
	end
