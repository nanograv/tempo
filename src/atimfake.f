c      $Id$
C	@(#)atimfake.f	9.19 2/3/94

	subroutine atimfake(oldpar,afmjd,nbin,nt,ncoord,alng,ipsr)
	implicit real*8 (a-h,o-z)
	parameter (n=31,pi=3.141592653589793d0)
	real*8 x(n)
	logical oldpar
        real*8 maxha

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
           nx = ntzrmjd/10
           fx = (ntzrmjd-10*nx)+ftzrmjd
	   write(50,1040)tzrsite,i,pname(1:8),tzrfrq,nx,fx
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
	fmjd1=afmjd+(wait-maxha)/24.+nspan/2880.	!Compute start time
        nmjd1 = int(fmjd1)      ! break into integer, fractional parts
        fmjd1 = fmjd1 - nmjd1
	fmjd1=nint(48*fmjd1)/48.d0 ! round to nearest half-hour?!
	nsets=(120*maxha+nspan-1)/nspan
 
	b=nspan/2 + 5
	a=-b
	bma=0.5*(b-a)
	bpa=0.5*(b+a)
	do 30 k=1,n
          x(k)=cos(pi*(k-0.5)/n)*bma+bpa
 30     continue
 
	i=0
	do 50 j=1,nsets				
          fmjd2 = fmjd1 + (j-1)*nspan/1440.d0
          ntmp1 = int(fmjd2*1.d8) ! round fmjd2 to 1.e-10 (messy since ints
          ntmp2 = int((fmjd2 *1.d8-ntmp1)*1.d2) ! don't have 10 significant
          fmjd2 = ntmp2*1.d-10 + ntmp1*1.d-8 !    digits)
          do 40 k=1,n
            i = i + 1
            if(i.gt.800) stop ' Nspan too small'
            ftmjd(i) = fmjd2 + x(k)/1440.d0
            ntmjd(i) = nmjd1
            if (ftmjd(i).ge.2.) then
              ftmjd(i) = ftmjd(i) - 2.
              ntmjd(i) = ntmjd(i) + 2
            elseif (ftmjd(i).ge.1.) then
              ftmjd(i) = ftmjd(i) - 1.
              ntmjd(i) = ntmjd(i) + 1
            elseif (ftmjd(i).lt.0.) then
              ftmjd(i) = ftmjd(i) + 1.
              ntmjd(i) = ntmjd(i) - 1
            endif
            tmin(i) = 1440.d0*((ntmjd(i)-ntmjd(1))+ftmjd(i)-ftmjd(1))
            nx = ntmjd(i)/10
            fx = (ntmjd(i)-10*nx)+ftmjd(i)
            write(50,1040) tzsite,i,pname(1:8),tzfreq,nx,fx
 40       continue
 50	continue
 1040   format(a1,i5,1x,a,f9.3,i4,f16.14)


	nt=i
	return
	end
