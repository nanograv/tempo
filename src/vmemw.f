c      $Id$
	subroutine vmemw(n,fctn,ct,dt,wgt,dn,terr,frq,fmjd,rfrq,npm)

C  Virtual memory handler.  Uses resid1.tmp only if "-r" was specified.

	implicit real*8 (a-h,o-z)
	include 'dim.h'
	include 'vcom.h'
	real*8 fctn(NPAP1)
	integer*2 npmsav
	common/vmem/ buf(NBUF),npmsav(NPTSMAX),ksav(NPTSMAX)
	save k

	if(n.eq.-1) return
	if(n.eq.1) k=0			!First TOA, set buf pointer
	ksav(n)=k
	npmsav(n)=npm			!Save value of npm
	do 10 i=1,npm			!Put data into buf
10	buf(i+k)=fctn(i)
	k=k+npm				!Update pointer
	buf(k+1)=ct
	buf(k+2)=dt
	buf(k+3)=wgt
	buf(k+4)=dn
	buf(k+5)=terr
	buf(k+6)=frq
	buf(k+7)=fmjd
	buf(k+8)=rfrq
	k=k+8				!Update pointer
	if(k.ge.NBUF-128) stop 'VMEM buffer overflow'	!Abort
	if(lresid1) then
	  if(n.eq.1) rewind 30
	  write(30) npm,(fctn(i),i=1,npm),ct,dt,wgt,dn,terr,frq,fmjd,rfrq
	  return
	endif
	end

	subroutine vmemr(n,fctn,ct,dt,wgt,dn,terr,frq,fmjd,rfrq,npmx)
	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*8 fctn(NPAP1)
	integer*2 npmsav
	common/vmem/ buf(NBUF),npmsav(NPTSMAX),ksav(NPTSMAX)

C###	if(n.eq.1) k=0			!Set pointer
	k=ksav(n)
	npm=npmx			!Silence Lahey warning
	npm=npmsav(n)
	do 20 i=1,npm
20	fctn(i)=buf(i+k)
	do 30 i=npm+1,NPAP1		!Zero the unused derivs
30	fctn(i)=0.d0
	k=k+npm				!Update pointer
	ct=buf(k+1)
	dt=buf(k+2)
	wgt=buf(k+3)
	dn=buf(k+4)
	terr=buf(k+5)
	frq=buf(k+6)
	fmjd=buf(k+7)
	rfrq=buf(k+8)
C###	k=k+8				!Update pointer

	return
	end
