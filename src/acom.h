c      $Id$
	logical gro,sim,xitoa,oldpar,psrframe
	integer parunit, nskip
	character psrname*12,obsflag*1,pardir*80,obskey*5
	real*8 jdmx,jst,jsp
	common pdec,pra,ba(3),bc(3),dm,dt,dt2,freq(NPAP1),
     +    ferr(NPAP1),fmin,hlt(36),hrd(36),jdmx,jsp,jst,wt,x(NPAP1),
     +    era,ec,erd,fmax,tmax,
     +    dither,xjdoff(2,NJUMP),dct(NJUMP),pmra,pmdec,pmrv,dt2sec,
     +    t0geo,nfit(NPAP1),mfit(NPAP1),n,nscan,nparam,nxoff,nprnt,
     +    nkeep,nfq,mjd1,mjd2,ncoord,gro,sim,xitoa,oldpar,psrframe,
     +	  parunit,nskip,ndmcalc
        common/acomch/psrname,pardir,obsflag,obskey(36) 












