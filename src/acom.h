c      $Id$

	logical gro,sim,xitoa,oldpar,psrframe,jumpout,eclcoord
	logical usestart, usefinish
        logical npulsein, npulseout
	integer parunit, nskip, iboot

	common pdec,pra,ba(3),bc(3),dm,dt,dt2,freq(NPAP1),
     +    ferr(NPAP1),fmin,hlt(36),hrd(36),wt,x(NPAP1),
     +    era,ec,erd,fmax,tmax,start,finish,amjd1,amjd2,posepoch,posep,
     +    dither,xjdoff(2,NJUMP),dct(NJUMP),pmra,pmdec,pmrv,dt2sec,
     +    t0geo,nfit(NPAP1),mfit(NPAP1),n,nscan,nparam,nxoff,nprnt,
     +    nkeep,nfq,ncoord,gro,sim,xitoa,oldpar,psrframe,jumpout,
     +	  eclcoord,usestart,usefinish,npulsein,npulseout,parunit,
     +    nskip,iboot,ndmcalc,
     +    nfcalc,gain,tres,ntoa


	character psrname*12,obsflag*1,pardir*80,obskey*5

        common/acomch/psrname,pardir,obsflag,obskey(36) 


        real*8 array(NPA,NPA) ! moved here from fit.f, djn, 8-Sep-98

        common/acomcov/array
