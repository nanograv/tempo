c      $Id$

	logical gro,sim,xitoa,oldpar,psrframe,jumpout,eclcoord
	logical usestart, usefinish
        logical usedmx
        logical npulsein, npulseout, infoout
	integer parunit, nskip, iboot
        integer infolen
        integer ssdmflag
        real*8 phimin

	common pdec,pra,ba(3),bc(3),dm,dt,dt2,freq(NPAP1),
     +    ferr(NPAP1),fmin,hlt(36),hrd(36),wt,x(NPAP1),era,ec,
     +    erd,fmax,emax,tmax,phimin,start,finish,amjd1,amjd2,posepoch,
     +    posep,dither,xjdoff(2,NJUMP),dct(NJUMP),
     +    dmx(NDMXMAX),dmxr1(NDMXMAX),dmxr2(NDMXMAX),dmxt,ndmx,usedmx,
     +    pmra,pmdec,pmrv,dt2sec,
     +    t0geo,nfit(NPAP1),mfit(NPAP1),n,nscan,nparam,nxoff,nprnt,
     +    nkeep,nfq,ncoord,gro,sim,xitoa,oldpar,psrframe,jumpout,
     +	  eclcoord,usestart,usefinish,npulsein,npulseout,
     +    parunit,nskip,iboot,ndmcalc,
     +    nfcalc,gain,tres,ntoa,nparam0,infolen,infoout,ssdmflag



	character psrname*12,obsflag*1,pardir*80,infotxt*160
        character obskey*5

        common/acomch/psrname,pardir,obsflag,infotxt,obskey(36) 


        real*8 array(NPA,NPA) ! moved here from fit.f, djn, 8-Sep-98

        common/acomcov/array
