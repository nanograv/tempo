c      $Id$

	logical gro,sim,xitoa,oldpar,psrframe,jumpout,eclcoord
	logical usestart, usefinish
        logical usedmx
        logical npulsein, npulseout, infoout, quiet, polystdout
        logical autotz
        logical useannorb   ! flag to use annual-orbital parallax
        logical usefixeddist! flag to use fixed input distance in ann-orb px
	integer parunit, nskip, iboot
        integer infolen
        integer ssdmflag
        real*8 phimin
        real*8 PAAscNode    ! position angle of ascending node
        real*8 solarn0      ! solar wind electron density at 1 AU (e-/cc)
        real*8 fixeddist    ! fixed input distance (instead of fit px) 
                            !    for ann-orb px

	common pdec,pra,ba(3),bc(3),dm,dt,dt2,freq(NPAP1),
     +    ferr(NPAP1),fmin,hlt(36),hrd(36),wt,x(NPAP1),era,ec,
     +    erd,fmax,emax,tmax,phimin,start,finish,amjd1,amjd2,posepoch,
     +    posep,dither,xjdoff(2,NJUMP),dct(NJUMP),
     +    dmx(NDMXMAX),dmxr1(NDMXMAX),dmxr2(NDMXMAX),dmxt,ndmx,usedmx,
     +    pmra,pmdec,pmrv,dt2sec,
     +    t0geo,gain,tres,
     +    PAAscNode,solarn0,fixeddist,
     +    nfit(NPAP1),mfit(NPAP1),n,nscan,nparam,nxoff,nprnt,
     +    nkeep,nfq,ncoord,gro,sim,xitoa,oldpar,psrframe,jumpout,
     +	  eclcoord,usestart,usefinish,npulsein,npulseout,
     +    parunit,nskip,iboot,ndmcalc,
     +    nfcalc,ntoa,nparam0,infolen,infoout,ssdmflag,
     +    quiet,polystdout,autotz,
     +    useannorb,usefixeddist



	character psrname*12,obsflag*1,pardir*80,infotxt*160
        character obskey*5

        common/acomch/psrname,pardir,obsflag,infotxt,obskey(36) 


        real*8 array(NPA,NPA) ! moved here from fit.f, djn, 8-Sep-98

        common/acomcov/array
