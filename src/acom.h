c      $Id$

	logical gro,sim,xitoa,oldpar,psrframe,jumpout,eclcoord
	logical usestart, usefinish
        logical usedmx, firstdmx, nonewdmx
        logical npulsein, npulseout, infoout, quiet, polystdout
        logical phisunout
        logical tz
        logical autotz
        logical useannorb   ! flag to use annual-orbital parallax
        logical usefixeddist! flag to use fixed input distance in ann-orb px
        logical jumpbarycenter ! flag to apply jumps in barycenter instead of
                               ! obs frame.  Jumps used to be done in
                               ! bary frame, but now (10-jul-04) default
                               ! is obs frame
        logical tdbif99fmt  ! true if tdb file is in format used
                            ! by alan irwin (reference if99)
        logical siteused      ! flags for which observing codes are used.
        logical nofitjump     ! blocks fitting for a given JUMP
        logical useglsfit     ! use generalized-least-squares fit
        logical usedmdata     ! use input DM measurements as "data"

	integer parunit, nskip, iboot
        integer infolen
        integer ssdmflag
        integer nflagjumps  ! number of tempo2-style flag-based jumps
        integer nflagefac   ! number of tempo2-style flag-based EFAC
        integer nflagequad  ! number of tempo2-style flag-based EQUAD
        integer nflagecorr  ! number of flag-based ecorr/jitter terms
        real*8 phimin
        real*8 PAAscNode    ! position angle of ascending node
        real*8 solarn0      ! solar wind electron density at 1 AU (e-/cc)
        real*8 solarn01     ! solar wind electron density time derivative
        real*8 fixeddist    ! fixed input distance (instead of fit px) 
                            !    for ann-orb px
        real*8 dmepoch, dmep

        real*8 ceecl, seecl ! cosine and sine of obliquity of the ecliptic

	common pdec,pra,ba(3),bc(3),dm,dt,dt2,freq(NPAP1),
     +    ferr(NPAP1),fmin,hlt(36),hrd(36),siteused(36),
     +    wt,x(NPAP1),era,ec,
     +    erd,fmax,emax,tmax,phimin,start,finish,amjd1,amjd2,posepoch,
     +    posep,dither,xjdoff(2,NJUMP),dct(NJUMP),nofitjump(NJUMP),
     +    dmepoch,dmep,
     +    dmx(NDMXMAX),dmxr1(NDMXMAX),dmxr2(NDMXMAX),dmxt,ndmx,usedmx,
     +    dmx1(NDMXMAX),dmxep(NDMXMAX),dmxf1(NDMXMAX),dmxf2(NDMXMAX),
     +    flagefac(NFLAGERR),flagequad(NFLAGERR),flagecorr(NFLAGERR),
     +    usedmx1,
     +    fdcof(NFDMAX),
     +    rnamp,rnidx,
     +    tcorr,
     +    pmra,pmdec,pmrv,dt2sec,
     +    t0geo,gain,tres,
     +    PAAscNode,solarn0,solarn01,fixeddist,
     +    ceecl, seecl,
     +    nfit(NPAP1),mfit(NPAP1),n,nscan,nparam,nxoff,nprnt,
     +    nkeep,nfq,ncoord,gro,sim,xitoa,oldpar,psrframe,jumpout,
     +	  eclcoord,usestart,usefinish,npulsein,npulseout,
     +    parunit,nskip,iboot,ndmcalc,nflagjumps,
     +    nflagefac,nflagequad,nflagecorr,
     +    nfcalc,ntoa,nparam0,infolen,infoout,phisunout,ssdmflag,
     +    quiet,polystdout,tz,autotz,firstdmx,nonewdmx,
     +    useannorb,usefixeddist,jumpbarycenter,useglsfit,
     +    usedmdata,
     +    tdbif99fmt



	character psrname*12,obsflag*1,pardir*80,infotxt*160
        character obskey*5
        character jumpflag*32, jumpflagval*32
        character infoflag*32
        character efacflag*32, efacflagval*32
        character equadflag*32, equadflagval*32
        character ecorrflag*32, ecorrflagval*32
        character eclcon*80
        character dcovfile*80

        common/acomch/psrname,pardir,obsflag,infotxt,obskey(36),
     +    eclcon,dcovfile,
     +    jumpflag(NJUMP),jumpflagval(NJUMP),infoflag,
     +    efacflag(NFLAGERR),efacflagval(NFLAGERR),
     +    equadflag(NFLAGERR),equadflagval(NFLAGERR),
     +    ecorrflag(NFLAGERR),ecorrflagval(NFLAGERR)


        real*8 array(NPA,NPA) ! moved here from fit.f, djn, 8-Sep-98

        common/acomcov/array
