c      $Id$
      subroutine arrtim(mode,xmean,sum,sumwt,dnpls,ddmch,ct2,
     +     alng,nsmax,nz,tz,nptsmax,nits,jits,
     +     buf,npmsav,ksav,nbuf,memerr,infofile,lw)

C  Reads pulse arrival times (at observatory, in UTC) from unit 4. Then
C  computes equivalent coordinate times at solar system barycenter,
C  and computes residuals from the assumed set of parameters.

C  RNM 27-Feb-91  Allow nsite=0 for geocentric data
C  DJN 18-Aug-92  Allow up to 36 sites

	implicit real*8 (a-h,o-z)
        save
	include 'dim.h'
	real*8 xmean(NPA),fctn(NPAP1),dnpls(*),ddmch(*),alng(36)
        real*8 buf(*)
        integer npmsav(*), ksav(*)
	real*4 gasdev
        integer ZEROWT(2000),idum
	integer ilen
        integer wflag
	character*80 card,card2,infile
	character asite*1,bsite*2,comment*8,aterr*9,afmjd*15
	logical first,offset,jdcatc,last,dithmsg,tz,track,search
	logical memerr
        character*160 infofile
	include 'acom.h'
	include 'bcom.h'
	include 'dp.h'
	include 'trnsfr.h'
	include 'orbit.h'
	include 'glitch.h'
	include 'tz.h'

	common /CONST/ PI,TWOPI,SECDAY,CONVD,CONVS,AULTSC,VELC,EMRAT,OBLQ,
     +              GAUSS,RSCHW,AULTVL
	common/crds/ rcb(6),rbe(6),rce(6),rse(3),rca(3),rsa(3),rba(3),
     +    rea(3),psid,epsd,pc,ps,tsid,prn(3,3),atut,ut1ut,etat,dtgr,
     +    tdis,bclt
	common/obsp/ site(3),pos(3),freqhz,bval,sitvel(3)
	data first/.true./,card2/'JUMP'/,idum/-1/

	offset=.true.
	jdcatc=.false.
	do 1 j=1,NPA
1	xmean(j)=0.

	rewind 50
        if (npulsein.or.npulseout) rewind 35
	do i=1,nskip
	  read(50,1002)
 1002	  format(a1)
	enddo

	dphase=0.d0
	deltat=0.d0
        pha1 = 0.d0
        pha2 = 0.d0
	equad  = 0.d0
	emin = 0.d0
	efac = 1.d0
	sigm = 0.d0
	amjd1=1000000.0
	amjd2=0.
	last=.false.
	dithmsg=.true.
	track=.false.
	search=.false.
	ct00=-1.d30
	fmin=0.
	fmax=1.d30
	emax=1.d30
	nblk=0
	ntrk=0
	dither=0.
        zawgt=0.
	mode=0
	nxoff=0
	n=0
	sum=0.
	sumsq=0.
	sumwt=0.
	modscrn=10
        nw=1
	if(first)then
	   call ztim(0,0.d0,nct,fct,wflag)
	   write (31,1011)
 1011	   format(/'    N     TOBS   WORD    NJMP    WGT',
     +     '     DTIME       TIME     DPHASE      PHASE'/)
	   nz=0
	   ntzref=0
	   first = .false.
	endif
	lu = 50

c       convert from mas/year (with cos(dec) term in) to radian/Julian century
c       in ra (without cos(dec) term), dec and rv.
	dpa = convd/3600.d0*pmra/10.d0/dcos(pdec)
	dpd = convd/3600.d0*pmdec/10.d0
	dpr = convd/3600.d0*pmrv/10.d0

c       parallax in radians from parallax in mas
	dpara = px*convd/3600.d3

	if(psrframe)then
C       A. Irwin stuff requires positions at epoch J2000. Compute space motion vectors 
C       and J2000 epoch position (Note all of this is in J2000 coordinates).
	   delt = (51545.d0 - posep)/36525.d0
	   pra2000 = pra + dpa*delt
	   pdec2000 = pdec + dpd*delt
c       radial velocity in au/Julian century
	   if(dpara.gt.0.d0)then
	      dvel = dpr/dpara
	   else
	      dvel=0.d0
	   endif
	   call space_motion(pra2000, pdec2000, dpa, dpd, dpara, dvel, 
     +   	aa_q, aa_m, 1)
c       aa_q/dpara is position vector in au at epoch J2000
c       aa_m/dpara is space velocity vector in au/Jcentury

	else
c       N. Wex transformations at pepoch
	   call pospm_conv(pra, pdec, dpa, dpd, aa_q, aa_m)
	endif

	comment='        '

C       The main loop starts here!

   10  	continue
   	read(lu,1010,end=40) card
1010	format(a80)
	if(card(1:1).lt.'A'.and.card(1:1).ne.'#') go to 50
	if((card(1:1).ge.'a').and.(card(1:1).le.'z')) go to 50
	if(card(1:4).eq.'END ') go to 45
	if(card(1:7).eq.'INCLUDE') then
	  lu=lu+1
	  j1 = 8
	  call citem(card,78,j1,infile,ilen)
	  write(31,1012) card(1:78)
1012	  format(1x,a78)
	  open(lu,file=infile(1:ilen),status='old',err=1013)
	  rewind lu
	else
	  call pcard(card,mode,zawgt,deltat,fmjd,dphase,sigm,offset,
     +     jdcatc,pha1,pha2,efac,emin,equad,jits,lu,track,trkmax,search,
     +         lw)
	endif
	go to 10

 1013	write(*,'(''Failed to open INCLUDE file: '',a)')infile
	STOP

50	continue
	if(card(1:35).eq.'                                   ') goto 10

c  default: Princeton format
	nfmt=0
c  first character blank: Parkes format
	if(card(1:1).eq.' ') nfmt=1
c  first character not blank, second character not blank:
c    could be either ITOA format, in which case first field is PSR name
c    or Princeton format, but with "scan number" starting at column 2
c    (Berkeley sometimes does this).  Differentiate between the 
c    two cases by searching for the "+" or "-" indicative of a pulsar name.
	if(card(1:1).ne.' ' .and. card(2:2).ne.' ') then
	  do i = 1, 80
	    if(card(i:i).eq." ") then
	      nfmt=0
	      goto 51
	    elseif (card(i:i).eq."+") then
	      nfmt=2
	      goto 51
	    elseif (card(i:i).eq."-") then
	      nfmt=2
	      goto 51
	    endif
	  enddo
	  nfmt = 2 ! shouldn't get here....but just in case
	endif
 51	continue

	if(nfmt.eq.0) then				! Princeton format
	  if(card(30:30).eq.'.') read(card,10500) asite,rfrq,
     +      nfmjd,ffmjd,aterr,ddm
10500	  format(a1,14x,f9.0,i5,f15.14,a9,15x,f10.0)
	  if(card(31:31).eq.'.') read(card,10501) asite,rfrq,
     +      nfmjd,ffmjd,aterr,ddm
10501	  format(a1,14x,f9.0,i6,f14.13,a9,15x,f10.0)
	  if(nfmjd.lt.30000)nfmjd=nfmjd+39126           ! Convert 1966 days to MJD

	  iz = 1
	  do i=1,9
 	    if(aterr(i:i).eq.' ') iz=i
	  end do
	  terr=0.
	  if (iz.le.9) read(aterr(iz:9),*,err=54,end=54) terr
54	  continue

	else if(nfmt.eq.1) then				! Parkes format
	  read(card,1051) rfrq,nfmjd,ffmjd,phs,terr,asite
1051	  format(25x,f9.0,i7,f14.13,f8.6,f8.1,8x,a1)
	  ffmjd=ffmjd+phs*p0/86400.d0
	
	else if(nfmt.eq.2) then				! ITOAF format
	  read(card,1052) nfmjd,ffmjd,terr,rfrq,ddm,bsite,comment
1052	  format(9x,i5,f14.13,f6.2,f11.4,f10.6,2x,a2,2x,a8)
	  asite=' '
	  do iobs=1,36
	     if(bsite.eq.obskey(iobs)(4:5))then
		asite=obskey(iobs)(1:1)
		goto 56
	     endif
	  enddo
	endif

 56	if(ffmjd.gt.1.d0) then
	  nfmjd=nfmjd+1
	  ffmjd=ffmjd-1.d0
	endif
	nsite = 99
	if (asite.ge.'0' .and. asite.le.'9') nsite = ichar(asite)-48
	if (asite.ge.'a' .and. asite.le.'z') nsite = ichar(asite)-87
	if (asite.eq.'@') nsite = -1

	fmjd=nfmjd+ffmjd
	if(rfrq.lt.fmin .or. rfrq.gt.fmax .or. 
     +      (terr .gt. emax) .or.
     +      (usestart .and. fmjd.lt.start) .or. 
     +      (usefinish .and. fmjd.gt.finish)) go to 10

C Arrival time
	n=n+1

        if (n.gt.nptsmax) then
	  if (tz) stop ' Too many points for tempo -z'
          memerr = .true.  !too many points.  continue through the file
	  goto 10          !to figure out how many points total we will have
        endif

C Store toa for tz reference phase
	if(ntzref.eq.0.and.fmjd.gt.pepoch)then
	   ntzref=n
	   ntzrmjd=nfmjd
	   ftzrmjd=ffmjd
	   tzrsite=asite
	   tzrfrq=rfrq
	endif

C  Get clock corrections
	if(nfmt.eq.2 .or. nsite.le.0)then ! no correction for...
	   clk2=deltat/86400.d0             !     ...ITOA, sites 0 and @
	else
	   call clockcor(fmjd,nsite,n,deltat,clk2)
	endif

	if(dither.ne.0.d0 .and. (.not.sim)) then
	  clk2=clk2+dither*gasdev(idum)/86400.d6
	  if(nits.ge.2.and.dithmsg) then
	    write(*,*) 'Iterating a solution with nonzero ',
     +        'DITHER is generally a very bad idea.'
	    dithmsg=.false.
	  endif
	endif
	ffmjd=ffmjd+clk2
	fmjd=fmjd+clk2
	go to 60  ! skip over end-of-file bits

c   Get here on end-of-file
40	if(lu.gt.50)then
	   close(lu)
	   lu=lu-1
	   go to 10
	endif

c   Get here on end-of-all-files
45	if(.not.gro) go to 100
	last=.true.
	nsite=0
	rfrq=0.
	t0geo=nint((amjd1+amjd2)/2.)
	fmjd=t0geo
	nfmjd=fmjd
        ffmjd=0.d0

c   Back to processing of all TOAs
 60	yrs=(fmjd-pepoch)/365.25d0

	if(nsite.gt.nsmax) write(*,1063) n,nsite,nsmax,nfmt,card
1063	format(i5,' *** Site',i3,' greater than',i3,
     +    ' read in. ***    nfmt:',i2/1x,a80)
	if(nsite.gt.0)then
	  site(1)=hrd(nsite)*dcos(hlt(nsite))*499.004786d0
	  site(2)=site(1)*dtan(hlt(nsite))
	  site(3)=alng(nsite)
	else
	  site(1)=0.d0
	  site(2)=0.d0
	  site(3)=0.d0
	endif

	freqhz=rfrq*1.0d6	!NB: rfrq=0 interpreted as infinite frequency
	ddmch(n)=ddm		!Save ddm
	if(nddm.eq.0) ddm=0	!Disable DM corrections if not wanted
	dmtot=dm+ddm
	if(ndmcalc.ge.2) then   
          fac = 1
	  do 61 i=1,ndmcalc-1
            fac = fac * yrs / real(i)
	    dmtot=dmtot + fac*dmcof(i) 
61          continue
	endif
        if (usedmx) then
          do i = 1, ndmx
            if (nfmjd+ffmjd.ge.min(dmxr1(i),dmxr2(i)-dmxt)
     +           .and.nfmjd+ffmjd.le.max(dmxr2(i),dmxr1(i)+dmxt)) then
              dmxr1(i)=min(dmxr1(i),nfmjd+ffmjd)
              dmxr2(i)=max(dmxr2(i),nfmjd+ffmjd)
              idmx = i
              goto 80
            endif
          end do
          if (ndmx+1.ge.NDMXMAX) then
            print *,"Error at TOA number ",n
            print *,"Too many dm offsets.  Maximum: NDMXMAX=",NDMXMAX
            print *,"  To change this, edit dim.h & recompile tempo"
            stop
          endif 
          ndmx = ndmx + 1
          dmxr1(ndmx) = nfmjd+ffmjd
          dmxr2(ndmx) = nfmjd+ffmjd
          idmx = ndmx
          nfit(NPAR6+ndmx)=1
          nparam=nparam+1
          mfit(nparam)=NPAR6+ndmx
 80       continue
        endif

        if (usedmx) then
          dmtot = dmtot + dmx(idmx)
        endif

	bval=dmtot/2.41d-16
	  
        wflag = 1
        if (nsite.ge.0) then
          call ztim(nfmjd,ffmjd,nct,fct,wflag)
        elseif (freqhz.lt.1) then ! barycenter, infinite frequency
          nct = nfmjd
          fct = ffmjd
          frq = freqhz * 1.e-6
        else                    ! barycenter, non-infinite frequency
          nct = nfmjd
          fct = ffmjd - bval/freqhz**2/SECDAY
          if (fct.lt.0.) then
            fct = fct + 1
            nct = nct - 1
          elseif (fct.ge.1.) then
            fct = fct - 1
            nct = nct + 1
          endif
          frq = freqhz * 1.e-6
        endif

	if(search.and.ct00.gt.0.d0.and.(abs(nct+fct-ct00).gt.trkmax)) 
     +       then
	  nblk=nblk+1
	  call pcard(card2,mode,zawgt,deltat,fmjd,dphase,sigm,offset,
     +     jdcatc,pha1,pha2,efac,emin,equad,jits,lu,track,trkmax,search)
	  if(nblk.ge.2) call pcard(card2,mode,zawgt,deltat,fmjd,
     +      dphase,sigm,offset,jdcatc,
     +      pha1,pha2,efac,emin,equad,jits,lu,track,trkmax,search)
	endif

	if(jdcatc) xjdoff(1,nxoff)=fmjd
	jdcatc=.false.
	if(nxoff.gt.0.and.x(NPAR2+nxoff).ne.0.d0)
     +    fct=fct+dct(nxoff)/86400.d0
	ct=nct+fct
	cp=p0+p1*(ct-pepoch)*86400.d0
	if (nsite.ge.0) then
	  call earth(cp,x(5),x(6),be,era,edc,pra,pdec,erd)
	endif
	if(ct.gt.ctmax) ctmax=ct
	wgt=1.0
	if(mode.eq.1.and.sigm.eq.0.d0) then
	  terr = sqrt(terr*terr+equad*equad)
	  terr=max(emin,efac*terr)
	  wgt=(1.d-6*terr/cp)**(-2)
	endif

	if(sigm.gt.0.d0) then
            terr=sigm
C  The following (added by JHT in version 7.03) is for PSR 1913+16 at 
C  Arecibo only.  NB: HA is abs(hour angle) in days
            if (zawgt .ne. 0.d0) then
              xlst = mod(tsid/6.283-0.1854+5.d0, 1.d0)
              ha = abs(xlst-0.8021)
              if (ha .gt. 0.0346) terr=sigm*(1.+1.8*(ha-0.0346)/0.0314)
            endif
          wgt=(1.d-6*terr/cp)**(-2)
        endif

        if (jits .eq. 0) then
            if (a1(1).ne.0.d0) then ! zero weight in case of orbital phase cut
              phi = dmod((ct-t0(1))*86400./pb(1) + 9000.,1.d0)
              if ((pha2.gt.pha1 .and. phi.gt.pha1 .and. phi.lt.pha2).or.
     +           (pha2.lt.pha1 .and.(phi.lt.pha2 .or. phi.gt.pha1)))then
                wgt = 0
                nz = nz + 1
                zerowt(nz) = n
              endif
            endif
            if (wflag.eq.0.and.wgt.ne.0) then ! zero weight in case of
                wgt = 0                       ! sun angle (phi) cut
                nz = nz + 1
                zerowt(nz) = n
              endif
          else 
            if (nw .le. nz) then
              if (n. eq. zerowt(nw)) then
                wgt = 0
                nw = nw + 1
              endif
            endif
        endif

	if(last) wgt=1.d-10*wgt

	wt=wgt

	if(xitoa) then				!Write the ITOA file
C  Write itoa file correctly, including observatory code.  (VMK, June94)
	  write(afmjd,1079) ffmjd
1079	  format(f15.13)
	  bsite='??'
	  do iobs=1,36
	     if (asite.eq.obskey(iobs)(1:1)) then
		bsite=obskey(iobs)(4:5)
		goto 70
	     endif
	  enddo

 70	  write(45,1080) psrname,nfmjd,afmjd(2:),terr,rfrq,
     +      ddmch(n),bsite,comment
1080	  format(a9,i5,a14,f6.2,f11.4,f10.6,2x,a2,2x,a8)
	endif

	call resid(nct,fct,dn,dphase,dnpls(n),nits,jits)

	if(track.and.n.gt.1) then
	  dt=dt+ntrk
 	  if (abs(ct-ct00).lt.trkmax) then
	    if(abs(dt+1.d0-dt00).lt.abs(dt-dt00)) then
	      dt=dt+1.d0
	      ntrk=ntrk+1
	    else if(abs(dt-1.d0-dt00).lt.abs(dt-dt00)) then
	      dt=dt-1.d0
	      ntrk=ntrk-1
	    endif
	  endif
	endif
	dt00=dt
	ct00=ct

c   write out tracking-corrected pulse number
	if (npulseout.and.jits.eq.0) then
	  write(35,fmt='(f14.0)') dn-ntrk
	endif

c   Write out "INFO" lines to info.tmp.  Write a line if:
c     -- this file has already been opened, in which case
c        there must be a line for every TOA, even if it
c        has zero length
c     -- this file has not been opened yet, but the current
c        TOA has a comment;  in this case, open the file
c        and write out lines of zero length for any previous
c        TOAs 
        if (infolen.ne.0 .or. infoout) then
          if (.not.infoout) then ! need to set up file,
            open(36,file=infofile)
            do i = 1, n-1
              write (36,fmt='()') ! write n-1 blank lines to the file
            end do        
            infoout = .true.
          endif
          write (36,fmt='(a)') infotxt(1:infolen)
        endif


C  Take a shortcut when called by TZ:
	if(tz) then
	  if(n.eq.2) ct2=ct
	  dnpls(n)=dt+dn
	  go to 10
	endif

C  DM-related partial derivatives

	x(16)=0.d0
        if (usedmx) then
          do i = 1, NDMXMAX
            x(NPAR6+i) = 0
          end do
        endif

	if(frq.gt.1.d0) then
          x(16)=f0*1.0d4/(2.41d0*frq**2)
          if (usedmx) x(NPAR6+idmx) = x(16)
        endif

	if(nfit(16).ge.2) then
          fac = 1.
	  do 89 i=1,nfit(16)-1
            fac = fac * yrs / real(i)
	    x(NPAR7+i)= fac * x(16) 
89        continue
	endif
	x(17)=f0*dtdpx
	x(36)=f0*dtdpmrv
	x(19)=f0*dtdppng
	do 90 j=1,nparam-1
90	fctn(j)=x(mfit(j+1))
	do 92 j=nparam,NPAP1
92	fctn(j)=0.
	do 93 j=1,nparam-1
93	xmean(j)=xmean(j)+wgt*fctn(j)
        sum=sum+wgt*dt
	sumsq=sumsq+wgt*dt*dt
	sumwt=sumwt+wgt
	call vmemw(n,fctn,ct,dt,wgt,dn,terr,frq,fmjd,rfrq,nparam-1,
     +     buf,npmsav,ksav,nbuf,memerr)

	nitsz=max0(nits,1)
	if(mod(n,20*modscrn).eq.1) write(*,1099) max0(nitsz,1)
1099	format(/'    N       MJD       Residual (p)  Residual (us)',
     +    '  Iteration (of',i2,')'/)
	if(n.gt.200) modscrn=100
	amjd1=min(amjd1,fmjd)
	amjd2=max(amjd2,fmjd)
	if(mod(n,modscrn).eq.1) write(*,1100)n,fmjd,dt,1d6*dt*p0,
     +    jits+1
1100	format(i6,f15.8,f11.6,f15.3,11x,i2)

        fmjdlast = fmjd

	if(.not.last) go to 10

c End of input file detected
100	continue

        if(mod(n,modscrn).ne.1) 
     +    write(*,1100)n,fmjdlast,dt,1d6*dt*p0,jits+1

	start=amjd1-1.d-3
	finish=amjd2+1.d-3

	if(.not.tz) then
	  sigma1=sqrt((sumsq-sum*sum/sumwt)/sumwt)*p0*1000.d0
	  m=-1
	  call vmemw(m,fctn,ct,dt,wgt,dn,terr,frq,fmjd,rfrq,nparam-1,
     +     buf,npmsav,ksav,nbuf,memerr)

	endif

        if (infoout) then
          close (36)
          infoout = .false.
        endif
9990	continue
	return
	end

