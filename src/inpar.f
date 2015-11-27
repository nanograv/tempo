c      $Id$
      subroutine zeropar(nits)

      implicit real*8 (A-H,O-Z)

      include 'dim.h'     
      include 'acom.h'
      include 'bcom.h'
      include 'trnsfr.h'
      include 'clocks.h'
      include 'dp.h'
      include 'orbit.h'
      include 'eph.h'
      include 'glitch.h'


      integer i

      do i=1,NPAP1
         nfit(i)=0
         x(i)=0.
      enddo
	
      pmra=0.
      pmdec=0.
      pmrv=0.
      p0=0.
      p1=0.
      p2=0.
      pepoch=50000.
      posepoch=0.
      dmepoch=0.
      f0=0.
      f1=0.
      f2=0.
      f3=0.
      do i=1,9
         f4(i)=0.
      enddo
      px=0.
      dm=0.
      do i=1,NDMCOFMAX
         dmcof(i)=0.
      enddo
      start=0.
      finish=0.
      ntoa=0.
      tres=0.

      do i=1,4
         a1(i)=0.
         e(i)=0.
         t0(i)=0.
         pb(i)=0.
         omz(i)=0.
      enddo
      t0asc=0.
      eps1=0.
      eps2=0.
      omdot=0.
      gamma=0.
      pbdot=0.
      si=0.
      shapmax=0. !!!! NEW in DDS
      shaphof=0.  ! NW: higher order Shapiro in DDS - scale factor
      cotchi0=0.  ! NW: higher order Shapiro in DDS - latitudinal time delay (RL06)
      am=0.
      am2=0.
      dr=0.
      dth=0.
      a0=0.
      b0=0.
      bp=0.
      bpp=0.
      xdot=0.
      edot=0.
      afac=0.
      om2dot=0.
      x2dot=0.
      eps1dot=0.
      eps2dot=0.
      dmvar1 = 0.
      dmvar2 = 999999.

      rnamp = 0.
      rnidx = 0.
      tcorr = 0.
      dcovfile = ""

      ntzrmjd = 0.
      ftzrmjd = 0.
      tzrfrq = 0.
      tzrsite = 0.

      nxoff=0
      nflagjumps=0
      do i=1,NJUMP
         xjdoff(1,i)=0.
         xjdoff(2,i)=0.
         dct(i)=0.
         nofitjump(i) = .false.
      enddo

      infoflag = ""
      
      ngl=0
      do i=1,NGLT
         glph(i)=0.
         glf0(i)=0.
         glf1(i)=0.
         glf0d(i)=0.
         gltd(i)=0.
         glepoch(i)=0.
      enddo

      nfbj = 0
      do i = 1, NFBJMAX
        tfbj(i) = 0.
        fbj(i) = 0.
      enddo

      ndmx = 0
      do i = 1, NDMXMAX
        dmx(i) = 0.
        dmxep(i) = 0.
        dmx1(i) = 0.
	dmxf1(i) = 0.
	dmxf2(i) = 0.
      enddo

      ndmcalc=0

      do i=1, NFDMAX
        fdcof(i) = 0.0
      enddo

      nbin=0
      nplanets=0
      nclk=0
      nephem=1
      nits=1
      ncoord=1
      nell1=0
      nshapho=0 ! NW: higher order Shapiro in DDS - flag to switch on higher oder corrections

      usestart = .false.
      usefinish = .false.

      usedmx = .false.
      firstdmx = .true.
      usedmx1 = 0
      ndmx = 0
      dmxt = 0
      
      usedmdata = .false.
      dmefac = 1.0

      eclcoord = .false.

      iboot = 0

      useannorb = .false.
      usefixeddist = .false.

      solarn0 = 10      
      solarn01 = 0      

      return
      end

c=======================================================================

      subroutine rdpar(nits)

C  Free format input
C  Line structure: "[fit] key value error/comment"
C  The error/comment is ignored by TEMPO

      implicit real*8 (A-H,O-Z)

      include 'dim.h'
      include 'acom.h'
      include 'bcom.h'
      include 'trnsfr.h'
      include 'clocks.h'
      include 'tz.h'
      include 'dp.h'
      include 'orbit.h'
      include 'eph.h'
      include 'glitch.h'

      character line*80, key*32, value*64, cfit*8, temp*80, cifit

      logical seteps            ! indicate when eps1 and/or eps2
                                ! had been set
      logical setepsdot         ! indicate when eps1dot and/or eps2dot
                                ! had been set
      logical set2dot           ! indicate when om2dot and/or x2dot
                                ! had been set
      logical setecl, setequ    ! indicate when some ecliptic or
	                        ! equatorial coordinate has been set
      logical setpb, setfb      ! indicate when orbital period/pdot or
                                ! orbital frequency/fdot has been set

      integer ffit(NFMAX)  ! local flags for which frequency derivatives f1 through f12 to fit
                           !  -1: not set at all
                           !   0 or 1: set as a fit parmaeter from a .par card
                           ! special meaning of the fit flag of F1:
                           !   if it has the value xx ("F1 1 xx") then set fit(xx)=1
                           ! after all parmaeters have been read:
                           !   --final fit flags are nfit(3), nfit(4), nfit(NPAR11+1) through nfit(NPA)
                           !   --any derivative f1 through f12 that has not been explicitly
                           !     set on or off bad has a higher derivative being fit will also
                           !     be fit itself

      integer i
 
      do i=1, NFMAX
        ffit(i) = -1
      enddo


      gain = 1.d0
      phimin = 0.d0

      seteps    = .false.
      setepsdot = .false.
      set2dot   = .false.
      setecl    = .false.
      setequ    = .false.
      setpb     = .false.
      setfb     = .false.

      ll=80

      eclcon = "DEFAULT"

      ijump  = 0  ! used for counting tempo2-style jump params
      iefac  = 0  ! used for counting tempo2-style efac params
      iequad = 0  ! used for counting tempo2-style equad params
      iecorr = 0  ! used for counting ecorr params

      nskip = 0  ! counts parameter lines, which are skipped if TOAs are
		 !    read from this file
	
      rewind(parunit)  ! probably not needed, but a good safety check

 10   read(parunit,'(a)',end=900)line
      nskip = nskip + 1

C  Get key, value and cfit
      jn=1
      call citem(line,ll,jn,key,lk)      
      if(lk.gt.32)then
        write(*,'('' Key overflow (32 max): '',a)')key(1:lk)
        stop
      endif
      if(key(1:1).eq.'#' .or. (key(1:1).eq.'C' .and. lk.eq.1))go to 10
      call upcase(key)
      call citem(line,ll,jn,value,lv)
      if(lv.gt.64)then
         write(*,'('' Value overflow (64 max): '',a)')value(1:lv)
         stop
      endif
      call citem(line,ll,jn,temp,lf)
      if(temp(1:1).eq.'#'.or.lf.eq.0)then
         cfit='0'
      else 
         cfit=temp
      endif
      ! cfit should never be a floating point number.  however, sometimes
      ! this happens when a par file is produced by the atnf catalogue
      ! and errors are printed in column 3.  Here we set cfit to '0'
      ! if it is determined to be floating point (checked for by
      ! the precense of a period)
      if(index(cfit,'.').ne.0) cfit = '0'
      ikey = keyidx(key)        ! extract xx in keys of form ssss_xx

C  Control parameters

      if(key(1:4).eq.'NPRN')then
         read(value,*)nprnt

      else if(key(1:4).eq.'NITS')then
         read(value,*)nits

      else if(key(1:4).eq.'IBOO')then
         read(value,*)iboot

      else if(key(1:4).eq.'NDDM')then
         read(value,*)nddm

      else if(key(1:6).eq.'DMVAR1')then
         read(value,*)dmvar1

      else if(key(1:6).eq.'DMVAR2')then
         read(value,*)dmvar2

      else if(key(1:6).eq.'PHIMIN')then
         read(value,*)phimin

      else if(key(1:4).eq.'DMX'.and.lk.eq.3)then
         usedmx = .true.
         read(value,*)dmxt

      else if(key(1:7).eq.'DMXFIX')then
         read(value,*)itmp
         if(itmp.gt.0)firstdmx=.false.

      else if(key(1:5).eq.'DMX1'.and.lk.eq.4)then
	 read(value,*)usedmx1

      else if(key(1:6).eq.'DMDATA')then
         read(value,*)itmp
         if(itmp.gt.0)usedmdata=.true.

      else if(key(1:6).eq.'DMEFAC')then
         read(value,*)dmefac

      else if(key(1:4).eq.'COOR')then
         if(value(1:5).eq.'B1950')ncoord=0

      else if(key(1:3).eq.'CLK')then
         call upcase(value)
         do i=0,nclkmax
            if(value(1:5).eq.clklbl(i)(1:5))then
               nclk=i
               go to 12
            endif
         enddo
         write(*,'(''Invalid CLK label: '',a)')value(1:5)
         stop
 12      continue

      else if(key(1:5).eq.'UNITS')then
         call upcase(value)
         if(value(1:3).ne.'TDB')then
            write(*,'(''Invalid UNITS: '',a)')value(1:5)
            stop
         endif

      else if(key(1:4).eq.'EPHE')then
         call upcase(value)
         do i=1,kephem
            if(value(1:ephnamel(i)).eq.ephname(i)(1:ephnamel(i))) then
		nephem=i
            	go to 14
	    endif
         enddo
         write(*,'(''Invalid EPHEM file name: '',a)')value(1:5)
         stop
 14      continue
      else if(key(1:6).eq.'TZRMJD')then
         itmp = index(value,'.')
         if (itmp.eq.0) then
           read(value,*)ntzrmjd
           ftzrmjd = 0
         else
           read(value(1:itmp-1),*)ntzrmjd
           read(value(itmp:64),*)ftzrmjd
         endif
         if (ntzrmjd.lt.30000) ntzrmjd = ntzrmjd + 39126  ! convert fut to mjd

      else if(key(1:6).eq.'TZRFRQ')then
         read(value,*)tzrfrq

      else if(key(1:7).eq.'TZRSITE')then
         tzrsite=value(1:1)
         
      else if(key(1:5).eq.'START')then
         read(value,*)start
         read(cfit,*)itmp
         if (itmp.gt.0) usestart=.true.

      else if(key(1:6).eq.'FINISH')then
         read(value,*)finish
         read(cfit,*)itmp
         if (itmp.gt.0) usefinish=.true.

C  Period/Frequency parameters

      else if(key(1:2).eq.'P0' .or. (key(1:1).eq.'P' .and. lk.eq.1))then
         read(value,*)p0
         read(cfit,*)nfit(2)

      else if(key(1:2).eq.'P1'.or.key(1:4).eq.'PDOT')then
         read(value,*)p1
         read(cfit,*)nfit(3)

      else if((key(1:1).eq.'F'.and..not.(key(1:2).eq.'FB')
     +		.and..not.(key(1:2).eq.'FD'))
     +		.or.(key(1:3).eq.'FB '.or.key(1:3).eq.'FD ')) then

        if (lk.eq.1) then 
          ifit = 0
        else 
          cifit=key(lk:lk)
          if(cifit.ge.'A')then
            call upcase(cifit)
            ifit = ichar(cifit)-55
          else
            read(key(2:lk),*) ifit
          endif
        endif

        if (ifit.eq.0) then
          read(value,*)f0
          read(cfit,*)nfit(2)
        else if (ifit.eq.1) then
          read(value,*)f1
          if(cfit(1:1).ge.'A')then
            call upcase(cfit)
            jfit = ichar(cfit(1:1))-55
          else
C  	    print *,"reading jfit from cfit which is",cfit
            read(cfit,*) jfit 
          endif
          if (jfit.eq.0) ffit(1) = 0
          if (jfit.gt.0) then
            ffit(1) = 1
            ffit(jfit) = 1  ! triggers fit F2 through jfit-1 as well (see note at top)
          endif
        else if (ifit.eq.2) then
          read(value,*)f2
          read(cfit,*)ffit(2)
        else if (ifit.eq.3) then
          read(value,*)f3
          read(cfit,*)ffit(3)
        else
          if (ifit.lt.4.or.ifit.gt.NFMAX) then
            print *,"Error: illegal freq derivative in following line"
            print *,line
            print *,"Derivative must be between 0 and ",NFMAX
            stop
          endif
          read(value,*)f4(ifit-3)
          read(cfit,*)ffit(ifit)
        endif

      else if(key(1:4).eq.'PEPO')then
         read(value,*)pepoch

C  Position parameters

      else if(key(1:3).eq.'PSR')then
         psrname=value(1:12)

      else if(key(1:3).eq.'DEC')then
         call decolon(value)
         read(value,*)pdec
         read(cfit,*)nfit(5)
	 setequ = .true.

      else if(key(1:2).eq.'RA')then
         call decolon(value)
         read(value,*)pra
         read(cfit,*)nfit(6)
	 setequ = .true.

      else if(key(1:4).eq.'PMDE')then
         read(value,*)pmdec
         read(cfit,*)nfit(7)
	 setequ = .true.

      else if(key(1:4).eq.'PMRA')then
         read(value,*)pmra
         read(cfit,*)nfit(8)
         setequ = .true.

      else if(key(1:4).eq.'BETA')then
         read(value,*)pdec
         read(cfit,*)nfit(5)
	 setecl = .true.

      else if(key(1:6).eq.'LAMBDA')then
         read(value,*)pra
         read(cfit,*)nfit(6)
	 setecl = .true.

      else if(key(1:6).eq.'PMBETA')then
         read(value,*)pmdec
         read(cfit,*)nfit(7)
	 setecl = .true.

      else if(key(1:8).eq.'PMLAMBDA')then
         read(value,*)pmra
         read(cfit,*)nfit(8)
         setecl = .true.

      else if(key(1:4).eq.'PMRV')then
         read(value,*)pmrv
         read(cfit,*)nfit(36)

      else if(key(1:2).eq.'PX')then
         read(value,*)px
         read(cfit,*)nfit(17)

      else if(key(1:5).eq.'POSEP')then
         read(value,*)posepoch

      else if(key(1:7).eq.'DMEPOCH')then
         read(value,*)dmepoch

      else if(key(1:4).eq.'DMX_') then
         if (ikey.gt.NDMXMAX) then
           write(*,'(''DMX key too high: '',a)')key
           stop
         endif
         ndmx = max(ndmx,ikey)
         read(value,*)dmx(ikey)
         read(cfit,*)nfit(NPAR6+2*ikey-1)

      else if(key(1:5).eq.'DMX1_') then
         if (ikey.gt.NDMXMAX) then
           write(*,'(''DMX key too high: '',a)')key
           stop
         endif
         ndmx = max(ndmx,ikey)
         read(value,*)dmx1(ikey)
         read(cfit,*)nfit(NPAR6+2*ikey)
C IHS force fit of DMX value if DMX1 is to be fit
         if(((nfit(NPAR6+2*ikey)).gt.0)
     +		.and.((nfit(NPAR6+2*ikey-1)).eq.0))
     +		write(*,'(''Fitting DMX value anyway in bin '',a)')key
         if((nfit(NPAR6+2*ikey)).gt.0) nfit(NPAR6+2*ikey-1)=1
C          if((nfit(NPAR6+2*ikey)).gt.0) write(*,'(''Fitting anyway'')')

      else if(key(1:6).eq.'DMXEP_') then
         if (ikey.gt.NDMXMAX) then
           write(*,'(''DMX key too high: '',a)')key
           stop
         endif
         ndmx = max(ndmx,ikey)
         read(value,*)dmxep(ikey)

      else if(key(1:6).eq.'DMXR1_') then
         if (ikey.gt.NDMXMAX) then
           write(*,'(''DMX key too high: '',a)')key
           stop
         endif
         ndmx = max(ndmx,ikey)
         read(value,*)dmxr1(ikey)
         
      else if(key(1:6).eq.'DMXR2_') then
         if (ikey.gt.NDMXMAX) then
           write(*,'(''DMX key too high: '',a)')key
           stop
         endif
         ndmx = max(ndmx,ikey)
         read(value,*)dmxr2(ikey)
         
      else if((key(1:3).eq.'DM0'.and.lk.eq.3).or.
     +         key(1:2).eq.'DM'.and.lk.eq.2)then
         read(value,*)dm
         if(cfit(1:1).le.'9')then
            itmp=ichar(cfit(1:1))-48
         else
            call upcase(cfit)
            itmp=ichar(cfit(1:1))-55
         endif
         nfit(16)=max(nfit(16),itmp)
         ndmcalc = max(ndmcalc,itmp)

      else if(key(1:2).eq.'DM'.and.
     +        key(3:3).ge.'0'.and.key(3:3).le.'9') then 
        read(key(3:lk),*)jj
        if (jj.lt.1.or.jj.gt.NDMCOFMAX) then
          write (*,'(''DM derivative '',i3,'' too high.)'')') jj
          stop
        endif
        read(value,*)dmcof(jj)
        if (cfit(1:1).gt.'0') nfit(16)=max(nfit(16),jj+1)
        ndmcalc=max(ndmcalc,jj+1)

      else if(key(1:2).eq.'FD'.and.
     +        key(3:3).ge.'1'.and.key(3:3).le.'9') then
        read(key(3:lk),*)jj
	if (jj.lt.1.or.jj.gt.NFDMAX) then
	  write (*,'(''FD coeff '',i3,'' too high.)'')') jj
	  stop
	endif
        read(value,*)fdcof(jj)
        read(cfit,*)nfit(NPAR10+jj)

C  Binary parameters

      else if(key(1:4).eq.'BINA')then
         call upcase(value)
         do i=1,NMODELS
            if(value(1:8).eq.bmodel(i)) goto 20
         enddo
         write(*,100) value(1:8)
 100     format(' WARNING: binary model - ',a,' - not recognized')
         goto 22         
c 20      nbin=i-1  ! ### Check this !!! (Works in Linux/Intel)
 20      nbin=i
         if(value(1:2).eq.'BT'.and.value(4:4).eq.'P')
     +      read(value,'(2x,i1)') nplanets
 22      continue

c need to check for this tempo2 keyword here
      else if(key(1:14).eq.'PLANET_SHAPIRO') then
	 continue

c next two lines by sets on 29 Aug 05
      else if(key(1:4).eq.'PLAN')then
         read(value,*)nplanets

      else if(key(1:4).eq.'A1_1'.or.(key(1:2).eq.'A1'.and.lk.eq.2))then
         read(value,*)a1(1)
         read(cfit,*)nfit(9)

      else if(key(1:3).eq.'E_1 '.or.(key(1:1).eq.'E'.and.lk.eq.1)
     +        .or. (key(1:3).eq.'ECC'.and.lk.eq.3)) then
         read(value,*)e(1)
         read(cfit,*)nfit(10)

      else if(key(1:4).eq.'T0_1'.or.(key(1:2).eq.'T0'.and.lk.eq.2))then
         read(value,*)t0(1)
         read(cfit,*)nfit(11)

      else if(key(1:4).eq.'PB_1'.or.(key(1:2).eq.'PB'.and.lk.eq.2))then
         read(value,*)pb(1)
         read(cfit,*)nfit(12)
         setpb = .true.

      else if(key(1:4).eq.'OM_1'.or.(key(1:2).eq.'OM'.and.lk.eq.2))then
         read(value,*)omz(1)
         read(cfit,*)nfit(13)

      else if(key(1:4).eq.'A1_2')then
         read(value,*)a1(2)
         read(cfit,*)nfit(26)

      else if(key(1:3).eq.'E_2 ')then
         read(value,*)e(2)
         read(cfit,*)nfit(27)

      else if(key(1:4).eq.'T0_2')then
         read(value,*)t0(2)
         read(cfit,*)nfit(28)

      else if(key(1:4).eq.'PB_2')then
         read(value,*)pb(2)
         read(cfit,*)nfit(29)

      else if(key(1:4).eq.'OM_2')then
         read(value,*)omz(2)
         read(cfit,*)nfit(30)

      else if(key(1:4).eq.'A1_3')then
         read(value,*)a1(3)
         read(cfit,*)nfit(31)

      else if(key(1:3).eq.'E_3 ')then
         read(value,*)e(3)
         read(cfit,*)nfit(32)

      else if(key(1:4).eq.'T0_3')then
         read(value,*)t0(3)
         read(cfit,*)nfit(33)

      else if(key(1:4).eq.'PB_3')then
         read(value,*)pb(3)
         read(cfit,*)nfit(34)

      else if(key(1:4).eq.'OM_3')then
         read(value,*)omz(3)
         read(cfit,*)nfit(35)

      else if(key(1:5).eq.'OMDOT')then
         if (lk.eq.5 .or. key(6:6).eq.'1') then
           read(value,*)omdot
           read(cfit,*)nfit(14)
         else if (key(6:6).ge.'2'.and.key(6:6).le.'9') then
           read(key(6:6),*)jj
           read(value,*)omdot2(jj)
           read(cfit,*)nfit(NPAR8+jj-1)
         endif

      else if(key(1:5).eq.'GAMMA')then
         read(value,*)gamma
         read(cfit,*)nfit(15)

      else if(key(1:5).eq.'PBDOT')then
         read(value,*)pbdot
         read(cfit,*)nfit(18)
         setpb = .true.

	! note: order is important: 'FBJ' should be before 'FB'
      else if(key(1:3).eq.'FBJ'.and.ikey.ge.1.and.ikey.le.NFBJMAX)then
         if (ikey.gt.nfbj) nfbj=ikey
         read(value,*) fbj(ikey)
         read(cfit,*)nfit(NPAR5+2*ikey)
            
	! note: order is important: 'FBJ' should be before 'FB'
      else if(key(1:2).eq.'FB')then
         if (lk.eq.2) then   ! "FB" -- treat it as "FB0"
           setfb = .true.
           read(value,*)fb(1)
           read(cfit,*)nfit(NPAR3+1)
         else if (key(3:3).ge.'0'.and.key(3:3).le.'9') then
           setfb = .true.
           read(key(3:lk),*)jj
           if (jj.eq.0.or.jj.eq.1) setfb = .true.
           read(value,*)fb(jj+1)
           read(cfit,*)nfit(NPAR3+jj+1)
         endif
 
      else if(key(1:4).eq.'TFBJ'.and.ikey.ge.1.and.ikey.le.NFBJMAX)then
         if (ikey.gt.nfbj) nfbj=ikey
         read(value,*) tfbj(ikey)
         read(cfit,*)nfit(NPAR5+2*ikey-1)

      else if(key(1:5).eq.'PPNGA')then
         read(value,*)nfit(19)

      else if(key(1:2).eq.'SI')then
         read(value,*)si
         read(cfit,*)nfit(20)

      else if(key(1:4).eq.'MTOT')then
         read(value,*)am
         read(cfit,*)nfit(21)

      else if(key(1:2).eq.'M2')then
         read(value,*)am2
         read(cfit,*)nfit(22)

      else if(key(1:5).eq.'DTHET')then
         read(value,*)dth
         read(cfit,*)nfit(23)

      else if(key(1:4).eq.'XDOT')then
         if (lk.eq.4 .or. key(5:5).eq.'1') then
           read(value,*)xdot
           read(cfit,*)nfit(24)
         else if (key(5:5).ge.'2'.and.key(5:5).le.'9') then
           read(key(5:5),*)jj
           read(value,*)xdot2(jj)
           read(cfit,*)nfit(NPAR4+jj-1)
         endif

      else if(key(1:4).eq.'EDOT')then
         if (lk.eq.4 .or. key(5:5).eq.'1') then
           read(value,*)edot
           read(cfit,*)nfit(25)
         else if (key(5:5).ge.'2'.and.key(5:5).le.'9') then
           read(key(5:5),*)jj
           read(value,*)edot2(jj)
           read(cfit,*)nfit(NPAR7+jj-1)
         endif

      else if(key(1:6).eq.'XOMDOT')then
         read(value,*)xomdot
         read(cfit,*)nfit(37)

      else if(key(1:6).eq.'XPBDOT')then
         read(value,*)xpbdot
         read(cfit,*)nfit(38)

      else if(key(1:6).eq.'OM2DOT')then
         read(value,*)om2dot
         read(cfit,*)nfit(39)
         set2dot=.true.

      else if(key(1:5).eq.'X2DOT')then
         read(value,*)x2dot
         read(cfit,*)nfit(40)
         set2dot=.true.

      else if(key(1:4).eq.'EPS1'.and.lk.eq.4)then
         read(value,*)eps1
         read(cfit,*)nfit(10)
         seteps=.true.

      else if(key(1:4).eq.'EPS2'.and.lk.eq.4)then
         read(value,*)eps2
         read(cfit,*)nfit(13)
         seteps=.true.

      else if(key(1:5).eq.'TASC'.and.lk.eq.4)then
         read(value,*)t0asc
         read(cfit,*)nfit(11)         
         
      else if(key(1:7).eq.'EPS1DOT'.and.lk.eq.7)then
         read(value,*)eps1dot
         read(cfit,*)nfit(39)
         setepsdot=.true.

      else if(key(1:7).eq.'EPS2DOT'.and.lk.eq.7)then
         read(value,*)eps2dot
         read(cfit,*)nfit(40)
         setepsdot=.true.

      else if(key(1:7).eq.'SHAPMAX')then
         read(value,*) shapmax
         read(cfit,*)nfit(20)

      else if(key(1:7).eq.'SHAPHOF')then ! NW: higher order Shapiro - scale factor
         read(value,*) shaphof
         read(cfit,*)nfit(39)
         nshapho=1

      else if(key(1:7).eq.'COTCHI0')then ! NW: higher order Shapiro - latitudinal bending delay
         read(value,*) cotchi0
         nshapho=1


C  Fixed binary parameters

      else if(key(1:9).eq.'PAASCNODE'.and.lk.eq.9) then !pos'n ang of ascending node
         read(value,*)PAAscNode
         useannorb = .true.

      else if(key(1:4).eq.'DIST'.and.lk.eq.4) then !fixed dist(kpc) for ann-orb px
         read(value,*)fixeddist
         usefixeddist = .true.


      else if(key(1:2).eq.'DR')then
         read(value,*)dr

      else if(key(1:2).eq.'A0')then
         read(value,*)a0

      else if(key(1:2).eq.'B0')then
         read(value,*)b0

      else if(key(1:2).eq.'BP')then
         read(value,*)bp

      else if(key(1:3).eq.'BPP')then
         read(value,*)bpp

      else if(key(1:4).eq.'AFAC')then
         read(value,*)afac

C  Glitches


      else if(key(1:4).eq.'GLEP'.and.ikey.ge.1.and.ikey.le.NGLT) then
        if (ikey.gt.ngl) ngl=ikey
        read(value,*) glepoch(ikey)
        
      else if(key(1:4).eq.'GLPH'.and.ikey.ge.1.and.ikey.le.NGLT) then
        if (ikey.gt.ngl) ngl=ikey
        read(value,*) glph(ikey)
        read(cfit,*) nfit(NPAR1+NGLP*(ikey-1)+1)

       else if(key(1:4).eq.'GLF1'.and.ikey.ge.1.and.ikey.le.NGLT) then
         if (ikey.gt.ngl) ngl=ikey
         read(value,*) glf1(ikey)
         read(cfit,*) nfit(NPAR1+NGLP*(ikey-1)+3)

       else if(key(1:5).eq.'GLF0D'.and.ikey.ge.1.and.ikey.le.NGLT) then
         if (ikey.gt.ngl) ngl=ikey
         read(value,*) glf0d(ikey)
         read(cfit,*) nfit(NPAR1+NGLP*(ikey-1)+4)

       else if(key(1:4).eq.'GLF0'.and.ikey.ge.1.and.ikey.le.NGLT) then
         if (ikey.gt.ngl) ngl=ikey
         read(value,*) glf0(ikey)
         read(cfit,*) nfit(NPAR1+NGLP*(ikey-1)+2)


       else if(key(1:4).eq.'GLTD'.and.ikey.ge.1.and.ikey.le.NGLT) then
         if (ikey.gt.ngl) ngl=ikey
         read(value,*) gltd(ikey)
         read(cfit,*) nfit(NPAR1+NGLP*(ikey-1)+5)

       else if(key(1:4).eq.'GAIN') then
         read (value,*) gain

       else if(key(1:4).eq.'JUMP'.and.ikey.ge.1.and.ikey.le.NJUMP) then
         read(value,*) dct(ikey)
         read(cfit,*) itmp
         if (itmp.eq.0) nofitjump(ikey) = .true.

       else if(key(1:4).eq.'JUMP'.and.ikey.eq.-1) then
	 ijump = ijump+1 ! could check that this does not exceed NJUMP..
	 nflagjumps = ijump
c flag-based jump. line should look like:
c JUMP -flag flag_value jump_value fitflag jump_err
	 if(value(1:1).eq.'-') then 
           jumpflag(ijump) = value
	   jumpflagval(ijump) = temp      ! this was read earlier
	   call citem(line,ll,jn,temp,lf) ! read jump val
	   ! If value not specified, default to value=0 and fit=1
	   if (lf.eq.0) then
	     dct(ijump) = 0.0
	     nofitjump(ijump) = .false.
           else ! Read value and fit flag
	     read(temp,*) dct(ijump)
	     temp = '0'                     ! default to no fit
	     call citem(line,ll,jn,temp,lf) ! read fit flag
	     read(temp,*) itmp
	     if (itmp.eq.0) nofitjump(ijump) = .true.
	   endif
         else
           print *,"Error: only flag-based TEMPO2-style JUMPs allowed"
	   stop
         endif

       else if (key(1:6).eq.'T2EFAC') then
         iefac = iefac+1
         nflagefac = iefac
         if (value(1:1).eq.'-') then
           efacflag(iefac) = value
           efacflagval(iefac) = temp
           call citem(line,ll,jn,temp,lf) ! read value from line
           if (lf.eq.0) then ! no value, default to 1.0 (?)
             flagefac(iefac) = 1.0
           else
             read(temp,*) flagefac(iefac)
           endif
         else
           print *, "Error: specify a flag/value pair for T2EFAC"
           stop
         endif

       else if (key(1:7).eq.'T2EQUAD') then
         iequad = iequad+1
         nflagequad = iequad
         if (value(1:1).eq.'-') then
           equadflag(iequad) = value
           equadflagval(iequad) = temp
           call citem(line,ll,jn,temp,lf) ! read value from line
           if (lf.eq.0) then ! no value, default to 0.0 (?)
             flagequad(iequad) = 0.0
           else
             read(temp,*) flagequad(iequad)
           endif
         else
           print *, "Error: specify a flag/value pair for T2EQUAD"
           stop
         endif

       else if (key(1:7).eq.'ECORR') then
         iecorr = iecorr+1
         nflagecorr = iecorr
         if (value(1:1).eq.'-') then
           ecorrflag(iecorr) = value
           ecorrflagval(iecorr) = temp
           call citem(line,ll,jn,temp,lf) ! read value from line
           if (lf.eq.0) then ! no value, default to 0.0 (?)
             flagecorr(iecorr) = 0.0
           else
             read(temp,*) flagecorr(iecorr)
           endif
         else
           print *, "Error: specify a flag/value pair for ECORR"
           stop
         endif

       else if(key(1:4).eq.'INFO') then
         read (value,*) infoflag

       else if(key(1:7).eq.'SOLARN0'.and.(lk.eq.7.or.ikey.eq.0)) then
         read (value,*) solarn0
       else if(key(1:7).eq.'SOLARN0'.and.ikey.eq.1) then
         read (value,*) solarn01

c Red-noise params for GLS (Cholesky) fit
       else if (key(1:5).eq.'RNAMP') then
         read (value,*) rnamp
       else if (key(1:5).eq.'RNIDX') then
         read (value,*) rnidx
       else if (key(1:5).eq.'TCORR') then
         read (value,*) tcorr

c GLS covariance matrix is contained in the given file
       else if (key(1:8).eq.'DCOVFILE') then
          read (value,*) dcovfile

c Choosing which convention (value) to use for obliquity of the eclptic:
       else if (key(1:3).eq.'ECL') then
         call upcase(value)
         read (value,*) eclcon

c Do nothing parameters
      else if(key(1:4).eq.'HEAD') then
      else if(key(1:4).eq.'TRES') then
      else if(key(1:4).eq.'NTOA') then
      else if(key(1:5).eq.'DMXF1') then
      else if(key(1:5).eq.'DMXF2') then
c Some more, these are from tempo2 .. could check the values ... 
      else if(key(1:7).eq.'TIMEEPH') then
      else if(key(1:9).eq.'T2CMETHOD') then
      else if(key(1:19).eq.'CORRECT_TROPOSPHERE') then
      else if(key(1:10).eq.'DILATEFREQ') then

      else if(key(1:3).eq.'TOA') then   ! end of parameter list 
	goto 900

      else 
         if(key.ne.'        ')
     +      write(*,'('' Unrecognized input key: '',a)')key

      endif

      goto 10

 900  continue

C  Set fit flags for frequency derivatives higher than 1
C     First set nfcalc to the highest frequency derivative that has been read in
C     and set nffit to the highest frequency derivative that is to be fit
      nfcalc = 0
      nffit  = 0
      do i = 1, NFMAX
        if (ffit(i).gt.-1) nfcalc = i
        if (ffit(i).gt.0 ) nffit  = i
      enddo
C     Now set the nfit values for F1 and higher derivatives.
C     Note that all nfit values have already been zeroed, so we need only re-set
C     those which are to be 1.
              ! F1 (nfit(3))
      if (ffit(1).gt.0 .or. ffit(1).eq.-1.and.nffit.ge.1) nfit(3) = 1
              ! F2 (nfit(4))
      if (ffit(2).gt.0 .or. ffit(2).eq.-1.and.nffit.ge.2) nfit(4) = 1
              ! F2 (nfit(4))
      do i = 3, nffit
        if (ffit(i).gt.0 .or. ffit(i).eq.-1.and.nffit.ge.2) 
     +                                          nfit(NPAR11+i-2) = 1
      enddo

C  Warnings

      if(nfit(3).lt.0.or.nfit(3).gt.NFMAX.or.nfcalc.lt.0.or.
     +    nfcalc.gt.NFMAX)
     +   write(*,'('' WARNING: Fit parameter for F1 out of range'')')

      if(nfit(16).lt.0.or.nfit(16).gt.NDMCOFMAX.or.ndmcalc.lt.0
     +   .or.ndmcalc.gt.NDMCOFMAX)
     +   write(*,'('' WARNING: Fit parameter for DM out of range'')')

C IHS addition to prevent pathological double-fit of DM0 in one case
C DJN modified to allow it if at least one DMX value is held fixed at zero

      if(ndmcalc.ge.2 .and. nfit(16).ge.1 .and. usedmx) then
        if (ndmx>0)  then
           do i = 1, ndmx
             if (nfit(NPAR6+2*i-1) .EQ. 0) goto 910
           enddo
        endif 
        write(*,'(''ERROR: Cannot fit for DM0 when combining'',
     +          '' DM polynomial and DMX offset fits.  Fit for the'',
     +          '' polynomial in the relevant section first,'',
     +          '' then freeze the DM0 to DMN coefficients '',
     +          '' while fitting for the DMX values.'')')
        stop
  910   continue
      endif

      if(setecl)then
	 if(setequ)then
	    write (*,'(''ERROR: cannot mix ecliptic and equatorial'',
     +              '' coordinates'')')
	    stop
	 else
	    eclcoord=.true.
	 endif
      endif

      if(setfb .and. setpb)then
        write (*,'(''ERROR: cannot mix orbital period/pbdot and'',
     +       '' orbital frequency/fbdot'')')
        stop
      endif

c     if binary frequencies are input but binary model requires 
c     binary periods, make the conversion if possible
      if (setfb .and. nbin.ne.10) then
        pb(1) = (1.d0/fb(1))/86400.d0
        pbdot = -fb(2)/(fb(1)*fb(1))*1.d12
        if (nfit(NPAR3+1).ne.0) then
          nfit(12) = nfit(NPAR3+1)
          nfit(NPAR3+1) = 0
        endif
        if (nfit(NPAR3+2).ne.0) then
          nfit(18) = nfit(NPAR3+2)
          nfit(NPAR3+2) = 0
        endif
        do i = 3, NFBMAX
          if (nfit(NPAR3+i).ne.0 .or. fb(i).ne.0) then
            write (*,'(''ERROR: cannot use orbital frequency '',
     +           '' derivative fb'',i2,'' with binary model '',a8)'),
     +           i, bmodel(nbin)
            stop
          endif
        enddo
      endif

c     in any case, make sure pb(1) is set, so orbital phase calculations
c     for tempo.lis, resid2.tmp can be done (albeit with limited accuracy, 
c     since this will be done with pre-fit orbital period).
      if (setfb .and. pb(1).eq.0) pb(1) = (1.d0/fb(1))/86400.d0


c     if binary periods are input but binary model requires 
c     binary frequencies, make the conversion 
      if(setpb .and. nbin.eq.10 .and. pb(1).ne.0) then
        fb(1) = (1.d0/pb(1))/86400.d0
        fb(2) = -(pbdot*1.d-12)/(pb(1)*86400.)**2
        if (nfit(12).ne.0) then
          nfit(NPAR3+1) = nfit(12)
          nfit(12) = 0
        endif
        if (nfit(18).ne.0) then
          nfit(NPAR3+2) = nfit(18)
          nfit(18) = 0
        endif
      endif

      call getecliptic


      if(nbin.eq.0.and.(nfit(9).ne.0.or.nfit(10).ne.0.or.nfit(11).ne.0
     +     .or.nfit(12).ne.0.or.nfit(13).ne.0.or.a1(1).ne.0.))then
         write(*,'('' ERROR: Binary model not defined'')')
         stop
      endif

      if(nbin.ne.8 .and. set2dot)then
         write(*,'('' WARNING: No OM2DOT or X2DOT in '',a,'' !!!'')') 
     +        bmodel(nbin)
      endif

      if(nbin.ne.9 .and. setepsdot)then
         write(*,'('' WARNING: No EPS1DOT or EPS2DOT in '',a,'' !!!'')') 
     +        bmodel(nbin)
      endif

      if(nbin.eq.4 .and. a0.ne.0) then
	write (*,'('' WARNING: Use AFAC, not A0, with model DDGR.'')')
      endif

      if(nbin.eq.9)then
         if(seteps .and. t0asc.eq.0.)then
            write(*,'('' WARNING: TASC not set, use TASC=T0 !!!'')')
            t0asc=t0(1)
            t0(1)=0.
         endif

         if(.not.seteps .and. t0(1).eq.0.)then
            write(*,'('' WARNING: T0 not set, use T0=TASC !!!'')')
            t0(1)=t0asc
            t0asc=0.
         endif

         if((nfit(14).ne.0.or.omdot.ne.0.) .and. setepsdot)then
            write(*,'('' WARNING: omdot is not used !!!'')')
            omdot=0.
            nfit(14)=0
         endif

         if((nfit(25).ne.0.or.edot.ne.0.) .and. setepsdot)then         
            write(*,'('' WARNING: edot is not used !!!'')')
            edot=0.
            nfit(25)=0
         endif
      endif

      if(t0asc.ne.0.and.nbin.ne.1.and.nbin.ne.9)then
         write(*,'('' WARNING: TASC only for BT and ELL1 models,'',
     +        '' use T0=TASC !!!'')')
         t0(1)=t0asc
         t0asc=0.         
      endif

      if(tz.and.(ntzrmjd.eq.0)) 
     +  write (*,'('' WARNING: TZ mode, reference TOA not set'')')


      return
      end

c=======================================================================

      subroutine decolon(w)

C  Remove ':' from line

      character w*(*), ww*80
      j=0
      ww=' '
      do i=1,len(w)
         if(w(i:i).ne.':')then
            j=j+1
            ww(j:j)=w(i:i)
         endif
      enddo
      w=ww
      return
      end

c=======================================================================

      subroutine upcase(w)

C  Converts string (up to first blank) to upper case.

      character*(*) w
      do 10 i=1,len(w)
         if(w(i:i).eq.' ') go to 20
         j=ichar(w(i:i))
         if(j.ge.97.and.j.le.122) w(i:i)=char(j-32)
 10   continue
 20   return
      end


c=======================================================================

      integer function keyidx(s)

C  Finds underscore in string, returns integer number after underscore
C  Example:  s="GLEP_10"  would return 10.
C  Returns -1 if no underscore is found  
C  Returns 0 if underscore but no number is found

      character*(*) s
      integer k

      k = index(s,'_')
      if (k.eq.0) then
        keyidx = -1
      else
        keyidx = 0
        do 10 i = k+1, len(s)
          if (s(i:i).ge.'0' .and. s(i:i).le.'9') then
            keyidx = 10*keyidx + ichar(s(i:i)) - 48
          else
            goto 20
          endif
 10     continue
 20     continue
      endif

      return
      end

