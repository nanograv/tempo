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
      f0=0.
      f1=0.
      f2=0.
      f3=0.
      do i=1,9
         f4(i)=0.
      enddo
      px=0.
      dm=0.
      do i=1,10
         dmcof(i)=0.
      enddo

      do i=1,4
         a1(i)=0.
         e(i)=0.
         t0(i)=0.
         pb(i)=0.
         omz(i)=0.
      enddo
      omdot=0.
      gamma=0.
      pbdot=0.
      si=0.
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

      nxoff=0
      do i=1,NJUMP
         xjdoff(1,i)=0.
         xjdoff(2,i)=0.
      enddo

      ngl=0
      do i=1,NGLT
         glph(i)=0.
         glf0p(i)=0.
         glf1p(i)=0.
         glf0d1(i)=0.
         gltd1(i)=0.
         glepoch(i)=0.
      enddo

      ndmcalc = 0

      nbin=0
      nclk=0
      nephem=1
      nits=1
      ncoord=1

      return
      end
C****************************************************************** ********
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

      character line*80, key*8, value*24, cfit*1, temp*80

      ll=80

      nskip = 0  ! counts parameter lines, which are skipped if TOAs are
		 !    read from this file
	
      rewind(parunit)  ! probably not needed, but a good safety check

 10   read(parunit,'(a)',end=900)line
      nskip = nskip + 1

C  Get key, value and cfit
      jn=1
      call citem(line,ll,jn,key,lk)
      if(key(1:1).eq.'#' .or. (key(1:1).eq.'C' .and. lk.eq.1))go to 10
      call upcase(key)
      call citem(line,ll,jn,value,lv)
      call citem(line,ll,jn,temp,lf)
      if(temp(1:1).eq.'#' .or. lf.eq.0)then
         cfit='0'
      else if(lf.eq.1)then
         cfit=temp(1:1)
      else
         write(*,'(''Illegal fit value: '',a)')temp(1:lf)
         stop
      endif

C  Control parameters

      if(key(1:4).eq.'NPRN')then
         read(value,*)nprnt

      else if(key(1:4).eq.'NITS')then
         read(value,*)nits

      else if(key(1:4).eq.'IBOO')then
         read(value,*)iboot

      else if(key(1:4).eq.'NDDM')then
         read(value,*)nddm

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

      else if(key(1:4).eq.'EPHE')then
         call upcase(value)
         do i=1,kephem
            if(value(1:5).eq.ephfile(i)(1:5)) then
		nephem=i
            	go to 14
	    endif
         enddo
         write(*,'(''Invalid EPHEM file name: '',a)')value(1:5)
         stop
 14      continue

      else if(key(1:6).eq.'TZRMJD')then
         read(value,*)tzrmjd

      else if(key(1:6).eq.'TZRFRQ')then
         read(value,*)tzrfrq

      else if(key(1:7).eq.'TZRSITE')then
         tzrsite=value(1:1)
         
C Period/Frequency parameters

      else if(key(1:2).eq.'P0' .or. (key(1:1).eq.'P' .and. lk.eq.1))then
         read(value,*)p0
         read(cfit,*)nfit(2)

      else if(key(1:2).eq.'P1'.or.key(1:4).eq.'PDOT')then
         read(value,*)p1
         read(cfit,*)nfit(3)

      else if(key(1:2).eq.'F0' .or. (key(1:1).eq.'F' .and. lk.eq.1))then
         read(value,*)f0
         read(cfit,*)nfit(2)

      else if(key(1:2).eq.'F1')then
         read(value,*)f1
         if(cfit.ge.'A')then
            call upcase(cfit)
            nfit(3)=ichar(cfit)-55
         else
            nfit(3)=ichar(cfit)-48
         endif

      else if(key(1:2).eq.'F2')then
         read(value,*)f2
         read(cfit,*)nfit(4)

      else if(key(1:2).eq.'F3')then
         read(value,*)f3
         read(cfit,*)ifit

      else if(key(1:1).eq.'F' .and.
     +           key(2:2).ge.'4' .and. key(2:2).le.'9') then
        read (key(2:2),*) jj
        read (value,*)f4(jj-3)
        read(cfit,*) ifit
        if (ifit.gt.0) nfit(3)=max(nfit(3),jj)

      else if(key(1:1).eq.'F' .and.
     +           key(2:2).ge.'A' .and. key(2:2).le.'C') then
        jj = ichar(key(2:2))-55  ! A=10, B=11, C=12
        read (value,*)f4(jj-3)
        read(cfit,*) ifit
        if (ifit.gt.0) nfit(3)=max(nfit(3),jj)

      else if(key(1:4).eq.'PEPO')then
         read(value,*)pepoch

C Position parameters

      else if(key(1:3).eq.'PSR')then
         psrname=value(1:12)

      else if(key(1:3).eq.'DEC')then
         call decolon(value)
         read(value,*)pdec
         read(cfit,*)nfit(5)

      else if(key(1:2).eq.'RA')then
         call decolon(value)
         read(value,*)pra
         read(cfit,*)nfit(6)

      else if(key(1:4).eq.'PMDE')then
         read(value,*)pmdec
         read(cfit,*)nfit(7)

      else if(key(1:4).eq.'PMRA')then
         read(value,*)pmra
         read(cfit,*)nfit(8)

      else if(key(1:4).eq.'PMRV')then
         read(value,*)pmrv
         read(cfit,*)nfit(36)

      else if(key(1:2).eq.'PX')then
         read(value,*)px
         read(cfit,*)nfit(17)

      else if(key(1:3).eq.'DM0'.or.key(1:2).eq.'DM'.and.lk.eq.2)then
         read(value,*)dm
         if(cfit.le.'9')then
            itmp=ichar(cfit)-48
         else
            call upcase(cfit)
            itmp=ichar(cfit)-55
         endif
         nfit(16)=max(nfit(16),itmp)

      else if(key(1:2).eq.'DM'.and.
     +        key(3:3).ge.'1'.and.key(3:3).le.'9') then 
        read (key(3:3),*) jj
        read (value,*) dmcof(jj)
        if (cfit.gt.'0') nfit(16)=max(nfit(16),jj+1)
        ndmcalc=max(ndmcalc,jj+1)

C Binary parameters

      else if(key(1:4).eq.'BINA')then
         call upcase(value)
         do i=1,10
            if(value(1:8).eq.bmodel(i))go to 20
         enddo
         write(*,100)value(1:8)
 100     format(' Warning: binary model - ',a,' - not recognized')
         go to 22
c 20      nbin=i-1               ! ### Check this !!! (Works in Linux/Intel)
 20      nbin=i
 22      continue

      else if(key(1:4).eq.'A1_1'.or.(key(1:2).eq.'A1'.and.lk.eq.2))then
         read(value,*)a1(1)
         read(cfit,*)nfit(9)

      else if(key(1:3).eq.'E_1 '.or.(key(1:1).eq.'E'.and.lk.eq.1))then
         read(value,*)e(1)
         read(cfit,*)nfit(10)

      else if(key(1:4).eq.'T0_1'.or.(key(1:2).eq.'T0'.and.lk.eq.2))then
         read(value,*)t0(1)
         read(cfit,*)nfit(11)

      else if(key(1:4).eq.'PB_1'.or.(key(1:2).eq.'PB'.and.lk.eq.2))then
         read(value,*)pb(1)
         read(cfit,*)nfit(12)

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
         read(value,*)omdot
         read(cfit,*)nfit(14)

      else if(key(1:5).eq.'GAMMA')then
         read(value,*)gamma
         read(cfit,*)nfit(15)

      else if(key(1:5).eq.'PBDOT')then
         read(value,*)pbdot
         read(cfit,*)nfit(18)

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
         read(value,*)xdot
         read(cfit,*)nfit(24)

      else if(key(1:4).eq.'EDOT')then
         read(value,*)edot
         read(cfit,*)nfit(25)

      else if(key(1:6).eq.'XOMDOT')then
         read(value,*)xomdot
         read(cfit,*)nfit(37)

      else if(key(1:6).eq.'XPBDOT')then
         read(value,*)xpbdot
         read(cfit,*)nfit(38)

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

C Glitches

      else if(key(1:6).eq.'GLEP_1')then
         read(value,*)glepoch(1)

      else if(key(1:6).eq.'GLPH_1')then
         read(value,*)glph(1)
         read(cfit,*)nfit(61)

      else if(key(1:6).eq.'GLF0_1')then
         read(value,*)glf0p(1)
         read(cfit,*)nfit(62)

      else if(key(1:6).eq.'GLF1_1')then
         read(value,*)glf1p(1)
         read(cfit,*)nfit(63)

      else if(key(1:7).eq.'GLF0D_1')then
         read(value,*)glf0d1(1)
         read(cfit,*)nfit(64)

      else if(key(1:6).eq.'GLTD_1')then
         call citem(line,ll,jn,cfit,ii)
         read(value,*)gltd1(1)
         read(cfit,*)nfit(65)

      else if(key(1:6).eq.'GLEP_2')then
         read(value,*)glepoch(2)

      else if(key(1:6).eq.'GLPH_2')then
         read(value,*)glph(2)
         read(cfit,*)nfit(66)

      else if(key(1:6).eq.'GLF0_2')then
         read(value,*)glf0p(2)
         read(cfit,*)nfit(67)

      else if(key(1:6).eq.'GLF1_2')then
         read(value,*)glf1p(2)
         read(cfit,*)nfit(68)

      else if(key(1:7).eq.'GLF0D_2')then
         read(value,*)glf0d1(2)
         read(cfit,*)nfit(69)

      else if(key(1:6).eq.'GLTD_2')then
         read(value,*)gltd1(2)
         read(cfit,*)nfit(70)

      else if(key(1:6).eq.'GLEP_3')then
         read(value,*)glepoch(3)

      else if(key(1:6).eq.'GLPH_3')then
         read(value,*)glph(3)
         read(cfit,*)nfit(71)

      else if(key(1:6).eq.'GLF0_3')then
         read(value,*)glf0p(3)
         read(cfit,*)nfit(72)

      else if(key(1:6).eq.'GLF1_3')then
         read(value,*)glf1p(3)
         read(cfit,*)nfit(73)

      else if(key(1:7).eq.'GLF0D_3')then
         read(value,*)glf0d1(3)
         read(cfit,*)nfit(74)

      else if(key(1:6).eq.'GLTD_3')then
         read(value,*)gltd1(3)
         read(cfit,*)nfit(75)

      else if(key(1:6).eq.'GLEP_4')then
         read(value,*)glepoch(4)

      else if(key(1:6).eq.'GLPH_4')then
         read(value,*)glph(4)
         read(cfit,*)nfit(76)

      else if(key(1:6).eq.'GLF0_4')then
         read(value,*)glf0p(4)
         read(cfit,*)nfit(77)

      else if(key(1:6).eq.'GLF1_4')then
         read(value,*)glf1p(4)
         read(cfit,*)nfit(78)

      else if(key(1:7).eq.'GLF0D_4')then
         read(value,*)glf0d1(4)
         read(cfit,*)nfit(79)

      else if(key(1:6).eq.'GLTD_4')then
         read(value,*)gltd1(4)
         read(cfit,*)nfit(80)

      else if(key(1:4).eq.'HEAD') then
c       (Do nothing) (DJN)

      else if(key(1:3).eq.'TOA') then   ! end of parameter list (DJN)
	goto 900

      else 
         if(key.ne.'        ')
     +        write(*,'('' Unrecognized input key: '',a)')key

      endif

      go to 10

 900  continue

      if(nfit(3).lt.0.or.nfit(3).gt.12)then
         write(*,'(''WARNING - Fit parameter for F1 out of range'')')
      endif

      if(nfit(16).lt.0.or.nfit(16).gt.10)then
         write(*,'(''WARNING - Fit parameter for DM out of range'')')
      endif

      if(nfit(61)+nfit(62)+nfit(63)+nfit(64)+nfit(65).ne.0)ngl=1
      if(nfit(66)+nfit(67)+nfit(68)+nfit(69)+nfit(70).ne.0)ngl=2
      if(nfit(71)+nfit(72)+nfit(73)+nfit(74)+nfit(75).ne.0)ngl=3
      if(nfit(76)+nfit(77)+nfit(78)+nfit(79)+nfit(80).ne.0)ngl=4

      if(nbin.eq.0.and.(nfit(9).ne.0.or.nfit(10).ne.0.or.nfit(11).ne.0
     :  .or.nfit(12).ne.0.or.nfit(13).ne.0))then
         write(*,'('' WARNING - Binary model not defined'')')
      endif

      return
      end

C***************************************************************************

      subroutine outpar(irh,irm,rsec,ers,decsgn,idd,idm,dsec,eds)

      implicit real*8 (A-H,O-Z)
      character decsgn*1, fit1*3

      include 'dim.h'     
      include 'acom.h'
      include 'bcom.h'
      include 'trnsfr.h'
      include 'clocks.h'
      include 'dp.h'
      include 'glitch.h'
      include 'eph.h'
      include 'tz.h'

      c=360.*3600./6.2831853
      cc=1.d-9*c*365.25*8.64d7
      fit1='  1'

      write(71,'(''PSR              '',a)')psrname

      irs=rsec
      rs=rsec-irs
      ids=dsec
      ds=dsec-ids
      if(nfit(6).gt.0)then
         write(71,1006)irh,irm,irs,rs,fit1,ers
      else
         write(71,1006)irh,irm,irs,rs
      endif
 1006 format('RA',i9.2,':',i2.2,':',i2.2,f9.8,a,f20.8)

      if(nfit(5).gt.0)then
         write(71,1005)decsgn,idd,idm,ids,ds,fit1,eds
      else
         write(71,1005)decsgn,idd,idm,ids,ds
      endif
 1005 format('DEC',6x,a,i2.2,':',i2.2,':',i2.2,f8.7,,a,f20.7)

      if(nfit(2).gt.0)then
         write(71,1002)f0,fit1,ferr(2)*1.d-9
      else
         write(71,1002)f0
      endif
 1002 format('F0',f24.16,a,f20.16)

      if(nfit(3).gt.0)then
         write(71,1003)f1,fit1,ferr(3)*1.d-18
      else
         write(71,1003)f1
      endif
 1003 format('F1',1p,d24.12,a,d20.12)

      if(f2.ne.0.)then
         if(nfit(4).gt.0)then
            write(71,1004)f2,fit1,ferr(4)*1.d-27
         else
            write(71,1004)f2
         endif
      endif
 1004 format('F2',1p,d24.12,a,d20.12)

      if(f3.ne.0.)then
         if(nfit(51).gt.0)then
            write(71,1051)f3,fit1,ferr(51)*1.d-36
         else
            write(71,1051)f3
         endif
      endif
 1051 format('F3',1p,d24.12,a,d20.12)

      do i = 1, 9
         if (f4(i).ne.0)then
            if(nfit(51+i).gt.0)then
               write(71,1052)i+3,f4(i),fit1,ferr(51+i)*(1.d-9)**(i+4)
            else
               write(71,1052)i+3,f4(i)
            endif
         endif
      enddo
 1052 format('F',z1,1p,d24.12,a,d20.12)

      write(71,'(''PEPOCH'',f20.6)')pepoch

      if(nfit(16).gt.0)then
         write(71,1016)dm,fit1,ferr(16)
      else
         write(71,1016)dm
      endif
 1016 format('DM',f24.6,a,f20.6)

      do i = 1, 9
         if(dmcof(i).ne.0)then
            if(nfit(40+i).gt.0)then
               write(71,1041)i,dmcof(i),fit1,ferr(40+i)
            else
               write(71,1041)i,dmcof(i)
            endif
         endif
      enddo
 1041 format('DM',z1,1p,d23.12,a,d20.12)

      if(pmra.ne.0.)then
         if(nfit(8).gt.0)then
            write(71,1008)pmra,fit1,ferr(8)*cc
         else
            write(71,1008)pmra
         endif
      endif
 1008 format('PMRA',f22.4,a,f20.4)

      if(pmdec.ne.0.)then
         if(nfit(7).gt.0)then
            write(71,1007)pmdec,fit1,ferr(7)*cc
         else
            write(71,1007)pmdec
         endif
      endif
 1007 format('PMDEC',f21.4,a,f20.4)

      if(pmrv.ne.0.)then
         if(nfit(36).gt.0)then
            write(71,1036)pmrv,fit1,10.*c*ferr(36)
         else
            write(71,1036)pmrv
         endif
      endif
 1036 format('PMRV',f22.4,a,f20.4)

      if(px.ne.0.)then
         if(nfit(17).gt.0)then
            write(71,1017)px,fit1,ferr(17)
         else
            write(71,1017)px
         endif
      endif
 1017 format('PX',f24.4,a,f20.4)

      if(ngl.gt.0)then
         do i=1,ngl
            ii=60+(i-1)*NGLP
            write(71,'(''GLEP_'',i1,f18.6)')i,glepoch(i)
            if(nfit(ii+1).gt.0)then
               write(71,1061)i,glph(i),fit1,ferr(ii+1)
            else
               write(71,1061)i,glph(i)
            endif
 1061       format('GLPH_',i1,f20.6,a,f20.6)
            if(nfit(ii+2).gt.0)then
               write(71,1062)i,glf0p(i),fit1,ferr(ii+2)
            else
               write(71,1062)i,glf0p(i)
            endif
 1062       format('GLF0_',i1,1p,d20.8,a,d20.8)
            if(nfit(ii+3).gt.0)then
               write(71,1063)i,glf1p(i),fit1,ferr(ii+3)
            else
               write(71,1063)i,glf1p(i)
            endif
 1063       format('GLF1_',i1,1p,d20.8,a,d20.8)
            if(nfit(ii+4).gt.0)then
               write(71,1064)i,glf0d1(i),fit1,ferr(ii+4)
            else
               write(71,1064)i,glf0d1(i)
            endif
 1064       format('GLF0D_',i1,1p,d19.8,a,d20.8)
            if(nfit(ii+5).gt.0)then
               write(71,1065)i,gltd1(1),fit1,ferr(ii+5)
            else
               write(71,1065)i,gltd1(1)
            endif
 1065       format('GLTD_',i1,f20.4,a,f20.4)
         enddo
      endif

      write(71,'(''EPHEM'',13x,,a)')ephfile(nephem)(1:5)
      write(71,'(''CLK'',15x,a)')clklbl(nclk)
      write(71,'(''TZRMJD '',f21.13)')tzrmjd
      write(71,'(''TZRFRQ '',f19.3)')tzrfrq
      write(71,'(''TZRSITE '',17x,a)')tzrsite
      if(nprnt.gt.0)write(71,'(''NPRNT'',i21)')nprnt
      if(nits.gt.0)write(71,'(''NITS'',i22)')nits
      if(iboot.gt.0)write(71,'(''IBOOT'',i21)')iboot
      if(nddm.gt.0)write(71,'(''NDDM'',i22)')nddm

      return
      end

C***************************************************************************

      subroutine outbinpar

      implicit real*8 (A-H,O-Z)

      include 'dim.h'     
      include 'acom.h'
      include 'bcom.h'
      include 'trnsfr.h'
      include 'dp.h'
      include 'orbit.h'

      character fit1*3

      fit1='  1'

      write(71,'(''BINARY'',12x,a)')bmodel(nbin)

      if(nfit(9).gt.0)then
         write(71,1009)a1(1),fit1,ferr(9)
      else
         write(71,1009)a1(1)
      endif
 1009 format('A1',f24.9,a,f20.9)

      if(nfit(10).gt.0)then
         write(71,1010)e(1),fit1,ferr(10)
      else
         write(71,1010)e(1)
      endif
 1010 format('E',f25.10,a,f20.10)

      if(nfit(11).gt.0)then
         write(71,1011)t0(1),fit1,ferr(11)
      else
         write(71,1011)t0(1)
      endif
 1011 format('T0',f24.9,a,f20.9)

      if(nfit(12).gt.0)then
         write(71,1012)pb(1),fit1,ferr(12)
      else
         write(71,1012)pb(1)
      endif
 1012 format('PB',f24.12,a,f20.12)

      if(nfit(13).gt.0)then
         write(71,1013)omz(1),fit1,ferr(13)*57.295
      else
         write(71,1013)omz(1)
      endif
 1013 format('OM',f24.6,a,f20.6)

      if(omdot.ne.0.)then
         if(nfit(14).gt.0)then
            write(71,1014)omdot,fit1,ferr(14)
         else
            write(71,1014)omdot
         endif
      endif
 1014 format('OMDOT',f21.7,a,f20.7)

      if(xomdot.ne.0.)then
         if(nfit(37).gt.0)then
            write(71,1037)xomdot,fit1,ferr(37)
         else
            write(71,1037)xomdot
         endif
      endif
 1037 format('XOMDOT',f20.7,a,f20.7)

      if(gamma.ne.0.)then
         if(nfit(15).gt.0)then
            write(71,1015)gamma,fit1,ferr(15)
         else
            write(71,1015)gamma
         endif
      endif
 1015 format('GAMMA',f21.7,a,f20.7)

      if(pbdot.ne.0.)then
         if(nfit(18).gt.0)then
            write(71,1018)pbdot*1.d12,fit1,ferr(18)*1.d6
         else
            write(71,1018)pbdot*1.d12
         endif
      endif
 1018 format('PBDOT',f21.7,a,f20.7)

      if(xpbdot.ne.0.)then
         if(nfit(38).gt.0)then
            write(71,1038)xpbdot*1.d12,fit1,ferr(38)*1.d6
         else
            write(71,1038)xpbdot*1.d12
         endif
      endif
 1038 format('XPBDOT',f20.7,a,f20.7)

      if(si.ne.0.)then
         if(nfit(20).gt.0)then
            write(71,1020)si,fit1,ferr(20)
         else
            write(71,1020)si
         endif
      endif
 1020 format('SINI',f22.6,a,f20.6)

      if(am.ne.0.)then
         if(nfit(21).gt.0)then
            write(71,1021)am,fit1,ferr(21)
         else
            write(71,1021)am
         endif
      endif
 1021 format('MTOT',f22.6,a,f20.6)

      if(am2.ne.0.)then
         if(nfit(22).gt.0)then
            write(71,1022)am2,fit1,ferr(22)
         else
            write(71,1022)am2
         endif
      endif
 1022 format('M2',f24.6,a,f20.6)

      if(dth.ne.0.)then
         if(nfit(23).gt.0)then
            write(71,1023)dth*1.d6,fit1,ferr(23)
         else
            write(71,1023)dth*1.d6
         endif
      endif
 1023 format('DTHETA',f20.6,a,f20.6)
      if(xdot.ne.0.)then
         if(nfit(24).gt.0)then
            write(71,1024)xdot*1.d12,fit1,ferr(24)*1.d12
         else
            write(71,1024)xdot*1.d12
         endif
      endif
 1024 format('XDOT',f22.6,a,f20.6)
      
      if(edot.ne.0.)then
         if(nfit(25).gt.0)then
            write(71,1025)edot*1.d12,fit1,ferr(25)*1.d12
         else
            write(71,1025)edot*1.d12
         endif
      endif
 1025 format('EDOT',f22.6,a,f20.6)

      if(afac.ne.0.)then
        write (71,10251)afac
      endif
10251 format('AFAC',f22.7)

      if(nbin.ge.9)then
         jbin=nbin-7
         do j=2,jbin
            jj=16+j*5
            if(nfit(jj).gt.0)then
               write(71,1026)j,a1(j),fit1,ferr(jj)
            else
               write(71,1026)j,a1(j)
            endif
 1026       format('A1_',i1,f22.9,a,f20.9)
            if(nfit(jj+1).gt.0)then
               write(71,1027)j,e(j),fit1,ferr(jj+1)
            else
               write(71,1027)j,e(j)
            endif
 1027       format('E_',i1,f23.9,a,f20.9)
            if(nfit(jj+2).gt.0)then
               write(71,1028)j,t0(j),fit1,ferr(jj+2)
            else
               write(71,1028)j,t0(j)
            endif
 1028       format('T0_',i1,f22.9,a,f20.9)
            if(nfit(jj+3).gt.0)then
               write(71,1029)j,pb(j),fit1,ferr(jj+3)
            else
               write(71,1029)j,pb(j)
            endif
 1029       format('PB_',i1,f22.12,a,f20.12)
            if(nfit(jj+4).gt.0)then
               write(71,1030)j,omz(j),fit1,ferr(jj+4)*57.295
            else
               write(71,1030)j,omz(j)
            endif
 1030       format('OM_',i1,f22.6,a,f20.6)
         enddo
      endif

      return
      end

C***************************************************************************
      subroutine decolon(w)

C  remove ':' from line

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

C*************************************************************************** 

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

