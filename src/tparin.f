      subroutine tparin(nostop,tz,lw,lpth,nparmax,nptsmax,
     +     version,npulsefile,infile,path,resfile1,hlpfile,parfile)

c     tparin -- Tempo PARamater INput
c     parses tempo command-line parameters and sets variables appropriately

      implicit real*8 (a-h,o-z)

      include 'dim.h'
      include 'acom.h'
      include 'vcom.h'

c     following variables are set in this routine:

      logical nostop
      logical lw
      logical tz
      integer lpth
      integer nparmax, nptsmax
      real*8 version
      character*160 npulsefile
      character*80 infile
      character*80 path
      character*80 resfile1
      character*80 hlpfile
      character*40 parfile

c     following variables are also set, but are in common blocks
c     and hence not formally declared here      
c      logical xitoa 
c      character*1 obsflag
c      logical gro
c      logical jumpout
c      logical lresid1
c      logical npulsein
c      logical npulseout
c      logical oldpar
c      logical psrframe
 

c     variables used internally in this routine:

      integer narg, iarg, ipar
      integer iargc
      integer i, ii
      character*160 s, s2

      character*160 getparm  ! external function


c     default values of parameters


      gro = .false.
      jumpout = .false.
      lresid1 = .false.
      nostop = .false.
      npulsein = .false.
      npulseout = .false.
      oldpar = .false.
      psrframe = .false.
      lw = .true.
      ssdmflag = 1
      tz = .false.
      xitoa = .false.
      nparmax = NPARDEF
      nptsmax = NPTSDEF
      parfile = 'def'

      infile = ''	

      call getenv ('TEMPO',path)
      lpth = index(path,' ')-1
      hlpfile = path(1:lpth)//'/tempo.hlp'	


      narg = iargc()
      iarg = 0
      ipar = 0
      if (narg.lt.1) goto 9999  ! error

      do while (iarg.lt.narg)

        iarg = iarg + 1
        call getarg (iarg,s)
        ii = index(s,' ') - 1

        if (s(1:1).eq.'-') then

          i = 2
          do while (i.le.ii)
            if (s(i:i).eq.'c') then
              nostop = .true.
            else if (s(i:i).eq.'d') then
              path = getparm(s,i,ii,iarg,narg)
              if (path.eq.'') goto 9999 ! error
              lpth = index(path,' ')-1
              hlpfile = path(1:lpth)//'/tempo.hlp'
            else if (s(i:i).eq.'f') then
              parfile = getparm(s,i,ii,iarg,narg)
              if (parfile.eq.'') goto 9999 ! error
            else if (s(i:i).eq.'g') then
              gro = .true.
              obsflag = 'P'
              if (s(i+1:i+1).ne.' ') obsflag = s(i+1:i+1)
              i = 9999
            else if (s(i:i).eq.'h') then
              goto 9999         ! error/help
            else if (s(i:i).eq.'j') then
              jumpout = .true.
            else if (s(i:i).eq.'l') then
              s2 = getparm(s,i,ii,iarg,narg)
              read (s2,*) nparmax
            else if (s(i:i).eq.'m') then
              s2 = getparm(s,i,ii,iarg,narg)
              read (s2,*) nptsmax
            else if (s(i:i).eq.'n') then
              i = i + 1
              if (s(i:i).eq.'i') then
                npulsein = .true.
                npulsefile = getparm(s,i,ii,iarg,narg)
              elseif (s(i:i).eq.'o') then
                npulseout = .true.
                npulsefile = getparm(s,i,ii,iarg,narg)
              else
                goto 9999       ! error
              endif
            else if (s(i:i).eq.'o') then
              oldpar = .true.
            else if (s(i:i).eq.'p') then
              psrframe = .true.
            else if (s(i:i).eq.'r') then
              lresid1 = .true.
              open(30,file=resfile1,form='unformatted',status='unknown')
            else if (s(i:i).eq.'s') then
              s2 = getparm(s,i,ii,iarg,narg)
              read (s2,*) ssdmflag
            else if (s(i:i).eq.'v') then
              write(*,1001) version
 1001         format (' Tempo v ',f6.3)
              stop
            else if (s(i:i).eq.'w') then
              lw = .false.
            else if (s(i:i).eq.'x') then
              xitoa = .true.
              open(45,file='itoa.out',status='unknown')
            else if (s(i:i).eq.'z') then
              tz = .true.
            else
              goto 9999  ! error
            endif
              
            i = i + 1
          enddo

        else

          ipar = ipar + 1  
          if (ipar.eq.1) then
            infile = s
          else
            goto 9999           ! error
          endif
          
        endif
        
      enddo

      if (ipar.lt.1 .and. .not. tz) goto 9999   ! error -- not enough parameters

      return

 9999 call system ('more '//hlpfile)
      stop
      end

C----------------------------------------------------------------

        character*160 function getparm(s,i,ii,iarg,narg)

        implicit none

        character*160 s
        integer i, ii, iarg, narg

        integer i0

        getparm = "" ! default (indicates error if this is read upon return)

        i0 = i+1
        if (i.eq.ii) then
          iarg = iarg + 1
          if (iarg.gt.narg) goto 9000
          call getarg(iarg,s)
          i0 = 1
          ii = index(s,' ')
        endif
        
        i = 9999


        getparm = s(i0:ii)

 9000   continue
        return
	end
