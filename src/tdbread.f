c      $Id$
      subroutine tdbread(JD1,JD2,ctatv)

c     tdbread -- read and interpolate TDB-TDT file

c     file should have been written by tdbgen, which uses the Fairhead
c     et al. approximate analytical formula to generate Chebyshev polynomial
c     coefficients, which are read and interpolated here.

      include 'tdbcom.h'

c     input variables: 
      real*8 jd1, jd2           ! jd1+jd2 is TDT:
                                ! jd1 is integer part,
                                ! jd2 is fractional part

      real*8 ctatv              ! output TDB-TDT

      real*8 buf(16)            ! max size of data buffer

      real*8 jda, jdb, t(2)

      integer i
      integer nr

      logical bigendian

      save buf                  ! For Linux

c following copied from JPL ephemeris routines -- is it also valid here?

c     Ephemeris always starts on an MJD that is some integer+0.5, and record
c     length is always an integer number of days.  So set jda to the nearest
c     value of int+0.5 below jd1+jd2, and set jdb to the remaining fraction
c     of a day.  Note:  assumes jd1 an integer and 0<=jd2<1.  


      if (jd2.ge.0.5d0) then
        jda = jd1 + 0.5d0
        jdb = jd2 - 0.5d0
      else
        jda = jd1 - 0.5d0
        jdb = jd2 + 0.5d0
      endif

      nr = int((jda-tdbd1)/tdbdt)+2 ! record number in file;  "+2" skips
                                ! skips over the hdr rec 

      if (nr.lt.1 .or. jd1+jd2.gt.tdbd2) then
        write (*,*) "Date ",jda+jdb," out of range of TDB-TDT table (",
     +       tdbd1,"-",tdbd2,")"
        stop
      endif

      if (nr.ne.tdbnrl) then
        tdbnrl = nr
        read (tdbunit, rec=nr) (buf(i),i=1,tdbncf)
	if (.not.bigendian()) call dbyterev(buf,tdbncf)
      endif

      t(1) = ((jda-((nr-2)*tdbdt+tdbd1))+jdb)/tdbdt ! fraction within record
      t(2) = 1.                 ! unused 

      call interp(buf,t,tdbncf,1,1,1,ctatv)

      return
      end

