c      $Id$
	subroutine setup(version,infile,obsyfile,alng,nsmax,parfile)

	implicit real*8 (a-h,o-z)
	parameter (TWOPI=6.28318530717958648d0)
	real*8 alng(36)
	include 'dim.h'
	include 'acom.h'
	integer time
	character timstr*24,obsnam*12,obsnum*35,damoyr*9,parfile*40
	character*80 infile,obsyfile
	data ault/499.004786d0/
	data obsnum/'123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

C###	call fdate(timstr)
	nsec=time()
	timstr=damoyr(40587+nsec/86400)
	nsec=mod(nsec,86400)
	nhr=nsec/3600
	nmin=(nsec-3600*nhr)/60
	nsec=mod(nsec,60)
	write(timstr(10:),1000) nhr,nmin,nsec
 1000	format(2x,i2.2,':',i2.2,':',i2.2,' UTC')

	write(31,30) version,timstr(1:23)
30	format(55(1h*)/1h*,53x,1h*/
     +    '*',17x,'TEMPO version',f7.3,16x,1h*/1h*,53x,1h*/,
     +    1h*,11x,'Analysis of pulsar timing data',12x,1h*/
     +    1h*,53x,1h*/
     +    1h*,14x,a24,15x,1h*/
     +    1h*,53x,1h*/55(1h*)///
     +    20x,'Observatory Coordinates'//
     +    5x,'Observatory     Geodetic lat   Geodetic long   Elevation'/
     +    22x,'ddmmss.ss       ddmmss.ss         m'/)

C  Read the geodetic coordinates of up to 36 observatories.
C  Icoord = 0 if geodetic; 1 = X, Y, Z geocentric coordinates, where
C  alat = X, along = Y, and elev = Z

	open(2,file=obsyfile,status='old',err=32)
	go to 34
 32	write(*,'(''Error opening Observatory file: '',a)')obsyfile
	STOP
 34	do 10 j=1,36
	 read(2,40,end=11) alat,along,elev,icoord,obsnam,obskey(j)
40	 format(3f15.0,2x,i1,2x,a12,8x,a5)
	 if(alat.ne.0.) write(31,50) obsnum(j:j),obsnam,alat,along,
     +        elev
50	 format(a1,3x,a12,f15.2,f16.2,f12.1)
	 if (icoord.eq.0) then
	   alat=ang(1,alat)
	   alng(j)=ang(1,along)
c		old approach is IAU 1964 power series.  tests show this is
c		good to 5.10^-10 for 1964 figure.  But series misapplied
c		below (erad series should be in alat, not hlt)
c		so actual error is 1.d-5.  Furthermore, want flexibility
c		of going to any figure to use full equation approach.
c		see AA K11.
c		IAU 1976 flattening f, equatorial radius a
		aa_f = 1.d0/298.257d0
		aa_a = 6378140.d0
		aa_c = 1.d0/sqrt(1.d0 + (-2.d0 + aa_f)*aa_f*sin(alat)*sin(alat))
		aa_arcf = (aa_a*aa_c + elev)*cos(alat)
		aa_arsf = (aa_a*(1.d0 - aa_f)*(1.d0 - aa_f)*aa_c + elev)*sin(alat)
		hlt(j) = datan2(aa_arsf,aa_arcf)
		erad = sqrt(aa_arcf*aa_arcf + aa_arsf*aa_arsf)
!	   delat=-692.7430d0*dsin(2.d0*alat)+1.1633d0*dsin(4.d0*alat)-
!     +       0.0026d0*dsin(6.d0*alat)
!	   hlt(j)=alat+delat*TWOPI/12.96d5
!	   erad=6378160.0d0*(0.998327073d0+0.001676438d0*
!     +       dcos(2.d0*hlt(j))-3.519d-6*dcos(4.d0*hlt(j))+
!     +       8.d-9*dcos(6.d0*hlt(j)))+elev
	   hrd(j)=erad/(2.99792458d8*ault)
	 else
	   erad=dsqrt(alat**2+along**2+elev**2)
	   hlt(j)=dasin(elev/erad)
	   alng(j)=datan2(-along,alat)
	   hrd(j)=erad/(2.99792458d8*ault)
	 end if
10	continue
11	nsmax=j-1
	if(nsmax.ne.36) then
	  do 12 i=nsmax+1,36
	   alng(i)=0.d0
	   hlt(i)=0.d0
12	   hrd(i)=0.d0
	end if
	close(2)
	len=index(infile,' ')-1
	if(oldpar.or.parunit.eq.50)then
	   write (31,60) infile(1:len)
 60	   format (/'Input data from ',a)
	else
	   len1=index(parfile,' ')-1
	   write (31,61) infile(1:len),parfile(1:len1)
 61	   format (/'Input data from ',a,',  Parameters from ',a)
	endif
	write(31,62)obsyfile
 62	format('Observatory data from ',a55)

	return
	end

