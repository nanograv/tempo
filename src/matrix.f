c      $Id$
	program matrix
	implicit real*8 (a-h,o-z)
	include 'dim.h'
	real*8 a(NPA,NPA),sig(NPA),gcor(NPA),fac(NPA)
	character*5 param(39)
	character fmt*32,fmt2*32,arg*32

	narg=iargc()
	arg=' '
	if(narg.ge.1) then
	  call getarg(1,arg)
          if(narg.gt.1.or.(arg(1:2).lt.'-0'.or.arg(1:2).gt.'-9')) then
	    print*,'Usage: matrix [-n]'
	    print*,
     +       '       n = number of decimal places [default 2]'
	    go to 999
	  endif
	endif
	ndec=2
	if(arg(1:1).eq.'-') read(arg(2:3),*,err=1) ndec
1	fmt='(f4.0,2x,a4,1h|,1x,35f 6.02)'
	write(fmt(23:27),1001) 4+ndec,ndec
1001    format(i2,'.',i2.2)

	open(72,file='matrix.tmp',status='old',form='unformatted',err=998)
	read(72) nn
	rewind 72

	do 10 n=1,nn
10	read(72) mm,j,param(j),gcor(j),sig(j),(a(j,k),k=1,mm)
	fmt2='(12x,25a'//fmt(23:24)//')'

	write(*,fmt2) (param(j),j=1,nn)
	write(*,1010) ('-',i=1,nn*(ndec+4)+2)
1010	format(10x,255a1)

	do 20 j=1,nn
20	write(*,fmt) float(j),param(j),(a(j,k),k=1,j)

	fmt2='(/7x,5hgcor:,'//fmt(20:)
	write(*,fmt2) (gcor(j),j=1,nn)
	fmt2='(12x,25a'//fmt(23:24)//')'
	write(*,1010) ('-',i=1,nn*(ndec+4)+2)
	write(*,fmt2) (param(j),j=1,nn)
	write(*,1010) ('-',i=1,nn*(ndec+4)+2)

	read(72,end=999) nboot
	write(*,1020) nboot
1020	format(/' Uncertainty factors from bootstrap method (',
     +    i4,' iterations):')

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,5h rms:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,5h  -2:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,5h  -1:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,5h  +1:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)

	read(72) mm,(fac(j),j=1,mm)
	fmt2='(7x,5h  +2:,'//fmt(20:)
	write(*,fmt2) (fac(j),j=1,mm)
	go to 999

998	print*,'Cannot open matrix.tmp'
999	end
