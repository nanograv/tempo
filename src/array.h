c  Adjustable arrays:
	integer dnplsoff,dnplsptr
	real*8 dnpls(1)
        integer ddmchoff,ddmchptr
        real*8 ddmch(1)
        integer bufoff,bufptr
        real*8 buf(1)
        integer npmsavoff,npmsavptr
        integer npmsav(1)
        integer ksavoff,ksavptr
        integer ksav(1)

        common /array/
     +    dnplsoff,ddmchoff,bufoff,npmsavoff,ksavoff,
     +    dnplsptr,ddmchptr,bufptr,npmsavptr,ksavptr,
     +    dnpls,   ddmch,   buf,   npmsav,   ksav

