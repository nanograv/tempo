c  Adjustable arrays:
	integer dnplsoff
	real*8 dnpls(1)
        integer ddmchoff
        real*8 ddmch(1)
        integer bufoff
        real*8 buf(1)
        integer npmsavoff
        integer npmsav(1)
        integer ksavoff
        integer ksav(1)

        common /array/
     +    dnpls,   ddmch,   buf,   npmsav,   ksav,
     +    dnplsoff,ddmchoff,bufoff,npmsavoff,ksavoff

