	subroutine tmalloc(nptsmax,nbuf)

        implicit none

	integer nptsmax, nbuf
        integer mallocx         ! external function in cutil.c

	include 'array.h'

        dnplsptr =  mallocx(dnpls,nptsmax,8,dnplsoff)
        ddmchptr =  mallocx(ddmch,nptsmax,8,ddmchoff)
        ksavptr =  mallocx(ksav,nptsmax,4,ksavoff)
        npmsavptr =  mallocx(npmsav,nptsmax,4,npmsavoff)
        bufptr =  mallocx(buf,nbuf,8,bufoff)

        return
        end
	

C----------------------------------------------------------------

	subroutine tfree()

        implicit none
        include 'array.h'

        call freex(dnplsptr)
        call freex(ddmchptr)
        call freex(ksavptr)
        call freex(npmsavptr)
        call freex(bufptr)

        return
        end
