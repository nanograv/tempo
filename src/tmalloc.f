	subroutine tmalloc(nptsmax,nbuf)

        implicit none

	integer nptsmax, nbuf

	include 'array.h'

        call mallocx(dnpls,nptsmax,8,dnplsoff)
        call mallocx(ddmch,nptsmax,8,ddmchoff)
        call mallocx(ksav,nptsmax,4,ksavoff)
        call mallocx(npmsav,nptsmax,4,npmsavoff)
        call mallocx(buf,nbuf,8,bufoff)

        return
        end
	

C----------------------------------------------------------------

	subroutine tfree()

        implicit none
        include 'array.h'

        call freex(dnpls)
        call freex(ddmch)
        call freex(ksav)
        call freex(npmsav)
        call freex(buf)

        return
        end
