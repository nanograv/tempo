c     $Id$
      subroutine equ2ecl(x)

c     rotate vector x(1..3) from equatorial to ecliptic coordinates
c     ==> rotate about "x" axis by angle epsilon

      implicit none

      real*8 x(3)

      real*8 tmpy, tmpz

      real*8 CE, SE

c     note:  
c     epsilon = 84381.412 arc sec at epoch of J2000.0 (IERS tech note 21, p. 19)
c     define:
c     ce = cos(epsilon), se = sin(epsilon)

      parameter (CE=0.91748213149438d0)
      parameter (SE=0.39777699580108d0)


c     n.b. x(1) remains unchanged

      tmpy = x(2)
      tmpz = x(3)

      x(2) = CE*tmpy + SE*tmpz
      x(3) = -SE*tmpy + CE*tmpz 

      return
      end

      

