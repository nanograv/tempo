      subroutine tasc2t0()

C  Converts time of ascending node [t0asc] to time of periastron [t0(1)],
C  (for BT model).

      implicit real*8(a-h,o-z)
      parameter(twopi=6.28318530717958648d0)
      include 'dim.h'
      include 'orbit.h'

      om = omz(1)*twopi/360.d0
      fe = DSQRT((1-e(1))/(1+e(1)))

      uasc=2.d0*DATAN(fe*DTAN(om/2.d0))
      dt=pb(1)/twopi*(uasc-e(1)*DSIN(uasc))

      t0(1)=t0asc+dt

      return
      end
