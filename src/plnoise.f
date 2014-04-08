c      $Id$

c Functions for computing covariance for power-law red noise
c spectra.

	real*8 function plnoise(dt,alpha,f0,amp,init)

c Compute cov for lag dt (in years).  Basically follows
c formula from van Haasteren et al (2011).  Set init to true
c the first time to set up coeffs, etc.  Input lag and freq in 
c years, output in in us^2.

        real*8 dt,alpha,f0,amp
        logical init

        real*8 fl
        data fl/0.01/

        real*8 pi
        data pi/3.141592653589793d0/

        parameter(NC=32)
        real*8 scale0, scale1
        real*8 log_coeffs(NC), signs(NC)
        real*8 ctmp,xtmp,lxtmp,ssum

c TODO replace alpha with "timing" spec idx?

        save pi, fl, scale0, scale1, signs, log_coeffs

        if (init) then
          scale0 = f0**(-2.0*alpha)/(4.0*pi**2)/(fl**(2.0-2.0*alpha))
          scale0 = scale0 * (365.24*86400.0)**2 * 1d12 ! us^2 output
          scale1 = gamma(-2.0+2.0*alpha) * cos(pi*alpha)
          do i=0,NC-1
            if (mod(i,2).eq.1) then
              signs(i+1) = -1.0
            else
              signs(i+1) = 1.0
            endif
            ctmp = 2.0d0*i + 2.0*alpha - 2.0
            if (ctmp.lt.0.0) then
              signs(i+1) = -1.0 * signs(i+1)
              ctmp = abs(ctmp)
            endif
            log_coeffs(i+1) = -log_gamma(2.0d0*i+1.0) - log(ctmp)
            !print *,"LC",i,log_coeffs(i+1)
          enddo
        endif

        if (dt.eq.0.0) then
          plnoise = (amp**2)*scale0*(-signs(1)*exp(log_coeffs(1)))
          return
        endif

        xtmp = abs(2.0d0*pi*fl*dt)
        lxtmp = log(xtmp)
        ssum = 0.0
        do i=0,NC-1
          ssum=ssum+signs(i+1)*exp(2.0d0*i*lxtmp + log_coeffs(i+1))
        enddo
        plnoise = (amp**2)*scale0*(-scale1*(xtmp**(2.0-2.0*alpha))-ssum)

        return
        end
