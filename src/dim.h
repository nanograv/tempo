c      $Id$

c       NPTSDEF  default maximum number of TOAs
c       NPT      maximum number of clock corrections

c       NGLT     maximum number of glitches
c       NGLP     number of parameters per glitch
c       NJUMP    maximum number of jumps
c       NFBMAX   maximum order of orbital frequency taylor expansion, orbit 1
c       NXDOTMAX maximum number of orbital size (x) derivatives, orbit 1
c       NFBJMAX  paximum number of orbital period jumps, orbit 1

c                parameters are in following order:
c       1       to NPAR1  sixty basic parameters
c       NPAR1+1 to NPAR2  glitch parameters (NGLT*NGLP)
c       NPAR2+1 to NPAR3  jump parameters (NJUMP)
c       NPAR3+1 to NPAR4  Orbital frequency & derivatives (1 to NFBMAX)
c       NPAR4+1 to NPAR5  Orbital size (ASINI) derivatives (3 to NXDOTMAX)
c       NPAR5+1 to NPAR6  Orbital period jump parameters
c       NPAR6+1 to NPA    DM offsets (1 to NDMXMAX)
c       NPA      total number of parameters (basic+glitch+jump)
c       NPAP1    total number of parameters plus one
c       NBUFDEF  default size of virtual memory buffer (why 35*NPTSMAX???)
c       NBOOTMAX max size of bootstrap Monte Carlo integraions
c       NCLKMAX  max number of clock correction files (obs to UTC, UTC to xx)
c       NEPHMAX  max number of ephemerides
c       NUTMAX   max number of ut1 corrections

	parameter (NPTSDEF=60000)
	parameter (NPT=10000)
	parameter (NGLT=9,NGLP=5)
	parameter (NJUMP=120)
	parameter (NFBMAX=20)
	parameter (NXDOTMAX=10)
        parameter (NFBJMAX=120)
        parameter (NDMXMAX=120)
        parameter (NPAR1=60)
        parameter (NPAR2=NPAR1+NGLT*NGLP)
	parameter (NPAR3=NPAR2+NJUMP)
	parameter (NPAR4=NPAR3+(NFBMAX))
        parameter (NPAR5=NPAR4+(NXDOTMAX-1))
	parameter (NPAR6=NPAR5+2*NFBJMAX)
        parameter (NPA=NPAR6+NDMXMAX)
        parameter (NPAP1=NPA+1)
	parameter (NPARDEF=28)
        parameter (NBOOTMAX=1024)
        parameter (NTZMAX=1000,NCLKMAX=5,NEPHMAX=5,NUTMAX=3000)
