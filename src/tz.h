c      $Id$
       character tzdir*80,tzfile*80,tztfile*80,tzsite*80
       character tzrsite*80
       character tzsitedef*1
       character polycofile*160
       character name(NTZMAX)*12,pname*12,params(6)*80
       integer ntzsite,ntzrsite
       integer nsp(NTZMAX),nco(NTZMAX)
       real*8 mxha(NTZMAX)
       integer ntmjd(NTZARR)
       real*8 ftmjd(NTZARR),tmin(NTZARR),tzof(NTZMAX)
       real*8 ftzrmjd
       integer ntzrmjd
       real*8 tzmjdstart
       integer lupolyco

c      non-character variables:
       common/tz/tzrfrq,ftzrmjd,ftmjd,tmin,tzof,mxha,
     +    ntzref,nsp,nco,nsets,params,nsite,
     +    ntzrmjd,ntmjd,
     +    tzmjdstart,lupolyco,ntzsite,ntzrsite
c      character variables:
       common/tzch/tzdir,tzfile,tztfile,tzrsite,tzsite,name,pname,
     +    tzsitedef,polycofile
