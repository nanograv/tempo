c      $Id$
       character tzdir*80,tzfile*80,tztfile*80,tzrsite*1,tzsite*1
       character name(NTZMAX)*12,pname*12,params(6)*80
       integer nsp(NTZMAX),nco(NTZMAX)
       real*8 mxha(NTZMAX)
       integer ntmjd(800)
       real*8 ftmjd(800),tmin(800),tzof(NTZMAX)
       real*8 ftzrmjd
       integer ntzrmjd

       common/tz/tzrfrq,ftzrmjd,ftmjd,tmin,tzof,mxha,
     +    ntzref,nsp,nco,nsets,params,nsite,ntzrmjd,ntmjd
       common/tzch/tzdir,tzfile,tztfile,tzrsite,tzsite,name,pname
