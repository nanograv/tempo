c      $Id$
       character tzdir*80,tzfile*80,tztfile*80,tzrsite*1,tzsite*1
       character name(NTZMAX)*12,pname*12,params(6)*80
       integer nsp(NTZMAX),nco(NTZMAX),mxha(NTZMAX)
       real*8 tfmjd(800),tmin(800),tzof(NTZMAX)

       common/tz/tzrfrq,tzrmjd,tfmjd,tmin,tzof,ntzref,nsp,nco,mxha,
     +    nsets,params,nsite
       common/tzch/tzdir,tzfile,tztfile,tzrsite,tzsite,name,pname
