c      $Id$
	subroutine tzinit(obsyfile,sitelng,num)

	implicit real*8 (a-h,o-z)
	character path*96,obsyfile*80
	character line*80,item*40
	real*8 maxhadef

	include 'dim.h'
	include 'tz.h'
        include 'acom.h'

	kd=index(tzdir,' ')-1
	path=tzdir(1:kd)//tzfile
	open(11,file=path,status='old',err=90)
	go to 92
 90	write(*,'(''Failed to open '',a)')path
	stop

 92	open(13,file='polyco.dat',status='unknown')

c Read default values of asite,maxha,nspan,ncoeff & freq from 1st line of tz.in
c Free format, but must be in order
	read(11,'(a)')line
	jn=1
	ll=80
	call citem(line,ll,jn,item,jl)                  ! asite
	if(jl.ne.1)then
	   write(*,'(''Invalid asite in tz.in: '',a)')item(1:jl)
	   stop
	endif
	tzsite=item(1:1)
	if(tzsite.ge.'0'.and.tzsite.le.'9')then
	   nsite=ichar(tzsite)-48
	else if(tzsite.ge.'a'.and.tzsite.le.'z')then
	   nsite = ichar(tzsite)-87
	else if (tzsite.eq.'@') then
	   nsite = -1
	else
	   write (*,*) 'Don''t understand observing site ''',tzsite,''''
	   stop
	endif

	call citem(line,ll,jn,item,jl) ! maxha
	read(item,*)maxhadef
	
	call citem(line,ll,jn,item,jl) ! nspan
	read(item,*)nspdef
	
	call citem(line,ll,jn,item,jl) ! ncoeff
	read(item,*)ncodef

	call citem(line,ll,jn,item,jl) ! freq
	if(jl.eq.0)then
	   freqdef=-1.		! No freq entered, assume tzref freq
	else 
	   read(item,*)freqdef
	endif
	
C Read pulsar list and nspan, ncoeff, maxha and freq overrides
	read(11,*)
	read(11,*)
        if (.not.quiet) then
  	  write(*,'(''TZ source list for site = '',a/)')tzsite
	  write(*,'(''    PSR        Nspan  Ncoeffs  Maxha    Freq'')')
	  write(*,'(''----------------------------------------------'')')
        endif
	
	do 310 i=1,ntzmax

	  read(11,'(a)',end=320)line

	  jn=1
	  name(i)=' '
	  call citem(line,ll,jn,item,jl)
	  if(jl.gt.12)jl=12
	  name(i)(1:jl)=item(1:jl)

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     nsp(i)=nspdef
	     nco(i)=ncodef
	     mxha(i)=maxhadef
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)nsp(i)
	  endif

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     nco(i)=ncodef
	     mxha(i)=maxhadef
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)nco(i)
	  endif

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     mxha(i)=maxhadef
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)mxha(i)
	  endif

	  call citem(line,ll,jn,item,jl)
	  if(jl.eq.0)then
	     tzof(i)=freqdef
	     go to 312
	  else
	     read(item,*)tzof(i)
	  endif

 312      continue
          if (.not.quiet) then
    	    if(tzof(i).gt.0.)then
              write(*,'(1x,a,2i8,x,f7.2,1x,f12.5)') name(i),nsp(i),
     +           nco(i), mxha(i),tzof(i)
            else
              write(*,'(1x,a,3i8,2x,''from tzref'')') 
     +           name(i),nsp(i),nco(i)
            endif
          endif
 310	continue
	print*,' Sorry, ',ntzmax,' pulsar limit'
 320	num=i-1
	close(11)

	open(11,file=obsyfile,status='old')
C  Allow for xyz observatory coordinates  djn 17 aug 92
	do 330 i=1,36
	  read(11,1020,end=340) alat, along, icoord
 1020	  format(2f13.0,21x,i1)
	  if (i.eq.nsite) then
	    if (icoord.eq.0) then
	      sitelng=ang(1,along)
	    else
	      sitelng = datan2(-along,alat)
	    endif
	  endif
 330	continue
 340	if(nsite.ge.i) then
	  print *,' Site',nsite,' not found in ',obsyfile
	  stop
	end if
	close(11)

	return
	end

