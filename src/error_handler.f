c      $Id$
	integer function error_handler (signal, code, sigcontext)
	integer signal, code, sigcontext(5)
	character*16 out
c	Give error handling possibilities for all FPE errors.
c	Leave it up to calling statement for enabling!
c	See manual page under f77_floatingpoint and related items.
	logical divide,over,invalid,under,inexact
	out = ''
	under  = loc(code).eq.204
	inexact = loc(code).eq.196
	over  = loc(code).eq.212
	divide = loc(code).eq.200
	if (under) then
		out = 'underflow'
	else if (inexact) then
		out = 'inexact'
	elseif (over) then
		out = 'overflow'
	else if (divide) then
		out = 'division'
	else
		invalid = .true.
		out = 'invalid'
	end if
c	write out for every intercepted floating point error, let setup in
c	main programme decide which to intercept.
c	Write your own error handling statements here.
c	'dump' is my own subroutine which prints out selected info.
c	abort is a system call that causes a core dump.
	write(0,
     +    '(''**** error: IEEE Exception '',a,'' at PC '',i12)')
     +    out,sigcontext(4)
!	write(0,
!	1 '(''**** error: IEEE Exception '',a,'' at PC '',i12)')
!	1 out,sigcontext(4)
!	call dump('error_handler',13)
!	call abort
	end  

	integer function loc(n)
C  Dummy function
	loc=n
	return
	end
