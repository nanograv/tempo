c      $Id$
c
      integer function sitea2n(s)

C     convert first character of ascii string s into a numerical
C     observatory code

      implicit none

      character*(*) s

      if(s(1:1).ge.'0'.and.s(1:1).le.'9')then
        sitea2n=ichar(s(1:1))-48
      else if(s(1:1).ge.'a'.and.s(1:1).le.'z')then
        sitea2n = ichar(s(1:1))-87
      else if (s(1:1).eq.'@') then
        sitea2n = -1
      else
        sitea2n = -2
      endif
      
      return
      
      end
