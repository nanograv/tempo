c      $Id$
      logical function bigendian()
      integer nbig
      data nbig/'big '/
      save

      if(nbig.eq.1651074848) then
         bigendian = .true.
      else if(nbig.eq.543648098) then
         bigendian=.false.
      else
         print*,'Unknown machine architecture, nbig=',nbig
      endif
      return
      end
