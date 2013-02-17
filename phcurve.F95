#include "gcdp.def"

subroutine phcurve(is,ic,ixp,iyp,izp,level)
      implicit none

      integer,intent(out) :: is,ic
      integer,intent(in) :: ixp,iyp,izp,level
! *** s: state (or kind) of curve segment
! *** c: index in each segment from 0 to 7
      integer icp
      integer cstate(0:7,0:11)
! *** define the order in each PH state in level 0 position ***
      data cstate(0:7,0)/0,1,2,3,4,5,6,7/ &
         ,cstate(0:7,1)/0,7,6,1,2,5,4,3/ &
         ,cstate(0:7,2)/0,3,4,7,6,5,2,1/ &
         ,cstate(0:7,3)/2,3,0,1,6,7,4,5/ &
         ,cstate(0:7,4)/6,5,2,1,0,3,4,7/ &
         ,cstate(0:7,5)/4,3,2,5,6,1,0,7/ &
         ,cstate(0:7,6)/2,5,4,3,0,7,6,1/ &
         ,cstate(0:7,7)/6,7,4,5,2,3,0,1/ &
         ,cstate(0:7,8)/4,7,0,3,2,1,6,5/ &
         ,cstate(0:7,9)/2,1,6,5,4,7,0,3/ &
         ,cstate(0:7,10)/6,1,0,7,4,3,2,5/ &
         ,cstate(0:7,11)/4,5,6,7,0,1,2,3/

! *** position in level 0 state ***
      if(iyp.eq.0) then
        if(ixp.eq.0) then
          if(izp.eq.0) then
            icp=0
          else
            icp=1
          endif
        else 
          if(izp.eq.0) then
            icp=3   
          else
            icp=2
          endif
        endif   
      else
        if(ixp.eq.0) then
          if(izp.eq.0) then
            icp=7
          else
            icp=6
          endif
        else 
          if(izp.eq.0) then
            icp=4   
          else
            icp=5
          endif
        endif   
      endif

      if(level.eq.1) then
        is=0
        ic=icp
      else
! *** decide state -> is ***
        if(is.eq.0) then
          if(ic.eq.0) then
            is=1
          else if(ic.eq.1.or.ic.eq.2) then
            is=2
          else if(ic.eq.3.or.ic.eq.4) then
            is=3
          else if(ic.eq.5.or.ic.eq.6) then
            is=4
          else
             is=5
          endif
        else if(is.eq.1) then
          if(ic.eq.0) then
            is=2
          else if(ic.eq.1.or.ic.eq.2) then
            is=0
          else if(ic.eq.3.or.ic.eq.4) then
            is=6
          else if(ic.eq.5.or.ic.eq.6) then
            is=7
          else
            is=8
          endif
        else if(is.eq.2) then
          if(ic.eq.0) then
            is=0
          else if(ic.eq.1.or.ic.eq.2) then
            is=1
          else if(ic.eq.3.or.ic.eq.4) then
            is=9
          else if(ic.eq.5.or.ic.eq.6) then
            is=10
          else
            is=11
          endif
        else if(is.eq.3) then
          if(ic.eq.0) then
            is=10
          else if(ic.eq.1.or.ic.eq.2) then
            is=8
          else if(ic.eq.3.or.ic.eq.4) then
            is=0
          else if(ic.eq.5.or.ic.eq.6) then
            is=9
          else
            is=6
          endif
        else if(is.eq.4) then
          if(ic.eq.0) then
            is=11
          else if(ic.eq.1.or.ic.eq.2) then
            is=6
          else if(ic.eq.3.or.ic.eq.4) then
            is=8
          else if(ic.eq.5.or.ic.eq.6) then
            is=5
          else
            is=0
          endif
        else if(is.eq.5) then
          if(ic.eq.0) then
            is=9
          else if(ic.eq.1.or.ic.eq.2) then
            is=7
          else if(ic.eq.3.or.ic.eq.4) then
            is=10
          else if(ic.eq.5.or.ic.eq.6) then
            is=0
          else
            is=4
          endif
        else if(is.eq.6) then
          if(ic.eq.0) then
            is=4
          else if(ic.eq.1.or.ic.eq.2) then
            is=11
          else if(ic.eq.3.or.ic.eq.4) then
            is=1
          else if(ic.eq.5.or.ic.eq.6) then
            is=3
          else
            is=9
          endif
        else if(is.eq.7) then
          if(ic.eq.0) then
            is=5
          else if(ic.eq.1.or.ic.eq.2) then
            is=9
          else if(ic.eq.3.or.ic.eq.4) then
            is=11
          else if(ic.eq.5.or.ic.eq.6) then
            is=8
          else
            is=1
          endif
        else if(is.eq.8) then
          if(ic.eq.0) then
            is=3
          else if(ic.eq.1.or.ic.eq.2) then
            is=10
          else if(ic.eq.3.or.ic.eq.4) then
            is=4
          else if(ic.eq.5.or.ic.eq.6) then
            is=1
          else
            is=7
          endif
        else if(is.eq.9) then
          if(ic.eq.0) then
            is=7
          else if(ic.eq.1.or.ic.eq.2) then
            is=5
          else if(ic.eq.3.or.ic.eq.4) then
            is=2
          else if(ic.eq.5.or.ic.eq.6) then
            is=6
          else
            is=3
          endif
        else if(is.eq.10) then
          if(ic.eq.0) then
            is=8
          else if(ic.eq.1.or.ic.eq.2) then
            is=3
          else if(ic.eq.3.or.ic.eq.4) then
            is=5
          else if(ic.eq.5.or.ic.eq.6) then
            is=11
          else
            is=2
          endif
        else if(is.eq.11) then
          if(ic.eq.0) then
            is=6
          else if(ic.eq.1.or.ic.eq.2) then
            is=4
          else if(ic.eq.3.or.ic.eq.4) then
            is=7
          else if(ic.eq.5.or.ic.eq.6) then
            is=2
          else
            is=10
          endif
        endif
! set ic
        ic=cstate(icp,is)
      endif

end subroutine


subroutine phixyzp(is,ic,ixp,iyp,izp,level)

      implicit none

      integer,intent(in) :: is,ic,level
      integer,intent(out) :: ixp,iyp,izp
! *** input is, ic and level and return ixp,iyp,iz
! *** s: state (or kind) of curve segment
! *** c: index in each segment from 0 to 7
      integer icp
      integer icstate(0:7,0:11)
! *** define the order in each PH state in level 0 position ***
      data icstate(0:7,0)/0,1,2,3,4,5,6,7/ &
        ,icstate(0:7,1)/0,3,4,7,6,5,2,1/ &
        ,icstate(0:7,2)/0,7,6,1,2,5,4,3/ &
        ,icstate(0:7,3)/2,3,0,1,6,7,4,5/ &
        ,icstate(0:7,4)/4,3,2,5,6,1,0,7/ &
        ,icstate(0:7,5)/6,5,2,1,0,3,4,7/ &
        ,icstate(0:7,6)/4,7,0,3,2,1,6,5/ &
        ,icstate(0:7,7)/6,7,4,5,2,3,0,1/ &
        ,icstate(0:7,8)/2,5,4,3,0,7,6,1/ &
        ,icstate(0:7,9)/6,1,0,7,4,3,2,5/ &
        ,icstate(0:7,10)/2,1,6,5,4,7,0,3/ &
        ,icstate(0:7,11)/4,5,6,7,0,1,2,3/
      integer ixind(0:7),iyind(0:7),izind(0:7)
! *** ixp,iyp,izp value as a function of icp
      data ixind/0,0,1,1,1,1,0,0/ &
       ,iyind/0,0,0,0,1,1,1,1/ &
       ,izind/0,1,1,0,0,1,1,0/

      if(level.eq.1) then
!        is=0
        icp=ic
      else
        icp=icstate(ic,is)
      endif
      ixp=ixind(icp)
      iyp=iyind(icp)
      izp=izind(icp)

end subroutine
