#include "gcdp.def"

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
