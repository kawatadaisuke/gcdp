#include "gcdp.def"

subroutine phcurven(id0,id1,is,ic,ixp,iyp,izp,level)
      implicit none

      integer,intent(in) :: id0,id1,level
      integer,intent(out) :: is(id0:id1+1),ic(id0:id1+1)
      integer,intent(in) :: ixp(id0:id1+1),iyp(id0:id1+1),izp(id0:id1+1)

! *** s: state (or kind) of curve segment
! *** c: index in each segment from 0 to 7
      integer i,icp(id0:id1+1)
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
      do i=id0,id1
        if(iyp(i).eq.0) then
          if(ixp(i).eq.0) then
            if(izp(i).eq.0) then
              icp(i)=0
            else
              icp(i)=1
            endif
          else 
            if(izp(i).eq.0) then
              icp(i)=3   
            else
              icp(i)=2
            endif
          endif   
        else
          if(ixp(i).eq.0) then
            if(izp(i).eq.0) then
              icp(i)=7
            else
              icp(i)=6
            endif
          else 
            if(izp(i).eq.0) then
              icp(i)=4   
            else
              icp(i)=5
            endif
          endif   
        endif
      enddo
      if(level.eq.1) then
        do i=id0,id1
          is(i)=0
          ic(i)=icp(i)
        enddo
      else
! *** decide state -> is ***
        do i=id0,id1
          if(is(i).eq.0) then
            if(ic(i).eq.0) then
              is(i)=1
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=2
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=3
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=4
            else
              is(i)=5
            endif
          else if(is(i).eq.1) then
            if(ic(i).eq.0) then
              is(i)=2
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=0
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=6
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=7
            else
              is(i)=8
            endif
          else if(is(i).eq.2) then
            if(ic(i).eq.0) then
              is(i)=0
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=1
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=9
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=10
            else
              is(i)=11
            endif
          else if(is(i).eq.3) then
            if(ic(i).eq.0) then
              is(i)=10
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=8
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=0
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=9
            else
              is(i)=6
            endif
          else if(is(i).eq.4) then
            if(ic(i).eq.0) then
              is(i)=11
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=6
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=8
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=5
            else
              is(i)=0
            endif
          else if(is(i).eq.5) then
            if(ic(i).eq.0) then
              is(i)=9
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=7
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=10
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=0
            else
              is(i)=4
            endif
          else if(is(i).eq.6) then
            if(ic(i).eq.0) then
              is(i)=4
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=11
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=1
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=3
            else
              is(i)=9
            endif
          else if(is(i).eq.7) then
            if(ic(i).eq.0) then
              is(i)=5
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=9
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=11
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=8
            else
              is(i)=1
            endif
          else if(is(i).eq.8) then
            if(ic(i).eq.0) then
              is(i)=3
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=10
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=4
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=1
            else
              is(i)=7
            endif
          else if(is(i).eq.9) then
            if(ic(i).eq.0) then
              is(i)=7
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=5
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=2
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=6
            else
              is(i)=3
            endif
          else if(is(i).eq.10) then
            if(ic(i).eq.0) then
              is(i)=8
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=3
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=5
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=11
            else
              is(i)=2
            endif
          else if(is(i).eq.11) then
            if(ic(i).eq.0) then
              is(i)=6
            else if(ic(i).eq.1.or.ic(i).eq.2) then
              is(i)=4
            else if(ic(i).eq.3.or.ic(i).eq.4) then
              is(i)=7
            else if(ic(i).eq.5.or.ic(i).eq.6) then
              is(i)=2
            else
              is(i)=10
            endif
          endif
! set ic(i)
          ic(i)=cstate(icp(i),is(i))
        enddo
      endif

end subroutine
