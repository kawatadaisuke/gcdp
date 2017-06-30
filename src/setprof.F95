#include "gcdp.def"
! *****************************************************
!    setprof.F95 for gcd+ 
!  22 June 2017  written by D.KAWATA
! ***************************************************** 


! *** exp disk profile as a function of x ***
#ifdef ADD_EXPDISK

! set M(<x) exp disk cummurative mass profile 
subroutine set_mxexpd(hd,xmin,xmax)
      use gcdp_const
      use gcdp_expdisk
      use gcdp_system

      implicit none
      include 'mpif.h'

      double precision,intent(in) :: hd,xmin,xmax
      integer i
      double precision mtot

      double precision mxexpdfunc
      external mxexpdfunc
  
      hd_expd=hd
      lxmin_expd=dlog(xmin)
      dlx_expd=(dlog(xmax)-lxmin_expd)/dble(NEXPDTAB-1)
      do i=0,NEXPDTAB-1
        lx_expd(i)=lxmin_expd+dlx_expd*dble(i)
        mx_expd(i)=mxexpdfunc(10.0d0**lx_expd(i),hd_expd)
      enddo
! normalise
      mtot=mxexpdfunc(xmax,hd_expd)
      do i=0,NEXPDTAB-1
        mx_expd(i)=mx_expd(i)/mtot
      enddo
! for check
      if(myrank.eq.0) then
        open(60,file='mxexpd.asc',status='unknown')
        do i=0,NEXPDTAB-1
          write(60,'(2(1pE13.5))') 10.0d0**(lx_expd(i)),mx_expd(i)
        enddo
        close(60)
      endif

end subroutine

!  M(<x)  assuming thin disk *****
double precision function mxexpdfunc(x,hd)
      implicit none
      double precision,intent(in) :: x,hd

      mxexpdfunc = (1.0d0-(1.0d0+x/hd)*dexp(-x/hd))
end function

#endif
