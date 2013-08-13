#include "gcdp.def"
! *********************************************
!  mesh_gen.F95 for GCD+
!  13  Aug. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
subroutine  mesh_gen()
      use gcdp_const
      use gcdp_system
      use gcdp_pm

      implicit none
!      include 'mpif.h'

      integer i,j,k
      double precision xp,yp,zp,xp0,yp0,zp0
!      character filen*60

! allocate memory space
      call allocate_pm()

! calculate grid space
      dx_m=SI_lbox(0)/dble(nx_m)
      dy_m=SI_lbox(1)/dble(ny_m)
      dz_m=SI_lbox(2)/dble(nz_m)
! first grid points
      xp0=0.5d0*SI_lbox(0)+0.5d0*dx_m
      yp0=0.5d0*SI_lbox(1)+0.5d0*dy_m
      zp0=0.5d0*SI_lbox(2)+0.5d0*dz_m

! mesh points
      do k=0,nz_m-1
        zp=zp0+dz_m*dble(k)
        do j=0,ny_m-1
          yp=yp0+dy_m*dble(j)
          do i=0,nx_m-1
            x_m(i,j,k)=xp0+dx_m*dble(i)
            y_m(i,j,k)=yp0
            z_m(i,j,k)=zp0
          enddo
        enddo
      enddo

end subroutine
#endif

