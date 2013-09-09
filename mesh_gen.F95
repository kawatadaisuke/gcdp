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
      include 'mpif.h'

      integer i,j,k
      double precision yp,zp
!      character filen*60

! allocate memory space
      call allocate_pm()

! calculate grid space
      dx_m=SI_lbox(0)/dble(nx_m)
      dy_m=SI_lbox(1)/dble(ny_m)
      dz_m=SI_lbox(2)/dble(nz_m)
! calculate grid space in k space, 2 pi/L seems better than 1/L
      dkx_m=2.0d0*M_PI/(SI_lbox(0))
      dky_m=2.0d0*M_PI/(SI_lbox(1))
      dkz_m=2.0d0*M_PI/(SI_lbox(2))

! first grid points
      x0_m=-0.5d0*SI_lbox(0)+0.5d0*dx_m
      y0_m=-0.5d0*SI_lbox(1)+0.5d0*dy_m
      z0_m=-0.5d0*SI_lbox(2)+0.5d0*dz_m
! half nx,ny,nz
      nx2_m=nx_m/2
      ny2_m=ny_m/2
      nz2_m=nz_m/2

      if(myrank.eq.0) then
        write(6,*) ' Mesh info'
        write(6,*) '   nx_m,ny_m,nz_m=',nx_m,ny_m,nz_m
        write(6,*) '   dx,dy,dz=',dx_m,dy_m,dz_m
        write(6,*) '   x0,y0,z0=',x0_m,y0_m,z0_m
      endif

! total number of mesh 
      nt_m=nx_m*ny_m*nz_m
! total box size
      vt_m=SI_lbox(0)*SI_lbox(1)*SI_lbox(2)
! mesh size
      dv_m=dx_m*dy_m*dz_m

! mesh points 1d
      do i=0,nx_m-1
        x_m1d(i)=x0_m+dx_m*dble(i)
      enddo
      do i=0,ny_m-1
        y_m1d(i)=y0_m+dy_m*dble(i)
      enddo
      do i=0,nz_m-1
        z_m1d(i)=z0_m+dz_m*dble(i)
      enddo
! mesh points 3d
      do k=0,nz_m-1
        zp=z_m1d(k)
        do j=0,ny_m-1
          yp=y_m1d(j)
          do i=0,nx_m-1
            x_m(i,j,k)=z_m1d(i)
            y_m(i,j,k)=yp
            z_m(i,j,k)=zp
          enddo
        enddo
      enddo

end subroutine
#endif

