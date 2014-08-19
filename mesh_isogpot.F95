#include "gcdp.def"
! *********************************************
!  mesh_isogpot.F95 for GCD+
!  31 Jul. 2014   written by D. Kawata
! *********************************************

! getting force at the mesh position

#ifdef TREEPM
#ifdef ISOPM
subroutine  mesh_isogpot()
      use gcdp_const
      use gcdp_system
      use gcdp_pm
#ifdef FFTW3
      use fftw3
#endif

      implicit none
      include 'mpif.h'

      integer i,j,k,i1
      integer ip,jp,kp,nxm2,nym2,nzm2
      integer(C_INTPTR_T) :: ifftw,jfftw,kfftw
! *** for test file *** 
      character fileo*60
      integer ierr

! initialization
      do k=1,lnz_m
        do j=1,ny_m
          do i=1,2*(nx_m/2+1)
            rsd_fftw(i,j,k)=0.0d0
          enddo
        enddo
      enddo
! set real space 1/r kernel (eq. 6-113) of Hockney & Eastwood
      nxm2=nx_m/2
      nym2=ny_m/2
      nzm2=nz_m/2
      do kfftw=1,lnz_m
        kp=kfftw+lzoff_mp(myrank)-1
        if(kp.gt.nzm2) then
          kp=nz_m-kp
        endif
        do jfftw=1,ny_m
          jp=jfftw-1
          if(jp.gt.nym2) then
            jp=ny_m-jp
          endif
          do ifftw=1,nx_m
            ip=ifftw-1
            if(ip.gt.nxm2) then
              ip=nx_m-ip
            endif
            if(ip.eq.0.and.(jp.eq.0.and.kp.eq.0)) then
! no difference in one point source test, but the latter makes more sense
!              rsdgpot_fftw(ifftw,jfftw,kfftw)=dv_m &
!                /dsqrt(dx_m**2+dy_m**2+dz_m**2)
              rsdgpot_fftw(ifftw,jfftw,kfftw)=dv_m
! should not be 1
!              rsdgpot_fftw(ifftw,jfftw,kfftw)=1.0d0
            else 
! factor dx_m*dy_m*dz_m (dv_m set in mesh_gen)for density -> mass
! factor 1.0d0/dsqrt(dx_m**2+dy_m**2+dz_m**2) for 1/r
              rsdgpot_fftw(ifftw,jfftw,kfftw)=dv_m &
               /(dsqrt((dble(ip)*dx_m)**2+(dble(jp)*dy_m)**2 &
               +(dble(kp)*dz_m)**2))
            endif
          enddo
        enddo
      enddo

! test output
!      write(fileo,'(a7,i3.3)') 'isogrho',myrank
!      open(60,file=fileo,status='unknown')
!      do k=0,lnz_m-1
!        do j=0,ny_m-1
!          do i=0,nx_m-1
!            write(60,'(3I10,(1pE13.5))') i,j,k+lzoff_m,rsdgpot_fftw(i+1,j+1,k+1)
!          enddo
!        enddo
!      enddo
!      close(60)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      stop      


! forward FFT
      call fftw_mpi_execute_dft_r2c(planfgpot_fftw,rsdgpot_fftw,ksdgpot_fftw)

! test output
!      write(fileo,'(a4,i3.3)') 'gpot',myrank
!      open(60,file=fileo,status='unknown')
!      do jfftw=1,lny_fftw
!        do kfftw=1,nz_fftw
!          do ifftw=1,nx_m/2+1
!            write(60,'(3I10,2(1pE13.5))') ifftw,kfftw,jfftw &
!              ,ksdgpot_fftw(ifftw,kfftw,jfftw)
!          enddo
!        enddo
!      enddo
!      close(60)

end subroutine
#endif
#endif





