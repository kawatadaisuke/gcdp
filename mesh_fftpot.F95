#include "gcdp.def"
! *********************************************
!  mesh_fftpot.F95 for GCD+
!  16 Aug. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
subroutine  mesh_fftpot()
      use gcdp_const
      use gcdp_system
      use gcdp_pm
#ifdef FFTW3
      use fftw3
#endif

      implicit none
      include 'mpif.h'

      integer i,j,k
! *** for test file *** 
      character fileo*60


! initialisation
      do k=1,nz_m
        do j=1,ny_m
          do i=1,2*(nx_m/2+1)
            rsd_fftw(i,j,k)=0.0d0
          enddo
        enddo
      enddo
! set real space data 
      do k=1,nz_m
        do j=1,ny_m
          do i=1,nx_m
            rsd_fftw(i,j,k)=rho_m(i-1,j-1,k-1)
          enddo
        enddo
      enddo

! forward FFT
      call fftw_execute_dft_r2c(planf_fftw,rsd_fftw,ksd_fftw)

! test output
      write(fileo,'(a3,i3.3)') 'ksd',myrank
      open(60,file=fileo,status='unknown')
      do k=1,nz_m
        do j=1,ny_m
          do i=1,nx_m/2+1
            write(60,'(3I10,2(1pE13.5))') i,j,k,ksd_fftw(i,j,k)*dv_m
          enddo
        enddo
      enddo
      close(60)

! backward FFT
      call fftw_execute_dft_c2r(planb_fftw,ksd_fftw,rsd_fftw)

! test output
      write(fileo,'(a3,i3.3)') 'rsd',myrank
      open(60,file=fileo,status='unknown')
      do k=1,nz_m
        do j=1,ny_m
          do i=1,nx_m
            write(60,'(3I10,(1pE13.5))') i,j,k,rsd_fftw(i,j,k)/dble(nt_m)
          enddo
        enddo
      enddo
      close(60)
      stop


end subroutine
#endif





