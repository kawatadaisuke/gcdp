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

      integer i,j,k,i1
      double precision kx,ky,kz,kabs2
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

! For normalisation see eq. A-9 and A-10 of Hockney & Eastwood

! test output
      write(fileo,'(a3,i3.3)') 'ksd',myrank
      open(60,file=fileo,status='unknown')

      do k=1,nz_m
        if(k.le.nz2_m+1) then
          kz=dkz_m*dble(k-1)
        else
          kz=-dkz_m*dble(nz_m-k+1)
        endif

!        write(6,*) ' k,kz=',k,kz,dkz_m,nz_m,nz2_m

        do j=1,ny_m
          if(j.le.ny2_m+1) then
            ky=dky_m*dble(j-1)
          else
            ky=-dky_m*dble(ny_m-j+1)
          endif
          if(k.eq.1.and.j.eq.1) then
! case of k=0, set power = 0 
! eq. (4a) of Efstathiou et al. eq. (4a)
            i1=2
            i=1
            ksd_fftw(i,j,k)=0.0d0*ksd_fftw(i,j,k)
          else
            i1=1
          endif
          do i=i1,nx_m/2+1
            kx=dkx_m*dble(i-1)
            kabs2=kx**2+ky**2+kz**2
! x Green function
            ksd_fftw(i,j,k)=4.0d0*M_PI*ksd_fftw(i,j,k)/kabs2            
          enddo

          do i=i1,nx_m/2+1
            kx=dkx_m*dble(i-1)
            kabs2=kx**2+ky**2+kz**2

            write(60,'(3I10,6(1pE13.5))') i,j,k,ksd_fftw(i,j,k)*dv_m &
            ,kx,ky,kz,kabs2
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
            write(60,'(5(1pE13.5))') x_m(i-1,j-1,k-1),y_m(i-1,j-1,k-1) &
             ,z_m(i-1,j-1,k-1) &
             ,rsd_fftw(i,j,k)/dble(nt_m) &
             ,dsqrt(x_m(i-1,j-1,k-1)**2+y_m(i-1,j-1,k-1)**2 &
              +z_m(i-1,j-1,k-1)**2)
          enddo
        enddo
      enddo
      close(60)
!      stop

! set real space potential data 
      do k=1,nz_m
        do j=1,ny_m
          do i=1,nx_m
            pot_m(i-1,j-1,k-1)=rsd_fftw(i,j,k)
          enddo
        enddo
      enddo

end subroutine
#endif





