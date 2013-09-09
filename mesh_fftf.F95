#include "gcdp.def"
! *********************************************
!  mesh_fftf.F95 for GCD+
!  9 Sep. 2013   written by D. Kawata
! *********************************************

! getting force at the mesh position

#ifdef TREEPM
subroutine  mesh_fftf()
      use gcdp_const
      use gcdp_system
      use gcdp_pm
#ifdef FFTW3
      use fftw3
#endif

      implicit none
      include 'mpif.h'

      integer i,j,k,i1
      integer ip2,ip1,im1,im2
      integer jp2,jp1,jm1,jm2
      integer kp2,kp1,km1,km2
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
! normalised (has to be done)
            pot_m(i-1,j-1,k-1)=rsd_fftw(i,j,k)/dble(nt_m)
          enddo
        enddo
      enddo

! set force
      do k=0,nz_m-1
        kp2=k+2
        if(kp2.gt.nz_m-1) then
          kp2=kp2-nz_m
        endif
        kp1=k+1
        if(kp1.gt.nz_m-1) then
          kp1=kp1-nz_m
        endif
        km2=k-2
        if(km2.lt.0) then
          km2=km2+nz_m
        endif
        km1=k-1
        if(km1.lt.0) then
          km1=km1+nz_m
        endif
        do j=0,ny_m-1
          jp2=j+2
          if(jp2.gt.ny_m-1) then
            jp2=jp2-ny_m
          endif
          jp1=j+1
          if(jp1.gt.ny_m-1) then
            jp1=jp1-ny_m
          endif
          jm2=j-2
          if(jm2.lt.0) then
            jm2=jm2+ny_m
          endif
          jm1=j-1
          if(jm1.lt.0) then
            jm1=jm1+ny_m
          endif
          do i=0,nx_m-1
            ip2=i+2
            if(ip2.gt.nx_m-1) then
              ip2=ip2-nx_m
            endif
            ip1=i+1
            if(ip1.gt.nx_m-1) then
              ip1=ip1-nx_m
            endif
            im2=i-2
            if(im2.lt.0) then
              im2=im2+nx_m
            endif
            im1=i-1
            if(im1.lt.0) then
              im1=im1+nx_m
            endif

            fx_m(i,j,k)=(1.0d0/dx_m) &
             *((2.0d0/3.0d0)*(pot_m(ip1,j,k)-pot_m(im1,j,k)) &
             -((1.0d0/12.0d0)*(pot_m(ip2,j,k)-pot_m(im2,j,k))))
            fy_m(i,j,k)=(1.0d0/dy_m) &
             *((2.0d0/3.0d0)*(pot_m(i,jp1,k)-pot_m(i,jm1,k)) &
             -((1.0d0/12.0d0)*(pot_m(i,jp2,k)-pot_m(i,jm2,k))))
            fz_m(i,j,k)=(1.0d0/dz_m) &
             *((2.0d0/3.0d0)*(pot_m(i,j,kp1)-pot_m(i,j,km1)) &
             -((1.0d0/12.0d0)*(pot_m(i,j,kp2)-pot_m(i,j,km2))))
     
          enddo
        enddo
      enddo

     write(fileo,'(a4,i3.3)') 'fpot',myrank
      open(60,file=fileo,status='unknown')
      do k=0,nz_m-1
        do j=0,ny_m-1
          do i=0,nx_m-1
            write(60,'(8(1pE13.5))') x_m(i,j,k),y_m(i,j,k),z_m(i,j,k) &
             ,pot_m(i,j,k)/dble(nt_m) &
!             ,pot_m(i,j,k) &
             ,fx_m(i,j,k),fy_m(i,j,k),fz_m(i,j,k) &
             ,dsqrt(x_m(i,j,k)**2+y_m(i,j,k)**2+z_m(i,j,k)**2)
          enddo
        enddo
      enddo
      close(60)


end subroutine
#endif





