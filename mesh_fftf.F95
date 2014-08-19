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
      integer kp2,kp1,km1,km2,kp
      double precision kx,ky,kz,kabs2
      double precision l2nx,l2ny,l2nz,kxl2nx,kyl2ny,kzl2nz,dnt_m
! *** for test file *** 
      character fileo*60
      integer ierr
#ifdef FFTW3_MPI
      integer ip,srank,rrank,koff,nsend
      integer istatus(MPI_STATUS_SIZE)
      integer npmp,ireqs,ireqr
      integer(C_INTPTR_T) :: ifftw,jfftw,kfftw
      double precision,allocatable :: tdvs(:),tdvr(:)
#endif

! initialisation
#ifdef FFTW3_MPI
      do k=1,lnz_m
#else
      do k=1,nz_m
#endif
        do j=1,ny_m
          do i=1,2*(nx_m/2+1)
            rsd_fftw(i,j,k)=0.0d0
          enddo
        enddo
      enddo
! for CIC deconvolution
      l2nx=0.5d0*dx_m
      l2ny=0.5d0*dy_m
      l2nz=0.5d0*dy_m

! set real space data 
#ifdef FFTW3_MPI
      do kfftw=1,lnz_fftw
        do jfftw=1,ny_fftw
          do ifftw=1,nx_fftw
            rsd_fftw(ifftw,jfftw,kfftw)=rho_m(ifftw-1,jfftw-1 &
              ,kfftw-1)
#else
      do k=1,nz_m
        do j=1,ny_m
          do i=1,nx_m
            rsd_fftw(i,j,k)=rho_m(i-1,j-1,k-1)
#endif
          enddo
        enddo
      enddo

#ifdef ISOPM
! test output
!      write(fileo,'(a7,i3.3)') 'ksdpfft',myrank
!      open(60,file=fileo,status='unknown')
!      do jfftw=1,lny_fftw
!        do kfftw=1,nz_fftw
!          do ifftw=i1,nx_m/2+1
!            write(60,'(3I10,2(1pE13.5))') ifftw,kfftw,jfftw &
!             ,ksdgpot_fftw(ifftw,kfftw,jfftw)
!          enddo
!        enddo
!      enddo
!      close(60)
#endif


! forward FFT
#ifdef FFTW3_MPI
      call fftw_mpi_execute_dft_r2c(planf_fftw,rsd_fftw,ksd_fftw)
#else
      call fftw_execute_dft_r2c(planf_fftw,rsd_fftw,ksd_fftw)
#endif

! For normalisation see eq. A-9 and A-10 of Hockney & Eastwood

! test output
!      write(fileo,'(a3,i3.3)') 'ksd',myrank
!      open(60,file=fileo,status='unknown')

#ifdef FFTW3_MPI

      do jfftw=1,lny_fftw
        j=jfftw+ljoff_fftw
        if(j.le.ny2_m+1) then
          ky=dky_m*dble(j-1)
        else
          ky=-dky_m*dble(ny_m-j+1)
        endif
        do kfftw=1,nz_fftw
          k=int(kfftw)
          if(k.le.nz2_m+1) then
            kz=dkz_m*dble(k-1)
          else
            kz=-dkz_m*dble(nz_m-k+1)
          endif

          if(k.eq.1.and.j.eq.1) then
! case of k=0, set power = 0 
! eq. (4a) of Efstathiou et al. eq. (4a)
            i1=2
            ifftw=1
            ksd_fftw(ifftw,kfftw,jfftw)=0.0d0
          else
            i1=1
          endif

          do ifftw=i1,nx_m/2+1
            i=int(ifftw)
            kx=dkx_m*dble(i-1)
            kabs2=kx**2+ky**2+kz**2

!            write(60,'(3I10,8(1pE13.5))') ifftw,kfftw,jfftw,ksd_fftw(ifftw,kfftw,jfftw) &
!            ,kx,ky,kz,kabs2,ksdgpot_fftw(ifftw,kfftw,jfftw)

! x Green function
#ifdef ISOPM
! convolve with Fourier component of 1/r kernel
            ksd_fftw(ifftw,kfftw,jfftw)=ksd_fftw(ifftw,kfftw,jfftw) &
             *ksdgpot_fftw(ifftw,kfftw,jfftw)
#else
            ksd_fftw(ifftw,kfftw,jfftw)=4.0d0*M_PI &
              *ksd_fftw(ifftw,kfftw,jfftw)/kabs2            
#endif


! PMTree long range force trunction
#ifndef PMONLY
            ksd_fftw(ifftw,kfftw,jfftw)=ksd_fftw(ifftw,kfftw,jfftw) &
              *dexp(-kabs2*SI_rsp2tpm)
#endif
! deconvolution for CIC
            if((kx.gt.0.0d0.and.ky.gt.0.0d0).and.kz.gt.0.0d0) then
              kxl2nx=kx*l2nx
              kyl2ny=ky*l2ny
              kzl2nz=kz*l2nz
              ksd_fftw(ifftw,kfftw,jfftw)=ksd_fftw(ifftw,kfftw,jfftw) &
               /(((dsin(kxl2nx)/(kxl2nx))**4) &
                *((dsin(kyl2ny)/(kyl2ny))**4) &
                *((dsin(kzl2nz)/(kzl2nz))**4)) 
            endif
#else
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
            ksd_fftw(i,j,k)=0.0d0
          else
            i1=1
          endif
          do i=i1,nx_m/2+1
            kx=dkx_m*dble(i-1)
            kabs2=kx**2+ky**2+kz**2

!            write(60,'(3I10,6(1pE13.5))') i,j,k,ksd_fftw(i,j,k) &
!            ,kx,ky,kz,kabs2

! x Green function
            ksd_fftw(i,j,k)=4.0d0*M_PI*ksd_fftw(i,j,k)/kabs2            
#ifndef PMONLY
! PMTree long range force trunction
            ksd_fftw(i,j,k)=ksd_fftw(i,j,k) &
             *dexp(-kabs2*SI_rsp2tpm)
#endif
! deconvolution for CIC
            if((kx.gt.0.0d0.and.ky.gt.0.0d0).and.kz.gt.0.0d0) then
              kxl2nx=kx*l2nx
              kyl2ny=ky*l2ny
              kzl2nz=kz*l2nz
              ksd_fftw(i,j,k)=ksd_fftw(i,j,k) &
               /(((dsin(kxl2nx)/(kxl2nx))**4) &
                *((dsin(kyl2ny)/(kyl2ny))**4) &
                *((dsin(kzl2nz)/(kzl2nz))**4)) 
            endif
#endif
          enddo

!          do i=i1,nx_m/2+1
!            kx=dkx_m*dble(i-1)
!            kabs2=kx**2+ky**2+kz**2
!
! note: ksd is complex, and two double output
! #ifdef FFTW3_MPI
!            write(60,'(3I10,6(1pE13.5))') i,j,k,ksd_fftw(i,kfftw,jfftw) &
! #else
!            write(60,'(3I10,6(1pE13.5))') i,j,k,ksd_fftw(i,j,k) &
! #endif
!            ,kx,ky,kz,kabs2
!          enddo

        enddo
      enddo

!      close(60)

! backward FFT
#ifdef FFTW3_MPI
      call fftw_mpi_execute_dft_c2r(planb_fftw,ksd_fftw,rsd_fftw)
#else
      call fftw_execute_dft_c2r(planb_fftw,ksd_fftw,rsd_fftw)
#endif

! test output
!      write(fileo,'(a3,i3.3)') 'rsd',myrank
!      open(60,file=fileo,status='unknown')
! #ifdef FFTW3_MPI
!       do k=1,lnz_fftw
! #else
!       do k=1,nz_m
! #endif
!         do j=1,ny_m
!           do i=1,nx_m
! #ifdef FFTW3_MPI
!             write(60,'(3I6,2(1pE13.5))') i,j,k &
!              ,rsd_fftw(i,j,k)/dble(nt_m) &
!              ,dsqrt(dble(i**2+j**2+k**2))
! #else
!             write(60,'(3I6,3(1pE13.5))') i,j,k &
!              ,rsd_fftw(i,j,k)/dble(nt_m) &
!              ,dsqrt(dble(i**2+j**2+k**2))
! #endif
!           enddo
!         enddo
!       enddo
!       close(60)

! set real space potential data 
#ifdef FFTW3_MPI

      do k=0,lnz_m+5
        do j=0,ny_m-1
          do i=0,nx_m-1
            pot_m(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      dnt_m=dble(nx_m)*dble(ny_m)*dble(nz_m)
      do kfftw=1,lnz_fftw
        do jfftw=1,ny_fftw
          do ifftw=1,nx_fftw
! normalised (has to be done)
! first 3 and last 3 arrays will be imported from the other cores -1+3=+2
            pot_m(ifftw-1,jfftw-1,kfftw+2)= &
#ifdef ISOPM
              rsd_fftw(ifftw,jfftw,kfftw)/(dnt_m)
#else
              rsd_fftw(ifftw,jfftw,kfftw)/dnt_m
#endif

#endif
          enddo
        enddo
      enddo

#ifdef FFTW3_MPI
      if(nprocs.gt.1) then
! send and receive potential of first 3 marginal grid 
! sending proc
        srank=myrank+1
        if(srank.ge.nprocs) then
          srank=srank-nprocs
        endif
! receiving proc
        rrank=myrank-1
        if(rrank.lt.0) then
         rrank=rrank+nprocs
        endif
! offset for sending pot_m lzoff(ips)-3-lzoff(myrank)+3
        koff=lzoff_mp(srank)-3
        if(koff.lt.0) then
          koff=koff+nz_m
        endif
        koff=koff-lzoff_m+3

        if(koff.lt.3.or.koff.ge.lnz_m+3) then
          write(6,*) ' Error in mesh_fftf(): setting first 3 marginal grid'
          write(6,*) ' myrank,srank,rrank,koff=',myrank,srank,rrank
          write(6,*) ' lzoff(srank,myrank)=',lzoff_mp(srank),lzoff_m
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif

        allocate(tdvs(0:nx_m*ny_m*3))
        allocate(tdvr(0:nx_m*ny_m*3))

        nsend=0
        do k=koff,koff+2
          do j=0,ny_m-1
            do i=0,nx_m-1
              tdvs(nsend)=pot_m(i,j,k)
              tdvr(nsend)=0.0d0
              nsend=nsend+1
            enddo
          enddo
        enddo
        if(nsend.ne.nx_m*ny_m*3) then
          write(6,*) 'Error in mesh_fftf():'
          write(6,*) ' setting sending data for first 3 merginal grid'
          write(6,*) ' nsend=',nsend,' should be',nx_m*ny_m*3
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif

        call MPI_ISEND(tdvs,nsend,MPI_DOUBLE_PRECISION,srank,1 &
         ,MPI_COMM_WORLD,ireqs,ierr)
        call MPI_IRECV(tdvr,nsend,MPI_DOUBLE_PRECISION,rrank,1 &
         ,MPI_COMM_WORLD,ireqr,ierr)
        call MPI_WAIT(ireqs,istatus,ierr)
        call MPI_WAIT(ireqr,istatus,ierr)

        nsend=0
        do k=0,2
          do j=0,ny_m-1
            do i=0,nx_m-1
              pot_m(i,j,k)=tdvr(nsend)
              nsend=nsend+1
            enddo
          enddo
        enddo

! send and receive potential of last 3 marginal grid 
! sending proc
        srank=myrank-1
        if(srank.lt.0) then
          srank=srank+nprocs
        endif
! receiving proc
        rrank=myrank+1
        if(rrank.ge.nprocs) then
          rrank=rrank-nprocs
        endif
! offset for sending pot_m lzoff(ips)+lnz(ips)-1+3-lzoff(myrank)+3-2
        koff=lzoff_mp(srank)+lnz_mp(srank)+2
        if(koff.ge.nz_m) then
          koff=koff-nz_m
        endif
        koff=koff-lzoff_m+3-2

        if(koff.lt.3.or.koff.ge.lnz_m+3) then
          write(6,*) ' Error in mesh_fftf(): setting last 3 marginal grid'
          write(6,*) ' myrank,srank,rrank,koff=',myrank,srank,rrank,koff
          write(6,*) ' lzoff(srank,myrank)=',lzoff_mp(srank),lzoff_m
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif

        nsend=0
        do k=koff,koff+2
          do j=0,ny_m-1
            do i=0,nx_m-1
              tdvs(nsend)=pot_m(i,j,k)
              tdvr(nsend)=0.0d0
              nsend=nsend+1
            enddo
          enddo
        enddo
        if(nsend.ne.nx_m*ny_m*3) then
          write(6,*) 'Error in mesh_fftf():'
          write(6,*) ' setting sending data for last 3 merginal grid'
          write(6,*) ' nsend=',nsend,' should be',nx_m*ny_m*3
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif

        call MPI_ISEND(tdvs,nsend,MPI_DOUBLE_PRECISION,srank,2 &
         ,MPI_COMM_WORLD,ireqs,ierr)
        call MPI_IRECV(tdvr,nsend,MPI_DOUBLE_PRECISION,rrank,2 &
         ,MPI_COMM_WORLD,ireqr,ierr)
        call MPI_WAIT(ireqs,istatus,ierr)
        call MPI_WAIT(ireqr,istatus,ierr)

        nsend=0
        do k=lnz_m+3,lnz_m+5
          do j=0,ny_m-1
            do i=0,nx_m-1
              pot_m(i,j,k)=tdvr(nsend)
              nsend=nsend+1
            enddo
          enddo
        enddo

        deallocate(tdvs)
        deallocate(tdvr)

     else

       do j=0,ny_m-1
         do i=0,nx_m-1
           pot_m(i,j,0)=pot_m(i,j,lnz_m-3+3)
           pot_m(i,j,1)=pot_m(i,j,lnz_m-2+3)
           pot_m(i,j,2)=pot_m(i,j,lnz_m-1+3)
           pot_m(i,j,lnz_m+3)=pot_m(i,j,3)
           pot_m(i,j,lnz_m+4)=pot_m(i,j,4)
           pot_m(i,j,lnz_m+5)=pot_m(i,j,5)
         enddo
       enddo

     endif

#endif

! set force
#ifdef FFTW3_MPI
      do k=0,lnz_m+1
! pot_m k+2 corresponds to k 
        kp=k+2
        kp2=k+4
        kp1=k+3
        km2=k
        km1=k+1
        if((kp2.lt.0.or.kp2.ge.lnz_m+6) &
         .or.(km2.lt.0.or.km2.ge.lnz_m+6)) then
          write(6,*) ' Error in mesh_fftf(): kp for force calculation'
          write(6,*) '  myrank,k,kp2,km2=',myrank,k,kp2,km2
        endif
#else
      do k=0,nz_m-1
        kp=k
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
#endif
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
             *((2.0d0/3.0d0)*(pot_m(ip1,j,kp)-pot_m(im1,j,kp)) &
             -((1.0d0/12.0d0)*(pot_m(ip2,j,kp)-pot_m(im2,j,kp))))
            fy_m(i,j,k)=(1.0d0/dy_m) &
             *((2.0d0/3.0d0)*(pot_m(i,jp1,kp)-pot_m(i,jm1,kp)) &
             -((1.0d0/12.0d0)*(pot_m(i,jp2,kp)-pot_m(i,jm2,kp))))
            fz_m(i,j,k)=(1.0d0/dz_m) &
             *((2.0d0/3.0d0)*(pot_m(i,j,kp1)-pot_m(i,j,km1)) &
             -((1.0d0/12.0d0)*(pot_m(i,j,kp2)-pot_m(i,j,km2))))
     
          enddo
        enddo
      enddo

!     write(fileo,'(a5,i3.3)') 'fmesh',myrank
!     open(60,file=fileo,status='unknown')

!#ifdef FFTW3_MPI
!      do k=0,lnz_m-1+2
!#else
!      do k=0,nz_m-1
!#endif
!        do j=0,ny_m-1
!          do i=0,nx_m-1
!#ifdef FFTW3_MPI
!            kp1=k+lzoff_m-1
!            if(kp1.lt.0) then
!              kp1=kp1+nz_m
!            endif
!            if(kp1.ge.nz_m) then
!              kp1=kp1-nz_m
!            endif
!#else
!            kp1=k
!#endif
!            write(60,'(3I10,3(1pE13.5))') i,j,kp1 &
!             ,fx_m(i,j,k),fy_m(i,j,k),fz_m(i,j,k)
!          enddo
!        enddo
!      enddo
!      close(60)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      stop


end subroutine
#endif





