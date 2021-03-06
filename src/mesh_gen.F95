#include "gcdp.def"
! *********************************************
!  mesh_gen.F95 for GCD+
!  13  Aug. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
subroutine  mesh_gen(flag)
      use gcdp_const
      use gcdp_system
      use gcdp_pm
      use gcdp_kernel
#ifdef FFTW3_MPI
      use fftw3
#endif

      implicit none
      include 'mpif.h'

      integer,intent(in) :: flag

      integer i,j,k
      double precision yp,zp
      double precision rrsh
!      character filen*60
#ifdef FFTW3_MPI
      integer ierr
      integer,allocatable :: tivs(:),tivr(:)
#endif
      double precision erfcc
      external erfcc

! allocate memory space
      if(flag.eq.0) then

#ifdef ISOPM
! set SI_lbox for isolated PM 
        SI_lbox(0)=2.0d0*ISOPMLMARGIN*MAXRLIM*dble(nx_m)/(dble(nx_m-5))
        SI_lbox(1)=SI_lbox(0)
        SI_lbox(2)=SI_lbox(0)
        if(myrank.eq.0) then
          write(6,*) ' PM box size (including inactive space)=' &
           ,SI_lbox(0),SI_lbox(1),SI_lbox(2)
        endif
#endif

#ifdef FFTW3_MPI
        lnz_m=int(lnz_fftw)
        lzoff_m=int(lkoff_fftw)

        if(myrank.eq.0.and.lzoff_m.ne.0) then
          write(6,*) ' Error in mesh_gen(): lzoff(0) should be zero'
          write(6,*) ' lzoff(0)=',lzoff_m
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
        if(myrank.eq.nprocs-1.and.lzoff_m+lnz_m-1.ne.nz_m-1) then
          write(6,*) ' Error in mesh_gen(): lzoff(nprocs-1) should be nz-1'
          write(6,*) ' lzoff,lnz(nprocs-1)=',lzoff_m,lnz_m
!          call MPI_ABORT(MPI_COMM_WORLD,ierr)
!          stop
        endif

!        write(6,*) ' myrank,lnz,lzoff=',myrank,lnz_m,lzoff_m

        allocate(lnz_mp(0:nprocs))
        allocate(lzoff_mp(0:nprocs))

        allocate(tivs(0:nprocs*2))
        allocate(tivr(0:nprocs*2))

        do i=0,2*nprocs-1
          tivs(i)=0
          tivr(i)=0
        enddo
        tivs(myrank)=lnz_m
        tivs(nprocs+myrank)=lzoff_m
        call MPI_ALLREDUCE(tivs,tivr,2*nprocs,MPI_INTEGER &
          ,MPI_SUM,MPI_COMM_WORLD,ierr)
        do i=0,nprocs-1
          lnz_mp(i)=tivr(i)
          lzoff_mp(i)=tivr(i+nprocs)
        enddo

! check
        do i=0,nprocs-1
!          if(myrank.eq.0) then
!            write(6,*) ' rank,lnz,lzoff=',i,lnz_mp(i),lzoff_mp(i)
!          endif
          if(i.ne.0.and.lzoff_mp(i-1).gt.lzoff_mp(i)-3) then
            if(myrank.eq.0) then
            write(6,*) &
             ' Error in mesh_gen(): need at least 3 grids in neighbour proc.'
            write(6,*) ' myrank,ip,lzoff(i-1,i),lnz=' &
             ,myrank,i,lzoff_mp(i-1),lzoff_mp(i),lnz_m
            endif
!            call MPI_ABORT(MPI_COMM_WORLD,ierr)
!            stop
          endif
          if(i.ne.nprocs-1.and.lzoff_mp(i)+lnz_mp(i)+2 &
            .gt.lzoff_mp(i+1)+lnz_mp(i+1)-1) then
            if(myrank.eq.0) then
            write(6,*) &
             ' Error in mesh_gen(): need at least 3 grids in neighbour proc.'
            write(6,*) ' myrank,ip,lzoff(i,i+1),lnz=' &
              ,myrank,i,lzoff_mp(i),lzoff_mp(i+1),lnz_m
            endif
!            call MPI_ABORT(MPI_COMM_WORLD,ierr)
!            stop
          endif
! proc i and procs-1 are not checked, but should be ok.
        enddo

        deallocate(tivs)
        deallocate(tivr)

        call allocate_pm(nx_m,ny_m,nz_m,lnz_m)
#else
        call allocate_pm(nx_m,ny_m,nz_m)
#endif

! half nx,ny,nz
        nx2_m=nx_m/2
        ny2_m=ny_m/2
        nz2_m=nz_m/2
! total number of mesh 
        nt_m=nx_m*ny_m*nz_m

      endif

! calculate grid space
      dx_m=SI_lbox(0)/dble(nx_m)
      dy_m=SI_lbox(1)/dble(ny_m)
      dz_m=SI_lbox(2)/dble(nz_m)
! calculate grid space in k space, 2 pi/L seems better than 1/L
      dkx_m=2.0d0*M_PI/(SI_lbox(0))
      dky_m=2.0d0*M_PI/(SI_lbox(1))
      dkz_m=2.0d0*M_PI/(SI_lbox(2))

! first grid points
#ifdef ISOPM
! set x=0 to -SI_lbox/4+0.54dx, for making SI_lbox/2-SI_lbox inactive.
      x0_m=-0.25d0*SI_lbox(0)+0.5d0*dx_m
      y0_m=-0.25d0*SI_lbox(1)+0.5d0*dy_m
      z0_m=-0.25d0*SI_lbox(2)+0.5d0*dz_m
#else
      x0_m=-0.5d0*SI_lbox(0)+0.5d0*dx_m
      y0_m=-0.5d0*SI_lbox(1)+0.5d0*dy_m
      z0_m=-0.5d0*SI_lbox(2)+0.5d0*dz_m
#endif

      if(flag.eq.0) then
        if(myrank.eq.0) then
          write(6,*) ' Mesh info'
          write(6,*) '   nx_m,ny_m,nz_m,lnz_m=',nx_m,ny_m,nz_m,lnz_m
          write(6,*) '   dx,dy,dz=',dx_m,dy_m,dz_m
          write(6,*) '   x0,y0,z0=',x0_m,y0_m,z0_m
        endif
      endif

! total box size
      vt_m=SI_lbox(0)*SI_lbox(1)*SI_lbox(2)
! mesh size
      dv_m=dx_m*dy_m*dz_m
! set TreePM scale length
! rs
      SI_rstpm=FRSTPM*(dv_m**THIRD)
! rs^2
      SI_rsp2tpm=(SI_rstpm**2)
! radius for short range force
      SI_rshtpm=FRSHTPM*SI_rstpm

      if(flag.eq.0) then
        if(myrank.eq.0) then
          write(6,*) ' rs for TreePM=',SI_rstpm
          write(6,*) '  rcut=',SI_rshtpm
        endif
      endif

! generate short force truncation factor table
      dnftpmtab=dble(NFTPMTAB)
!      if(myrank.eq.0) then
!        open(60,file='fshfac.dat',status='unknown')
!      endif
      do i=0,NFTPMTAB
        rstpm_tb(i)=dble(i)/dnftpmtab
        rrsh=(dble(i)/dnftpmtab)*SI_rshtpm
        fshtpm_tb(i)=erfcc(rrsh/(2.0d0*SI_rstpm)) &
           +(rrsh/(SI_rstpm*dsqrt(M_PI))) &
           *dexp(-(rrsh**2)/(4.0d0*(SI_rstpm**2)))
!        if(myrank.eq.0.and.i.ne.0) then
!          write(60,'(3(1pE13.5))') rrsh/SI_rshtpm,fshtpm_tb(i) &
!           ,(fshtpm_tb(i)-fshtpm_tb(i-1))/fshtpm_tb(i)
!        endif
      enddo
!      if(myrank.eq.0) then
!        close(60)  
!      endif

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

end subroutine
#endif

