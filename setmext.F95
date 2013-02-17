#include "gcdp.def"
! *********************************************
!  setmext.F95 for GCD+  ver. f03.0
!  04  Feb. 2013   written by D. Kawata
! *********************************************

#ifdef MEXT
subroutine  setmext()
      use gcdp_const
      use mext
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer i,nval,ierr
      character filen*60
! for work
      double precision,allocatable :: tdvr(:)

      if(myrank.eq.0) then
        open(50,file='./ini/mext.dat',status='old',form='unformatted')
        read(50) SI_nmext
        read(50) SI_dlr,SI_lri,SI_lro
        do i=0,SI_nmext-1
          read(50) rmext(i),mextr(i)
        enddo
        close(50)
! *** sending numbers to the other procs ***
      endif
      call MPI_BCAST(SI_nmext,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! *** send SIs ***  

      allocate(tdvr(0:2))

      if(myrank.eq.0) then
        tdvr(0)=SI_dlr
        tdvr(1)=SI_lri
        tdvr(2)=SI_lro
      endif
      nval=3
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        SI_dlr=tdvr(0)
        SI_lri=tdvr(1)
        SI_lro=tdvr(2)
      endif

      deallocate(tdvr)

! *** send rmext ***
!      write(6,*) ' myrank,nmext=',myrank,SI_nmext,SI_dlr
      nval=SI_nmext
      call MPI_BCAST(rmext,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! *** send mextr ***
      call MPI_BCAST(mextr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! for test
!      write(filen,'(a4,i3.3,a4)') 'mext',myrank,'.dat'
!      open(60,file=filen,status='unknown')
!      do i=0,SI_nmext-1
!        write(60,*) rmext(i),mextr(i)
!      enddo
!      close(60)

end subroutine
#endif

