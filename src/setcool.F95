#include "gcdp.def"
! *********************************************
!  setcool.F95 for GCD+ ver. f03.0
!  03 Feb. 2013   written by D. Kawata
! *********************************************

#if defined(COOL) || defined(BBAT)
subroutine  setcool()
      use gcdp_const
      use gcdp_cool
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer i,j,k,ierr
      integer iz,inh,imet,it,nval,nc
      integer lmet10,dmet10,lnh10,dnh10
      integer z100,met10,nh10
      character czred*5,cmet*5,cnh*5,filen*60
      double precision crad,hrad
! for work
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:)

      allocate(tivr(0:3))

      if(myrank.eq.0) then
        open(50,file='./cool12feb/setcool.dat',status='old')
        read(50,'(I6)') nz_crtb
        read(50,150) nnh_crtb,lnh10,dnh10
        read(50,150) nmet_crtb,lmet10,dmet10
        read(50,'(I6)') nt_crtb
        read(50,'(F5.2)') SI_zeor
  150   format(3I6)
        close(50)

        open(51,file='./cool12feb/setcoolzred.dat',status='old')
        do iz=0,nz_crtb-1
          read(51,*) z100
          z_crtb(iz)=dble(z100)/100.0d0
          write(czred,'(a1,i4.4)') 'z',z100
          do inh=0,nnh_crtb-1
            nh10=lnh10+inh*dnh10
            nh_crtb(inh)=dble(nh10)/10.0d0
            if(nh10.ge.0) then
              write(cnh,'(a2,i3.3)') 'hp',nh10
            else
              write(cnh,'(a2,i3.3)') 'hm',-nh10
            endif
            do imet=0,nmet_crtb-1
              met10=lmet10+imet*dmet10
              met_crtb(imet)=dble(met10)/10.0d0
              if(met10.ge.0) then
                write(cmet,'(a2,i3.3)') 'Zp',met10
              else
                write(cmet,'(a2,i3.3)') 'Zm',-met10
              endif
! *** set filename ***
              write(filen,'(a18,a5,a5,a5,a4)') &
!               12345678901234567890
               './cool12feb/chrate',czred,cmet,cnh,'.dat'
              open(50,file=filen,status='old')
              do i=1,4
                read(50,*)
              enddo
              do i=0,nt_crtb-1
                read(50,151) t_crtb(i),cr_crtb(i,imet,inh,iz) &
                 ,hr_crtb(i,imet,inh,iz),myu_crtb(i,imet,inh,iz) &
                 ,ne_crtb(i,imet,inh,iz)
  151           format(f5.3,f8.3,f8.3,f6.3,f8.3) 
              enddo
              close(50)
            enddo
          enddo
        enddo
        close(51)
! *** sending numbers to the other procs ***
        tivr(0)=nz_crtb
        tivr(1)=nnh_crtb
        tivr(2)=nmet_crtb
        tivr(3)=nt_crtb
        write(6,*) ' cooling table highest z=',z_crtb(nz_crtb-1)
        write(6,*) ' EoR z=',SI_zeor
      endif
      call MPI_BCAST(tivr,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nz_crtb=tivr(0)
        nnh_crtb=tivr(1)
        nmet_crtb=tivr(2)
        nt_crtb=tivr(3)
      endif

      deallocate(tivr)


! *** send zeor ***
      call MPI_BCAST(SI_zeor,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! *** send zred grid ***
      nval=nz_crtb

      allocate(tdvr(0:nval))

      if(myrank.eq.0) then
        do iz=0,nz_crtb-1
          tdvr(iz)=z_crtb(iz)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do iz=0,nz_crtb-1
          z_crtb(iz)=tdvr(iz)
        enddo
      endif
! *** send nh grid ***
      nval=nnh_crtb

      deallocate(tdvr)
      allocate(tdvr(0:nval))

      if(myrank.eq.0) then
        do inh=0,nnh_crtb-1
          tdvr(inh)=nh_crtb(inh)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do inh=0,nnh_crtb-1
          nh_crtb(inh)=tdvr(inh)
        enddo
      endif
! *** send metal grid ***
      nval=nmet_crtb

      deallocate(tdvr)
      allocate(tdvr(0:nval))

      if(myrank.eq.0) then
        do imet=0,nmet_crtb-1
          tdvr(imet)=met_crtb(imet)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do imet=0,nmet_crtb-1
          met_crtb(imet)=tdvr(imet)
        enddo
      endif
! *** send T grid ***
      nval=nt_crtb

      deallocate(tdvr)
      allocate(tdvr(0:nval))

      if(myrank.eq.0) then
        do i=0,nt_crtb-1
          tdvr(i)=t_crtb(i)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        do i=0,nt_crtb-1
          t_crtb(i)=tdvr(i)
        enddo
      endif
! *** send cooling rate ***
      nval=nt_crtb*nmet_crtb*nnh_crtb*nz_crtb

      deallocate(tdvr)
      allocate(tdvr(0:nval))

      if(myrank.eq.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                tdvr(nc)=cr_crtb(i,imet,inh,iz)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                cr_crtb(i,imet,inh,iz)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
! *** send heating rate ***
      if(myrank.eq.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                tdvr(nc)=hr_crtb(i,imet,inh,iz)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                hr_crtb(i,imet,inh,iz)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
! *** send myu ***
      if(myrank.eq.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                tdvr(nc)=myu_crtb(i,imet,inh,iz)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                myu_crtb(i,imet,inh,iz)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif

! *** send ne ***
      if(myrank.eq.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                tdvr(nc)=ne_crtb(i,imet,inh,iz)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do iz=0,nz_crtb-1
          do inh=0,nnh_crtb-1
            do imet=0,nmet_crtb-1
              do i=0,nt_crtb-1
                ne_crtb(i,imet,inh,iz)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          enddo
        enddo
      endif

      deallocate(tdvr)

end subroutine
#endif

