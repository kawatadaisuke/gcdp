#include "gcdp.def"
!*********************************************
!   setyield.F for GCD+ ver f03.0
!   03 Feb. 2013    written by D.Kawata  
! *********************************************

#if defined(SF_EFD) || defined(SF_ZFD)
subroutine  setyield()
      use gcdp_const
      use gcdp_system
      use gcdp_yields

      implicit none
      include 'mpif.h'

      integer i,j,k
      integer nzg
      integer nval,nc,ierr
! *** metallicity and age ***
      double precision llz,luz,dlz,llt,lut,dlt
      double precision x,mu,ml,mlsnii
      double precision,allocatable :: tdvr(:)
      character chead*22,cread*200
      character fileo*60

      if(myrank.eq.0) then
        open(50,file='./ini/tzy.cmp.dat',status='old')
        read(50,'(4(1pE13.5),I10)') x,mu,ml,mlsnii,SI_nsp
        read(50,'(a22)') chead
        if(myrank.eq.0) then
          write(6,*) 'IMF: x,mu,ml=',x,mu,ml
          write(6,*) ' SNe II mass limit, N mass group=',mlsnii,SI_nsp
          write(6,*) chead
        endif
      endif
      call MPI_BCAST(SI_nsp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(myrank.eq.0) then
        do j=0,NYTZ+1
          i=0
          t_ytb(0) =0.0d0
          mej_ytb(i,j)=0.0d0
          nsn_ytb(i,j)=0.0d0
          mzHe_ytb(i,j)=0.0d0
          mzC_ytb(i,j)=0.0d0
          mzN_ytb(i,j)=0.0d0
          mzO_ytb(i,j)=0.0d0
          mzNe_ytb(i,j)=0.0d0
          mzMg_ytb(i,j)=0.0d0
          mzSi_ytb(i,j)=0.0d0
          mzFe_ytb(i,j)=0.0d0
          mzZ_ytb(i,j)=0.0d0
          nsp_ytb(j)=0
! *** Unit Gyr ***
          do i=1,NYTT+1
            read(50,'(a200)') cread
            if(cread(1:1).eq.'#') then
              read(cread(2:14),'(1pE13.5)') tspi_ytb(nsp_ytb(j),j)
              nsp_ytb(j)=nsp_ytb(j)+1
              read(50,'(a200)') cread
            endif
            read(cread(1:195),150) t_ytb(i),z_ytb(j),mzHe_ytb(i,j) &
             ,mzC_ytb(i,j),mzN_ytb(i,j),mzO_ytb(i,j),mzNe_ytb(i,j) &
             ,mzMg_ytb(i,j),mzSi_ytb(i,j),mzFe_ytb(i,j),mzZ_ytb(i,j) &
             ,mej_ytb(i,j),nsn_ytb(i,j),nsw_ytb(i,j)
  150       format(14(1pE13.5))
! *** set initial number of SW stars ***
            if(i.eq.1) then
              nsw_ytb(0,j)=nsw_ytb(i,j)
            endif
          enddo
! *** set final mass group time ***
          tspi_ytb(nsp_ytb(j),j)=t_ytb(NYTT+1)
          nsp_ytb(j)=nsp_ytb(j)+1
        enddo
! *** read mgr group SNII end 
        read(50,'(a9,I10)') chead,SI_snii
        write(6,*) ' mass group with SNeII<',SI_snii
        close(50)
      endif

      call MPI_BCAST(SI_snii,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! sending the data to the other ranks
      nval=12

      allocate(tdvr(0:nval*(NYTZ+2)*(NYTT+2)-1))

      nc=0      
      if(myrank.eq.0) then
        do j=0,NYTZ+1
          do i=0,NYTT+1
            tdvr(nc)=mzHe_ytb(i,j)
            tdvr(nc+1)=mzC_ytb(i,j)
            tdvr(nc+2)=mzN_ytb(i,j)
            tdvr(nc+3)=mzO_ytb(i,j)
            tdvr(nc+4)=mzNe_ytb(i,j)
            tdvr(nc+5)=mzMg_ytb(i,j)
            tdvr(nc+6)=mzSi_ytb(i,j)
            tdvr(nc+7)=mzFe_ytb(i,j)
            tdvr(nc+8)=mzZ_ytb(i,j)
            tdvr(nc+9)=mej_ytb(i,j)
            tdvr(nc+10)=nsn_ytb(i,j)
            tdvr(nc+11)=nsw_ytb(i,j)
            nc=nc+nval
          enddo
        enddo
!        write(6,*) ' setyield nc check=',nc,nval*(NYTZ+2)*(NYTT+2)
      endif
      nc=nval*(NYTZ+2)*(NYTT+2)
      call MPI_BCAST(tdvr,nc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        nc=0
        do j=0,NYTZ+1
          do i=0,NYTT+1
            mzHe_ytb(i,j)=tdvr(nc)
            mzC_ytb(i,j)=tdvr(nc+1)
            mzN_ytb(i,j)=tdvr(nc+2)
            mzO_ytb(i,j)=tdvr(nc+3)
            mzNe_ytb(i,j)=tdvr(nc+4)
            mzMg_ytb(i,j)=tdvr(nc+5)
            mzSi_ytb(i,j)=tdvr(nc+6)
            mzFe_ytb(i,j)=tdvr(nc+7)
            mzZ_ytb(i,j)=tdvr(nc+8)
            mej_ytb(i,j)=tdvr(nc+9)
            nsn_ytb(i,j)=tdvr(nc+10)
            nsw_ytb(i,j)=tdvr(nc+11)
            nc=nc+nval
          enddo
        enddo
      endif

      deallocate(tdvr)

! set mass group number and time
      call MPI_BCAST(nsp_ytb,NYTZ+2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      do j=0,NYTZ+1

        allocate(tdvr(0:nsp_ytb(j)-1))

        if(myrank.eq.0) then
          do i=0,nsp_ytb(j)-1
            tdvr(i)=tspi_ytb(i,j)
          enddo
        endif
        call MPI_BCAST(tdvr,nsp_ytb(j),MPI_DOUBLE_PRECISION &
         ,0,MPI_COMM_WORLD,ierr)
        if(myrank.ne.0) then
          do i=0,nsp_ytb(j)-1
            tspi_ytb(i,j)=tdvr(i)
          enddo
        endif             

        deallocate(tdvr)

      enddo     

! *** reset Metallicity ***
      llz = -6.0d0
      luz = 0.0d0
      dlz = (luz-llz)/dble(NYTZ)
      z_ytb(0)=0.0d0
      do  i=1,NYTZ+1
        z_ytb(i)=10.0d0**(llz+dlz*dble(i-1))
      enddo
! *** reset Age ***
      llt = -3.0d0
      lut = 1.2d0
      dlt = (lut-llt)/dble(NYTT)

      do  i=1,NYTT+1
        t_ytb(i)=10.0d0**(llt+dlt*dble(i-1))
      enddo

! for test
!     write(fileo,'(a3,i3.3)') 'tzy',myrank
!     open(60,file=fileo,status='unknown')
!     do  j=0,NYTZ+1
!       write(60,'(a1,I8)') '#',nsp_ytb(j)
!       do i=0,nsp_ytb(j)-1
!         write(60,'(a1,1pE13.5)') '#',tspi_ytb(i,j)
!       enddo
!       do i=0,NYTT+1
!         write(60,150) t_ytb(i),z_ytb(j),mzHe_ytb(i,j) &
!         ,mzC_ytb(i,j),mzN_ytb(i,j),mzO_ytb(i,j),mzNe_ytb(i,j) &
!          ,mzMg_ytb(i,j),mzSi_ytb(i,j),mzFe_ytb(i,j),mzZ_ytb(i,j) &
!          ,mej_ytb(i,j),nsn_ytb(i,j),nsw_ytb(i,j)
!       enddo
!     enddo
!     close(60)

end subroutine
#endif

