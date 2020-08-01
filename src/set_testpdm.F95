#include "gcdp.def"
! ***************************************
!      set_testpdm.F95 for GCD+ 
!  9 September 2016    produced by D.KAWATA
! *************************************** 

#if defined(TESTP) || defined(TESTP_HYBRID)
subroutine set_testpdm(ndm,ndmt,nadm,nadmt)
      use gcdp_const
      use gcdp_system
#ifdef DM
      use gcdp_dm
#endif
     
      implicit none
      include 'mpif.h'

      integer,intent(inout) :: ndm,ndmt,nadm,nadmt
      integer i,ip,nc,pn
      integer nval,nprp
      integer ipsta,ipend
! work
      integer jsta,jend,ierr,rrank,i1,i2,i3,i4
      integer,allocatable :: jstap(:),jendp(:)
      integer,allocatable :: istatus(:)
      integer,allocatable :: ireqs(:),ireqr(:)
      double precision,allocatable :: tdvr(:)

! allocate array for MPI
      allocate(istatus(MPI_STATUS_SIZE))
      allocate(ireqs(0:nprocs-1))
      allocate(ireqr(0:nprocs-1))
      allocate(jstap(0:nprocs-1))
      allocate(jendp(0:nprocs-1))

! reading test particles and add them to the end of dm particle data
      if(myrank.eq.0) then
! *** read DM data only with rank 0 
        write(6,*) ' reading testp.dat'
        open(51,file='ini/testp.dat',status='old',form='unformatted')
        read(51) SI_nttpdm
        do ip=0,nprocs-1
          call para_range(0,SI_nttpdm-1,nprocs,ip,jstap(ip),jendp(ip))
        enddo
      endif
      call MPI_BCAST(SI_nttpdm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call para_range(0,SI_nttpdm-1,nprocs,myrank,jsta,jend)
! number of particles in each proc
      SI_ntpdm=jend-jsta+1

      nc=ndm
      nval=6
! *** Read and  Set Position ***
      rrank=0
      if(myrank.eq.0) then
        ipsta=0
        ipend=nprocs-1
      else
        ipsta=myrank
        ipend=myrank
      endif
      do ip=ipsta,ipend
        if(myrank.eq.0) then
! number of baryon particle read for rank ip
          nprp=jendp(ip)-jstap(ip)+1

          allocate(tdvr(0:nprp*nval-1))

        else

          allocate(tdvr(0:SI_ntpdm*nval-1))
  
        endif

        if(myrank.eq.0) then
! read the data only by rank 0
          do i=0,nprp-1
            read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
              ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5)
          enddo
          if(ip.ne.myrank) then
            call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
             ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        else
          rrank=0
          call MPI_IRECV(tdvr,SI_ntpdm*nval,MPI_DOUBLE_PRECISION &
           ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
          call MPI_WAIT(ireqr(ip),istatus,ierr)
        endif
        if(myrank.eq.ip) then
          do i=0,SI_ntpdm-1
            x_dm(nc)=tdvr(i*nval)
            y_dm(nc)=tdvr(i*nval+1)
            z_dm(nc)=tdvr(i*nval+2)
            vnx_dm(nc)=tdvr(i*nval+3)
            vny_dm(nc)=tdvr(i*nval+4)
            vnz_dm(nc)=tdvr(i*nval+5)
! set m and rho same as id=0
            m_dm(nc)=m_dm(0)
            rho_dm(nc)=rho_dm(0)
            id_dm(nc)=ndmt+jsta+i
            nc=nc+1
          enddo
        endif

        deallocate(tdvr)
   
      enddo
! change ndm
      ndm=nc
      if(myrank.eq.0) then
        close(51)
      endif

! set initla value of test particle info
      do i=0,ndm-1
        rperi_dm(i)=INF
        rapo_dm(i)=0.0d0
        zmax_dm(i)=0.0d0
        xc_dm(i)=x_dm(i)
        yc_dm(i)=y_dm(i)
        zc_dm(i)=z_dm(i)
#ifdef SIM1D
        h_dm(i)=(1.0d0+ERRH)*ETAHDM*((m_dm(i)/rho_dm(i)))
#elif defined(SIM2D)
        h_dm(i)=(1.0d0+ERRH)*ETAHDM*((m_dm(i)/rho_dm(i))**0.5d0)
#else
        h_dm(i)=(1.0d0+ERRH)*ETAHDM*((m_dm(i)/rho_dm(i))**THIRD)
#endif
        omgh_dm(i)=1.0d0
        zetah_dm(i)=0.0d0
        list_adm(i)=i
        flagt_dm(i)=0
        dtr_dm(i)=INF
        dt_dm(i)=0.0d0
        dtmnb_dm(i)=INF
      enddo

      deallocate(istatus)
      deallocate(ireqs)
      deallocate(ireqr)
      deallocate(jstap)
      deallocate(jendp)

! make new list_adm, test particle first and then the other dm
      nc=0
      do i=ndm-SI_ntpdm,ndm-1
        list_adm(nc)=i
        nc=nc+1
      enddo
      do i=0,ndm-SI_ntpdm-1
        list_adm(nc)=i
        nc=nc+1
     enddo
#ifdef TESTP_HYBRID
! ndm is set to be nc (see above)
     nadm=ndm
#else     
     nadm=SI_ntpdm
#endif
! update the number of DM particles
     ndmt=ndmt+SI_nttpdm
#ifdef TESTP_HYBRID     
     nadmt=ndmt
#else
     nadmt=SI_nttpdm
#endif
     SI_ndm1t=ndmt
     SI_ndm1=ndm
     if(myrank.eq.0) then
        write(6,*) ' new ndm and ndmt=',ndm,ndmt
        write(6,*) ' Ntottal test particle and in myrank 0 =',SI_nttpdm,SI_ntpdm
     endif

! test output
!      open(60,file='fixdm.asc',status='unknown')
!      do i=nadm,ndm-1
!        pn=list_adm(i)
!        write(60,'(3(1pE13.5),I10)') x_dm(pn),y_dm(pn),z_dm(pn),id_dm(pn)
!      enddo
!      close(60)
!      open(60,file='testp.asc',status='unknown')
!      do i=0,nadm-1
!        pn=list_adm(i)
!        write(60,'(3(1pE13.5),I10)') x_dm(pn),y_dm(pn),z_dm(pn),id_dm(pn)
!      enddo
!      close(60)

end subroutine
#endif
