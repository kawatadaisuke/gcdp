#include "gcdp.def"
! *****************************************************
!    del_dmpart.F95 for gcd+ 
!  28 June 2016  written by D.KAWATA
! ***************************************************** 

#ifdef DM
#ifdef DEL_DMPART
! add or delete particles
! delete only 
subroutine del_dmpart(ndm,ndmt)
      use gcdp_const
      use gcdp_dm
      use gcdp_system
 
      implicit none
      include 'mpif.h'

      integer,intent(inout) :: ndm,ndmt
      integer i,j,ip
      integer ndc,ndmdel,ndm0
      integer pn,pnt,pntsta,pndc,pnswp
      integer ndmtdctg,dndmtdc,ndmdelt
! for work
      integer ierr
      integer,allocatable :: ndmdcproc(:),tivr(:)
      integer,allocatable :: pntdc(:)
      integer,allocatable :: listdc(:),listdmdel(:)
      integer,allocatable :: pnold(:),flagdel(:)
      double precision,allocatable :: tdvr(:)
      real ran1
      external ran1

      allocate(ndmdcproc(0:nprocs-1))

! original number of DM particle
      ndm0=ndm
! target number at the next step
      ndmtdctg=int((SI_mtdmdc0-SI_deldmdt*(TM_tot+TM_dt))/SI_mdmp0)
      dndmtdc=SI_ndmtdc-ndmtdctg
! use the value at myrank=0
      call MPI_BCAST(dndmtdc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(dndmtdc.lt.0) then
        if(myrank.eq.0) then
          write(6,*) 'Error in del_dmpart: ndelp is negative=',dndmtdc
          write(6,*) 'ndm target, current=',ndmtdctg,SI_ndmtdc
        endif
        stop
      endif
      if(myrank.eq.0) then
        write(6,*) 'N delete particles=',dndmtdc
      endif

      if(dndmtdc.gt.0) then
!  number of particles for the component deleting the particle in each proc
        do ip=0,nprocs-1
          ndmdcproc(ip)=0
        enddo

        allocate(listdc(0:ndm-1))

        ip=myrank
        ndc=0
        do i=0,ndm-1
          if(id_dm(i).ge.SI_dmdelci(0).and.id_dm(i).le.SI_dmdelci(1)) then
            ndmdcproc(ip)=ndmdcproc(ip)+1
! store the list of del component
            listdc(ndc)=i
            ndc=ndc+1
          endif
        enddo

        allocate(tivr(0:nprocs-1))
     
        call MPI_ALLREDUCE(ndmdcproc,tivr,nprocs,MPI_INTEGER &
          ,MPI_SUM,MPI_COMM_WORLD,ierr)
        do ip=0,nprocs-1
          ndmdcproc(ip)=tivr(ip)
        enddo

        deallocate(tivr)

! identify the particle deleting

        allocate(pntdc(0:dndmtdc-1))

        if(myrank.eq.0) then
! do random search only in myrank=0
          do i=0,dndmtdc-1
! poition within the component
 90         pnt=int(ran1(idum)*real(SI_ndmtdc))
            do j=0,i-1
! check for hitting the same particle
              if(pnt.eq.pntdc(j)) then
                goto 90
               endif
            enddo
            pntdc(i)=pnt
          enddo
        endif
        call MPI_BCAST(pntdc,dndmtdc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! set starting pnt
        pntsta=0
        do ip=0,myrank-1
          pntsta=pntsta+ndmdcproc(ip)
        enddo
! identifying deleting particles

        allocate(listdmdel(0:dndmtdc-1))

        ndmdel=0
        do i=0,dndmtdc-1  
          pnt=pntdc(i)
          if(pnt.ge.pntsta.and.pnt.lt.pntsta+ndc) then
            pndc=listdc(pnt-pntsta)
            listdmdel(ndmdel)=pndc
            ndmdel=ndmdel+1
          endif
        enddo

        deallocate(listdc)
        deallocate(pntdc)

! check total deleted particles
        ndmdelt=0
        call MPI_ALLREDUCE(ndmdel,ndmdelt,1,MPI_INTEGER,MPI_SUM &
          ,MPI_COMM_WORLD,ierr)
        if(ndmdelt.ne.dndmtdc) then
          if(myrank.eq.0) then
            write(6,*) 'Error in del_dmpart: N del error.'
            write(6,*) ' ndmdelt, dndmtdc=',ndmdelt,dndmtdc
          endif
          stop
        endif

        if(ndmdel.gt.0) then
! new list of particles from old list pn
          ndm=ndm0-ndmdel
     
          allocate(pnold(0:ndm-1))
          allocate(flagdel(0:ndm0-1))

! flag for deleted particles
          do i=0,ndm0-1
            flagdel(i)=0
          enddo        
          do i=0,ndmdel-1
            pndc=listdmdel(i)
            flagdel(pndc)=1
          enddo
          ndm=0
          do i=0,ndm0-1
            if(flagdel(i).eq.0) then
! survived particles
              pnold(ndm)=i
              ndm=ndm+1
            endif
          enddo
          if(ndm.ne.ndm0-ndmdel) then
            write(6,*) 'Error in del_dmpart: ndm inconsistent after delete'
            write(6,*) 'myrank,ndm,ndm0,ndmdel=',myrank,ndm,ndm0,ndmdel 
            stop
          endif

          deallocate(flagdel)

! reallocate the DM data using pnold
! integer
  
          allocate(tivr(0:ndm-1))

! flagt
          do i=0,ndm-1
            tivr(i)=flagt_dm(pnold(i))
          enddo
          do i=0,ndm-1
            flagt_dm(i)=tivr(i)
          enddo
! id
          do i=0,ndm-1
            tivr(i)=id_dm(pnold(i))
          enddo
          do i=0,ndm-1
            id_dm(i)=tivr(i)
          enddo

          deallocate(tivr)
          allocate(tdvr(0:ndm-1))

! double
! xc
          do i=0,ndm-1
            tdvr(i)=xc_dm(pnold(i))
          enddo
          do i=0,ndm-1
            xc_dm(i)=tdvr(i)
          enddo
! yc
          do i=0,ndm-1
            tdvr(i)=yc_dm(pnold(i))
          enddo
          do i=0,ndm-1
            yc_dm(i)=tdvr(i)
          enddo
! zc
          do i=0,ndm-1
            tdvr(i)=zc_dm(pnold(i))
          enddo
          do i=0,ndm-1
            zc_dm(i)=tdvr(i)
          enddo
! vnx
          do i=0,ndm-1
            tdvr(i)=vnx_dm(pnold(i))
          enddo
          do i=0,ndm-1
            vnx_dm(i)=tdvr(i)
          enddo
! vny
          do i=0,ndm-1
            tdvr(i)=vny_dm(pnold(i))
          enddo
          do i=0,ndm-1
            vny_dm(i)=tdvr(i)
          enddo
! vnz
          do i=0,ndm-1
            tdvr(i)=vnz_dm(pnold(i))
          enddo
          do i=0,ndm-1
            vnz_dm(i)=tdvr(i)
          enddo
! m
          do i=0,ndm-1
            tdvr(i)=m_dm(pnold(i))
          enddo
          do i=0,ndm-1
            m_dm(i)=tdvr(i)
          enddo
! h
          do i=0,ndm-1
            tdvr(i)=h_dm(pnold(i))
          enddo
          do i=0,ndm-1
            h_dm(i)=tdvr(i)
          enddo
! dvx
          do i=0,ndm-1
            tdvr(i)=dvx_dm(pnold(i))
          enddo
          do i=0,ndm-1
            dvx_dm(i)=tdvr(i)
          enddo
! dvy
          do i=0,ndm-1
            tdvr(i)=dvy_dm(pnold(i))
          enddo
          do i=0,ndm-1
            dvy_dm(i)=tdvr(i)
          enddo
! dvz
          do i=0,ndm-1
            tdvr(i)=dvz_dm(pnold(i))
          enddo
          do i=0,ndm-1
            dvz_dm(i)=tdvr(i)
          enddo
! dt
          do i=0,ndm-1
            tdvr(i)=dt_dm(pnold(i))
          enddo
          do i=0,ndm-1
            dt_dm(i)=tdvr(i)
          enddo
! vdt
          do i=0,ndm-1
            tdvr(i)=vdt_dm(pnold(i))
          enddo
          do i=0,ndm-1
            vdt_dm(i)=tdvr(i)
          enddo
! lt
          do i=0,ndm-1
            tdvr(i)=lt_dm(pnold(i))
          enddo
          do i=0,ndm-1
            lt_dm(i)=tdvr(i)
          enddo
! rho
          do i=0,ndm-1
            tdvr(i)=rho_dm(pnold(i))
          enddo
          do i=0,ndm-1
            rho_dm(i)=tdvr(i)
          enddo
! zetah
          do i=0,ndm-1
            tdvr(i)=zetah_dm(pnold(i))
          enddo
          do i=0,ndm-1
            zetah_dm(i)=tdvr(i)
          enddo
! omgh
          do i=0,ndm-1
            tdvr(i)=omgh_dm(pnold(i))
          enddo
          do i=0,ndm-1
            omgh_dm(i)=tdvr(i)
          enddo
! div_v_dm
          do i=0,ndm-1
            tdvr(i)=div_v_dm(pnold(i))
          enddo
          do i=0,ndm-1
            div_v_dm(i)=tdvr(i)
          enddo
! dtr
          do i=0,ndm-1
            tdvr(i)=dtr_dm(pnold(i))
          enddo
          do i=0,ndm-1
            dtr_dm(i)=tdvr(i)
          enddo
#if defined(COSM) && defined(BOXSIM)
! asc
          do i=0,ndm-1
            tdvr(i)=asc_dm(pnold(i))
          enddo
          do i=0,ndm-1
            asc_dm(i)=tdvr(i)
          enddo
! tc
          do i=0,ndm-1
            tdvr(i)=tc_dm(pnold(i))
          enddo
          do i=0,ndm-1
            tc_dm(i)=tdvr(i)
          enddo
#endif

          deallocate(tdvr)
          deallocate(pnold)

! list_adm will be reset at ddecdm
          do i=0,ndm-1
            list_adm(i)=i
          enddo   

          deallocate(listdmdel)

        endif

        SI_ndm1=ndm
        SI_ndmtdc=SI_ndmtdc-dndmtdc
        ndmt=ndmt-dndmtdc
        SI_ndm1t=ndmt

      endif

      deallocate(ndmdcproc)

end subroutine

#endif
#endif
