#include "gcdp.def"
! *****************************************************
!    add_dmpart.F95 for gcd+ 
!  26 June 2017  written by D.KAWATA
! ***************************************************** 

#ifdef DM
#ifdef ADD_DMPART
! add particles
subroutine add_dmpart(ndm,ndmt)
      use gcdp_const
      use gcdp_dm
      use gcdp_system
 
      implicit none
      include 'mpif.h'

      integer,intent(inout) :: ndm,ndmt
      integer i,j,ip
      integer ndmadd,ndm0
      integer dndmtadd,dnaddp
      double precision,allocatable :: xadd0(:),yadd0(:),zadd0(:) &
        ,xaddp(:),yaddp(:),zaddp(:)
      double precision pfr,rp,ph
      double precision raddp,vsigrp,vsigphp,vsigzp,vrotp,vradp,rp
! for work
      integer nval,rrank,nc
      integer jsta,jend,ipsta,ipend,nprp
      integer,allocatable :: jstap(:),jendp(:)
      integer,allocatable :: istatus(:)
      integer,allocatable :: ireqs(:),ireqr(:)
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:)

      real ran1,gasdev
      external ran1,gasdev

! original number of DM particle
      ndm0=ndm
! target number at the next step
      dndmtadd=int(SI_dmdtadd*(TM_tot+TM_dt))/SI_mdmp0)-SI_ndmtadded
! use the value at myrank=0
      call MPI_BCAST(dndmtadd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(dndmtadd.lt.0) then
        if(myrank.eq.0) then
          write(6,*) 'Error in add_dmpart: naddp is negative=',dndmtadd
          write(6,*) 'ndm target, already added=' &
            ,int(SI_dmdtadd*(TM_tot+TM_dt))/SI_mdmp0),SI_ndmtadded
        endif
        stop
      endif
      if(myrank.eq.0) then
        write(6,*) 'N particles added =',dndmtadd
      endif

      if(dndmtadd.gt.0) then

! allocate array for MPI
        allocate(istatus(MPI_STATUS_SIZE))
        allocate(ireqs(0:nprocs))
        allocate(ireqr(0:nprocs))
        allocate(jstap(0:nprocs-1))
        allocate(jendp(0:nprocs-1))


!  number of particles for the component adding the particle in each proc
        do ip=0,nprocs-1
          ndmdcproc(ip)=0
        enddo

! calculate the position of added particles only in myrank=-
        if(myrank.eq.0) then

          allocate(xadd0(0:dndmtadd-1))
          allocate(yadd0(0:dndmtadd-1))
          allocate(zadd0(0:dndmtadd-1))

#ifdef ADD_EXPDISK
          do i=0,dndmtadd-1
! calculate the position of 
   69       pfr=ran1(idum)
            do j=1,NEXPDTAB-1
              if(mx_expd(j).gt.pfr) then
                goto 71
              endif
            enddo
            goto 69
 71         if(j.eq.0) then
              j=1
            endif
            rp=10.0d0**(lx_expd(j-1)+(lx_expd(j)-lx_expd(j-1)) &
              /(dlog10(mx_expd(j))-dlog10(mx_expd(j-1))))
            ph=2.0d0*M_PI*dble(ran1(idum))
            xadd0(i)=rp*dcos(ph)
            yadd0(i)=rp*dsin(ph)
            zadd0(i)=0.0d0
          enddo          
#endif

! number of new particles allocated in each proc
          do ip=0,nprocs-1
            call para_range(0,dndmtadd-1,nprocs,ip,jstap(ip),jendp(ip))
          enddo 

        endif
! distribute new particles position to the other proc
        call para_range(0,dndmtadd-1,nprocs,ip,jsta,jend)
! number of particles added
        dnaddp=jend-jsta+1

! allocate xadd
        if(dnaddp.gt.0) the
          allocate(xadd(0:dnaddp-1))
          allocate(yadd(0:dnaddp-1))
          allocate(zadd(0:dnaddp-1))
        endif

! distribute new particle position from myrank 0 to the other proc
        nval=3
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
! number of particles sending from myrank to ip
            nprp=jendp(ip)-jstap(ip)+1

            if(nprp.gt.0) then

              allocate(tdvr(0:nprp*nval-1)

              do i=0,nprp-1
                tdvr(i)=xadd0(jstap(ip)+i)
                tdvr(i+nprp)=yadd0(jstap(ip)+i)
                tdvr(i+nprp*2)=zadd0(jstap(ip)+i)
              enddo
              if(ip.ne.myrank) then
                call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
                 ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            endif
          else if(dnaddp.gt.0) then
      
            allocate(tdvr(0:dnaddp*nval-1)    

            call MPI_IRECV(tdvr,ndaddp*nval,MPI_DOUBLE_PRECISION &
               ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(dnaddp.gt.0.and.myrank.eq.ip) then
            nc=0
            do i=0,dnaddp-1
              xaddp(nc)=tdvr(i)
              yaddp(nc)=tdvr(i+nval)
              zaddp(nc)=tdvr(i+nval*2)
              nc=nc+1
            enddo
          endif

          if((myrank.eq.0.and.nprp.gt.0).or.dnaddp.gt.0) then

            deallocate(tdvr)            

          endif
        enddo

        if(myrank.eq.0) then
    
          deallocate(xadd0)
          deallocate(yadd0)
          deallocate(zadd0)

        endif

        if(dnaddp.gt.0) then
! new list of particles from old list pn
          ndm=ndm0+ndmdel

! add new particle. 
! integer 
          nval=2

          allocate(tivr(0:ndm0*nval-1))         

! save the original ones
          do i=0,ndm0-1
            tivr(i)=flagt_dm(i)
            tivr(i+ndm0)=id_dm(i)
          enddo
  
! reallocate DM variables
#ifdef TESTPDMF
          call reallocate_dm_int(ndm+SI_ntpdm)
#else
          call reallocate_dm_int(ndm)
#endif

          do i=0,ndm0-1
            flagt_dm(i)=tivr(i)
            id_dm(i)=tivr(i+ndm0)
          enddo            
! new valuse
          do i=ndm0,ndm-1
! active particle (similar to RANAZRLIM in update.F95
            flagt_dm(i)=1
            id_dm(i)=ndmt+jsta+i
          enddo

          deallocate(tivr)
          
! xc?_dm, vn?_dm
          nval=6

          allocate(tdvr(0:ndm0*nval-1))

          do i=0,ndm0-1
            tdvr(i)=xc_dm(i)
            tdvr(i+ndm0)=yc_dm(i)
            tdvr(i+ndm0*2)=zc_dm(i)
            tdvr(i+ndm0*3)=vnx_dm(i)
            tdvr(i+ndm0*4)=vny_dm(i)
            tdvr(i+ndm0*5)=vnz_dm(i)
          enddo

! reallocate DM variables
#ifdef TESTPDMF
          call reallocate_dm_d1(ndm+SI_ntpdm)
#else
          call reallocate_dm_d1(ndm)
#endif

           do i=0,ndm0-1
            xc_dm(i)=tdvr(i)
            yc_dm(i)=tdvr(i+ndm0)
            zc_dm(i)=tdvr(i+ndm0*2)
            vnx_dm(i)=tdvr(i+ndm0*3)
            vny_dm(i)=tdvr(i+ndm0*4)
            vnz_dm(i)=tdvr(i+ndm0*5)
           enddo
! new values
           nc=0
           do i=ndm0,ndm-1
             xc_dm(i)=xaddp(nc)
             yc_dm(i)=yaddp(nc)
             zc_dm(i)=zaddp(nc)
             raddp=dsqrt(xc_dm(i)**2+y_dm(i)**2+z_dm(i))
             vsigrp20.0d0/VUKMS
             vsigphp=20.0d0/VUKMS
             vrotp=220.0/VUKMS+gasdev(idum)*vsigphp
             vradp=gasdev(idum)*vsigrp
             vsigzp=20.0d0/VUKMS
! positive = clock-wise 
             vnx_dm(i)=vrotp*yc_dm(i)/raddp+vradp*xc_dm(i)/raddp
             vny_dm(i)=-vrotp*xc_dm(i)/raddp+vradp*yc_dm(i)/raddp
             vnz_dm(i)=dble(gasdev)*vsigzp
           enddo

! put the reallocate the DM data using pnold
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

! allocate array for MPI
        deallocate(istatus)
        deallocate(ireqs)
        deallocate(ireqr)
        deallocate(jstap)
        deallocate(jendp)


      endif

      deallocate(ndmdcproc)

end subroutine

#endif
#endif
