#include "gcdp.def"
! *****************************************************
!    add_dmpart.F95 for gcd+ 
!  3 Oct. 2017  written by D.KAWATA
! ***************************************************** 

#ifdef DM
#ifdef ADD_DMPART
! add particles
subroutine add_dmpart(ndm,ndmt)
      use gcdp_const
      use gcdp_dm
      use gcdp_system
      use gcdp_expdisk
 
      implicit none
      include 'mpif.h'

      integer,intent(inout) :: ndm,ndmt
      integer i,j,ip
      integer ndmadd,ndm0
      integer dndmtadd,dnaddp
      double precision mtotaddtn
      double precision,allocatable :: xadd0(:),yadd0(:),zadd0(:) &
        ,xaddp(:),yaddp(:),zaddp(:),radd0(:)
      double precision pfr,rp,ph
      double precision rxyp,vrotp,vradp
! for work
      integer nval,rrank,nc
      integer jsta,jend,ipsta,ipend,nprp,ierr
      integer,allocatable :: jstap(:),jendp(:)
      integer,allocatable :: istatus(:)
      integer,allocatable :: ireqs(:),ireqr(:)
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:)

      real ran1,gasdev
      external ran1,gasdev

! original number of DM particle
      ndm0=ndm
! total mass added by the next step
      mtotaddtn=SI_sfr0dmadd*(SI_tsfrdmadd/exp(SI_tinidmadd/SI_tsfrdmadd) &
       -(SI_tsfrdmadd/exp((TM_tot+TM_dt)/SI_tsfrdmadd)))
! target number at the next step
      dndmtadd=int(mtotaddtn/SI_mdmp0)-SI_ndmtadded
!      if(myrank.eq.0) then
!        write(6,*) ' dn,SI_dmdtadd,T,SI_mdmp0=',dndmtadd,SI_dmdtadd &
!         ,TM_tot+TM_dt,SI_mdmp0,SI_ndmtadded
!      endif
! use the value at myrank=0
      call MPI_BCAST(dndmtadd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(dndmtadd.lt.0) then
        if(myrank.eq.0) then
          write(6,*) 'Error in add_dmpart: naddp is negative=',dndmtadd
          write(6,*) 'ndm target, already added=' &
            ,int(SI_dmdtadd*(TM_tot+TM_dt)/SI_mdmp0),SI_ndmtadded
        endif
        stop
      endif
      if(myrank.eq.0) then
        write(6,*) 'N particles added =',dndmtadd,' SFR=',SI_sfr0dmadd &
         *dexp(-(TM_tot+TM_dt)/SI_tsfrdmadd)
      endif

      if(dndmtadd.gt.0) then

! allocate array for MPI
        allocate(istatus(MPI_STATUS_SIZE))
        allocate(ireqs(0:nprocs))
        allocate(ireqr(0:nprocs))
        allocate(jstap(0:nprocs-1))
        allocate(jendp(0:nprocs-1))

! calculate the position of added particles only in myrank=-
        if(myrank.eq.0) then

          allocate(xadd0(0:dndmtadd-1))
          allocate(yadd0(0:dndmtadd-1))
          allocate(zadd0(0:dndmtadd-1))
          allocate(radd0(0:dndmtadd-1))

#ifdef ADD_EXPDISK
          do i=0,dndmtadd-1
! calculate the radius
   69       pfr=ran1(idum)
            do j=1,NEXPDTAB-1
              if(mx_expd(j).gt.pfr) then
                goto 71
              endif
            enddo
            write(6,*) ' Warining in add_dmpart: Radd not found. pfr=',pfr
            goto 69
 71         if(j.eq.0) then
              j=1
            endif
            radd0(i)=10.0d0**(lx_expd(j-1)+(lx_expd(j)-lx_expd(j-1)) &
              *(dlog10(pfr)-lmx_expd(j-1))/(lmx_expd(j)-lmx_expd(j-1)))
          enddo
          do i=0,dndmtadd-1
! calculate x,y,z position
            ph=2.0d0*M_PI*dble(ran1(idum))
            xadd0(i)=radd0(i)*dcos(ph)
            yadd0(i)=radd0(i)*dsin(ph)
            zadd0(i)=0.0d0
          enddo          

          deallocate(radd0)

#endif
! number of new particles allocated in each proc
          do ip=0,nprocs-1
            call para_range(0,dndmtadd-1,nprocs,ip,jstap(ip),jendp(ip))
          enddo 

        endif
! distribute new particles position to the other proc
        call para_range(0,dndmtadd-1,nprocs,myrank,jsta,jend)
! number of particles added
        dnaddp=jend-jsta+1

! allocate xadd
        if(dnaddp.gt.0) then
          allocate(xaddp(0:dnaddp-1))
          allocate(yaddp(0:dnaddp-1))
          allocate(zaddp(0:dnaddp-1))
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

              allocate(tdvr(0:nprp*nval-1))

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
      
            allocate(tdvr(0:dnaddp*nval-1))

            call MPI_IRECV(tdvr,dnaddp*nval,MPI_DOUBLE_PRECISION &
               ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(dnaddp.gt.0.and.myrank.eq.ip) then
            nc=0
            do i=0,dnaddp-1
              xaddp(nc)=tdvr(i)
              yaddp(nc)=tdvr(i+dnaddp)
              zaddp(nc)=tdvr(i+dnaddp*2)
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
          ndm=ndm0+dnaddp

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
            list_adm(i)=i
          enddo            
! new valuse
          nc=0
          do i=ndm0,ndm-1
! active particle (similar to RANAZRLIM in update.F95
            flagt_dm(i)=1
            id_dm(i)=ndmt+jsta+nc
            list_adm(i)=i
            nc=nc+1
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
            vrotp=220.0/VUKMS+gasdev(idum)*SI_sigph0add
            vradp=gasdev(idum)*SI_sigr0add
            rxyp=dsqrt(xc_dm(i)**2+yc_dm(i)**2)
! positive = clock-wise 
            vnx_dm(i)=vrotp*yc_dm(i)/rxyp+vradp*xc_dm(i)/rxyp
            vny_dm(i)=-vrotp*xc_dm(i)/rxyp+vradp*yc_dm(i)/rxyp
            vnz_dm(i)=dble(gasdev(idum))*SI_sigz0add
            nc=nc+1
          enddo

          deallocate(tdvr)

! m, h, dvx, dt, lt
          nval=7

          allocate(tdvr(0:ndm0*nval-1))

          do i=0,ndm0-1
            tdvr(i)=m_dm(i)
            tdvr(i+ndm0)=h_dm(i)
            tdvr(i+ndm0*2)=dvx_dm(i)
            tdvr(i+ndm0*3)=dvy_dm(i)
            tdvr(i+ndm0*4)=dvz_dm(i)
            tdvr(i+ndm0*5)=dt_dm(i)
            tdvr(i+ndm0*6)=lt_dm(i)
          enddo

! reallocate DM variables 2
#ifdef TESTP
          call reallocate_dm_d2(ndm+SI_ntpdm)
#else
          call reallocate_dm_d2(ndm)
#endif

          do i=0,ndm0-1
            m_dm(i)=tdvr(i)
            h_dm(i)=tdvr(i+ndm0)
            dvx_dm(i)=tdvr(i+ndm0*2)
            dvy_dm(i)=tdvr(i+ndm0*3)
            dvz_dm(i)=tdvr(i+ndm0*4)
            dt_dm(i)=tdvr(i+ndm0*5)
            lt_dm(i)=tdvr(i+ndm0*6)
          enddo
! new particles
          do i=ndm0,ndm-1
! set the same mass with the other prticles
            m_dm(i)=SI_mdmp0
            h_dm(i)=SI_hdmp0
            dvx_dm(i)=0.0d0
            dvy_dm(i)=0.0d0
            dvz_dm(i)=0.0d0
! this routine is called before updating TM_lt
            dt_dm(i)=TM_dt
            lt_dm(i)=TM_lt
          enddo

          deallocate(tdvr)

! rho,zeta,omgh,div_v, dtr
! no cosmo or box sim
          nval=6

          allocate(tdvr(0:ndm0*nval-1))

          do i=0,ndm0-1
            tdvr(i)=rho_dm(i)
            tdvr(i+ndm0)=zetah_dm(i)
            tdvr(i+ndm0*2)=omgh_dm(i)
            tdvr(i+ndm0*3)=div_v_dm(i)
            tdvr(i+ndm0*4)=dtr_dm(i)
            tdvr(i+ndm0*5)=tadd_dm(i)
          enddo

! reallocate DM variables 3
#ifdef TESTPDMF
      call reallocate_dm_d3(ndm+SI_ntpdm)
#else
      call reallocate_dm_d3(ndm)
#endif

          do i=0,ndm0-1
            rho_dm(i)=tdvr(i)
            zetah_dm(i)=tdvr(i+ndm0)
            omgh_dm(i)=tdvr(i+ndm0*2)
            div_v_dm(i)=tdvr(i+ndm0*3)
            dtr_dm(i)=tdvr(i+ndm0*4)
            tadd_dm(i)=tdvr(i+ndm0*5)
          enddo
! new particles
          if(ndm0.gt.0) then
            do i=ndm0,ndm-1
              rho_dm(i)=rho_dm(0)
            enddo
          else
            do i=ndm0,ndm-1
! set an arbitrary number 
              rho_dm(i)=0.1d0
            enddo
          endif
          do i=ndm0,ndm-1
            zetah_dm(i)=0.0d0
            omgh_dm(i)=1.0d0
            div_v_dm(i)=0.0d0
            dtr_dm(i)=TM_dt
            tadd_dm(i)=TM_tot
          enddo

          deallocate(tdvr)

#ifdef TESTP
          nval=3

          allocate(tdvr(0:ndm0*nval-1))

          do i=0,ndm0-1
            tdvr(i)=rperi_dm(i)
            tdvr(i+ndm0)=rapo_dm(i)
            tdvr(i+ndm0*2)=zmax_dm(i)
          enddo

! reallocate DM variables
#ifdef TESTPDMF
          call reallocate_dm_testp(ndm+SI_ntpdm)
#else
          call reallocate_dm_testp(ndm)
#endif

          do i=0,ndm0-1
            rperi_dm(i)=tdvr(i)
            rapo_dm(i)=tdvr(i+ndm0)
            zmax_dm(i)=tdvr(i+ndm0*2)
          enddo
! for new particles
          do i=ndm0,ndm-1
            rxyp=dsqrt(xc_dm(i)**2+y_dm(i)**2)
            rperi_dm(i)=rxyp
            rapo_dm(i)=rxyp
            zmax_dm(i)=dabs(zc_dm(i))
          enddo
   
          deallocate(tdvr)

#endif
        endif

        SI_ndmtadded=SI_ndmtadded+dndmtadd
        SI_ndm1=ndm
        ndmt=ndmt+dndmtadd
        SI_ndm1t=ndmt

!        open(60,file='addp.dat',status='unknown')
!        do i=ndm0,ndm-1
!          write(60,'(6(1pE13.5),I10)') xc_dm(i),yc_dm(i),zc_dm(i) &
!            ,vnx_dm(i),vny_dm(i),vnz_dm(i),id_dm(i)
!        enddo
!        close(60)

! allocate array for MPI
        deallocate(istatus)
        deallocate(ireqs)
        deallocate(ireqr)
        deallocate(jstap)
        deallocate(jendp)

      endif

end subroutine

#endif
#endif
