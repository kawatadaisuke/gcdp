#include "gcdp.def"
! ***************************************
!      set_ini.F95 for GCD+ ver. f03.0
!  26 Jun., 2014    produced by D.KAWATA
! *************************************** 

subroutine set_ini(npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
      ,step,totstep,lstep,flagcont)
      use gcdp_const
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#ifdef METAL
      use gcdp_metal
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
      use gcdp_yields
#endif
#endif
#ifdef DM
      use gcdp_dm
#endif
#ifdef MHD
      use gcdp_mhd
#endif
      use gcdp_system
#ifdef TREEPM
      use gcdp_pm
#endif
     
      implicit none
      include 'mpif.h'

      integer,intent(out) :: npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
       ,step,totstep,lstep,flagcont
      integer MNL
      parameter (MNL=3)
      character fin*20,tcha*1,fileo*60,ctmp*60,mname*20
      integer i,j,il,ip,nc,pn,itmp
      integer ipsta,ipend
      integer flagrf
      integer nif,idf,ndmf,ngf,nsf,ngt,nst
      integer flagfdm,flagfg,flagfs
      integer nval,nival,nprp
      double precision gamf,gkmm
#ifdef COSM
! *** nd0, nd1 for number of particles for level 0 and 1 ***
      integer nd(0:MNL),ndl,ndt(0:MNL),nplev
#endif
! * om0 (km s^-1 kpc^-1) *      
      double precision Ini_OMGg,dxi(0:MNL)
! * (10^12 M@), (100 kpc),(1.0e4 K),(4.71e8 yr) * /     
      integer flagtime
      double precision Ini_T0,Ini_ze,Ini_zi
! *** to read recfast.out **
      double precision zrec(2),xrec(4,2),xe(4),hnf,henf
      double precision Crec(6)
      integer Ini_nd
! *** physical value ***
      double precision rhob,rhob0,hubi,iai,rhoba
! * Mass of SPH, DM, total *      
      double precision m,mdm
      double precision kmmp
      double precision hmin,hmint,etahmax
! *** Initial parameter ***
! *** Ini_hyd: 0: N-body only, 1: include gas in high reso region ***
! *** Ini_mul: 0: only level0, 1: level 0 and 1 ***
      integer Ini_hyd,Ini_mul
! work
      integer jsta,jend,ierr,rrank
      integer,allocatable :: jstap(:),jendp(:)
      integer,allocatable :: istatus(:)
      integer,allocatable :: ireqs(:),ireqr(:)
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:)
#if defined(SF_EFD) || defined(SF_ZFD)
! work
      double precision,allocatable :: mfd(:)
! * for feedback  common *
! * Total *
      double precision,allocatable :: tmej_sn(:),tnsn_sn(:),tmzHe_sn(:) &
       ,tmzC_sn(:),tmzN_sn(:),tmzO_sn(:),tmzNe_sn(:),tmzMg_sn(:) &
       ,tmzSi_sn(:),tmzFe_sn(:),tmzZ_sn(:),nsw_sn(:)
#endif

! allocate array for MPI
      allocate(istatus(MPI_STATUS_SIZE))
      allocate(ireqs(0:nprocs))
      allocate(ireqr(0:nprocs))
      allocate(jstap(0:nprocs-1))
      allocate(jendp(0:nprocs-1))

#ifdef COSM
      if(myrank.eq.0) then
        write(6,*) ' H0,HUB0 = ',H0_1,HUB0
      endif
#endif

#ifdef COOL
      SI_flagrad=0
#else
      SI_flagrad=1
#endif

      if(myrank.eq.0) then
!*****   Open Initial Data File ./ini/input.dat
        open(50,file='ini/input.dat',status='old')      
! *** Read Initial Parameter ***
! * Number of Particle etc. *  
        read(50,*) Ini_mul
! * Interval of Output Data, Distribution *
        read(50,*) Ini_nd
        read(50,*) SI_nof
! * Initial Values unit is our Unit *
#ifdef COSM
        read(50,*) Ini_OMGg
        read(50,*) Ini_ze
#else
        read(50,*) Ini_zi
        read(50,*) flagtime
        read(50,*) TM_tot,TM_t0
        read(50,*) flagfg,flagfdm,flagfs
#endif
#ifdef TREEPM
        read(50,*) nx_m,ny_m,nz_m
        if((mod(nx_m,2).ne.0)) then
          write(6,*) ' Error in set_ini(): nx_m,ny_m,nz_m=',nx_m,ny_m,nz_m
          write(6,*) ' nx_m should be able to be devided by 2.'
          call MPI_FINALIZE(ierr)
          stop
        endif
#endif
#ifndef MULTI_LEV
        if(Ini_mul.gt.0) then
          write(6,*) ' Error: MULTI_LEV is off in gcdp.def, but Ini_mul=' &
           ,Ini_mul
          call MPI_FINALIZE(ierr)
          stop
        endif
#endif
        close(50)
      endif

! send the input parameters
! integer 
#ifdef COSM
      nval=3
#else
      nval=7
#endif
#ifdef TREEPM
      nval=nval+3
#endif

      allocate(tivr(0:nval-1))
 
      if(myrank.eq.0) then
        tivr(0)=Ini_mul
        tivr(1)=Ini_nd
        tivr(2)=SI_nof
#ifdef COSM
        i=3
#else
        tivr(3)=flagtime
        tivr(4)=flagfg
        tivr(5)=flagfdm
        tivr(6)=flagfs
        i=7
#endif
#ifdef TREEPM
        tivr(i)=nx_m
        tivr(i+1)=ny_m
        tivr(i+2)=nz_m
#endif  
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Ini_mul=tivr(0)
      Ini_nd=tivr(1)
      SI_nof=tivr(2)
#ifdef COSM
      i=3
#else
      flagtime=tivr(3)
      flagfg=tivr(4)
      flagfdm=tivr(5)
      flagfs=tivr(6)
      i=7
#endif
#ifdef TREEPM
      nx_m=tivr(i)
      ny_m=tivr(i+1)
      nz_m=tivr(i+2)
#endif  

      deallocate(tivr)

! double
#ifdef COSM
      nval=2
#else
      nval=3
#endif

      allocate(tdvr(0:nval-1))

      if(myrank.eq.0) then
#ifdef COSM
        tdvr(0)=Ini_OMGg
        tdvr(1)=Ini_ze
#else
        tdvr(0)=Ini_zi
        tdvr(1)=TM_tot
        tdvr(2)=TM_t0
#endif
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#ifdef COSM
      Ini_OMGg=tdvr(0)
      Ini_ze=tdvr(1)
#else
      Ini_zi=tdvr(0)
      TM_tot=tdvr(1)
      TM_t0=tdvr(2)
#endif
      deallocate(tdvr)

! *** initialization for number of particles ***   
      npt=0
      ndmt=0
      nst=0
      ngt=0
      np=0
      ng=0
      ndm=0
      ns=0
#ifndef COSM
      Ini_ze=Ini_zi
#endif

#ifdef COSM
! *****   cosmological simulation data input *****
#ifdef STAR
      SI_tnsf=0
#endif
! *** cosmological simulation ***
      SI_omgb=Ini_OMGg
! *** Check Input Parameter ***
      if(Ini_OMGg.lt.0.0d0.or.Ini_OMGg.gt.1.0d0) then
        if(myrank.eq.0) then
          write(6,*) ' Error in set_ini() : OMGg = ',Ini_OMGg
        endif
        stop
      endif
! *** reading initial condision for cosmological simulation
! *** Open Initial Data File for level 1 ***
      ndmt = 0
      ndm = 0
#ifdef DM
! get ndm for each rank
      do il=Ini_mul,0,-1
        if(myrank.eq.0) then
! read all the data by rank 0
          write(fin,'(a9,i1,a5)') 'ini/level',il,'d.dat'
          open(51,file=fin,status='old',form='unformatted')
! *** read the data ***
          read(51) nif,idf,ndl,nd(il),flagrf
          write(6,*) il,'level ndm=',ndl
          close(51)
        endif
        call MPI_BCAST(ndl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! count the number DM for each core for each level
        call para_range(0,ndl-1,nprocs,myrank,jsta,jend)
        ndm=ndm+jend-jsta+1
      enddo 

! allocate dm variables
      call allocate_dm_all(ndm)      

      ndm=0
      do il=Ini_mul,0,-1
        nval=5
      
        allocate(tivr(0:nval-1))

        if(myrank.eq.0) then
          write(fin,'(a9,i1,a5)') 'ini/level',il,'d.dat'
          open(51,file=fin,status='old',form='unformatted')
          read(51) nif,idf,ndl,nd(il),flagrf

! store the range of the particle for each rank
          do ip=0,nprocs-1
            call para_range(0,ndl-1,nprocs,ip,jstap(ip),jendp(ip))
          enddo

          tivr(0)=nif
          tivr(1)=idf
          tivr(2)=ndl
          tivr(3)=nd(il)
          tivr(4)=flagrf
        endif
        call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        nif=tivr(0)
        idf=tivr(1)
        ndl=tivr(2)
        nd(il)=tivr(3)
        flagrf=tivr(4)

        deallocate(tivr)

! count the number DM for each core for each level
        call para_range(0,ndl-1,nprocs,myrank,jsta,jend)
! number of particle in each level
        nplev=jend-jsta+1

! for double data
        nval=5

        allocate(tdvr(0:nval-1))

        if(myrank.eq.0) then
! *** flagrf: 0: read only position and velocity, otherwise read also mass
! *** dx (Mpc, Comoving), zstart,omegam,lambda ***
          read(51) dxi(il),SI_a,SI_omg0,SI_lam0,SI_h0
          tdvr(0)=dxi(il)
          tdvr(1)=SI_a
          tdvr(2)=SI_omg0
          tdvr(3)=SI_lam0
          tdvr(4)=SI_h0
        endif
        call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        dxi(il)=tdvr(0)
        SI_a=tdvr(1)
        SI_omg0=tdvr(2)
        SI_lam0=tdvr(3)
        SI_h0=tdvr(4)
        dxi(il) = dxi(il)*10.0d0

        deallocate(tdvr)

        ndt(il)=0
        nd(il)=nplev
        if(flagrf.eq.0) then
          nval=6
        else
          nval=7
        endif

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
! number DM particle read for rank ip
            nprp=jendp(ip)-jstap(ip)+1

            allocate(tdvr(0:nprp*nval-1))

          else

            allocate(tdvr(0:nplev*nval-1))
 
          endif

          if(myrank.eq.0) then
            if(flagrf.eq.0) then
              do i=0,nprp-1
                read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                 ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5)
              enddo
            else
              do i=0,nprp-1
                read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                 ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5),tdvr(i*nval+6)
              enddo
            endif
            if(ip.ne.myrank) then
              call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
               ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif
          else
            call MPI_IRECV(tdvr,nplev*nval,MPI_DOUBLE_PRECISION &
             ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(myrank.eq.ip) then
            do i=0,nplev-1
              x_dm(ndm)=tdvr(i*nval)
              y_dm(ndm)=tdvr(i*nval+1)
              z_dm(ndm)=tdvr(i*nval+2)     
              vx_dm(ndm)=tdvr(i*nval+3)
              vy_dm(ndm)=tdvr(i*nval+4)
              vz_dm(ndm)=tdvr(i*nval+5)     
              if(flagrf.ne.0) then
                m_dm(ndm)=tdvr(i*nval+6)
              endif
              id_dm(ndm)=ndmt+jsta+i
              ndm=ndm+1
            enddo
          endif

          deallocate(tdvr)

        enddo
        if(myrank.eq.0) then
#ifdef BOXSIM
          if(il.eq.0) then
            read(51) SI_lbox0(0),SI_lbox0(1),SI_lbox0(2)
! Mpc -> 100 kpc
            SI_lbox0(0)=SI_lbox0(0)*10.0d0
            SI_lbox0(1)=SI_lbox0(1)*10.0d0
            SI_lbox0(2)=SI_lbox0(2)*10.0d0
          endif
#endif
          close(51)
        endif
#ifdef BOXSIM
        if(il.eq.0) then
          nval=3
          call MPI_BCAST(SI_lbox0,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! set box size at the initial redshift
          SI_lbox(0)=SI_a*SI_lbox0(0)
          SI_lbox(1)=SI_a*SI_lbox0(1)
          SI_lbox(2)=SI_a*SI_lbox0(2)
          if(myrank.eq.0) then
            write(6,*) ' box size (comoving)=',SI_lbox0(0),SI_lbox0(1) &
             ,SI_lbox0(2)
            write(6,*) ' box size (initial) =',SI_lbox(0),SI_lbox(1) &
            ,SI_lbox(2)
          endif
        endif
#endif

! update ndmt       
        ndmt=ndmt+ndl
        ndt(il)=ndl

      enddo
!endif DM
#endif

! *** for Gas particles
! *** Open Initial Data File for level 1 ***
      npt = 0
      np = 0
#ifdef GAS
      il=Ini_mul
      nval=5

      allocate(tivr(0:nval-1))

      if(myrank.eq.0) then
        write(fin,'(a9,i1,a5)') 'ini/level',il,'p.dat'
        open(51,file=fin,status='old',form='unformatted')
        read(51) nif,idf,npt,itmp,flagrf
! store the range of the particle for each rank
        do ip=0,nprocs-1
          call para_range(0,npt-1,nprocs,ip,jstap(ip),jendp(ip))
        enddo
        tivr(0)=nif
        tivr(1)=idf
        tivr(2)=npt
        tivr(3)=itmp
        tivr(4)=flagrf
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nif=tivr(0)
      idf=tivr(1)
      npt=tivr(2)
      itmp=tivr(3)
      flagrf=tivr(4)
! *** read the data ***
      call para_range(0,npt-1,nprocs,myrank,jsta,jend)
      nplev=jend-jsta+1

      deallocate(tivr)

      call allocate_baryon_all(nplev)
#ifdef MHD
      call allocate_mhd(nplev)
#endif

      nval=5

      allocate(tdvr(0:nval-1))

      if(myrank.eq.0) then
! *** dx (Mpc, Comoving), zstart,omegam,lambda ***
        read(51) dxi(il),SI_a,SI_omg0,SI_lam0,SI_h0
        tdvr(0)=dxi(il)
        tdvr(1)=SI_a
        tdvr(2)=SI_omg0
        tdvr(3)=SI_lam0
        tdvr(4)=SI_h0        
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      dxi(il)=tdvr(0)
      SI_a=tdvr(1)
      SI_omg0=tdvr(2)
      SI_lam0=tdvr(3)
      SI_h0=tdvr(4)
      dxi(il) = dxi(il)*10.0d0

      deallocate(tdvr)

      np=0
      if(flagrf.eq.0) then
        nval=7
      else
        nval=8
      endif
      rrank=0
      if(myrank.eq.0) then
        ipsta=0
        ipend=nprocs-1
      else
        ipsta=myrank
        ipend=myrank
      endif

! *** Read and  Set Position ***
      do ip=ipsta,ipend
        if(myrank.eq.0) then
! number of baryon particle read for rank ip
          nprp=jendp(ip)-jstap(ip)+1
       
          allocate(tdvr(0:nprp*nval-1))

        else

          allocate(tdvr(0:nplev*nval-1))
  
        endif

        if(myrank.eq.0) then
! read the data only by rank 0
! * x, y, z, a(t)b(t)p(q) (Mpc) *		           

          if(flagrf.eq.0) then
            do i=0,nprp-1
              read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5),tdvr(i*nval+6)
            enddo
          else
            do i=0,nprp-1
              read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                 ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5) &
                 ,tdvr(i*nval+6),tdvr(i*nval+7)
            enddo
          endif
          if(ip.ne.myrank) then
            call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
             ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        else
          call MPI_IRECV(tdvr,nplev*nval,MPI_DOUBLE_PRECISION &
           ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
          call MPI_WAIT(ireqr(ip),istatus,ierr)
        endif
        if(myrank.eq.ip) then
          do i=0,nplev-1
            xc_p(np)=tdvr(i*nval)
            yc_p(np)=tdvr(i*nval+1)
            zc_p(np)=tdvr(i*nval+2)
            vnx_p(np)=tdvr(i*nval+3)
            vny_p(np)=tdvr(i*nval+4)
            vnz_p(np)=tdvr(i*nval+5)
            rho_p(np)=tdvr(i*nval+6)
            if(flagrf.ne.0) then
              m_p(np)=tdvr(i*nval+7)
            endif
#if defined(SF_EFD) || defined(SF_ZFD)
            flagfd_p(np)=0
            flagrfd_p(np)=0
            ts_p(np)=0.0d0
#endif
            id_p(np)=jsta+i
            np = np+1
          enddo
        endif

        deallocate(tdvr)

      enddo
      if(myrank.eq.0) then
        close(51)
      endif

#ifdef MHD

      if(myrank.eq.0) then
        open(57,file='ini/gas-mhd.dat',status='old',form='unformatted')
        write(6,*) 'reading ini/gas-mhd.dat'
      endif

! read MHD data
      nc=0
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
! number of baryon particle read for rank ip
          nprp=jendp(ip)-jstap(ip)+1
       
          allocate(tdvr(0:nprp*nval-1))

        else

          allocate(tdvr(0:nplev*nval-1))
  
        endif
        if(myrank.eq.0) then
! read the data only by rank 0
! bx, by, bz
          do i=0,nprp-1
            read(57) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2)
          enddo
          if(ip.ne.myrank) then
            call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
             ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        else
          call MPI_IRECV(tdvr,nplev*nval,MPI_DOUBLE_PRECISION &
           ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
          call MPI_WAIT(ireqr(ip),istatus,ierr)
        endif
        if(myrank.eq.ip) then
          do i=0,nplev-1
            bx_p(nc)=tdvr(i*nval)
            by_p(nc)=tdvr(i*nval+1)
            bz_p(nc)=tdvr(i*nval+2)
            pbx_p(nc)=bx_p(nc)
            pby_p(nc)=by_p(nc)
            pbz_p(nc)=bz_p(nc)
            dbx_p(nc)=0.0d0
            dby_p(nc)=0.0d0
            dbz_p(nc)=0.0d0
            pdbx_p(nc)=0.0d0
            pdby_p(nc)=0.0d0
            pdbz_p(nc)=0.0d0
            divb_p(nc)=0.0d0
            arotb_p(nc)=0.0d0
            agradb_p(nc)=0.0d0
#ifdef TIMEDIS
            alpb_p(nc)=ALPB_MIN
#else
            alpb_p(nc)=ALPHA_CONST
#endif
#ifdef DCLEAN
            psi_p(nc)=0.0d0
            ppsi_p(nc)=0.0d0
            dpsi_p(nc)=0.0d0
            pdpsi_p(nc)=0.0d0
#endif
            nc=nc+1
          enddo
        endif

        deallocate(tdvr)

      enddo
      if(myrank.eq.0) then
        close(57)
      endif
!endif MHD
#endif
! endif Gas
#endif

      SI_h0 = SI_h0/100.0d0
#ifdef GAS
      Ini_OMGg = Ini_OMGg
#else
      Ini_OMGg = 0.0d0
#endif
!*********************************
!  Make Time tabele and Set TM_t0
! ******************************** 
!      write(6,*) ' Ini_ze = ',Ini_ze
      call makettable(1.0d0/(1.0d0+Ini_ze))
! *****   Set Initial Parameter   *****
! *** Initial Time, Set TM_tot ***
      call tfroma(SI_a,TM_tot)
      Ini_zi = 1.0d0/SI_a-1.0d0
      if(myrank.eq.0) then
        write(6,*) ' cosmology: Omg0,lam0,h,omb=',SI_omg0,SI_lam0,SI_h0,SI_omgb
      endif
#ifdef GAS
! read only with rank 0
      if(myrank.eq.0) then
        open(50,file='./ini/recfast.out',status='old')
        read(50,*)
        read(50,'(a1,6(1pE13.5))') tcha,Crec(1),Crec(2),Crec(3) &
         ,Crec(4),Crec(5),Crec(6)
        read(50,*)
        write(6,*) ' Ob, Odm, lam, H0, T0, Yp'
        write(6,*) Crec(1),Crec(2),Crec(3),Crec(4),Crec(5),Crec(6)      
! *** check ***
        if(dabs(Crec(1)-Ini_OMGg).gt.0.001d0) then
          write(6,*) ' Error in reading recfast.out: Omgb'
          write(6,*) ' input OMGg and in recfast=',Ini_OMGg,Crec(1)
          call MPI_FINALIZE(ierr)
          stop
        endif
        if(dabs((Crec(2)-(SI_omg0-Ini_OMGg))).gt.0.001d0) then
          write(6,*) ' Error in reading recfast.out: Omgdm'
          write(6,*) ' Omgdm (sim)',SI_omg0-Ini_OMGg
          call MPI_FINALIZE(ierr)
          stop
        endif
        if(dabs(Crec(3)-SI_lam0).gt.0.001d0) then
          write(6,*) ' Error in reading recfast.out: lam0'
          call MPI_FINALIZE(ierr)
          stop
        endif
        if(dabs(Crec(4)-100.0d0*SI_h0).gt.0.001d0) then
          write(6,*) ' Error in reading recfast.out: h0'
          call MPI_FINALIZE(ierr)
          stop
        endif
        zrec(1) = 10000.0d0
        do i=1,4
          xrec(i,1)=0.0d0
        enddo
        do i=1,2000
          do j=1,4
            xrec(j,2)=xrec(j,1)
          enddo  
          zrec(2) = zrec(1)
! *** read z xe xh xhe Tmat ***
          read(50,'(5(1pE13.5))') zrec(1),xrec(1,1),xrec(2,1),xrec(3,1),xrec(4,1)
          if(zrec(1).lt.Ini_zi) then
            goto 71
          endif
        enddo
        if(myrank.eq.0) then
          write(6,*) ' Error in reading recfast: z(sim)=',Ini_zi 
        endif
   71   if(myrank.eq.0) then
          write(6,*) ' zregrid(1,2), zi=',zrec(1),zrec(2),Ini_zi
        endif
        close(50)
        do i=1,4
          xe(i)=xrec(i,2)+(xrec(i,1)-xrec(i,2))*(Ini_zi-zrec(2)) &
           /(zrec(1)-zrec(2))
        enddo
! *** set initial temperature ***
        Ini_T0 = xe(4)/TUK
        if(myrank.eq.0) then
          write(6,*) ' Temperature   = ',Ini_T0,' (10^4 K)'
          write(6,*) ' xe, xH, xHe = ',xe(1),xe(2),xe(3)
        endif
      endif
! send Ini_T0 and Crec(6)
      call MPI_BCAST(Ini_T0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Crec(6),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
! *** Set Coordinate and Data of particle **
      nag=np
      ng=np
      ns=0	
! *** Set Hubble constant ***
      iai = 1.0d0/SI_a
      hubi = SI_h0*HUB0*dsqrt(SI_omg0*iai*iai*iai &
       -(SI_omg0+SI_lam0-1.0)*iai*iai+SI_lam0)
#ifdef BOXSIM
      SI_hub=hubi
#endif
! *** Calculate Initial Value ***
      rhob0=3.0d0*HUB0*SI_h0*HUB0*SI_h0/(8.0d0*M_PI*G)
      rhoba=rhob0*iai*iai*iai
! *** Calculate Initial Baryon Density ***
      rhob=Ini_OMGg*rhob0*iai*iai*iai
! *** Lambda ***
      SI_lam0 = SI_lam0*HUB0*HUB0*SI_h0*SI_h0

! **********************************************************
!    Set Coordinate and Velocity and Density of particle
! **********************************************************
! *** Invert to this program Unit ***
#ifdef DM
      do i=0,ndm-1
        x_dm(i)=x_dm(i)*(SI_a*10.0d0)
        y_dm(i)=y_dm(i)*(SI_a*10.0d0)
        z_dm(i)=z_dm(i)*(SI_a*10.0d0)
#ifdef BOXSIM
        asc_dm(i)=SI_a
        tc_dm(i)=TM_tot
#endif
      enddo
! *** mass of DM for level 1 particles ***
      nc=0
      do il=Ini_mul,0,-1
        if(il.eq.Ini_mul) then
          mdm = (SI_omg0-Ini_OMGg)*rhob0*(dxi(il)**3)
        else
          mdm = SI_omg0*rhob0*(dxi(il)**3)
        endif
! *****   Set Initial Mass for DM *****
        do i=nc,nc+nd(il)-1
          if(flagrf.eq.0) then
            m_dm(i)=mdm
          endif
! assume density is uniform because tough to get the info
          if(il.eq.Ini_mul) then
            rho_dm(i)=(SI_omg0-Ini_OMGg)*rhoba
          else
            rho_dm(i)=SI_omg0*rhoba
          endif 
        enddo        
        nc = nc+nd(il)
        if(myrank.eq.0) then

          if(flagrf.eq.0) then
            write(6,*) ' mdm (level ',il,') = ',mdm
          else
            write(6,*) ' mdm (level ',il,') = ',m_dm(nc-1)
          endif
        endif
      enddo
! *** Set Initial Velocity ***
! *** Local velocity + Hubble Flow ***
! *** for DM Particles ***
      do i=0,ndm-1      
#ifdef BOXSIM
! v = a^2 vpec, tx = a x vpec
        vnx_dm(i)=(SI_a)*vx_dm(i)/VUKMS
        vny_dm(i)=(SI_a)*vy_dm(i)/VUKMS
        vnz_dm(i)=(SI_a)*vz_dm(i)/VUKMS
#else
        vnx_dm(i)=hubi*x_dm(i)+vx_dm(i)/VUKMS
        vny_dm(i)=hubi*y_dm(i)+vy_dm(i)/VUKMS
        vnz_dm(i)=hubi*z_dm(i)+vz_dm(i)/VUKMS
#endif
      enddo
#endif
#ifdef GAS
! *** mass of Gas ***
      m = Ini_OMGg*rhob0*(dxi(Ini_mul)**3)
! * Calculate physical Value *
      kmmp=(KCGS/(MP*MYU))/K_MU

      if(myrank.eq.0) then      
        write(6,*) ' Rho_b,b0   =',rhob,rhob0
        write(6,*) ' Hub_i = ',hubi
        if(flagrf.eq.0) then
          write(6,*) ' mgas       =',m
        else
          write(6,*) ' mgas       =',m_p(0)
        endif
        write(6,*) ' dxi(0)     =',dxi(0)*SI_a
        write(6,*) ' a (zi)     = ',SI_a
      endif

! *** for gas particle ***
      do i=0,np-1
! *** positions ***
        xc_p(i)=xc_p(i)*(SI_a*10.0d0)
        yc_p(i)=yc_p(i)*(SI_a*10.0d0)
        zc_p(i)=zc_p(i)*(SI_a*10.0d0)
! *** velocities ***
#ifdef BOXSIM
        vnx_p(i)=(SI_a)*vnx_p(i)/VUKMS
        vny_p(i)=(SI_a)*vny_p(i)/VUKMS
        vnz_p(i)=(SI_a)*vnz_p(i)/VUKMS
        asc_p(i)=SI_a
        tc_p(i)=TM_tot
        ascg_p(i)=SI_a
        tcg_p(i)=TM_tot
#else
        vnx_p(i)=hubi*xc_p(i)+vnx_p(i)/VUKMS
        vny_p(i)=hubi*yc_p(i)+vny_p(i)/VUKMS
        vnz_p(i)=hubi*zc_p(i)+vnz_p(i)/VUKMS
#endif
      enddo
! *****   Set Initial Mass and Pressure
!   and Sound Velocity  and Thermal Energy*****
      do i=0,np-1
        if(flagrf.eq.0) then
          m_p(i)=m
        endif
        rho_p(i)=rhob*(1.0d0+rho_p(i))
        p_p(i)=kmmp*rho_p(i)*Ini_T0
! just for test
!        p_p(i)=kmmp*rho_p(i)*(Ini_T0+10.0d0**(dble(i)*1.5e-3-1.0d0))
        cs_p(i)=dsqrt(GAM*p_p(i)/rho_p(i))
        u_p(i)=p_p(i)/((GAM-1.0d0)*rho_p(i))
        myu_p(i)=MYU
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**THIRD)
! * Initialize for star *
#ifdef SF_EFD
        Gsn_p(i) = 0.0d0
#endif
#ifdef STAR
        ts_p(i)=0.0d0
#endif
        list_ap(i)=i
      enddo

! *** Set Initial Metallicity ***
! *** H and He mass fraction w.r.t. w.r.t. the total number of p *
      hnf = 1.0d0-Crec(6)
      henf = Crec(6)
      do i=0,np-1
#ifdef METAL
        mzHe_p(i)=henf*m_p(i)*MUSM
        mzC_p(i)=0.0d0
        mzN_p(i)=0.0d0
        mzO_p(i)=0.0d0
        mzNe_p(i)=0.0d0
        mzMg_p(i)=0.0d0
        mzSi_p(i)=0.0d0
        mzFe_p(i)=0.0d0
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** total ***
        tnsn_p(i)=0.0d0
        tmej_p(i)=0.0d0
        tmzHe_p(i)=0.0d0
        tmzC_p(i)=0.0d0
        tmzN_p(i)=0.0d0
        tmzO_p(i)=0.0d0
        tmzNe_p(i)=0.0d0
        tmzMg_p(i)=0.0d0
        tmzSi_p(i)=0.0d0
        tmzFe_p(i)=0.0d0
        tmzZ_p(i)=0.0d0
#endif
      enddo
! end for GAS
#endif
! ***** end of input for cosmological simulation *****
!else COSM
#else
      Ini_mul=-1
! ***** non-cosmological simulation *****
      SI_a = 1.0d0/(Ini_zi+1.0d0)
! ***  Time bin ***
      if(flagtime.ne.0) then
        TM_tot = TM_tot/TMUGYR
        TM_t0 = TM_t0/TMUGYR
      endif
! * Set Time Parameter *      
      TM_dt = 0.0d0
      Ini_OMGg=0.0d0
#ifdef DM
      if(flagfdm.eq.0) then
        if(myrank.eq.0) then
! *** read DM data only with rank 0 
          write(6,*) ' reading dm.dat'
          open(51,file='ini/dm.dat',status='old',form='unformatted')
          read(51) nif,idf,ndmt,ndm
          if(myrank.eq.0) then
            write(6,*) ' ndmt,ndm=',ndmt,ndm
          endif
! store the range of the particle for each rank
          do ip=0,nprocs-1
            call para_range(0,ndmt-1,nprocs,ip,jstap(ip),jendp(ip))
          enddo
        endif
        call MPI_BCAST(ndmt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call para_range(0,ndmt-1,nprocs,myrank,jsta,jend)
! allocate DM variables
        ndm=jend-jsta+1
        call allocate_dm_all(ndm)

        nc=0
        nval=8
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

            allocate(tdvr(0:ndm*nval-1))
  
          endif

          if(myrank.eq.0) then
! read the data only by rank 0
            do i=0,nprp-1
              read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5) &
                ,tdvr(i*nval+6),tdvr(i*nval+7)
            enddo
            if(ip.ne.myrank) then
              call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
               ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif
          else
            rrank=0
            call MPI_IRECV(tdvr,ndm*nval,MPI_DOUBLE_PRECISION &
             ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(myrank.eq.ip) then
            do i=0,ndm-1
              x_dm(nc)=tdvr(i*nval)
              y_dm(nc)=tdvr(i*nval+1)
              z_dm(nc)=tdvr(i*nval+2)
              vnx_dm(nc)=tdvr(i*nval+3)
              vny_dm(nc)=tdvr(i*nval+4)
              vnz_dm(nc)=tdvr(i*nval+5)
              m_dm(nc)=tdvr(i*nval+6)
              rho_dm(nc)=tdvr(i*nval+7)
              id_dm(nc)=jsta+i
              nc=nc+1
            enddo
          endif

          deallocate(tdvr)
   
        enddo
        if(myrank.eq.0) then
          read(51) mname
          write(6,*) mname
        endif
#ifdef BOXSIM
        if(myrank.eq.0) then
          read(51) SI_lbox(0),SI_lbox(1),SI_lbox(2)
          write(6,*) ' box size=',SI_lbox(0),SI_lbox(1),SI_lbox(2)
        endif
        nval=3
        call MPI_BCAST(SI_lbox,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
        if(myrank.eq.0) then
          close(51)
          write(6,*) ' mdm(0)=',m_dm(0)
        endif
      else if(myrank.eq.0) then
        write(6,*) '##### Warning in set_ini()#####'
        write(6,*) ' DM data not read, but DM is on in gcdp.def.'
        write(6,*) ' This likely causes a problem.'
      endif
#endif

! *** for Gas particles
! *** Open Initial Data File for level 1 ***
      ng=0
      ns=0
      np=0
#ifdef GAS
! get ng for each rank
      if(flagfg.le.0) then
        if(myrank.eq.0) then
          write(6,*) ' reading gas.dat. flagfg=',flagfg
          open(51,file='ini/gas.dat',status='old',form='unformatted')
          read(51) nif,idf,ngt,ng,flagrf
          close(51) 
          write(6,*) ' nif,idf,ngt,ng,flagrf=',nif,idf,ngt,ng,flagrf
        endif
        call MPI_BCAST(ngt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call para_range(0,ngt-1,nprocs,myrank,jsta,jend)
        ng=jend-jsta+1
      endif
#endif
! get ns for each rank
#ifdef STAR
      SI_tnsf=0
      if(flagfs.eq.0) then
        if(myrank.eq.0) then
          open(51,file='ini/star.dat',status='old',form='unformatted')
          write(6,*) ' reading star.dat. flagfs=',flagfs
          read(51) nif,idf,nst,ns
          close(51)
          write(6,*) ' nif,idf,nst,ns=',nif,idf,nst,ns
        endif      
        call MPI_BCAST(nst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call para_range(0,nst-1,nprocs,myrank,jsta,jend)
        ns=jend-jsta+1
      endif
#endif
#if defined(GAS) || defined(STAR)
      np=ng+ns

! initial allocate memory space
! all baryon
      call allocate_baryon_all(np)
#ifdef MHD
      call allocate_mhd(np)
#endif
#endif
#ifdef GAS
! *** Read and Set Position and basic values ***
      if(flagfg.le.0) then
        nval=5

        allocate(tivr(0:nval-1))

        if(myrank.eq.0) then
          open(51,file='ini/gas.dat',status='old',form='unformatted')
          read(51) nif,idf,ngt,ng,flagrf
! store the range of the particle for each rank
          do ip=0,nprocs-1
            call para_range(0,ngt-1,nprocs,ip,jstap(ip),jendp(ip))
          enddo
          tivr(0)=nif
          tivr(1)=idf
          tivr(2)=ngt
          tivr(3)=ng
          tivr(4)=flagrf
        endif
        call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        nif=tivr(0)
        idf=tivr(1)
        ngt=tivr(2)
        ng=tivr(3)
        flagrf=tivr(4)
        call para_range(0,ngt-1,nprocs,myrank,jsta,jend)
        ng=jend-jsta+1

        deallocate(tivr)

! flagrf is not 0: read flagfd_p info from gas.dat
! since pv31.19
#ifndef STAR
        if(myrank.eq.0.and.flagrf.ne.0) then
          write(6,*) ' set flagrf=0'
        endif    
        flagrf=0
#endif    
#ifdef STAR
        if(flagrf.eq.0) then
          nval=9
          nival=0
        else
          nval=10
          nival=1
        endif
#else
        nval=9
#endif
! *** Read and  Set Position ***
        nc=0
! set do loop range
        if(myrank.eq.0) then
          ipsta=0
          ipend=nprocs-1
        else
          ipsta=myrank
          ipend=myrank
        endif
        rrank=0
        do ip=ipsta,ipend
          if(myrank.eq.0) then
! number of baryon particle read for rank ip
            nprp=jendp(ip)-jstap(ip)+1

            allocate(tdvr(0:nprp*nval-1))
#ifdef STAR
            if(flagrf.ne.0) then
              allocate(tivr(0:nprp*nival-1))
            endif
#endif

          else

            allocate(tdvr(0:ng*nval-1))
#ifdef STAR
            if(flagrf.ne.0) then
              allocate(tivr(0:ng*nival-1))
            endif
#endif
  
          endif
          if(myrank.eq.0) then
! read the data only by rank 0
#ifdef STAR
            if(flagrf.eq.0) then 
#endif
              do i=0,nprp-1
                read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                  ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5) &
                  ,tdvr(i*nval+6),tdvr(i*nval+7),tdvr(i*nval+8)
              enddo
#ifdef STAR
            else
              do i=0,nprp-1
                read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                  ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5) &
                  ,tdvr(i*nval+6),tdvr(i*nval+7),tdvr(i*nval+8) &
                  ,tdvr(i*nval+9),tivr(i*nival)
              enddo
            endif 
#endif
            if(ip.ne.myrank) then
              call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
               ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
#ifdef STAR
              if(flagrf.ne.0) then 
                call MPI_ISEND(tivr,nprp*nival,MPI_INTEGER &
                 ,ip,1,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
#endif
            endif
          else
            call MPI_IRECV(tdvr,ng*nval,MPI_DOUBLE_PRECISION &
             ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
#ifdef STAR
            if(flagrf.ne.0) then 
              call MPI_IRECV(tivr,ng*nival,MPI_INTEGER &
               ,rrank,1,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
#endif
          endif
          if(myrank.eq.ip) then     
            do i=0,ng-1
              xc_p(nc)=tdvr(i*nval)
              yc_p(nc)=tdvr(i*nval+1)
              zc_p(nc)=tdvr(i*nval+2)
              vnx_p(nc)=tdvr(i*nval+3)
              vny_p(nc)=tdvr(i*nval+4)
              vnz_p(nc)=tdvr(i*nval+5)
              m_p(nc)=tdvr(i*nval+6)
              rho_p(nc)=tdvr(i*nval+7)
              u_p(nc)=tdvr(i*nval+8)
              id_p(nc)=jsta+i
#if defined(SF_EFD) || defined(SF_ZFD)
              if(flagrf.eq.0) then
                flagfd_p(nc)=0
                flagrfd_p(nc)=0
                ts_p(nc)=0.0d0
              else
                ts_p(nc)=tdvr(i*nval+9)
                flagfd_p(nc)=tivr(i*nival)
                if(flagfd_p(nc).lt.0) then
                  flagrfd_p(nc)=1
                else
                  flagrfd_p(nc)=0
                endif
              endif
#endif
              nc=nc+1
            enddo
          endif

          deallocate(tdvr)
#ifdef STAR
          if(flagrf.ne.0) then
            deallocate(tivr)
          endif
#endif  
        enddo
        if(nc.ne.ng) then
          write(6,*) ' Error in reading gas.dat at rank',myrank
          write(6,*) ' ng in gas.dat gas.dat=',ng,nc
          call MPI_FINALIZE(ierr)
          stop
        endif

        if(myrank.eq.0) then
          read(51) mname
          if(myrank.eq.0) then
            write(6,*) mname
          endif
          read(51) gamf
          if(dabs(GAM-gamf)/GAM.gt.0.1d0) then
            if(myrank.eq.0) then
              write(6,*) 
              if(myrank.eq.0) then
                write(6,*) ' Error(): gamma is incosnsitent'
                write(6,*) ' ini/input.dat, ini/gas.dat',GAM,gamf
              endif
              call MPI_FINALIZE()
              stop
            endif
          endif    
        endif
#ifdef BOXSIM
        if(myrank.eq.0) then
          read(51) SI_lbox(0),SI_lbox(1),SI_lbox(2)
          write(6,*) ' box size=',SI_lbox(0),SI_lbox(1),SI_lbox(2)
        endif
        nval=3
        call MPI_BCAST(SI_lbox,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

#ifdef BCX_FIX
#ifdef MHD
        nval=16
#else
        nval=6
#endif

        allocate(tdvr(0:nval-1))

        if(myrank.eq.0) then
          read(51) BCX_ndx
          tivr(0)=BCX_ndx
        endif
        call MPI_BCAST(BCX_ndx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then
#ifdef MHD
          read(51) BCX_rhox(0),BCX_px(0),BCX_vx(0),BCX_vy(0),BCX_vz(0)
          read(51) BCX_rhox(1),BCX_px(1),BCX_vx(1),BCX_vy(1),BCX_vz(1)
          read(51) BCX_bx(0),BCX_by(0),BCX_bz(0)
          read(51) BCX_bx(1),BCX_by(1),BCX_bz(1)
          tdvr(0)=BCX_rhox(0)
          tdvr(1)=BCX_px(0)
          tdvr(2)=BCX_vx(0)
          tdvr(3)=BCX_vy(0)
          tdvr(4)=BCX_vz(0)
          tdvr(5)=BCX_rhox(1)
          tdvr(6)=BCX_px(1)
          tdvr(7)=BCX_vx(1)
          tdvr(8)=BCX_vy(1)
          tdvr(9)=BCX_vz(1)
          tdvr(10)=BCX_bx(0)
          tdvr(11)=BCX_by(0)
          tdvr(12)=BCX_bz(0)
          tdvr(13)=BCX_bx(1)
          tdvr(14)=BCX_by(1)
          tdvr(15)=BCX_bz(1)
#else
          read(51) BCX_rhox(0),BCX_px(0),BCX_vx(0)
          read(51) BCX_rhox(1),BCX_px(1),BCX_vx(1)
          tdvr(0)=BCX_rhox(0)
          tdvr(1)=BCX_px(0)
          tdvr(2)=BCX_vx(0)
          tdvr(3)=BCX_rhox(1)
          tdvr(4)=BCX_px(1)
          tdvr(5)=BCX_vx(1)
#endif
        endif
        call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#ifdef MHD
        BCX_rhox(0)=tdvr(0)
        BCX_px(0)=tdvr(1)
        BCX_vx(0)=tdvr(2)
        BCX_vy(0)=tdvr(3)
        BCX_vz(0)=tdvr(4)
        BCX_rhox(1)=tdvr(5)
        BCX_px(1)=tdvr(6)
        BCX_vx(1)=tdvr(7)
        BCX_vy(1)=tdvr(8)
        BCX_vz(1)=tdvr(9)
        BCX_bx(0)=tdvr(10)
        BCX_by(0)=tdvr(11)
        BCX_bz(0)=tdvr(12)
        BCX_bx(1)=tdvr(13)
        BCX_by(1)=tdvr(14)
        BCX_bz(1)=tdvr(15)
        BCX_ndx=5
#else
        BCX_rhox(0)=tdvr(0)
        BCX_px(0)=tdvr(1)
        BCX_vx(0)=tdvr(2)
        BCX_rhox(1)=tdvr(3)
        BCX_px(1)=tdvr(4)
        BCX_vx(1)=tdvr(5)
#endif
        if(myrank.eq.0) then
          write(6,*) ' boundary is fixed from the edge+-dx=',BCX_ndx
        endif
!ifdef BCX_FIX
#endif
!ifdef BOX_SIM
#endif
        if(myrank.eq.0) then
          close(51)
        endif
        if(myrank.eq.0) then
          write(6,*) ' mgas(0)=',m_p(0),myrank
        endif
      else
        if(myrank.eq.0) then
          write(6,*) &
           ' Error: not reading gas.dat, but GAS is defined.'
          write(6,*) ' input flagfg=',flagfg
        endif
        call MPI_FINALIZE(ierr)
        stop
      endif
#ifdef METAL
! *** reading metal data ***
      if(flagfg.eq.-1.or.flagfg.le.-3) then
        if(myrank.eq.0) then
          write(6,*) ' reading ini/gas-metal.dat'
          open(51,file='ini/gas-metal.dat',status='old',form='unformatted')
        endif
        nval=9
! *** Read metal info ***
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

            allocate(tdvr(0:ng*nval-1))
  
          endif
          if(myrank.eq.0) then
! read the data only by rank 0
            do i=0,nprp-1
              read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5) &
                ,tdvr(i*nval+6),tdvr(i*nval+7),tdvr(i*nval+8)
            enddo
            if(ip.ne.myrank) then
              call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
               ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif
          else
            call MPI_IRECV(tdvr,ng*nval,MPI_DOUBLE_PRECISION &
             ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(myrank.eq.ip) then
            do i=0,ng-1
              mzHe_p(i)=tdvr(i*nval)
              mzC_p(i)=tdvr(i*nval+1)
              mzN_p(i)=tdvr(i*nval+2)
              mzO_p(i)=tdvr(i*nval+3)
              mzNe_p(i)=tdvr(i*nval+4)
              mzMg_p(i)=tdvr(i*nval+5)
              mzSi_p(i)=tdvr(i*nval+6)
              mzFe_p(i)=tdvr(i*nval+7)
              mzZ_p(i)=tdvr(i*nval+8)
            enddo
          endif

          deallocate(tdvr)

        enddo
        if(myrank.eq.0) then
          close(51)
        endif
      endif


#else
      if(flagfg.eq.-1.or.flagfg.le.-3) then
        if(myrank.eq.0) then
          write(6,*) ' Error: flagfg is inconsistent=',flagfg
          write(6,*) ' METAL is NOT defined, but trying to read gas-metal.dat'
        endif
        call MPI_FINALIZE(ierr)
        stop
      endif
#endif
      np=ng  
#ifdef METAL
      if(flagfg.ne.-1.and.flagfg.gt.-3) then
! *** assume primordial abundance ***
        do i=0,ng-1
          mzHe_p(i)=0.24d0*MUSM*m_p(i)
          mzC_p(i)=0.0d0
          mzN_p(i)=0.0d0
          mzO_p(i)=0.0d0
          mzNe_p(i)=0.0d0
          mzMg_p(i)=0.0d0
          mzSi_p(i)=0.0d0
          mzFe_p(i)=0.0d0
          mzZ_p(i)=0.0d0
        enddo
      endif
#endif

#ifdef MHD
      if(myrank.eq.0) then
! *** reading B data ***
        open(51,file='ini/gas-mhd.dat',status='old',form='unformatted')
        write(6,*) 'reading ini/gas-mhd.dat'
      endif 
      nval=3
! *** Read metal info ***
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

          allocate(tdvr(0:ng*nval-1))

        endif
        if(myrank.eq.0) then
! read the data only by rank 0
          do i=0,nprp-1
            read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2)
          enddo
          if(ip.ne.myrank) then
            call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
             ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        else
          call MPI_IRECV(tdvr,ng*nval,MPI_DOUBLE_PRECISION &
           ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
          call MPI_WAIT(ireqr(ip),istatus,ierr)
        endif
        if(myrank.eq.ip) then
          do i=0,ng-1
            bx_p(i)=tdvr(i*nval)
            by_p(i)=tdvr(i*nval+1)
            by_p(i)=tdvr(i*nval+2)
            pbx_p(i)=bx_p(i)
            pby_p(i)=by_p(i)
            pbz_p(i)=bz_p(i)
            dbx_p(i)=0.0d0
            dby_p(i)=0.0d0
            dbz_p(i)=0.0d0
            pdbx_p(i)=0.0d0
            pdby_p(i)=0.0d0
            pdbz_p(i)=0.0d0
            divb_p(i)=0.0d0
            arotb_p(i)=0.0d0
#ifdef TIMEDIS
            alpb_p(i)=ALPB_MIN
#else
            alpb_p(i)=ALPHA_CONST
#endif
#ifdef DCLEAN
            psi_p(i)=0.0d0
            ppsi_p(i)=0.0d0
            dpsi_p(i)=0.0d0
            pdpsi_p(i)=0.0d0
#endif
          enddo
        endif
    
        deallocate(tdvr)
     
      enddo
      if(myrank.eq.0) then
        close(51)
      endif
!endif MHD
#endif

! ***** end of reading files *****
! *** the other properties **
      do i=0,ng-1
#ifdef SF_EFD
        Gsn_p(i) = 0.0d0
#endif
        flagc_p(i) = 0
        p_p(i) = (GAM-1.0d0)*rho_p(i)*u_p(i)
        cs_p(i)=dsqrt(GAM*p_p(i)/rho_p(i))
        myu_p(i)=MYU
#ifdef STAR
        ts_p(i)=0.0d0
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** total ***
        tnsn_p(i)=0.0d0
        tmej_p(i)=0.0d0
        tmzHe_p(i)=0.0d0
        tmzC_p(i)=0.0d0
        tmzN_p(i)=0.0d0
        tmzO_p(i)=0.0d0
        tmzNe_p(i)=0.0d0
        tmzMg_p(i)=0.0d0
        tmzSi_p(i)=0.0d0
        tmzFe_p(i)=0.0d0
        tmzZ_p(i)=0.0d0
#endif
! *** particle list ***
        list_ap(i) = i
      enddo
#ifdef STAR
      SI_tnsf=0
      if(flagfs.eq.0) then
        nval=4

        allocate(tivr(0:nval-1))

        if(myrank.eq.0) then
          write(6,*) ' reading star.dat. flagfs=',flagfs
          open(51,file='ini/star.dat',status='old',form='unformatted')
          read(51) nif,idf,nst,ns
          write(6,*) ' nif,idf,nst,ns=',nif,idf,nst,ns
! store the range of the particle for each rank
          do ip=0,nprocs-1
            call para_range(0,nst-1,nprocs,ip,jstap(ip),jendp(ip))
          enddo
          tivr(0)=nif
          tivr(1)=idf
          tivr(2)=nst
          tivr(3)=ns
        endif
        call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        nif=tivr(0)
        idf=tivr(1)
        nst=tivr(2)
        ns=tivr(3)
        call para_range(0,nst-1,nprocs,myrank,jsta,jend)
        ns=jend-jsta+1
! keep the total number of formed star 
        SI_tnsf=nst

        deallocate(tivr)

        nval=10
        nival=1
! *** Read and  Set Position ***
        nc=ng
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
            allocate(tivr(0:nprp*nival-1))

          else

            allocate(tdvr(0:ns*nval-1))
            allocate(tivr(0:ns*nival-1))
  
          endif
          if(myrank.eq.0) then
! read the data only by rank 0
            do i=0,nprp-1
              read(51) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5) &
                ,tdvr(i*nval+6),tdvr(i*nval+7),tdvr(i*nval+8) &
                ,tdvr(i*nval+9),tivr(i*nival)
            enddo
            if(ip.ne.myrank) then
              call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
               ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
              call MPI_ISEND(tivr,nprp*nival,MPI_INTEGER &
               ,ip,1,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif
          else
            call MPI_IRECV(tdvr,ns*nval,MPI_DOUBLE_PRECISION &
             ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
            call MPI_IRECV(tivr,ns*nival,MPI_INTEGER &
             ,rrank,1,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(myrank.eq.ip) then
            do i=0,ns-1
              xc_p(nc)=tdvr(i*nval)
              yc_p(nc)=tdvr(i*nval+1)
              zc_p(nc)=tdvr(i*nval+2)
              vnx_p(nc)=tdvr(i*nval+3)
              vny_p(nc)=tdvr(i*nval+4)
              vnz_p(nc)=tdvr(i*nval+5)
              m_p(nc)=tdvr(i*nval+6)
! tdvr(i*nval+7):
! in case of stellar mass needs two values, like initial and current
! in current version, stellar mass does not change, therefore this is not used
              rho_p(nc)=tdvr(i*nval+8)
! *** age -> formation time ***
              ts_p(nc)=TM_tot-tdvr(i*nval+9)
#if defined(SF_EFD) || defined(SF_ZFD)
              flagfd_p(nc)=tivr(i*nival)
#endif
              id_p(nc)=ngt+jsta+i
              nc=nc+1
            enddo
          endif

          deallocate(tdvr)
          deallocate(tivr)

        enddo
        if(myrank.eq.0) then
          read(51) mname
          close(51)
          write(6,*) ' mstar(0)=',m_p(ng)
        endif

#ifdef METAL
        if(myrank.eq.0) then
          open(52,file='ini/star-metal.dat',status='old',form='unformatted')
        endif
        nval=9
! *** Read metal info ***
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

            allocate(tdvr(0:ns*nval-1))
  
          endif
          if(myrank.eq.0) then
! read the data only by rank 0
            do i=0,nprp-1
              read(52) tdvr(i*nval),tdvr(i*nval+1),tdvr(i*nval+2) &
                ,tdvr(i*nval+3),tdvr(i*nval+4),tdvr(i*nval+5) &
                ,tdvr(i*nval+6),tdvr(i*nval+7),tdvr(i*nval+8)
            enddo
            if(ip.ne.myrank) then
              call MPI_ISEND(tdvr,nprp*nval,MPI_DOUBLE_PRECISION &
               ,ip,0,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif
          else
            call MPI_IRECV(tdvr,ns*nval,MPI_DOUBLE_PRECISION &
             ,rrank,0,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(myrank.eq.ip) then
            do i=0,ns-1
              mzHe_p(i+ng)=tdvr(i*nval)
              mzC_p(i+ng)=tdvr(i*nval+1)
              mzN_p(i+ng)=tdvr(i*nval+2)
              mzO_p(i+ng)=tdvr(i*nval+3)
              mzNe_p(i+ng)=tdvr(i*nval+4)
              mzMg_p(i+ng)=tdvr(i*nval+5)
              mzSi_p(i+ng)=tdvr(i*nval+6)
              mzFe_p(i+ng)=tdvr(i*nval+7)
              mzZ_p(i+ng)=tdvr(i*nval+8)
            enddo
          endif

          deallocate(tdvr)

        enddo
        if(myrank.eq.0) then
          close(52)
        endif
#endif
        np=ng+ns

        do i=ng,np-1
! *** the other properties **
#ifdef SF_EFD
          Gsn_p(i) = 0.0d0
#endif
          flagc_p(i) = 1
          flagt_p(i) = 0
! *** set zero for unused data ***
          dtr_p(i)=INF
          p_p(i) = 0.0d0
          cs_p(i) = 0.0d0
          u_p(i) = 0.0d0
          myu_p(i)=0.0d0
! *** particle list ***
          list_ap(i) = i
#if defined(SF_EFD) || defined(SF_ZFD)
          flagrfd_p(i)=0
#endif
        enddo
      endif

#ifdef METAL
! initialise all the metal 
      do i=0,np-1
        mzHe0_p(i)=mzHe_p(i)
        mzC0_p(i)=mzC_p(i)
        mzN0_p(i)=mzN_p(i)
        mzO0_p(i)=mzO_p(i)
        mzNe0_p(i)=mzNe_p(i)
        mzMg0_p(i)=mzMg_p(i)
        mzSi0_p(i)=mzSi_p(i)
        mzFe0_p(i)=mzFe_p(i)
        mzZ0_p(i)=mzZ_p(i)
      enddo
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** get yield even if there is no star, but for feedback gas ***
      call setyield()

      allocate(tmej_sn(0:np))
      allocate(tnsn_sn(0:np))
      allocate(tmzHe_sn(0:np))
      allocate(tmzC_sn(0:np))
      allocate(tmzN_sn(0:np))
      allocate(tmzO_sn(0:np))
      allocate(tmzNe_sn(0:np))
      allocate(tmzMg_sn(0:np))
      allocate(tmzSi_sn(0:np))
      allocate(tmzFe_sn(0:np))
      allocate(tmzZ_sn(0:np))
      allocate(nsw_sn(0:np))

      TM_dt=0.0d0
      TMsf_dt=0.0d0

      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)
! *** set tmz*_p for both star and gas ***
      do i = 0,ng+ns-1
        if((flagfd_p(i).ge.1.and.flagfd_p(pn).le.SI_snii) &
         .or.(flagfd_p(i).gt.SI_snii.and.flagrfd_p(i).eq.1)) then
          flagrfd_p(i)=0
! since pv32.9
!          flagfd_p(i)=-flagfd_p(i)
        endif
! * update tm*_p() *
! * Unit MUSM *
        tmej_p(i)=tmej_p(i)
        tnsn_p(i)=tnsn_sn(i)
! * Unit Solar Mass *
        tmzHe_p(i)=tmzHe_sn(i)
        tmzC_p(i)=tmzC_sn(i)
        tmzN_p(i)=tmzN_sn(i)
        tmzO_p(i)=tmzO_sn(i)
        tmzNe_p(i)=tmzNe_sn(i)
        tmzMg_p(i)=tmzMg_sn(i)
        tmzSi_p(i)=tmzSi_sn(i)
        tmzFe_p(i)=tmzFe_sn(i)
        tmzZ_p(i)=tmzZ_sn(i)
! * check mass *
      enddo
! check feedback energy
#ifdef SF_EFD
      TMsf_dt=DTSF
      TM_dt=TMsf_dt
      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)

      allocate(mfd(0:np))

      do i=0,ng+ns-1
        pn=i
        if(flagfd_p(pn).lt.0) then
! *** number of supernovae ***
          mfd(pn)=tnsn_sn(pn)-tnsn_p(pn)
          if(mfd(pn).lt.0.0d0.and.dabs(mfd(pn)).lt.MININF*tnsn_p(pn)) then
            mfd(pn)=0.0d0
          endif
          if(nsw_sn(pn).lt.0.0d0) then
            nsw_sn(pn)=0.0d0
          endif
          Gsn_p(pn)=mfd(pn)*SNEU+SWEU*nsw_sn(pn)
#ifdef CONSTFDE
! *** erg/s unit 5.752e43
          Gsn_p(pn)=(1.739e-6)*TMsf_dt
#endif
          Gsn_p(pn)=Gsn_p(pn)/m_p(pn)
        endif
      enddo

      deallocate(mfd)

#endif

      deallocate(tmej_sn)
      deallocate(tnsn_sn)
      deallocate(tmzHe_sn)
      deallocate(tmzC_sn)
      deallocate(tmzN_sn)
      deallocate(tmzO_sn)
      deallocate(tmzNe_sn)
      deallocate(tmzMg_sn)
      deallocate(tmzSi_sn)
      deallocate(tmzFe_sn)
      deallocate(tmzZ_sn)
      deallocate(nsw_sn)

#endif
      TM_dt=0.0d0
! *** end for STAR
#endif
! *** end for GAS
#endif
      npt=ngt+nst
! *** end of non-cosmological simulation input ***
#endif
! *** end of data setting ***

! *** Set the data in common for cosmological and non-cosmological mode ***
   92 nag=ng
      nadm=ndm
      nas=ns
      nagravg=ng
      np=ng+ns
      step=0
      lstep=0
      totstep = 0
! ***  Time bin ***
! * Set Time Parameter. This is used at ddecb() for setting tm*_p *
      TM_dt=0.0d0 
! ***   about Output   ***
      OF_odt = (TM_t0-TM_tot)/dble(Ini_nd)
      OF_todt = OF_odt+TM_tot

      if(myrank.eq.0) then
        write(6,*) '#####   Initial Condition Information   #####'
#ifdef COSM
        write(6,*) ' h, Omg0, lamda0 = ',SI_h0,SI_omg0,SI_lam0
        write(6,*) ' OMGg        = ',Ini_OMGg
#endif
        write(6,*) ' total Nb Ng, Ns Ndm = ',npt,ngt,nst,ndmt
#ifdef MULTI_LEV
        do i=0,Ini_mul
          write(6,*) i,' level n,dx  = ',ndt(i),dxi(i)
        enddo
#endif
        write(6,*) ' gamma       = ',GAM
        write(6,*) ' simulate z  = ',1.0d0/SI_a-1.0d0,'->',Ini_ze
        write(6,*) '          t  = ',TM_tot,'->',TM_t0
        write(6,*) '  output dt  = ',OF_odt
        write(6,*) ' # of output = ',Ini_nd
        write(6,*) ' # of output file for each output=',SI_nof
        if(SI_nof.gt.nprocs) then
          write(6,*) ' Error: It has to be less than nprocs=',nprocs
          stop
        endif
        if(Ini_mul.eq.0) then
          write(6,*) '   level 0 Only'
        endif
      endif        

! *** setting the other parameters ***
#ifdef STAR
! *** for star ***	
      do i=ng,ng+ns-1
        x_p(i)=xc_p(i)
        y_p(i)=yc_p(i)
        z_p(i)=zc_p(i)
        vx_p(i)=vnx_p(i)
        vy_p(i)=vny_p(i)
        vz_p(i)=vnz_p(i)
#if defined(COSM) && defined(BOXSIM)
        vvnx_p(i)=vnx_p(i)/SI_a
        vvny_p(i)=vny_p(i)/SI_a
        vvnz_p(i)=vnz_p(i)/SI_a
#else
        vvnx_p(i)=vnx_p(i)
        vvny_p(i)=vny_p(i)
        vvnz_p(i)=vnz_p(i)
#endif
#ifdef SIM1D
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i)))
#elif defined(SIM2D)
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**0.5d0)
#else
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**THIRD)
#endif
        omgh_p(i)=1.0d0
        zetah_p(i)=0.0d0
        flagt_p(i)=0
        dtr_p(i)=INF
        dt_p(i)=INF
        dtmnb_p(i)=INF
        dtg_p(i)=INF
        ax_p(i)=0.0d0
        ay_p(i)=0.0d0
        az_p(i)=0.0d0
      enddo
#endif
#ifdef DM
! *** for DM ***	
      do i=0,ndm-1
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
#endif
#ifdef GAS
      do i=0,ng-1
        x_p(i)=xc_p(i)
        y_p(i)=yc_p(i)
        z_p(i)=zc_p(i)
        vx_p(i)=vnx_p(i)
        vy_p(i)=vny_p(i)
        vz_p(i) = vnz_p(i)
#if defined(COSM) && defined(BOXSIM)
        vvnx_p(i)=vnx_p(i)/SI_a
        vvny_p(i)=vny_p(i)/SI_a
        vvnz_p(i)=vnz_p(i)/SI_a
#else
        vvnx_p(i)=vnx_p(i)
        vvny_p(i)=vny_p(i)
        vvnz_p(i)=vnz_p(i)
#endif

#ifdef SIM1D
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i)))
#elif defined(SIM2D)
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**0.5d0)
#else
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**THIRD)
#endif
        omgh_p(i)=1.0d0
        zetah_p(i)=0.0d0
        as_p(i)=p_p(i)/(rho_p(i)**GAM)
      enddo
! *** For Gas, Set flagc ***
      do i=0,ng-1
        flagc_p(i)=0
        flagt_p(i)=0
        dtr_p(i)=INF
        dt_p(i)=INF
        dtmnb_p(i)=INF
!        alpv_p(i)=V_ALPHAMAX
        alpv_p(i)=V_ALPHAINI
!        alpu_p(i)=U_ALPHAMAX
        alpu_p(i)=0.0d0
        hvsigdt_p(i)=INF
        dtg_p(i)=INF
        ax_p(i)=0.0d0
        ay_p(i)=0.0d0
        az_p(i)=0.0d0
#ifdef METALDIFF
        vsig_p(i)=0.0d0
#endif
      enddo
#endif      

#ifdef METAL
! initialise all the metal 
      do i=0,np-1
        mzHe0_p(i)=mzHe_p(i)
        mzC0_p(i)=mzC_p(i)
        mzN0_p(i)=mzN_p(i)
        mzO0_p(i)=mzO_p(i)
        mzNe0_p(i)=mzNe_p(i)
        mzMg0_p(i)=mzMg_p(i)
        mzSi0_p(i)=mzSi_p(i)
        mzFe0_p(i)=mzFe_p(i)
        mzZ0_p(i)=mzZ_p(i)
      enddo
#endif

#ifdef COSM
      if(Ini_mul.gt.0) then
        SI_ndm1=nd(Ini_mul)
        SI_ndm1t=ndt(Ini_mul)
      else
#endif
        SI_ndm1t=ndmt
        SI_ndm1=ndm
#ifdef COSM
      endif
#endif
      flagcont=0
      SI_flagout = 0
      SI_mms = 0.0d0
#ifdef STAR
      SI_esns = 0.0d0
      SI_zsns = 0.0d0
! *** star formation time ***
      TMsf_t=TM_tot
#endif
      idum = -111      

#ifdef SGRAV
! *** set softening length limit ***
#ifdef SF
      SI_eps=(XHSOL*(DU/MP)/NSTH)**(1.0d0/3.0d0)*ETAH
#else
      SI_eps=(XHSOL*(DU/MP)/NSTH)**(1.0d0/3.0d0)*ETAH
#endif
      if(myrank.eq.0) then
        write(6,*) ' softening limit = mp^1/3 x',SI_eps
#if defined(GAS) || defined(STAR)
        write(6,*) '  for gas and/or star ',SI_eps*(m_p(0)**THIRD)
        write(6,*) '  corresponding nh(Zsolar)=' &
         ,0.706d0*(DU/MP)*((ETAH/SI_eps)**3)
#endif
#if defined(GEPSLIM_G) && defined(GAS)
        write(6,*) '  minimu h=',FH*SI_eps*(m_p(0)**THIRD)
#endif
#ifdef DM
        write(6,*) '  for DM ',SI_eps*(m_dm(0)**THIRD)
#endif
      endif
! *** check the minimum h ***
      hmin=INF
#if defined(GAS) || defined(STAR)
      do i=0,np-1
        if(h_p(i)/(m_p(i)**THIRD).lt.hmin) then 
          hmin=h_p(i)/(m_p(i)**THIRD)
        endif
      enddo
#endif
#ifdef DM
      do i=0,ndm-1
        if(h_dm(i)/(m_dm(i)**THIRD).lt.hmin) then 
          hmin=h_dm(i)/(m_dm(i)**THIRD)
       endif
      enddo
#endif
      hmint=INF
      call MPI_ALLREDUCE(hmin,hmint,1,MPI_DOUBLE_PRECISION &
       ,MPI_MIN,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0) then
        write(6,*) ' minimum h/(m^1/3)=',hmin
      endif
#ifdef COSM
      SI_eps0=SI_eps
! *** set redshift when SI_eps becomes effective ***
! *** assume mass of DM is always bigger ***
#ifdef DM
      etahmax=ETAHDM
! *** try slightly higher density ***
      SI_epsa=SI_eps*2.0d0*(m_dm(0)**THIRD)/(dxi(Ini_mul)*etahmax)
#else
      etahmax=ETAH
      SI_epsa=SI_eps*2.0d0*(m_p(0)**THIRD)/(dxi(Ini_mul)*etahmax)
#endif
      if(myrank.eq.0) then
        write(6,*) ' SI_epsa=',SI_epsa
      endif
      if(SI_a.lt.SI_epsa) then
        SI_eps=(1.0d0-SI_epsa+SI_a)*SI_eps
      endif
      if(hmint.lt.SI_eps) then
        if(myrank.eq.0) then
          write(6,*) ' Warning in set_ini(): SI_eps may be too big'
          write(6,*) ' hmint,SI_eps=',hmint,SI_eps
#ifdef GAS
          write(6,*) ' m/rho**THIRD=',(m_p(0)/rho_p(0))**THIRD
#else
          write(6,*) ' m/rho**THIRD=',(m_dm(1)/rho_dm(1))**THIRD
#endif
        endif
!        call MPI_ABORT(MPI_COMM_WORLD,ierr)
!        stop
      endif
#endif
! *** fixed softening case.
#ifdef DM
#ifdef FIXEDDMEPS
! *** set to minimum eps 
      do i=0,ndm-1
        h_dm(i)=SI_eps*(m_dm(i)**THIRD)
      enddo
#endif
#endif
#ifdef STAR
! *** set softening for stars
      do i=ng,ng+ns-1
        pn=list_ap(i)
        h_p(pn)=SI_eps*(m_p(pn)**THIRD)
        omgh_p(pn)=1.0d0
        zetah_p(pn)=0.0d0
        myu_p(pn)=MYU
      enddo
#endif
#else
      SI_eps=0.0d0
#endif

#ifdef FIXEDP
      open(50,file='ini/fixed.dat',status='old')      
      read(50,*) SI_fid(0),SI_fid(1)
      if(myrank.eq.0) then
        write(6,*) ' Fixed particles id between ',SI_fid(0),SI_fid(1)
      endif
      if(SI_fid(0).gt.SI_fid(1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1.'
        endif
        stop
      endif
! *** change nag and nagravg
      nag=0
      nc=0
      do i=0,ng-1
        pn=list_ap(i)
        if(id_p(pn).ge.SI_fid(0).and.id_p(pn).le.SI_fid(1)) then
          nalist(nc)=pn
          nc=nc+1
! *** set dvx and etc.
          dvx_p(pn)=0.0d0
          dvy_p(pn)=0.0d0
          dvz_p(pn)=0.0d0
          ndu_p(pn)=0.0d0
          pdu_p(pn)=0.0d0
        else
          list(nag)=pn
          nag=nag+1
        endif
      enddo
      nagravg=nag
! *** change list_ap putting fixedp at the end
      do i=0,nag-1
        list_ap(i)=list(i)
      enddo
      do i=0,nc-1
        list_ap(nag+i)=nalist(i)
      enddo
#endif

#ifdef FIXEDDMP
      open(50,file='ini/fixeddm.dat',status='old')      
      read(50,*) SI_fiddm(0),SI_fiddm(1)
      if(myrank.eq.0) then
        write(6,*) ' Fixed DM particles id between ',SI_fiddm(0),SI_fiddm(1)
      endif
      if(SI_fiddm(0).gt.SI_fiddm(1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1 for DM.'
        endif
        call MPI_FINALIZE()
        stop
      endif
! *** change nadm
      nadm=0
      nc=0
      do i=0,ndm-1
        if(id_dm(i).ge.SI_fiddm(0).and.id_dm(i).le.SI_fiddm(1)) then
          nalist(nc)=i
          nc=nc+1
! *** set dvx and etc.
          dvx_dm(i)=0.0d0
          dvy_dm(i)=0.0d0
          dvz_dm(i)=0.0d0
        else
          list(nadm)=i
          nadm=nadm+1
        endif
      enddo
! *** change list_ap putting fixedp at the end
      do i=0,nadm-1
        list_adm(i)=list(i)
      enddo
      do i=0,nc-1
        list_adm(nadm+i)=nalist(i)
      enddo
#endif

#ifdef RANAZ
      open(50,file='ini/ranaz.dat',status='old')      
      read(50,*) SI_ranazi(0),SI_ranazi(1),SI_ranazi(2)
      read(50,*) SI_ranazt
      if(myrank.eq.0) then
        if(SI_ranazi(0).eq.0) then
          write(6,*) ' RANAZ baryon particles id between ' &
           ,SI_ranazi(1),SI_ranazi(2)
       else
          write(6,*) ' RANAZ DM particles id between ' &
           ,SI_ranazi(1),SI_ranazi(2)
        endif
        write(6,*) ' until t (Gyr)=',SI_ranazt
      endif
      SI_ranazt=SI_ranazt/TMUGYR
      if(SI_ranazi(1).gt.SI_ranazi(2)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1 RANAZ.'
        endif
        call MPI_FINALIZE()
        stop
      endif
#endif

#ifdef ARTFDM
      open(50,file='ini/artfdm.dat',status='old')      
      read(50,*) SI_artfdm(0),SI_artfdm(1)
      if(myrank.eq.0) then
        write(6,*) ' ARTF DM particles id between ' &
         ,SI_artfdm(0),SI_artfdm(1)
      endif
      if(SI_artfdm(0).gt.SI_artfdm(1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1 for DM.'
        endif
        call MPI_FINALIZE()
        stop
      endif
#endif

! *** output initial parameters ***
      if(myrank.eq.0) then
! ****** Initialzation for output file *****
        call system("rm ./output/ana/*")
        call system("rm ./output/data/*")           
        call system("rm ./output/*anim*")
        open(60,file='./output/ana/inipara.dat',status='unknown')
        write(60,'(a34,8I10)') 'npt,ndmt,np,ng,ns,ndm,ndm1,ndm1t=' &
         ,npt,ndmt,np,ng,ns,ndm,SI_ndm1,SI_ndm1t
#ifdef COSM
        write(60,'(a19,4(1pE13.5))') 'h,Omg0,lamda0,omgb=' &
         ,SI_h0,SI_omg0,SI_lam0,Ini_OMGg
#endif
!        write(60,'(a5,1pE13.5)') 'OFdt=',OF_odt
! since pv30.4
        write(60,'(a5,1pE25.17)') 'OFdt=',OF_odt
#ifdef COSM
!                                      1234567890123456789012345678901234
        write(60,'(a14,2(1pE25.17))') ' SI_eps0,epsa=' &
         ,SI_eps0,SI_epsa
#else
        write(60,'(a8,2(1pE25.17))') ' SI_eps=',SI_eps
#endif
        write(60,'(a6,2(1pE13.5))') 'zi,ze=',1.0d0/SI_a-1.0d0,Ini_ze
        write(60,'(a6,2(1pE12.5))') 'ti,te=',TM_tot,TM_t0
        write(60,'(a4,1pE13.5)') 'GAM=',GAM
#ifdef GAS
        write(60,'(a3,(1pE13.5))') 'mb=',m_p(0)
#endif
#ifdef COSM
        do i=0,Ini_mul
          write(60,'(i2,a12,I10,1pE13.5)') i,' level n,dx=',ndt(i),dxi(i)
        enddo
        nc = 0
#ifdef DM
        do il=Ini_mul,0,-1
          nc = nc+nd(il)
          write(60,'(a4,i2,a2,(1pE13.5))') 'mdm(l',il,')=',m_dm(nc-1)
        enddo
#endif
#endif
        write(60,'(a21,2(1pE13.5))') ' MGROOT,THETA,ALPHAT=',MGROOT,THETA,ALPHAT
        write(60,'(a10,I10)') ' NATTABLE=',NATTABLE
        write(60,'(a24,4(1pE13.5))')  &
!          12345678901234567890123456789012345
         ' CCFL,CDYN,CGRAV,FITLIM=',CCFL,CDYN,CGRAV,FITLIM
        write(60,'(a35,7(1pE13.5))')  &
!          12345678901234567890123456789012345
          ' AV ALPHAINI,MAX,MIN,BETA,SHARE,U,C' &
          ,V_ALPHAINI,V_ALPHAMAX,V_ALPHAMIN,V_BETA,SHAREV,U_ALPHAMAX,C_AV
        write(60,'(a18,3(1pE13.5))')  &
         ' ETAH,ETAHDM,ERRH=',ETAH,ETAHDM,ERRH
        write(60,'(a8,3(1pE13.5))')  &
         ' DHFLIM=',DHFLIM
        write(60,'(a38,5(1pE13.5))')  &
!          123456789012345678901234567890123456789012345678901
         ' MGTL,MGTU,CCOOLDT,CCOOLEDT,CCOOLRIDT=' &
         ,MGTL,MGTU,CCOOLDT,CCOOLEDT,CCOOLRIDT
        write(60,'(a36,5(1pE13.5))')  &
!         123456789012345678901234567890123456789012345678901
         ' TLLIMIT,THLIMIT,CLIMIT,EPSC,DU_UPU=' &
         ,TLLIMIT,THLIMIT,CLIMIT,EPSC,DU_UPU
        write(60,'(a22,2(1pE13.5))') &
!         12345678901234567890123456789012345
         ' DTSF,MINTMSFD,DLDUFD=',DTSF,MINTMSFD,LDUFD 
        write(60,'(a30,9(1pE13.5))') &
         ' NSTH,NSTHLIM,CSEFF,SNEU,SWEU=',NSTH,NSTHLIM,CSEFF,SNEU,SWEU
!         12345678901234567890123456789012345
#ifdef BOXSIM
        write(60,'(a10,3(1pE13.5))') ' box size=',SI_lbox(0) &
         ,SI_lbox(1),SI_lbox(2)
#ifdef COSM
!                                     12345678901234567890123456789012345
        write(60,'(a22,3(1pE13.5))') ' box size (comoving) =',SI_lbox0(0) &
         ,SI_lbox0(1),SI_lbox0(2)
#endif
#endif
#ifdef TREEPM
        write(60,'(a15,3I10)') ' Mesh nx,ny,nz=',nx_m,ny_m,nz_m
#endif
#ifdef BCX_FIX
        write(60,'(a9,I10)') ' BCX_ndx=',BCX_ndx
#ifdef MHD
        write(60,'(a25,5(1pE13.5))') 'Left BCX_rho,p,vx,vy,vz=' &
          ,BCX_rhox(0),BCX_px(0),BCX_vx(0),BCX_vy(0),BCX_vz(0)
        write(60,'(a13,3(1pE13.5))')'BCX_bx,by,bz=' &
          ,BCX_bx(0),BCX_by(0),BCX_bz(0)
        write(60,'(a26,5(1pE13.5))') 'Right BCX_rho,p,vx,vy,vz=' &
            ,BCX_rhox(1),BCX_px(1),BCX_vx(1),BCX_vy(1),BCX_vz(1)
        write(60,'(a13,3(1pE13.5))')'BCX_bx,by,bz=' &
          ,BCX_bx(1),BCX_by(1),BCX_bz(1)
#else
        write(60,'(a18,3(1pE13.5))') 'Left BCX_rho,p,vx=',BCX_rhox(0) &
         ,BCX_px(0),BCX_vx(0)
        write(60,'(a19,3(1pE13.5))') 'Right BCX_rho,p,vx=' &
         ,BCX_rhox(1),BCX_px(1),BCX_vx(1)
#endif
#endif
#ifdef MAXR
        write(60,'(a9,1pE13.5)') ' MAXRLIM=',MAXRLIM
#endif
#ifdef METALDIFF
        write(60,'(a9,1pE13.5)') ' ZDFCFAC=',ZDFCFAC
#endif
#ifdef UTHLIM_JM
        write(60,'(a9,(1pE13.5))') ' ULJMFAC=',ULJMFAC
#endif
#ifdef COSM
!                           12345678901234567890123456789012345
        write(60,'(a33)') '###  cosmological simulation  ###'
#endif
#ifdef MULTI_LEV
        write(60,'(a33)') '###  multi level  ###            '
#endif
#ifdef GAS
        write(60,'(a33)') '###  including GAS  ###          '
#endif
#ifdef DM
        write(60,'(a33)') '###  including DM  ###           '
#endif
#ifdef SGRAV
        write(60,'(a33)') '###  include self-gravity  ###   '
#endif
#ifdef COOL
        write(60,'(a33)') '###  include cooling  ###        '
#ifdef NOUVB
        write(60,'(a33)') '###  no UVB heating  ###         '
#endif
#endif
#ifdef STAR
        write(60,'(a22)') '###  influde star ### '
#endif
#ifdef SF
        write(60,'(a33)') '###  influde star formation  ### '
#endif
#ifdef SF_EFD
        write(60,'(a33)') '   ###  include energy feedback  ###'
#ifdef DTFD
        write(60,'(a18)') '   ### DTFD on ###'
#endif
#endif
#ifdef SF_ZFD
        write(60,'(a33)') '   ###  include metal feedback  ###'
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
        write(60,'(a34,I8,I8)') ' star mass group: SI_nsp, SI_snii=' &
          ,SI_nsp,SI_snii
#endif
#ifdef COREDW
        write(60,'(a33)') '   ###  core dw/ds for P grad   ###'
#endif
#ifdef GEPSLIM_G
        write(60,'(a31,1pE13.5)') &
         '###  gas softening limit on FH=',FH
#endif
#ifdef GEPSLIM_VH
        write(60,'(a36)') '   ###  gas no lim for smoothing ###'
#endif
#ifdef UTHLIM_JM
        write(60,'(a33)') '   ###  UTHLIM (JM) on  ###'
#endif
#ifdef METALDIFF
        write(60,'(a33)') '   ###  METALDIFF on ###'
#endif
#ifdef MAXR
        write(60,'(a33)') '   ###  MAXR on ###'
#endif
#ifdef MEXT
        write(60,'(a33)') '   ###  MEXT on ###'
#endif
#ifdef FIXEDDMEPS
!                          12345678901234567890123456789012345
        write(60,'(a30)') '   ###  Fixed DM softening ###'
#endif
#ifdef RANAZ
!                          1234567890123456789012345678901234567890
        write(60,'(a38)') '   ###  move particles azimuthally ###'
        write(60,'(a10,3I10,1pE13.5)') '   i012,t=',SI_ranazi(0) &
         ,SI_ranazi(1),SI_ranazi(2),SI_ranazt*TMUGYR
#endif
#ifdef QPTREEF
        write(60,'(a24)') '   ###  QP TREE F on ###'
#endif
#ifdef FASTG
        write(60,'(a20)') '   ###  FASTG on ###'
#endif
#ifdef DTDADT
!                          123456789012345678901234567890123456
        write(60,'(a21)') '   ###  DTDADT on ###'
#endif
#ifdef MINDT_SFD
!                          123456789012345678901234567890123456
        write(60,'(a24)') '   ###  MINDT_SFD on ###'
#endif
#ifdef NOAC
        write(60,'(a18)') '   ###  AC off ###'
#endif
#ifdef NOAVSWITCH
        write(60,'(a25)') '   ###  AV switch off ###'
#endif
#ifdef RHOINI
        write(60,'(a27)') '   ###  use initial RHO ###'
#endif
#ifdef PINI
        write(60,'(a25)') '   ###  use initial P ###'
#endif
#ifdef GLASS
!                          123456789012345678901234567890123456
        write(60,'(a24,F7.4,a4)') '   ###  GLASS on FGLASS=' &
        ,FGLASS,' ###'
#endif
#ifdef FIXEDP
        write(60,'(a28)') '   ###  fixed particlese ###'
        write(60,'(a12,2I10)') '   ID range=',SI_fid(0),SI_fid(1)
#endif
#ifdef FIXEDDMP
!                          123456789012345678901234567890123456
        write(60,'(a31)') '   ###  fixed DM particlese ###'
        write(60,'(a12,2I10)') '   ID range=',SI_fiddm(0),SI_fiddm(1)
#endif
#ifdef ARTFDM
!                          123456789012345678901234567890123456
        write(60,'(a30)') '   ###  artF DM particlese ###'
        write(60,'(a12,2I10)') '   ID range=',SI_artfdm(0),SI_artfdm(1)
#endif
#ifdef NOLOWMSF
!                          123456789012345678901234567890123456
        write(60,'(a30)') '   ###  NOWMSF on ###'
        write(60,'(a30)') '   ###  Do not creat star ID>SI_snii ###'
#endif

        close(60)
        write(6,*) ' end writing inipara.dat'
      endif

!#ifdef GAS
!      write(fileo,'(a4,i3.3)') 'inip',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ng-1
!        pn=list_ap(i) 
!        write(60,'(10(1pE13.5),2I10)') xc_p(pn),yc_p(pn),zc_p(pn) &
!         ,vnx_p(pn),vny_p(pn),vnz_p(pn),m_p(pn),rho_p(pn) &
!#ifdef METAL
!         ,h_p(pn),mzZ_p(pn),id_p(pn),flagc_p(pn)
!#else
!         ,h_p(pn),u_p(pn),id_p(pn),flagc_p(pn)
!#endif
!      enddo
!      close(60)
!#ifdef STAR
!      write(fileo,'(a4,i3.3)') 'inis',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ns-1
!        pn=list_ap(i+ng) 
!        write(60,'(10(1pE13.5),2I10)') xc_p(pn),yc_p(pn),zc_p(pn) &
!         ,vnx_p(pn),vny_p(pn),vnz_p(pn),m_p(pn),rho_p(pn) &
!#ifdef METAL
!         ,h_p(pn),mzZ_p(pn),id_p(pn),flagfd_p(pn)
!#else
!         ,h_p(pn),ts_p(pn),id_p(pn),flagfd_p(pn)
!#endif
!      enddo
!      close(60)
!#endif
!#endif
!#ifdef DM
!      write(fileo,'(a4,i3.3)') 'inid',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ndm-1
!        pn=list_adm(i)
!        write(60,'(8(1pE13.5),I10)') xc_dm(pn),yc_dm(pn),zc_dm(pn) &
!         ,vnx_dm(pn),vny_dm(pn),vnz_dm(pn),m_dm(pn),h_dm(pn),id_dm(pn)
!      enddo
!      close(60)
!#endif

      deallocate(istatus)
      deallocate(ireqs)
      deallocate(ireqr)
      deallocate(jstap)
      deallocate(jendp)

end subroutine
