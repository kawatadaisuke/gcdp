#include "gcdp.def"
! *************************************
!    main.f for GCD+ ver.f03.0
!  24  Jan. 2013    written by D.KAWATA
! *************************************
!70*******************************************************************

program gcdp
      use gcdp_const
      use gcdp_system
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#ifdef METAL
      use gcdp_metal
#endif
#endif
#ifdef DM
      use gcdp_dm
#endif
#ifdef FFTW3_MPI  
      use fftw3
#endif

      implicit none
      include 'mpif.h'
! * Total Number of Particles, Active *   
      integer npt,np,pn,nst,ngt
! * Total Number of Gas Particles, Active *       
      integer ng,nag,nagravg
! * Total Number of DM Particles, Active *        
      integer ndmt,ndm,nadm
! * Number of Star Particles, Active, formed stars*	  
      integer ns,nas
      integer nagt,nagravgt,nast,nadmt
! * Number of system steps *	  
      integer step
! * Local step in system step *	  
      integer lstep
! * Total Time Steps *	  
      integer totstep
! * Flag for Kin.eq , Ene.eq, dt from gravity, decomposition *
      integer flagv,flagu,flagdtg
      integer i
! * flag for output *
      integer flagend
      integer flagcont
      integer flagt
! * for continue *	
      integer tstep
      integer ierr
! *** for test ***
      character fileo*60
      double precision rhonth
!      double precision cs10,cs100
! * for time *
      double precision ectime,dtiddec
      double precision tout,tddec,tcout
      double precision itime,ftime,ntime,pstime
      double precision JOBTIME,MARGIN
! for work
      integer tivs(0:3),tivr(0:3)
! *****   parameter values   *****
! 1 week
!      parameter (JOBTIME=604800.0d0)
! 60 hours
      parameter (JOBTIME=216000.0d0)
! 36 hours
!      parameter (JOBTIME=129600.0d0)
! 12 hours
!      parameter (JOBTIME=43200.0d0)
! 8 hours
!      parameter (JOBTIME=28800.0d0)
! 6 hours
!      parameter (JOBTIME=21600.0d0)

      Parameter (MARGIN=1800.0d0)
      integer OUTSTEP
      parameter (OUTSTEP=10000)

! ***** MPI Initialization *****
      call MPI_INIT(ierr)
#ifdef FFTW3_MPI
      call fftw_mpi_init()
#endif
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)


! *** print version ***
      if(myrank.eq.0) then
        write(6,*) ' GCD+ ver. f03.0 15/01/14'
        print *, "Process ",myrank, " of ", nprocs, " is alive"
      endif
      if(nprocs.gt.NCPU) then
        if(myrank.eq.0) then 
          write(6,*) ' Error: increase NCPU=',NCPU
        endif
        call MPI_FINALIZE()
        stop
      endif

! ***** For time *****      
      itime=MPI_WTIME()
      ftime = JOBTIME-MARGIN
      pstime=itime
! **********   Initialization   **********
      flagt = 0
      flagend = 0
      flagdtg=0

! *** setting time ***
      WTM_ini=0.0d0
      WTM_out=0.0d0
      WTM_ddecb=0.0d0
      WTM_ddecdm=0.0d0
      WTM_cont=0.0d0
      WTM_setval=0.0d0
      WTM_calcvu=0.0d0
      WTM_tree=0.0d0
      WTM_starf=0.0d0
      WTM_starfd=0.0d0
      WTM_setvaldm=0.0d0
      WTM_time=0.0d0
      WTM_setvals=0.0d0
      WTM_setvalfd=0.0d0
      WTM_feedb=0.0d0
      WTM_bf=0.0d0
      WTM_dmf=0.0d0
      WTM_update=0.0d0
      WTM_updateu=0.0d0
#ifdef TREEPM
      WTM_pm=0.0d0
#endif

! *************************
!  Set Initial Condition
! *************************
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()      
#endif
      call set_ini(npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
       ,step,totstep,lstep,flagcont)
! *** get total number ofactive particles ***
      tivs(0)=nag
      tivs(1)=nadm
      tivs(2)=nas
      tivs(3)=nagravg
      call MPI_ALLREDUCE(tivs,tivr,4,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      nagt=tivr(0)
      nadmt=tivr(1)
      nast=tivr(2)
      nagravgt=tivr(3)

      if(flagcont.ne.0) then
! *** use pax_p read in set_ini ***
        flagdtg=1
      endif

! initial set up
! setting kernel for SPH and softening
      call setkernel()
#ifdef TREEPM
#ifdef FFTW3
      call fftw3_init()
#endif
! initial set up for TreePM
      call mesh_gen(0)
#endif
#ifdef MEXT
! set up the external potential
      call setmext()
#endif

#ifdef GAS
! estimate of dynamical time in the softening limit
      if(myrank.eq.0) then
        pn=0
#ifdef METAL
        rhonth=(NSTH/(DU/MP))/((m_p(pn)-((mzZ_p(pn)+mzHe_p(pn)) &
         /MUSM))/m_p(pn))
#else
        rhonth=(NSTH/(DU/MP))/((m_p(pn)/MUSM)/m_p(pn))
#endif
        write(6,*) ' rho_th, tdyn=',rhonth &
         ,dsqrt(3.0d0*M_PI/(16.0d0*G*rhonth))*TMUGYR

!        cs10=dsqrt(GAM*TPRHO*(MYU/1.25d0)*0.001d0)
!        cs100=dsqrt(GAM*TPRHO*(MYU/1.25d0)*0.01d0)
!        write(6,*) ' cs at 10 and 100 K=',cs10*VUKMS,cs100*VUKMS
!        write(6,*) ' mjeans for 100 and 10 K (MYU=1.25)='
!     &   ,(M_PI**(5.0d0/2.0d0))*(cs100**3)
!     &   /(6.0d0*(G**(3.0d0/2.0d0))*dsqrt(rhonth))
!     &   ,(M_PI**(5.0d0/2.0d0))*(cs10**3)
!     &   /(6.0d0*(G**(3.0d0/2.0d0))*dsqrt(rhonth))

      endif
#endif

      tstep = totstep+int(OUTSTEP)
      if(npt.gt.0) then
#if defined(GAS) || defined(STAR)
#if defined(SF_EFD) || defined(SF_ZFD)
        call setyield()
#endif
#endif
#if defined(COOL) || defined(BBAT)
        call setcool()
        if(1.0d0/SI_a-1.0d0.lt.SI_zeor) then
          SI_flagrad=-NITEOR
        endif
#endif
      endif

      if(flagcont.ne.0) then
        tivs(0)=nag
        tivs(1)=nadm
        tivs(2)=nas
        tivs(3)=nagravg
        tivr(0)=0
        tivr(1)=0
        tivr(2)=0
        tivr(3)=0
        call MPI_ALLREDUCE(tivs,tivr,4,MPI_INTEGER &
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(myrank.eq.0) then         
          write(6,*) ' k=',TM_k,' kn=',TM_kn,' systime=',TM_sdt &
           ,' initime=',TM_ldtmin
          write(6,'(a15,4i8)') &
!            12345678901234567890
           ' tnag,nadm,nas,nagravg=',tivr(0),tivr(1),tivr(2),tivr(3)
        endif          
      else
        TM_kn = 0
        TM_k = 0
      endif
      nst=0
      ngt=0
#ifdef GAS
      tivs(0)=ng
      tivs(1)=ns
      call MPI_ALLREDUCE(tivs,tivr,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      ngt=tivr(0)
      nst=tivr(1)
#endif
      if(myrank.eq.0) then      
        write(6,*) '**********************************************'
        write(6,*) ' Total Number of GAS&STAR and DM =',npt,ndmt
        write(6,*) ' Total Number of GAS and STAR    =',ngt,nst
        write(6,*) ' number of particles in rank 0 np=',np
        write(6,*) '   GAS  =',ng,' active = ',nag
        write(6,*) '   STAR =',ns,' active = ',nas
        write(6,*) '   DM   =',ndm,' active = ',nadm
        write(6,*) '   DM (high resolution)  =',SI_ndm1
        write(6,*) '  Specific heat      =',GAM
#ifdef COSM
         write(6,*) '   ###  cosmological simulation  ###'
#endif
#ifdef MULTI_LEV
         write(6,*) '   ###  multi level  ###'
#endif
#ifdef GAS
         write(6,*) '   ###  including GAS  ###'
#endif
#ifdef DM
         write(6,*) '   ###  including DM  ###'
#endif
#ifdef MHD
         write(6,*) '   ###  including MHD  ###'
#endif
#ifdef TIMEDIS
         write(6,*) '   ###  include time dependent resistivity  ###'
#endif
#ifdef BBAT
         write(6,*) '   ###  include Biermann battery source  ###'
#endif
#ifdef SGRAV
         write(6,*) '   ###  include self-gravity  ###'
#endif
#ifdef COOL
        write(6,*) ' ###  include cooling  ###'
#endif
#ifdef STAR
        write(6,*) ' ###  include star ###'
#endif
#ifdef SF
        write(6,*) ' ###  include star formation  ###'
#endif
#ifdef SF_EFD
        write(6,*) '   ###  include energy feedback  ###'
#endif
#ifdef SF_ZFD
        write(6,*) '   ###  include metal feedback  ###'
#endif
#ifdef BCX_FIX
        write(6,*) ' Fix boundary on X-axis'
#endif
#ifdef BCX_PERIOD
        write(6,*) ' Periodic boundary on X-axis'
#endif
#ifdef BCY_PERIOD
        write(6,*) ' Periodic boundary on Y-axis'
#endif
#ifdef BCZ_PERIOD
        write(6,*) ' Periodic boundary on Z-axis'
#endif
#ifdef SIM1D
        write(6,*) ' 1 D simulation'
#elif defined(SIM2D)
        write(6,*) ' 2 D simulation'
#endif
        write(6,*) '**********************************************'
      endif

! *** check for self-gravity in 1D or 2D ***
#if defined(SIM1D) || defined(SIM2D)
#ifdef SGRAV
        write(6,*) ' Error: self-gravity cannot include in 1D or 2D sim.'
        call MPI_FINALIZE()
        stop
#endif
#endif

! *** for time keeping (moved here since pv28.4)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_ini=WTM_ini+ectime-ntime
#endif

! **********************
!   Domain Decomposition
!  **********************
!      if(nprocs.gt.1.and.flagcont.eq.0) then
      if(flagcont.eq.0) then
! *** for baryon ***
        if(step.eq.0.and.flagcont.eq.0) then
#ifdef RHOINI
! *** keep the initial density information ****
          if(1.eq.2) then
#endif
#ifdef GAS
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ntime=MPI_WTIME()
#endif
          call ddecb(npt,np,ng,ns,nag,nas,nagravg,-1,flagcont)
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ectime=MPI_WTIME()
          WTM_ddecb=WTM_ddecb+ectime-ntime
#endif
#endif
#ifdef DM
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ntime=MPI_WTIME()
#endif
          call ddecdm(ndmt,ndm,nadm,-1)
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ectime=MPI_WTIME()
          WTM_ddecdm=WTM_ddecdm+ectime-ntime
#endif
#endif
#ifdef RHOINI
          endif
#endif
        else
#ifdef GAS
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ntime=MPI_WTIME()
#endif
          call ddecb(npt,np,ng,ns,nag,nas,nagravg,totstep,flagcont)
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ectime=MPI_WTIME()
          WTM_ddecb=WTM_ddecb+ectime-ntime
#endif
#endif
#ifdef DM
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ntime=MPI_WTIME()
#endif
          call ddecdm(ndmt,ndm,nadm,totstep)
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ectime=MPI_WTIME()
          WTM_ddecdm=WTM_ddecdm+ectime-ntime
#endif
#endif
        endif
#ifdef WTIME
        if(myrank.eq.0) then
          write(6,*) ' initial ddec time for baryon,DM=',WTM_ddecb,WTM_ddecdm
        endif
        dtiddec=WTM_ddecb+WTM_ddecdm
#endif
      endif

! ***********************************
!   Gas Tree Building for set_value
! ***********************************
#ifdef GAS
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
! *** for star it is needed also in continue mode ***
      call gtreebuild(ng,ns)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_tree=WTM_tree+ectime-ntime
#endif
#endif
      if(flagcont.le.2) then
#ifdef SGRAV
#ifdef DM
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
        call dmtreebuild(SI_ndm1)
#ifdef MULTI_LEV
        call ldmtreebuild(ndm)
#endif
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_tree=WTM_tree+ectime-ntime
#endif
#endif
#endif
     
! ****************************************
!   Set Value of rho, |div(v)|, |rot(v)|
! ****************************************
! * 1 <- Not Calculate Density *		
#ifdef GAS
! do not use flag=1, because ddecb does not send rho_p
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
#ifdef RHOINI
        call set_value(ng,nag,ns,1)
#else
        call set_value(ng,nag,ns,2)
#endif
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_setval=WTM_setval+ectime-ntime
#endif

! *** for gas properties ***
        if(flagcont.eq.0) then
          do i=0,ng-1
            pn=list_ap(i)
#ifdef PINI
            u_p(pn)=p_p(pn)/((GAM-1.0d0)*rho_p(pn))
#else
            p_p(pn)=(GAM-1.0d0)*rho_p(pn)*u_p(pn)
#endif
            cs_p(pn)=dsqrt(GAM*p_p(pn)/rho_p(pn))
            as_p(pn)=p_p(pn)/(rho_p(pn)**GAM)
          enddo
        endif
#endif

! **************************************************
!  set values for DM ( Initial ) flag=1 use vnx_dm
! **************************************************
#ifdef DM
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
        call set_valuedm(ndm,nadm,1)
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_setvaldm=WTM_setvaldm+ectime-ntime
#endif
#endif

! **********************************************
!   Tree Building & Compute Mass after set h_p
! **********************************************
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
#ifdef GAS
! *** need to be done after h is changed ***
        call gtreebuild(ng,ns)
#endif
#ifdef SGRAV
#if defined(GAS) || defined(STAR)
        call treebuild(np)
#endif
#ifdef DM
#ifndef FIXEDDMEPS
        call dmtreebuild(SI_ndm1)
#ifdef MULTI_LEV
        call ldmtreebuild(ndm)
#endif
#endif
#endif
#endif
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_tree=WTM_tree+ectime-ntime
#endif

! *** initialise flagrfd_p ***
#if defined(GAS) || defined(STAR)
#if defined(SF_EFD) || defined(SF_ZFD)
        do i=0,ng+ns-1        
          flagrfd_p(list_ap(i))=0
        enddo
#endif
#endif


! *************************************
!  Calculate dv/dt & du/dt ( Initial )
! *************************************
        flagv = 0
        flagu = 0

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        if(myrank.eq.0) then
          write(6,*) ' begin tree force',ectime-ntime
        endif

        call calc_dv_du(np,ng,ndm,nag,nadm,nas,nagravg &
         ,nagt,nadmt,nast,nagravgt,flagv,flagu,1)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        if(myrank.eq.0) then
          write(6,*) ' finish tree force',ectime-ntime
        endif

#ifdef TREEPM
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
! TreePM long-range force
! ***** adding PM force *****
! assign density
        call mesh_setrho(np,ndm)
! calculate potential
        call mesh_fftf()
! assign force to particles
!        call mesh_setfp(ng,nag,nas,nadm)

#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_pm=WTM_pm+ectime-ntime
      if(myrank.eq.0) then
        write(6,*) ' mesh force end',ectime-ntime
      endif
#endif


#ifdef DM
! test output
      write(fileo,'(a4,i3.3)') 'fdmp',myrank
      open(60,file=fileo,status='unknown')
      do i=0,nadm-1
        if(x_dm(i)**2+y_dm(i)**2+z_dm(i)**2.gt.0.0d0) then
        write(60,'(13(1pE13.5),2I8)') x_dm(i),y_dm(i),z_dm(i) &
         ,dvx_dm(i),dvy_dm(i),dvz_dm(i) &
         ,dsqrt(x_dm(i)**2+y_dm(i)**2+z_dm(i)**2) &
         ,-(dvx_dm(i)*x_dm(i)+dvy_dm(i)*y_dm(i)+dvz_dm(i)*z_dm(i)) &
         /dsqrt(x_dm(i)**2+y_dm(i)**2+z_dm(i)**2) &
         ,omgh_dm(i),zetah_dm(i),dvz_dm(i),rho_dm(i),h_dm(i),nnb_dm(i),i
        else
         write(60,'(13(1pE13.5),2I8)') x_dm(i),y_dm(i),z_dm(i) &
        ,dvx_dm(i),dvy_dm(i),dvz_dm(i) &
         ,dsqrt(x_dm(i)**2+y_dm(i)**2+z_dm(i)**2) &
         ,-(dvx_dm(i)*x_dm(i)+dvy_dm(i)*y_dm(i)+dvz_dm(i)*z_dm(i)) &
         ,omgh_dm(i),zetah_dm(i),dvz_dm(i),rho_dm(i),h_dm(i),nnb_dm(i),i
        endif
      enddo
      close(60)
#endif
#if defined(GAS) || defined(STAR)
      write(fileo,'(a3,i3.3)') 'fbp',myrank
      open(60,file=fileo,status='unknown')
      do i=0,nag-1
        if(x_p(i)**2+y_p(i)**2+z_p(i)**2.gt.0.0d0) then
        write(60,'(13(1pE13.5),2I8)') x_p(i),y_p(i),z_p(i) &
         ,ax_p(i),ay_p(i),az_p(i) &
         ,dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) &
         ,-(ax_p(i)*x_p(i)+ay_p(i)*y_p(i)+az_p(i)*z_p(i)) &
         /dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) &
         ,omgh_p(i),zetah_p(i),az_p(i),rho_p(i),h_p(i),nnb_p(i),i
        else
         write(60,'(13(1pE13.5),2I8)') x_p(i),y_p(i),z_p(i) &
        ,ax_p(i),ay_p(i),az_p(i) &
         ,dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) &
         ,-(ax_p(i)*x_p(i)+ay_p(i)*y_p(i)+az_p(i)*z_p(i)) &
         ,omgh_p(i),zetah_p(i),az_p(i),rho_p(i),h_p(i),nnb_p(i),i
        endif
      enddo
      close(60)
#endif

      goto 9999

#endif

      endif
! *****   >>>  Start Iteration   <<<   *****
   77 if(flagcont.eq.0) then	  
! **********************************
!  Analize various System Value
! **********************************
        if(TM_tot.ge.TM_t0-0.0000001d0) then
          SI_flagout = 1
          flagend=1
        endif
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
        call output(npt,ndmt,np,ng,ndm,ns,step,totstep)
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_out=WTM_out+ectime-ntime
#endif        
!****************************
!   Check condition of end
! ****************************
        if(flagend.ne.0) then
          goto 999
        endif
      else
        flagcont = 0
      endif

      
! **********************************
!  Calculate Individual time step
! **********************************
      if(TM_kn.eq.TM_k) then
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
        call systime(np,ng,ndm,ns,step,flagdtg)
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_time=WTM_time+ectime-ntime
#endif
#ifdef STAR
        TMsf_t=TM_tot
        TMsf_dt=0.0d0
#endif         
        if(myrank.eq.0) then        
          write(6,*) step,',',totstep,'step, k=',TM_k, &
           ' sdt=',TM_sdt,' idt=',TM_dt,', Ttot=',TM_tot &
#ifdef COOL
           ,' SI_a=',SI_a,'flagrad=',SI_flagrad
#else
           ,' SI_a=',SI_a 
#endif
        endif          
! * Number of local steps *
        lstep = 0
! * Local time in system time *		
        TM_lt = 0.0
      endif

! **********   Iteration in system time step   **********      
  777 if(flagt.ne.0) then  
        if(myrank.eq.0) then
          write(6,*) ' This JOB is Stopped for TIME LIMIT.' 
        endif
        goto 9999
      endif


! ***********************
!  Calculate local time
!  ***********************
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
      call lctime(np,ng,ndm,ns,nag,nadm,nas,lstep,nagravg,flagdtg)
      if(flagdtg.ne.0) then
! *** use pax_p read in set_ini ***
        flagdtg=1
      endif


#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_time=WTM_time+ectime-ntime
#endif      

#ifdef RANAZ
! ***************************************
!  move particles azimuthally at randum
! ***************************************
      call ranazimuth(ng,ndm,ns,nag,nadm,nas,nagravg)
#endif

! *******************************
!   Star formation and feedback
! *******************************
#ifdef STAR
! *** star foramtion every time step ***
#ifdef SF
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
      call starf(ng,ns,nag,nas,nagravg)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_starf=WTM_starf+ectime-ntime
#endif
#endif

! *** feedback energy every minimum time step ***
      TMsf_dt=TM_dt
      TMsf_t=TM_tot
#if defined(SF_EFD) || defined(SF_ZFD)

#ifdef SF_EFD
      do i=0,ng-1
        Gsn_p(list_ap(i)) = 0.0d0
      enddo
#endif
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
      call starfd(ng,ns,nag,nas,nagravg)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_starfd=WTM_starfd+ectime-ntime
#endif
      TMsf_t=TM_tot+TM_dt
#endif
#endif

! *** this is necessary ***
      tivs(0)=nag
      tivs(1)=nadm
      tivs(2)=nas
      tivs(3)=nagravg
      call MPI_ALLREDUCE(tivs,tivr,4,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      nagt=tivr(0)
      nadmt=tivr(1)
      nast=tivr(2)
      nagravgt=tivr(3)

      if(0.eq.0) then
        tivs(0)=ng
        tivs(1)=ns
        call MPI_ALLREDUCE(tivs,tivr,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
        if(myrank.eq.0) then
!          write(6,*) 'after starf and starfd'
          write(6,'(a41,5i10,4(1pE11.3))') &
!            12345678901234567890123456789012345678901
           ' tstep,tnag,nadm,nas,nagravg,dt,lt,sdt,t=' &
           ,totstep,nagt,nadmt,nast,nagravgt &
           ,TM_dt,TM_lt,TM_sdt,TM_tot
          write(6,'(a19,2i10,1pE11.3)') &
!            123456789012345678901234
           ' total ng,ns,wtime=',tivr(0),tivr(1),ntime-pstime
        endif
        pstime=ntime
      endif

! ************************************************************
!  Update Artificial Viscosity and Thermal Conductivity
!  parameters, Alpha_v and Alpha_u
!  using the values at step n
! ************************************************************
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
#ifdef GAS
      call updatealpvu(nag)
#endif
      
! ****************************
!   Update Smoothing Length
! ****************************
      call updatesl(ng,nag,ns,nas,nadm)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_update=WTM_update+ectime-ntime
#endif

#ifdef WTIME
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            ntime=MPI_WTIME()
#endif
! ***************************
!  Update Velocity V(n+1/2)
!  ***************************
      call updatev(np,ndm)

! ***************************************************************
!  Update Scale Factor (SI_a) at TM_tot+TM_dt
! ***************************************************************
#ifdef COSM
      call afromt(TM_tot+TM_dt,SI_a)
      if(SI_a.lt.SI_epsa) then
        SI_eps=(1.0d0-SI_epsa+SI_a)*SI_eps0
      else
        SI_eps=SI_eps0
      endif
#ifdef BOXSIM
! a^n+1 is also used in updatev_pos
! update SI_lbox at n+1 step
      SI_lbox(0)=SI_a*SI_lbox0(0)
      SI_lbox(1)=SI_a*SI_lbox0(1)
      SI_lbox(2)=SI_a*SI_lbox0(2)
! get Hubble parameter at n+1 step
      call hubfroma(SI_a,SI_hub)
#endif
#endif

! **************************
!  Update Virtual Position 
! **************************
      call updatev_pos(np,ndm)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_update=WTM_update+ectime-ntime
#endif


! **********************
!   Domain Decomposition
!  **********************
!      if(nprocs.gt.1.and.totstep.ne.0) then
      if(totstep.ne.0) then
#ifdef COOL
!        if(SI_flagrad.ge.-1.or.SI_flagrad.eq.-10) then
#endif
          if(mod(totstep,NDSTEPDDEC).eq.0) then
#ifdef GAS
! *** for baryon ***
#ifdef WTIME
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            ntime=MPI_WTIME()
#endif

            call ddecb(npt,np,ng,ns,nag,nas,nagravg,totstep,flagcont)

#ifdef WTIME
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            ectime=MPI_WTIME()
            WTM_ddecb=WTM_ddecb+ectime-ntime
#endif
#endif
#ifdef DM
! *** for DM ***
#ifdef WTIME
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            ntime=MPI_WTIME()
#endif
            call ddecdm(ndmt,ndm,nadm,totstep)
#ifdef WTIME
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            ectime=MPI_WTIME()
            WTM_ddecdm=WTM_ddecdm+ectime-ntime
#endif
#endif
          endif


! ***********************************************
!  Update Virtual Position after ddecb and ddecdm
! ***********************************************
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ntime=MPI_WTIME()
#endif
! ***************************
!  Update Velocity V(n+1/2)
!  ***************************
          call updatev(np,ndm)
! **************************
!  Update Virtual Position 
! **************************
          call updatev_pos(np,ndm)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_update=WTM_update+ectime-ntime
#endif
      endif

! ************************
!  Update Time and Step
! ************************
      lstep=lstep+1
      totstep=totstep+1
      TM_tot=TM_tot+TM_dt
      TM_lt=TM_lt+TM_dt

! *********************************
!  Update Virtual Velocity V(n+1)
!  *********************************
#ifdef GAS
      call updatev_vn(ng+ns)
#endif
! ********************
!   Update Position 
! ********************
#ifdef BCX_FIX
      call boundary_xfix(np,npt)
#endif
      call updatepos(ng,nag,nadm,nas)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      ectime=MPI_WTIME()
      WTM_update=WTM_update+ectime-ntime
#endif
! *********************************
!   Gas and DM Tree Building & Compute Mass
! *********************************
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
#ifdef GAS
      if(nagt.gt.0) then
        call gtreebuild(ng,ns)
      endif
#endif
#ifdef SGRAV
#ifdef DM
      if(nagravgt.gt.0.or.nadmt.gt.0) then
        call dmtreebuild(SI_ndm1)
#ifdef MULTI_LEV
        call ldmtreebuild(ndm)
#endif
      endif
#endif
#endif
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_tree=WTM_tree+ectime-ntime
#endif

! *************************************************
!  Set New Value of rhon+1,|div(v)|n+1,|rho(v)|n+1
! *************************************************
#ifdef GAS
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif

!      write(6,*) ' set_value start',myrank

      call set_value(ng,nag,ns,2)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_setval=WTM_setval+ectime-ntime
#endif
#endif

! ********************
!  set values for DM 
! ********************
#ifdef DM
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
      call set_valuedm(ndm,nadm,0)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_setvaldm=WTM_setvaldm+ectime-ntime
#endif
#endif

! **********************************************
!   Tree Building & Compute Mass after set h_p
! **********************************************
#ifdef GAS
! *** need to be done after h is changed ***
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif

!      write(6,*) ' gtree2 start',myrank

      call gtreebuild(ng,ns)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_tree=WTM_tree+ectime-ntime
#endif
#endif
#ifdef SGRAV
#if defined(GAS) || defined (STAR)

!      write(6,*) 'myrank,nagravgt=',myrank,nagravgt

!      if(nagravgt.gt.0.or.nadmt.gt.0) then
        call treebuild(np)
!      endif
#endif
#ifdef DM
#ifndef FIXEDDMEPS
      if(nadmt.gt.0.or.nagravgt.gt.0) then
        call dmtreebuild(SI_ndm1)
#ifdef MULTI_LEV
        call ldmtreebuild(ndm)
#endif
      endif
#endif
#endif
#endif

#ifdef GAS
! **********************************
!  update pressure & sound velocity
! **********************************
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif 

!      write(6,*) ' vpcs start',myrank

      call virtualp_cs(np,nag)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_updateu=WTM_updateu+ectime-ntime
#endif

#ifdef MHD
! ************************
!  Update B, predictor
! ************************
      call updatevb(ng)
#endif

! ************************
!  Calculate (du/dt)n+1 
! ************************
      flagv =  1
      flagu = 2

!      write(6,*) ' calcdvdu 1 start',myrank
!      call system('rm *calvu*')

      call calc_dv_du(np,ng,ndm,nag,nadm,nas,nagravg,nagt,nadmt,nast &
       ,nagravgt,flagv,flagu,0)
! *************************************
!  Update Internal Energy  First!
! ************************************
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif

!      write(6,*) ' updateu start',myrank

      call updateu(np,nag)
! **********************************
!  Update Pressure & Suound Velocity
! **********************************
      call updatep_cs(nag)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_updateu=WTM_updateu+ectime-ntime
#endif
#endif      
#if defined(GAS) || defined(STAR)
#if defined(SF_EFD) || defined(SF_ZFD)
! *** initialise flagrfd_p ***
        do i=0,ng+ns-1        
          flagrfd_p(list_ap(i))=0
        enddo
#endif
#endif

!***************************
! update B (with corrector)
!***************************
#ifdef MHD
! update B for active particles
      call updatebn(nag)
#endif

! ************************
!  Calculate (du/dt)n+1 
! ************************
      flagv =  0
      flagu = 0

!      write(6,*) ' calcdvdu 2 start',myrank

      call calc_dv_du(np,ng,ndm,nag,nadm,nas,nagravg,nagt,nadmt &
       ,nast,nagravgt,flagv,flagu,0)

!      write(fileo,'(a6,i3.3)') 'calvu2',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ng-1
!        pn=list_ap(i)
!        write(60,'(21(1pE13.5),3I10)') x_p(pn),y_p(pn),z_p(pn) &
!        ,ndu_p(pn),dvx_p(pn),dvy_p(pn),dvz_p(pn),div_v_p(pn),arot_v_p(pn) &
!        ,d2u_p(pn),p_p(pn),u_p(pn),ax_p(pn),ay_p(pn),az_p(pn) &
!        ,rho_p(pn),omgh_p(pn),hvsigdt_p(pn),alpv_p(pn),ram_p(pn) &
!        ,myu_p(pn),id_p(pn),nnb_p(pn),flagfd_p(pn)
!      enddo
!      close(60)
!      stop


#ifdef TREEPM
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
! ***** adding PM force *****
! generate new mesh
        call mesh_gen(1)
! assign density
        call mesh_setrho(np,ndm)
! calculate potential
        call mesh_fftf()
! assign force to particles
        call mesh_setfp(ng,nag,nas,nadm)
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_pm=WTM_pm+ectime-ntime
#endif
#endif

!      write(fileo,'(a4,i3.3,i4.4)') 'test',myrank,totstep
!      open(60,file=fileo,status='unknown')
!      do i=0,ndm-1
!        write(60,'(12(1pE13.5))') x_dm(i),y_dm(i),z_dm(i) &
!         ,dvx_dm(i),dvy_dm(i),dvz_dm(i) &
!         ,vnx_dm(i),vnx_dm(i),vnx_dm(i) &
!         ,xc_dm(i),yc_dm(i),zc_dm(i)
!      enddo

! *************************
!  Update Velocity V(n+1)
! **************************
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ntime=MPI_WTIME()
#endif
#ifdef BCX_FIX
      call boundary_xfix(np,npt)
#endif
      call updatevn(ng,nag,nadm,nas,nagravg)
#ifdef METALDIFF
! ******************************
!  Metal diffusion
! if continue from cont*.dat
! ****************************** 
      call metal_diff(nag)
#endif
#ifdef WTIME
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ectime=MPI_WTIME()
      WTM_update=WTM_update+ectime-ntime
#endif

! ***************************
!   Make file for continuing
! ***************************
      ntime=MPI_WTIME()
      if(ntime-itime.gt.ftime) then
        flagt=1
      endif
      call MPI_BCAST(flagt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(totstep.eq.tstep.or.flagt.eq.1) then
        tstep = tstep+OUTSTEP
!        flagt=1        
        if(myrank.eq.0) then
          write(6,*) ' making continue files as totstep=',totstep
          write(6,*) ' step,k,kn=',step,TM_k,TM_kn
        endif
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
        if(TM_kn.ne.TM_k)	then
          call forcont(npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
           ,step,totstep,lstep)
        else
          call forcont(npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
           ,step+1,totstep,lstep)
        endif
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_cont=WTM_cont+ectime-ntime
#endif
      endif

! *** set flagfd_p for feedback gas particle ***
! *** negative flagfd_p for gas (new feedback gas) is used in cont.F
#if defined(SF_EFD) || defined(SF_ZFD)      
      do i=0,nag-1
        pn=list_ap(i)
        if(flagfd_p(pn).gt.0) then
          flagfd_p(pn)=-flagfd_p(pn)
        endif
      enddo
#endif

      if(TM_kn.ge.TM_k) then
        goto 99
      endif

      goto 777
! ***** end local time step *****
   99 step=step+1
      goto 77
! ****** end system time step *****
  999 if(myrank.eq.0) then
         write(6,*) ' *** The simulation ended  ***'
         write(6,*) '  t,t_end=',TM_tot,TM_t0
      endif
      call forcont(npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
       ,step,totstep,lstep)
 9999 ntime=MPI_WTIME()
      if(myrank.eq.0) then
        write(6,*) ' Calculation time summary'
        write(6,*) ' Total=',ntime-itime
#ifdef WTIME
        write(6,*) ' 1:input=',WTM_ini
        write(6,*) ' 2:output=',WTM_out
        write(6,*) ' 3:output for continue files=',WTM_cont
        write(6,*) ' 4:initial domain decomposition=',dtiddec
        write(6,*) ' domain decomposition b, DM=',WTM_ddecb,WTM_ddecdm
        write(6,*) ' total-(1-4 time)=',ntime-itime-WTM_ini-WTM_out &
         -WTM_cont-dtiddec
        write(6,*) ' set_value()=',WTM_setval
        write(6,*) ' calc_dv_du()=',WTM_calcvu
        write(6,*) ' treebuild() =',WTM_tree
        write(6,*) ' treef() for baryon=',WTM_bf
        write(6,*) ' treef() for DM=',WTM_dmf
        write(6,*) ' set_valuedm()=',WTM_setvaldm
        write(6,*) ' set_values()=',WTM_setvals
        write(6,*) ' update()=',WTM_update
        write(6,*) ' updateu()=',WTM_updateu
        write(6,*) ' sys,lctime()=',WTM_time
        write(6,*) ' star formation=',WTM_starf
        write(6,*) ' feedback=',WTM_starfd
        write(6,*) '  setvalfd=',WTM_setvalfd
        write(6,*) '  feedb=',WTM_feedb
#ifdef TREEPM
        write(6,*) ' PM force=',WTM_pm
#endif
#endif
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

#ifdef FFTW3
#ifdef FFTW3_MPI
      call fftw_mpi_cleanup()
#endif
#endif

      call MPI_FINALIZE(ierr)

end program gcdp
