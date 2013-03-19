#include "gcdp.def"
! ***************************************
!    stime.F for gcd+ ver.f03.0
!  27 Jan. 2013   written by D.KAWATA
! *************************************** 
! **********************************************
!    Definition of Function for Time Interval
! **********************************************

! *****   Definition of systime()   *****
!   This function calculate individual time step.
!     ->  dt_p[]
!   And Calculate minimum time step. = TM->dt
!   return k   (dtmax = dtmin x 2~k)
      
subroutine systime(np,ng,ndm,ns,step,flagdtg)
      use gcdp_const
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif
#ifdef DM
      use gcdp_dm
#endif
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np,ng,ndm,ns,step,flagdtg
      integer i,pn,ierr
      double precision zn,tn,dtmin,dtmaxeor,tdvs
! * Max time step *      
      double precision dtmax,dtu,dtc,dtzd
! * minimum limit due to output step
      double precision dtolim
! *** for work
      character fileo*60
      integer knp
      double precision dtcr1,dtcr2,adv
#if defined(GAS) || defined(STAR)
      double precision,allocatable :: dti(:)
#endif
#ifdef DM
      double precision,allocatable :: dtdmi(:)
#endif
#ifdef DTDADT
      double precision,allocatable :: av
#endif
#ifdef COSM
#ifdef COOL
      integer,allocatable :: listb(:)
      double precision,allocatable :: ram(:)
#endif
#endif

! **** Initialization ***
      dtmin = dble(INF)
      dtmax = 0.0d0
      dtolim=DTLIMIT*10.0d0
      dtu=INF
      dtc=INF
#if defined(GAS) || defined(STAR)
      allocate(dti(0:np))
#endif
#ifdef GAS

!      write(fileo,'(a4,i3.3)') 'stmg',myrank
!      open(60,file=fileo,status='unknown')

! *** for SPH Particles ***
      do i=0,ng-1
        pn = list_ap(i)
#ifdef FIXEDP
        if(id_p(pn).lt.SI_fid(0).or.id_p(pn).gt.SI_fid(1)) then
#endif
! *** Courant-Friedrich-Levy time criterion ***
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).eq.0.and.hvsigdt_p(pn).gt.0.0d0) then
#else
        if(hvsigdt_p(pn).gt.0.0d0) then
#endif
          dtcr1=CCFL*hvsigdt_p(pn)
        else
          dtcr1=INF
        endif          
! *** The force should not change too much ***
        adv=(dvx_p(pn)+ax_p(pn))**2+(dvy_p(pn)+ay_p(pn))**2 &
         +(dvz_p(pn)+az_p(pn))**2
        adv=dsqrt(adv)
        if(adv.ne.0.0d0) then
          dtcr2=CDYN*dsqrt(h_p(pn)*0.5d0/adv)
        else
          dtcr2=INF
        endif            
        dti(pn)=dtcr1
        if(dtcr1.gt.dtcr2) then
          dti(pn)=dtcr2
        endif
! *** time constraints for metal mixing 
! #ifdef METALDIFF
!        if(vsig_p(pn).gt.0.0d0) then
!          dtzd=CCFL*0.50d0*h_p(pn)/vsig_p(pn)
!          if(dtzd.lt.dti(pn)) then
!            dti(pn)=dtzd
!          endif
!        endif
! #endif
! *** time constraints for feedback gas
#ifdef SF_EFD
#ifdef DTFD
        if(flagfd_p(pn).ne.0.and.Gsn_p(pn).gt.0.0d0) then
          dtu=LDUFD*u_p(pn)/(Gsn_p(pn)/TMsf_dt)
          if(dtu.lt.dti(pn)) then
            dti(pn)=dtu
          endif
        endif
#endif
#endif
! *** check
        if(dti(pn).le.DTLIMIT*2.0d0) then
          write(6,*) ' Error in systime() : pn,id',pn,id_p(pn)
          write(6,*) ' gas lt_p,dt_p,TM_lt',lt_p(pn),dt_p(pn),TM_lt
          write(6,*) '  dt1,dt2',dtcr1,dtcr2,dti(pn)
          write(6,*) '  dvx,h,hvsig',dvx_p(pn),dvy_p(pn),dvz_p(pn),h_p(pn) &
           ,hvsigdt_p(pn) 
          write(6,*) '  ax',ax_p(pn),ay_p(pn),az_p(pn)
          write(6,*) ' rho,nnb=',rho_p(pn),nnb_p(pn)
          write(6,*) ' u,du=',u_p(pn),ndu_p(pn)
#ifdef COOL
          write(6,*) ' cooling time is=',dabs(u_p(pn)/ram_p(pn))
#endif
#ifdef SF_EFD
          write(6,*) ' Gsn,Tsf_dt,dtu,flagfd=',Gsn_p(pn),TMsf_dt,dtu,flagfd_p(pn)
#endif
! #ifdef METALDIFF
!          write(6,*) ' dt from metal diff,vsig=',dtzd,vsig_p(pn)
! #endif
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
#ifdef FIXEDP
        else
          dti(pn)=INF
       endif
#endif

!       write(60,'(16(1pE13.5),2I6)') x_p(pn),y_p(pn),dtcr1,dtcr2,dti(pn) &
!        ,dtu,u_p(pn),Gsn_p(pn) &
!        ,h_p(pn),hvsigdt_p(pn),adv,dvx_p(pn),dvy_p(pn),dvz_p(pn) &
!        ,ax_p(pn),ay_p(pn),id_p(pn),nnb_p(pn)

      enddo

!      close(60)
!      if(step.eq.1) then 
!        stop
!      endif

#ifdef COSM
#ifdef COOL
      if(SI_flagrad.le.-1.and.SI_flagrad.gt.-NITEOR) then
        allocate(listb(0:ng))
! ram can be over 0:np-1
        allocate(ram(0:np))
! initialise
        do i=0,np-1
          ram(i)=0.0d0
        enddo
        do i=0,ng-1
          listb(i)=list_ap(i)
        enddo
        zn = 1.0d0/SI_a-1.0d0
        call cool(ng,listb,zn,np,ram)
        deallocate(listb)
        do i=0,ng-1
          pn = list_ap(i)
! *** check cooling time ***
          if(dabs(ndu_p(pn)-ram(pn)).gt.0.0d0) then
            if(dti(pn).gt.CCOOLRIDT*dabs(u_p(pn)/(ndu_p(pn)-ram(pn)))) then
              dti(pn)=CCOOLRIDT*dabs(u_p(pn)/(ndu_p(pn)-ram(pn)))
            endif
          endif
        enddo
        deallocate(ram)
      else
        do i=0,ng-1
          pn = list_ap(i)
! *** check cooling time ***
#ifdef SF_EFD
          if(flagfd_p(pn).eq.0.and.dabs(ndu_p(pn)-ram_p(pn)).gt.0.0d0) then
#else
          if(dabs(ndu_p(pn)-ram_p(pn)).gt.0.0d0) then
#endif
            dtc=CCOOLDT*dabs(u_p(pn)/(ndu_p(pn)-ram_p(pn)))
            if(dti(pn).gt.dtc) then
              dti(pn)=dtc
            endif
          endif
        enddo
      endif
#endif
#endif
#ifdef STAR
! *** for Star Particles ***
      if(ns.gt.0) then
        do i=ng,np-1
          pn = list_ap(i)           
! *** The force should not change too much ***
          adv=dvx_p(pn)*dvx_p(pn)+dvy_p(pn)*dvy_p(pn)+dvz_p(pn)*dvz_p(pn)
          adv=dsqrt(adv)
          if(adv.ne.0.0d0) then
            dti(pn)=CGRAV*dsqrt(0.5d0*h_p(pn)/adv)
          else
            dti(pn) = INF
          endif 
! *** higher order ***
#ifdef DTDADT
          if(flagdtg.ne.0) then
! *** The force should not change too much ***
            av=(pax_p(pn)-ax_p(pn))**2 &
             +(pay_p(pn)-ay_p(pn))**2+(paz_p(pn)-az_p(pn))**2
            av=dsqrt(av(pn))/dt_p(pn)
            if(av.ne.0.0d0) then
              dtcr2=CGRAV*(adv/av)
            else
              dtcr2=INF
            endif 
            if(dti(pn).gt.dtcr2) then
              dti(pn)=dtcr2
            endif
          endif
#endif

! *** check delta t ***
          if(dti(pn).le.DTLIMIT*2.0d0) then
            write(6,*) ' Error in systime() : pn,id',pn,id_p(pn)
            write(6,*) ' star lt_p,dt_p,TM_lt',lt_p(pn),dt_p(pn),TM_lt
            write(6,*) '  dt=',dtcr2,dti(i)
            write(6,*) '  dvx,h',dvx_p(pn),dvy_p(pn),dvz_p(pn),h_p(pn)
            write(6,*) '  ax',ax_p(pn),ay_p(pn),az_p(pn)
            write(6,*) ' rho,nnb=',rho_p(pn),nnb_p(pn)
            write(6,*) ' u,du=',u_p(pn),ndu_p(pn)
#ifdef SF_EFD
            write(6,*) ' Gsn,Tsf_dt,dtu,flagfd=' &
             ,Gsn_p(pn),TMsf_dt,dtu,flagfd_p(pn)
#endif
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
            stop
          endif
        enddo
      endif
#endif

! *** find minimum and max time step ***
      do i=0,np-1
        if(dti(i).lt.dtmin) then
          dtmin=dti(i)
! for check
          pn = i
!
        endif
      enddo

      do i=0,ng-1
        pn=list_ap(i)           
        if(dti(pn).ne.INF) then
          if(dti(pn).gt.dtmax) then
            dtmax=dti(pn)
          endif
        endif
      enddo
      if(ng.le.0) then
        do i=ng,ng+ns-1
          pn=list_ap(i)           
          if(dti(pn).ne.INF) then
            if(dti(pn).gt.dtmax) then
              dtmax=dti(pn)
            endif
          endif
        enddo
      endif
! endif for GAS
#endif

#ifdef DM
! *** for DM Particles ***
      if(ndm.gt.0) then
        allocate(dtdmi(0:ndm))
! * Softning length *	           
        do i=0,ndm-1
#ifdef FIXEDDMP
          if(id_dm(i).lt.SI_fiddm(0).or.id_dm(i).gt.SI_fiddm(1)) then
#endif

! *** The force should not change too much ***
          adv=dvx_dm(i)*dvx_dm(i)+dvy_dm(i)*dvy_dm(i)+dvz_dm(i)*dvz_dm(i)
          adv=dsqrt(adv)
          if(adv.ne.0.0d0) then
            dtdmi(i)=CGRAV*dsqrt(h_dm(i)*0.5d0/adv)
          else
            dtdmi(i)=INF
          endif
#ifdef DTDADT
! *** higher order ***
          if(flagdtg.ne.0) then
! *** The force should not change too much ***
            av=(pdvx_dm(i)-dvx_dm(i))**2 &
             +(pdvy_dm(i)-dvy_dm(i))**2+(pdvz_dm(i)-dvz_dm(i))**2
            av=dsqrt(av)/dt_dm(i)
            if(av.ne.0.0d0) then
              dtcr2(i)=CGRAV*(adv/av)
            else
              dtcr2(i) = INF
            endif 
            if(dtdmi(i).gt.dtcr2(i)) then
              dtdmi(i)=dtcr2(i)
            endif
          endif
#endif
#ifdef FIXEDDMP
          else
            dtdmi(i)=INF
          endif
#endif
        enddo
        do i=0,ndm-1
          if(dtdmi(i).lt.dtmin) then
            dtmin = dtdmi(i)
            pn=i
          endif          
        enddo
      endif
      if(np.eq.0) then
        do i=0,ndm-1
          if(dtdmi(i).gt.dtmax.and.dtdmi(i).ne.INF) then
           dtmax = dtdmi(i)
          endif 
        enddo
      endif
#endif
! *** get dtmax and dtmin for the whole node ***
      if(nprocs.gt.1) then
        tdvs=dtmin
        call MPI_ALLREDUCE(tdvs,dtmin,1,MPI_DOUBLE_PRECISION &
         ,MPI_MIN,MPI_COMM_WORLD,ierr)
        tdvs=dtmax      
        call MPI_ALLREDUCE(tdvs,dtmax,1,MPI_DOUBLE_PRECISION &
         ,MPI_MAX,MPI_COMM_WORLD,ierr)
      endif

#if defined(SF) || defined(SF_EFD) || defined(SF_ZFD)
! *** minimum time step for SF and feedback
      if(dtmin.gt.MINTMSFD) then
        dtmin=MINTMSFD
      endif
#endif

! *** define time steps ***
      TM_k=int(dlog(dtmax/dtmin)/dlog(2.0d0))
      if(step.le.5) then
        TM_k=0
! #ifdef MULTI_LEV
      else if(TM_k.gt.8) then
        TM_k = 8 
      endif
! #else
! *** individual time step off
#ifdef NOIDT
      TM_k=0       
#endif
! *** limit by 4 seems better: since pv.29.5 (6/11 2008) ***
!      else if(TM_k.gt.4) then
!        TM_k = 4 
!      endif
! #endif
#ifdef COSM
#ifdef COOL
      if(SI_flagrad.le.-1.and.SI_flagrad.gt.-NITEOR) then
        SI_flagrad=SI_flagrad-1
        TM_k=0
      endif
#endif
#endif
! for test
!      TM_k=0
! *** define largest time bin ***
      dtmax=2.0d0**(dble(TM_k))*dtmin

! *** Check EoR ***
#ifdef COSM
#ifdef COOL
      dtmaxeor=0.0d0
      if(SI_flagrad.eq.0) then
        call afromt(TM_tot+dtmax,zn)  
        zn = 1.0d0/zn-1.0d0
        if(zn.lt.SI_zeor) then
          SI_flagrad = -1
          zn = 1.0d0/(SI_zeor+1.0d0)
          call tfroma(zn,tn)
          dtmax = tn-TM_tot
          if(dtmin.lt.dtmax) then
            TM_k = int(dlog(dtmax/dtmin)/dlog(2.0d0))+1
            dtmin = dtmax/(2.0d0**(dble(TM_k)))         
          else
            TM_k = 0
            dtmin = dtmax
          endif 
          dtmaxeor = dtmax
        endif
      endif
#endif
#endif

! ***  Check Output time ***
      if(dtmin+TM_tot.gt.OF_todt) then
        dtmin = OF_todt-TM_tot
        dtmax = dtmin
        TM_k=0
        SI_flagout = 1
      else if(dtmax+TM_tot+dtolim.gt.OF_todt) then
        dtmax = OF_todt-TM_tot
        TM_k = int(dlog(dtmax/dtmin)/dlog(2.0d0))+1
        dtmin = dtmax/(2.0d0**(dble(TM_k)))
        SI_flagout = 1
      endif
! ***  Check end time ***
      if(dtmin+TM_tot.gt.TM_t0) then
        dtmin=TM_t0-TM_tot
        dtmax = dtmin
        TM_k = 0
        SI_flagout = 1        
      else if(dtmax+TM_tot.gt.TM_t0) then
        dtmax=TM_t0-TM_tot
        TM_k=int(dlog(dtmax/dtmin)/dlog(2.0d0))+1
        dtmin=dtmax/(2.0d0**(dble(TM_k)))
        SI_flagout = 1        
      endif
#ifdef COOL
      if(SI_flagrad.eq.-1.and.dtmax.lt.dtmaxeor) then
        SI_flagrad = 0
      endif
#endif

! *** Update local Time Interval ***
      TM_dt = dtmin
      TM_ldtmin = dtmin
      TM_sdt = dtmax

#if defined(SF) || defined(SF_EFD) || defined(SF_ZFD)
! *** set TM_mindtsfd ***
      i=int(dlog(MINTMSFD/dtmin)/dlog(2.0d0))
      TM_mindtsfd=dtmin*(2.0d0**dble(i))
   90 if(TM_mindtsfd.gt.MINTMSFD) then
        TM_mindtsfd=TM_mindtsfd*0.5d0
      endif
      if(TM_mindtsfd.gt.MINTMSFD) then 
        i=i-1
        if(myrank.eq.0) then
          write(6,*) ' TM_mindtsfd devided by 2',i
        endif
        goto 90
      endif
#endif

#ifdef GAS
! *** Define Individual time step ***
! *** Gas and Star Particles ***
      do i=0,np-1
#ifdef FIXEDP
        if(id_p(i).lt.SI_fid(0).or.id_p(i).gt.SI_fid(1)) then
#endif
        if(dti(i).lt.dtmax) then
          knp=int(dlog(dti(i)/dtmin)/dlog(2.0d0))
          dt_p(i)=dtmin*2.0d0**dble(knp)
        else
          dt_p(i)=dtmax
        endif
#ifdef FIXEDP
        else
          dt_p(i)=INF
          vdt_p(i)=0.0d0
        endif
#endif
      enddo
! *** Set Individual Time ***
      do i=0,np-1
        lt_p(i) = 0.0d0
        ltg_p(i) = 0.0d0
        flagt_p(i)=0
        dtr_p(i)=0.0d0
        dtmnb_p(i)=INF
      enddo
      deallocate(dti)
#endif

#ifdef DM
! *** DM Particles ***
      do i=0,ndm-1
#ifdef FIXEDDMP
        if(id_dm(i).lt.SI_fiddm(0).or.id_dm(i).gt.SI_fiddm(1)) then
#endif
        if(dtdmi(i).lt.dtmax) then
          knp=int(dlog(dtdmi(i)/dtmin)/dlog(2.0d0))
          dt_dm(i)=dtmin*2.0d0**dble(knp)
        else
          dt_dm(i)=dtmax
        endif          
#ifdef FIXEDDMP
        else
          dt_dm(i)=INF
          vdt_dm(i)=0.0d0
        endif
#endif
      enddo
! *** Set Individual Time ***
      do i=0,ndm-1
        lt_dm(i) = 0.0d0
        flagt_dm(i)=0
        dtr_dm(i)=0.0d0
        dtmnb_dm(i)=INF
      enddo        
      deallocate(dtdmi)
#endif

end subroutine
