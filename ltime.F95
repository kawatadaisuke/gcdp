#include "gcdp.def"
! ***************************************
!    ltime.c for gcd+ ver. f03.0
!   27 Jan. 2013    written by D.KAWATA
! *************************************** 

! **********************************************
!    Definition of Function for Time Interval
! **********************************************

! *************************************************
!    Definition lctime()
! Finding active particle.    -> list_ap[].
!  Except local step = 0, finding minimum time step
! and update individual time step of active particle
! **************************************************
subroutine lctime(np,ng,ndm,ns,nag,nadm,nas,lstep,nagravg,flagdtg)
      use gcdp_const
#if defined(GAS) || defined (STAR)
      use gcdp_baryon
#endif
#ifdef DM
      use gcdp_dm
#endif
      use gcdp_system

      implicit none
      include 'mpif.h'
      
      integer,intent(in) :: np,ng,ndm,ns,lstep,flagdtg
      integer,intent(inout) :: nag,nadm,nas,nagravg
      integer i,j,pn,ierr
! * Number of Non-active particles, active, active also gravity *	  
      integer nanp,na
      double precision ttmp,dtu,dtc
      integer kp
      double precision upt,lwt
      double precision tlim,dtp
! *** to check minimum time step from cooling ***
      integer pnt
      double precision dtmin,dtcmin,dtumin,dt1min,dt2min
! for work
      integer knp
      double precision adv,av,dtcr1,dtcr2,tdvs
#ifdef DTDADT
      double precision av
#endif
      integer,allocatable :: list(:),nalist(:),talist(:)
      double precision,allocatable :: dti(:)

      dtu=INF
      dtc=INF

! *** get dt from gravity only for gas ***
#ifdef GAS
      if(lstep.eq.0) then

        allocate(dti(0:ng))

        do i=0,ng-1
          pn=list_ap(i)
#ifdef FIXEDP
          if(id_p(pn).lt.SI_fid(0).or.id_p(pn).gt.SI_fid(1)) then
#endif
! *** The gravitational force should not change too much ***
          adv=ax_p(pn)**2+ay_p(pn)**2+az_p(pn)**2
          adv=dsqrt(adv)
          if(adv.ne.0.0d0) then
#if defined(GEPSLIM_VH) || defined(GEPSLIM_G)
            if(h_p(pn).lt.SI_eps*(m_p(pn)**THIRD)) then
              dti(i)=CDYN*dsqrt(0.5d0*SI_eps*(m_p(pn)**THIRD)/adv)
            else
#endif
              dti(i)=CDYN*dsqrt(h_p(pn)*0.5d0/adv)
#if defined(GEPSLIM_VH) || defined(GEPSLIM_G)
            endif
#endif
          else
            dti(i)=INF
          endif
#ifdef DTDADT
          if(flagdtg.ne.0) then
! *** d a/dt ***
            av=(pax_p(pn)-ax_p(pn))**2+(pay_p(pn)-ay_p(pn))**2 &
             +(paz_p(pn)-az_p(pn))**2
            av=dsqrt(av)/dt_p(pn) 
            if(av.ne.0.0d0) then
              dtcr2=CDYN*(adv/av)
            else
              dtcr2=INF
            endif
            if(dti(i).gt.dtcr2) then
              dti(i)=dtcr2
            endif
          endif
#endif
! *** block time step ***
          knp=int(dlog(dti(i)/dt_p(pn))/dlog(2.0d0))
          if(knp.gt.FASTGLIMIT) then
            knp=FASTGLIMIT
          else if(knp.lt.0) then
            knp=0
          endif
          dtg_p(pn)=dt_p(pn)*(2.0d0**dble(knp))
          if(dtg_p(pn)+ltg_p(pn).gt.TM_sdt) then
            dtg_p(pn)=TM_sdt-ltg_p(pn)
          endif
#ifdef FIXEDP
          else
            dtg_p(pn)=INF
          endif
#endif
        enddo

        deallocate(dti)

! *** set this for initial kick later ***
        nagravg=ng
      endif
#else 
      nag=0
      nagravg=0
#endif

! *** for active particle ***		 
      if(lstep.ne.0) then
! *** update Individual time ***
#ifdef GAS
! *** Gas Particles ***
        do i=0,nag-1
          pn=list_ap(i)
          lt_p(pn)=TM_lt
          if(flagt_p(pn).ne.0) then
            if(flagt_p(pn).eq.2) then
! *** also update dtg_p for feedback particles made in starfd.F
              dtg_p(pn)=dtr_p(pn)
            endif
! *** correction for dt ***
            dt_p(pn)=dtr_p(pn)
            flagt_p(pn)=0
            dtr_p(pn)=0.0d0
          endif
        enddo
#ifdef STAR
! *** Star Particles ***
        do i=ng,ng+nas-1
          pn=list_ap(i)
          lt_p(pn)=TM_lt
          if(flagt_p(pn).ne.0) then
! *** correction for dt ***
            dt_p(pn)=dtr_p(pn)
            flagt_p(pn)=0
            dtr_p(pn)=0.0d0
          endif
        enddo		  
#endif
#endif
#ifdef DM
! *** DM Particles ***
        do i=0,nadm-1
          pn=list_adm(i)
          lt_dm(pn)=TM_lt
          if(flagt_dm(pn).ne.0) then
! *** correction for dt ***
            dt_dm(pn)=dtr_dm(pn)
            flagt_dm(pn)=0
            dtr_dm(pn)=0.0d0
          endif
        enddo
#endif
! *** update Individual time step ***
! *** calculate time interval of active particle ***
        dtmin=INF
        pn=0      
#ifdef GAS
! *** for Gas Particles ***
        if(nag.gt.0) then

          allocate(dti(0:nag))

          do i=0,nag-1
            pn=list_ap(i)			 
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
#if defined(SF_EFD) || defined(SF_ZFD)
#ifdef GEPSLIM_VH
              if(flagfd_p(pn).ne.0.and.h_p(pn).lt.SI_eps*(m_p(pn)**THIRD)) then
                dtcr2=CDYN*dsqrt(0.5d0*SI_eps*(m_p(pn)**THIRD)/adv)
              else
#endif
#endif
                dtcr2=CDYN*dsqrt(h_p(pn)*0.5d0/adv)
#if defined(SF_EFD) || defined(SF_ZFD)
#ifdef GEPSLIM_VH
              endif
#endif
#endif
            else
              dtcr2=INF
            endif			  
            dti(i)=dtcr1
            if(dtcr1.gt.dtcr2) then
              dti(i)=dtcr2
            endif
! *** time constraints for feedback gas
            dtu=INF
#ifdef SF_EFD
#ifdef DTFD
            if(flagfd_p(pn).ne.0.and.Gsn_p(pn).gt.0.0d0) then
              dtu=LDUFD*u_p(pn)/(Gsn_p(pn)/TMsf_dt)
              if(dtu.lt.dti(i)) then
                dti(i)=dtu
              endif
            endif
#endif
#endif
! *** time constraints from cooling ***
            dtc=INF
#ifdef COOL
! *** check cooling time
#ifdef SF_EFD
            if(flagfd_p(pn).eq.0.and.dabs(ndu_p(pn)-ram_p(pn)).gt.0.0d0)then
#else
            if(dabs(ndu_p(pn)-ram_p(pn)).gt.0.0d0) then
#endif
              dtc=CCOOLDT *dabs(u_p(pn)/(ndu_p(pn)-ram_p(pn)))
              if(dtc.lt.dti(i)) then
                dti(i)=dtc
              endif
            endif
#endif

! *** check minimum dt
            if(dti(i).lt.dtmin) then
              pnt=pn
              dtmin=dti(i)
              dtcmin=dtc
              dtumin=dtu
              dt1min=dtcr1
              dt2min=dtcr2
            endif

! #ifdef COOL
! c *** if cooling time is short, keep the currnt delta t 01/04/09 ***
!            if(dabs(ram_p(pn)).gt.0.0d0) then
!              if(dti(i).gt.0.10d0*dabs(u_p(pn)/ram_p(pn))) then
!                dti(i)=0.1d0*dabs(u_p(pn)/ram_p(pn))
!              endif
!            endif
! #endif
! *** cooling time ***
#ifdef COSM
            if(SI_flagrad.le.-NITEOR.or.SI_flagrad.gt.-1) then
#endif
              if(dti(i).le.DTLIMIT*2.0d0) then
              write(6,*) ' Error in localtime() : pn,id',pn,id_p(pn)
              write(6,*) '  lt_p,dt_p,TM_lt',lt_p(pn),dt_p(pn),TM_lt
              write(6,*) '  dt1,dt2',dtcr1,dtcr2,dti(i)
              write(6,*) '  dvx,h,hvsig',dvx_p(pn),dvy_p(pn),dvz_p(pn) &
               ,h_p(pn),hvsigdt_p(pn)
              write(6,*) '  ax',ax_p(pn),ay_p(pn),az_p(pn)
              write(6,*) ' rho,nnb=',rho_p(pn),nnb_p(pn)
              write(6,*) ' u,du=',u_p(pn),ndu_p(pn)
#ifdef COOL
              write(6,*) ' cooling time is=',dabs(u_p(pn)/ram_p(pn))
#endif
#ifdef SF_EFD
              write(6,*) ' Gsn,Tsf_dt,dtu,flagfd=' &
               ,Gsn_p(pn),TMsf_dt,dtu,flagfd_p(pn)
#endif
              call MPI_ABORT(MPI_COMM_WORLD,ierr)
              stop
            endif
#ifdef COSM
            endif
#endif

          enddo			

!*** update individual time step to discreat time step ***
          do i=0,nag-1
            pn = list_ap(i)			 
            if(dti(i).lt.0.25d0*dt_p(pn)) then
              dt_p(pn)=0.125d0*dt_p(pn)
            else if(dti(i).lt.0.5d0*dt_p(pn)) then
              dt_p(pn)=0.25d0*dt_p(pn)
            else if(dti(i).lt.dt_p(pn)) then
              dt_p(pn) = 0.5d0*dt_p(pn)
            else if(dti(i).gt.8.0d0*dt_p(pn)) then
              if(dmod(lt_p(pn)+0.5d0*dt_p(pn),8.0d0*dt_p(pn)).lt.dt_p(pn)) then
#ifdef FASTG
              if(dmod(ltg_p(pn)+dtg_p(pn)-lt_p(pn)+0.5*dt_p(pn) &
               ,8.0d0*dt_p(pn)).lt.dt_p(pn)) then
#endif
                dt_p(pn)=8.0d0*dt_p(pn)
#ifdef FASTG
              endif
#endif
              endif				
            else if(dti(i).gt.4.0d0*dt_p(pn)) then
              if(dmod(lt_p(pn)+0.5d0*dt_p(pn),4.0d0*dt_p(pn)).lt.dt_p(pn)) then
#ifdef FASTG
              if(dmod(ltg_p(pn)+dtg_p(pn)-lt_p(pn)+0.5*dt_p(pn) &
               ,4.0d0*dt_p(pn)).lt.dt_p(pn)) then
#endif
                dt_p(pn)=4.0d0*dt_p(pn)
#ifdef FASTG
              endif
#endif
              endif				
            else if(dti(i).gt.2.0d0*dt_p(pn)) then
              if(dmod(lt_p(pn)+0.5d0*dt_p(pn),2.0d0*dt_p(pn)) &
               .lt.dt_p(pn)) then
#ifdef FASTG
              if(dmod(ltg_p(pn)+dtg_p(pn)-lt_p(pn)+0.5*dt_p(pn) &
               ,2.0d0*dt_p(pn)).lt.dt_p(pn)) then
#endif
! this also ensure synchronize with dtg_p 
                dt_p(pn)=2.0d0*dt_p(pn)
#ifdef FASTG
!               else
!                 write(6,*) ' Warning in ltime(), myrank,dt,dtg,lt,ltg=' &
!                 ,myrank,dt_p(pn),dtg_p(pn),lt_p(pn),ltg_p(pn)
              endif
#endif
              endif				
            endif				
            if(dt_p(pn)+lt_p(pn).gt.TM_sdt) then
              dt_p(pn)=TM_sdt-lt_p(pn)
            endif
          enddo

          do i=0,nagravg-1
            pn=list_ap(i)
            ltg_p(pn)=TM_lt
! *** The gravitational force should not change too much ***
            adv=ax_p(pn)**2+ay_p(pn)**2+az_p(pn)**2
            adv=dsqrt(adv)
            if(adv.ne.0.0d0) then
#ifdef GEPSLIM_VH
              if(h_p(pn).lt.SI_eps*(m_p(pn)**THIRD)) then
                dti(i)=CDYN*dsqrt(0.5d0*SI_eps*(m_p(pn)**THIRD)/adv)
              else
#endif
                dti(i)=CDYN*dsqrt(h_p(pn)*0.5d0/adv)
#ifdef GEPSLIM_VH
              endif
#endif
            else
              dti(i)=INF
            endif
#ifdef DTDADT
! *** d a/dt ***
            av=(pax_p(pn)-ax_p(pn))**2+(pay_p(pn)-ay_p(pn))**2 &
             +(paz_p(pn)-az_p(pn))**2
            av=dsqrt(av)/dt_p(pn)
            if(av.ne.0.0d0) then
              dtcr2=CDYN*(adv/av)
            else
              dtcr2=INF
            endif
            if(dti(i).gt.dtcr2) then
              dti(i)=dtcr2
            endif
#endif
! *** update dtg ***
            if(dti(i).lt.0.25d0*dtg_p(pn)) then
              dtg_p(pn)=0.125d0*dtg_p(pn)
            else if(dti(i).lt.0.5d0*dtg_p(pn)) then
              dtg_p(pn)=0.25d0*dtg_p(pn)
            else if(dti(i).lt.dtg_p(pn)) then
              dtg_p(pn) = 0.5d0*dtg_p(pn)
            else if(dti(i).gt.2.0d0*dtg_p(pn)) then
              if(dmod(ltg_p(pn)+0.5d0*dtg_p(pn),2.0d0*dtg_p(pn)) &
               .lt.dtg_p(pn)) then
                dtg_p(pn)=2.0d0*dtg_p(pn)
              endif				
            endif	
            if(dtg_p(pn)+ltg_p(pn).gt.TM_sdt) then
              dtg_p(pn)=TM_sdt-ltg_p(pn)
            endif
          enddo

          deallocate(dti)

        endif
#endif

#ifdef DM
! ***** for DM Particles *****
        if(nadm.gt.0) then

          allocate(dti(0:nadm))

          do i=0,nadm-1
            pn=list_adm(i)			 
! *** The force should not change too much ***
            adv=dvx_dm(pn)*dvx_dm(pn)+dvy_dm(pn)*dvy_dm(pn) &
             +dvz_dm(pn)*dvz_dm(pn)
            adv=dsqrt(adv)
            if(adv.ne.0.0d0) then
              dti(i)=CGRAV*dsqrt(0.5d0*h_dm(pn)/adv)
            else
              dti(i)=INF
            endif
#ifdef DTDADT
! *** higher order ***
            av=(pdvx_dm(pn)-dvx_dm(pn))**2+(pdvy_dm(pn)-dvy_dm(pn))**2 &
             +(pdvz_dm(pn)-dvz_dm(pn))**2
            av=dsqrt(av)/dt_dm(pn)
            if(av.gt.0.0d0) then
              dtcr2=CGRAV*(adv/av)
            else
              dtcr2=INF
            endif		  
            if(dti(i).gt.dtcr2) then
              dti(i)=dtcr2
            endif
#endif
          enddo
! *** update individual time step to discreat time step ***
          do i=0,nadm-1
            pn=list_adm(i)
            if(dti(i).lt.0.25d0*dt_dm(pn)) then
              dt_dm(pn)=0.125d0*dt_dm(pn)
            else if(dti(i).lt.0.5d0*dt_dm(pn)) then
              dt_dm(pn)=0.25d0*dt_dm(pn)
            else if(dti(i).lt.dt_dm(pn)) then
              dt_dm(pn)=0.5d0*dt_dm(pn)
            else if(dti(i).gt.2.0d0*dt_dm(pn)) then
              if(dmod(lt_dm(pn)+0.5d0*dt_dm(pn),2.0d0*dt_dm(pn)) &
               .lt.dt_dm(pn)) then
                dt_dm(pn)=2.0d0*dt_dm(pn)
              endif
            endif
            if(dt_dm(pn)+lt_dm(pn).gt.TM_sdt) then
              dt_dm(pn)=TM_sdt-lt_dm(pn)
            endif
          enddo

          deallocate(dti)

        endif
#endif
#ifdef STAR
! ***** for Star Particles *****
        if(nas.gt.0) then

          allocate(dti(ng:ng+nas))

          do i=ng,nas+ng-1
            pn = list_ap(i)
! *** The force should not change too much ***
            adv=dvx_p(pn)*dvx_p(pn)+dvy_p(pn)*dvy_p(pn)+dvz_p(pn)*dvz_p(pn)
            adv=dsqrt(adv)
            if(adv.ne.0.0d0) then
              dti(i)=CGRAV*dsqrt(0.5d0*h_p(pn)/adv)
            else
              dti(i)=INF
            endif			  
#ifdef DTDADT
! *** higher order ***
            av=(pax_p(pn)-ax_p(pn))**2+(pay_p(pn)-ay_p(pn))**2 &
             +(paz_p(pn)-az_p(pn))**2
            av=dsqrt(av)/dt_p(pn)
            if(av.ne.0.0d0) then
              dtcr2=CGRAV*(adv/av)
            else
              dtcr2=INF
            endif			  
            if(dti(i).gt.dtcr2) then
              dti(i)=dtcr2
            endif
#endif

! *** check delta t ***
            if(dti(i).le.DTLIMIT*2.0d0) then
              write(6,*) ' Error in localtime() : pn,id',pn,id_p(pn)
              write(6,*) '  lt_p,dt_p,TM_lt',lt_p(pn),dt_p(pn),TM_lt
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
! *** update individual time step to discreat time step ***
          do i=ng,nas+ng-1
            pn=list_ap(i) 
            if(dti(i).lt.0.25d0*dt_p(pn)) then
              dt_p(pn)=0.125d0*dt_p(pn)
            else if(dti(i).lt.0.5d0*dt_p(pn)) then
              dt_p(pn)=0.25d0*dt_p(pn)
            else if(dti(i).lt.dt_p(pn)) then
              dt_p(pn)=0.5d0*dt_p(pn)
            else if(dti(i).gt.2.0d0*dt_p(pn)) then
              if(dmod(lt_p(pn)+0.5d0*dt_p(pn),2.0d0*dt_p(pn)).lt.dt_p(pn)) then 
                dt_p(pn) = 2.0d0*dt_p(pn)
              endif
            endif
            if(dt_p(pn)+lt_p(pn).gt.TM_sdt) then
              dt_p(pn)=TM_sdt-lt_p(pn)
            endif			  
          enddo

          deallocate(dti)

        endif
#endif
        tlim=TM_lt
#ifdef GAS
! *** check for the time step limiter for baryon 
        do i=0,np-1
#ifdef FIXEDP
          if(id_p(i).lt.SI_fid(0).or.id_p(i).gt.SI_fid(1)) then         
#endif
          dtmnb_p(i)=dtmnb_p(i)*FITLIM
          if((flagt_p(i).eq.0.and.dt_p(i).gt.dtmnb_p(i)) &
           .or.(flagt_p(i).ne.0.and.dtmnb_p(i).lt.dtr_p(i))) then
            kp=int((TM_sdt-tlim)/dtmnb_p(i))-1
            dtp=(TM_sdt-dtmnb_p(i)*dble(kp))-lt_p(i)
            if(dtp.lt.dt_p(i)) then
              dt_p(i)=dtp
              dtr_p(i)=dtmnb_p(i) 
              flagt_p(i)=1
            endif
          endif
#ifdef FIXEDP
          endif
#endif
        enddo
#endif
#ifdef DM
        do i=0,ndm-1
#ifdef FIXEDDMP
          if(id_dm(i).lt.SI_fiddm(0).or.id_dm(i).gt.SI_fiddm(1)) then
#endif
          dtmnb_dm(i)=dtmnb_dm(i)*FITLIM
! *** 4 x dtj should be 1/2**n times ***
          if((flagt_dm(i).eq.0.and.dt_dm(i).gt.dtmnb_dm(i)) &
           .or.(flagt_dm(i).ne.0.and.dtmnb_dm(i).lt.dtr_dm(i))) then
            kp=int((TM_sdt-tlim)/dtmnb_dm(i))-1
            dtp=(TM_sdt-dtmnb_dm(i)*dble(kp))-lt_dm(i)
            if(dtp.lt.dt_dm(i)) then
              dt_dm(i)=dtp
              dtr_dm(i)=dtmnb_dm(i)
              flagt_dm(i)=1
            endif  
          endif
#ifdef FIXEDDMP
          endif
#endif
        enddo
#endif
! *** find minimum time step ***
! *** for Gas and Star Particle ***
        TM_dt=INF
#ifdef GAS   
        if(np.gt.0) then   
          TM_dt=lt_p(0)+dt_p(0)
#ifdef DM
        else
          TM_dt=lt_dm(0)+dt_dm(0)
#endif
        endif
#else
        TM_dt=lt_dm(0)+dt_dm(0)
#endif
#ifdef GAS
        do i=1,np-1
          ttmp=lt_p(i)+dt_p(i)
          if(ttmp.lt.TM_dt) then
            TM_dt = ttmp
          endif
        enddo

#endif
#ifdef DM
! *** for DM Particle ***
        do i=0,ndm-1
          ttmp=lt_dm(i)+dt_dm(i)
          if(ttmp.lt.TM_dt) then
            TM_dt=ttmp
#ifdef COSM
            if(SI_flagrad.le.-NITEOR.or.SI_flagrad.gt.-1) then
#endif
            if(TM_dt.le.TM_lt+DTLIMIT) then
              write(6,*) ' Error in localtime(): tlim=',tlim
              write(6,*) '  lt_dm,dt_dm,TM_lt,TM_sdt',lt_dm(i),dt_dm(i) &
               ,TM_lt,TM_sdt
              write(6,*) '  dvx,h',dvx_dm(i),dvy_dm(i),dvz_dm(i),h_dm(i)
              write(6,*) ' rho,nnb=',rho_dm(i),nnb_dm(i)
            adv=dvx_dm(i)*dvx_dm(i)+dvy_dm(i)*dvy_dm(i)+dvz_dm(i)*dvz_dm(i)
            adv=dsqrt(adv)
            if(adv.ne.0.0d0) then
              dtcr2=CGRAV*dsqrt(h_dm(i)*0.5d0/adv)
            endif
              write(6,*) ' dtcr2=',dtcr2
              write(6,*) '  flagt,dtr,dtmnb',flagt_dm(i),dtr_dm(i),dtmnb_dm(i) 
              write(6,*) ' sdt-tlim,/dtmnb=',TM_sdt-tlim &
               ,(TM_sdt-tlim)/dtmnb_dm(i)
              write(6,*) ' x=',x_dm(i),y_dm(i),z_dm(i)
              do j=0,nadm
                if(i.eq.list_adm(j)) then
                  goto 90
                endif
              enddo
   90         write(6,*) myrank,' id=',i,id_dm(i),' in list_adm,nadm',j,nadm
              call MPI_ABORT(MPI_COMM_WORLD,ierr)
              stop
            endif
#ifdef COSM
            endif
#endif

          endif
        enddo
! * check *	  
#ifdef COSM
        if(SI_flagrad.le.-NITEOR.or.SI_flagrad.gt.-1) then
#endif

        if(TM_dt.le.TM_lt+DTLIMIT) then
          write(6,*) ' Error in localtime() : '
          write(6,*) '  local timestep is zero at DM',TM_dt,TM_lt
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
#ifdef COSM
        endif
#endif
#endif
! *** get minimum time step among nodes ***
        tdvs=TM_dt
        call MPI_ALLREDUCE(tdvs,TM_dt,1,MPI_DOUBLE_PRECISION &
         ,MPI_MIN,MPI_COMM_WORLD,ierr)
! * define local time bin *	  
        TM_dt=TM_dt-TM_lt
      endif
#if defined(SF) || defined(SF_EFD) || defined(SF_ZFD)
! *** check SF and FD time step
      if(TM_dt.gt.TM_mindtsfd) then
        TM_dt=TM_mindtsfd
      endif
#endif
! *** check maximum local time bin ***
! *** ensure the block time step ***

      if(TM_dt.gt.TM_ldtmin) then
        TM_kn=int(dlog(MGTU*TM_dt/TM_ldtmin)/dlog(2.0d0))
      else
        TM_kn=int(dlog(MGTL*TM_dt/TM_ldtmin)/dlog(2.0d0))
      endif
      TM_dt=TM_ldtmin*(2.0d0**int(TM_kn))
      TM_kn=int(dlog((MGTU*TM_dt+TM_lt)/TM_ldtmin)/dlog(2.0d0))
      if(TM_k.eq.TM_kn) then
        TM_dt=TM_sdt-TM_lt
      else if(TM_kn.gt.TM_k) then
        write(6,*) ' Error in localtime(): local time is overshoot'
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
        stop
      endif

! *** find active particle ***
      lwt=TM_lt+MGTL*TM_dt
      upt=TM_lt+MGTU*TM_dt
#if defined(GAS) || defined(STAR)

      allocate(talist(0:np))
      allocate(list(0:np))
      allocate(nalist(0:np))

      do i=0,np-1
        talist(i)=list_ap(i)
      enddo
#endif
#ifdef GAS
! *** for Gas Particles ***
! *** compare with hydro dt ***
      do i=0,nagravg-1
        pn=talist(i)
        ltg_p(pn)=TM_lt
        if(dtg_p(pn).lt.dt_p(pn)) then
          dtg_p(pn)=dt_p(pn)
        else
          knp=int(dlog((dtg_p(pn)+TM_dt*0.1d0)/dt_p(pn))/dlog(2.0d0))
          if(knp.gt.FASTGLIMIT) then
            knp=FASTGLIMIT
            dtg_p(pn)=dt_p(pn)*(2.0d0**dble(knp))
          endif
        endif
! *** initial kick ***
        vnx_p(pn)=vnx_p(pn)+0.5d0*dtg_p(pn)*ax_p(pn)
        vny_p(pn)=vny_p(pn)+0.5d0*dtg_p(pn)*ay_p(pn)
        vnz_p(pn)=vnz_p(pn)+0.5d0*dtg_p(pn)*az_p(pn)
      enddo

      na=0
      nagravg=0
      nanp=0
      do i=0,ng-1
        pn=talist(i)		 
        if(lwt.lt.ltg_p(pn)+dtg_p(pn).and.ltg_p(pn)+dtg_p(pn).lt.upt) then
          list_ap(nagravg)=pn
          nagravg=nagravg+1
        else if(lwt.lt.lt_p(pn)+dt_p(pn).and.lt_p(pn)+dt_p(pn).lt.upt) then
          list(na)=pn
          na=na+1
        else
          nalist(nanp)=talist(i)
          nanp=nanp+1
        endif
      enddo
! *** Make list_ap[] ***
      nag=na+nagravg
      j=nagravg
      do i=0,na-1
        list_ap(j)=list(i)
        j=j+1
      enddo
      do i=0,nanp-1
        list_ap(j)=nalist(i)
        j=j+1
      enddo
#endif
#ifdef STAR	
! *** for Star Particles ***
      if(ns.gt.0) then
        na=0
        nanp=0
        do i=ng,np-1
          pn=talist(i)		   
! *** set dtg_p and ltg_p
          dtg_p(pn)=dt_p(pn)
          ltg_p(pn)=lt_p(pn)
          if(lwt.lt.lt_p(pn)+dt_p(pn).and. lt_p(pn)+dt_p(pn).lt.upt) then
            list_ap(ng+na)=pn
            na=na+1
          else
            nalist(nanp)=pn
            nanp=nanp+1
          endif
        enddo
! *** Make list_ap[] ***
        j=ng+na
        nas = na
        do i=0,nanp-1
          list_ap(j)=nalist(i)
          j=j+1
        enddo
      endif
#endif	

#if defined(GAS) || defined(STAR)
      deallocate(talist)
      deallocate(list)
      deallocate(nalist)
#endif

#ifdef DM
! *** for DM Particles ***
      na = 0
      nanp=0
      if(ndm.gt.0) then

        allocate(nalist(0:ndm))

        do i=0,ndm-1
          if(lwt.lt.lt_dm(i)+dt_dm(i).and.lt_dm(i)+dt_dm(i).lt.upt) then
            list_adm(na)=i
            na=na+1
          else
            nalist(nanp)=i
            nanp=nanp+1
          endif
        enddo
        nadm=na
        j=na
        do i=0,nanp-1
          list_adm(j)=nalist(i)
          j=j+1
        enddo

        deallocate(nalist)

      endif
#endif      
! *** setup virtual timestep ***
#ifdef GAS
! *** for Gas and Star Particles ***
      do i=0,np-1
#ifdef FIXEDP
        if(id_p(i).lt.SI_fid(0).or.id_p(i).gt.SI_fid(1)) then
#endif
        vdt_p(i)=TM_lt+TM_dt-lt_p(i)
#ifdef FIXEDP
        else
          vdt_p(i)=0.0d0
        endif
#endif
      enddo		

! *** correction individual time bin of active particle ***
      do i=0,nag-1
        dt_p(list_ap(i))=vdt_p(list_ap(i))
      enddo
      do i=0,nagravg-1
        pn=list_ap(i)
        dtc=TM_lt+TM_dt-ltg_p(pn)
        dtu=0.5d0*(dtc-dtg_p(pn))
! *** correction for velocity due to the difference between
!  initial and final dt***
        vnx_p(pn)=vnx_p(pn)+ax_p(pn)*dtu
        vny_p(pn)=vny_p(pn)+ay_p(pn)*dtu
        vnz_p(pn)=vnz_p(pn)+az_p(pn)*dtu
        dtg_p(pn)=dtc
      enddo

#ifdef STAR
! *** correction individual time bin of active particle ***
      do i=ng,ng+nas-1
        dt_p(list_ap(i))=vdt_p(list_ap(i))
        dtg_p(list_ap(i))=vdt_p(list_ap(i))
      enddo		
#endif
#endif

#ifdef DM
! *** for DM Particles ***
      do i=0,ndm-1
#ifdef FIXEDDMP
        if(id_dm(i).lt.SI_fiddm(0).or.id_dm(i).gt.SI_fiddm(1)) then
#endif

        vdt_dm(i)=TM_lt+TM_dt-lt_dm(i)
#ifdef FIXEDDMP
        else
          vdt_dm(i)=0.0d0
        endif
#endif
      enddo		
! *** correction individual time bin of active particle ***
      do i=0,nadm-1
        dt_dm(list_adm(i))=vdt_dm(list_adm(i))
      enddo
#endif

end subroutine
