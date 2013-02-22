#include "gcdp.def"
! *********************************************************
!  esolve.F95  solving energy equation for GCD+ ver. f03.0
!  3  Feb. 2013   written by D.KAWATA
! *********************************************************
!
! assuming as T increase, heating -> cooling
!          no heating -> cooling -> heating 
!  use dt_p for individual time step
!
! common valubes used for work... dti,tx,du0,ulim,hu,lu,fx,adv,av
!   pn_nfp,node
!
! for subtime step, update asp (and p_p), assume rho,Gsn_p and etc.,
!  and UVB are constant
!
! nnb_p should be updated.   
!

#ifdef COOL
subroutine energy_solver(ngs,nlist0,plist,current_z,flagdu)
      use gcdp_const
      use gcdp_baryon
#ifdef METAL
      use gcdp_metal
#endif
      use gcdp_system

      implicit none

      integer,intent(in) :: nlist0
      integer,intent(in) :: ngs,flagdu
      integer,intent(in) :: plist(0:nlist0)
      double precision,intent(in) :: current_z

      integer it,maxnt,nlist
      double precision gkmm,dudt,dtc
      double precision tcool,tdyn     
! * for iteration *
! * Number of Notfinished particle, temp *      
      integer i,pn,tnlist,ie,mite
      integer nptc
      character filen*60
! work allocatable
      integer,allocatable :: flagu(:),palist(:),flagu1(:),flagtc(:),node(:)
      integer,allocatable :: pn_nfp(:),clist(:)
      double precision,allocatable :: asp(:),uhlim(:),rhogam(:),gamrhogam1(:)
      double precision,allocatable :: tx(:),adv(:),av(:),ulim(:),dti(:),du0(:) &
       ,ram(:),hu(:),lu(:),fhx(:),fx(:)

      mite=50000
      maxnt=1000000000
      gkmm=(GAM-1.0d0)/((KCGS/(MP*MYU))/K_MU)

      allocate(ulim(0:ngs))
      allocate(uhlim(0:ngs))
      allocate(tx(0:ngs))
      allocate(adv(0:ngs))
      allocate(av(0:ngs))
      allocate(asp(0:ngs))
      allocate(flagtc(0:ngs))
      allocate(rhogam(0:ngs))
      allocate(gamrhogam1(0:ngs))
      allocate(dti(0:ngs))
      allocate(du0(0:ngs))
      allocate(ram(0:ngs))
      allocate(hu(0:ngs))
      allocate(lu(0:ngs))
      allocate(fhx(0:ngs))
      allocate(fx(0:ngs))
      allocate(flagu(0:ngs))
      allocate(flagu1(0:ngs))

      allocate(pn_nfp(0:nlist0))
      allocate(palist(0:nlist0))

      nptc=nlist0
      do i=0,nlist0-1
        pn = plist(i)
        pn_nfp(i)=pn
        ulim(pn)=(TLLIMIT/gkmm)*(MYU/myu_p(pn))
        ulim(pn)=ulim(pn)*(GAM-1.0d0)/(rho_p(pn)**(GAM-1.0d0))
! *** highest temperature limit ***
        uhlim(pn)=(THLIMIT/gkmm)*(MYU/myu_p(pn))
! entropy
        uhlim(pn)=uhlim(pn)*(GAM-1.0d0)/(rho_p(pn)**(GAM-1.0d0))
! for check
        tx(pn)=myu_p(pn)
! initial ram
        adv(pn)=ram_p(pn)
        av(pn)=0.0d0
! initial entropy
        asp(pn)=as_p(pn)
! flag for tc
        flagtc(pn)=1
! parameters for u vs. entropy conversion
        gamrhogam1(pn)=(GAM-1.0d0)/(rho_p(pn)**(GAM-1.0d0))
        rhogam(pn)=rho_p(pn)**GAM
      enddo
#ifdef SF_EFD
      if(flagdu.eq.1) then
        do i=0,nlist0-1
          pn = plist(i)
! *** only flagdu.eq.1
! *** check cooling time-scale vs. dynamical time scale for SNe II particle
          if(flagc_p(pn).eq.-1) then
            if(ram_p(pn).gt.0.0d0) then
              tcool=u_p(pn)/ram_p(pn)
              tdyn=dsqrt(3.0d0*M_PI/(16.0d0*G*rho_p(pn)))
              if(tcool.gt.tdyn.and.div_v_p(pn).gt.0.0d0) then
! *** turn on cooling if tcool > tdyn
                flagc_p(pn)=0
                ram_p(pn)=0.0d0
                adv(pn)=0.0d0
              endif
            else
              flagc_p(pn)=0
            endif
          endif
        enddo
      endif
#endif

      it=0
   71 it=it+1
      if(it.eq.maxnt) then
        write(6,*) ' Error in esolve() myrank=',myrank
        write(6,*) ' it > maxnt, nptc=',nptc
      endif

!      write(6,*) 'it,nptc=',it,nptc

! *** check cooling time vs. dt_p and set dti ***
      do i=0,nptc-1
        pn=pn_nfp(i)
! *** calculate total du/dt
        if(flagdu.eq.0) then
          dudt=0.5d0*(pdu_p(pn)+ndu_p(pn))-adv(pn)
        else
          dudt=ndu_p(pn)-adv(pn)
        endif
        dudt=dabs(dudt)
#ifdef SF_EFD
! *** adiabatic if flagc_p=-1
        if(flagc_p(pn).eq.0.and.dudt.gt.0.0d0) then
#else
        if(dudt.gt.0.0d0) then
#endif
          dtc=CCOOLEDT*u_p(pn)/dudt
          dti(pn)=dtc
          if(av(pn)+dti(pn).gt.dt_p(pn)) then
            dti(pn)=dt_p(pn)-av(pn)
            flagtc(pn)=0              
          endif 
        else
          dti(pn)=dt_p(pn)-av(pn)
          flagtc(pn)=0
        endif
      enddo

      nlist = 0
! *** make a list for iteration ***
      do i=0,nptc-1
        pn=pn_nfp(i)
        if(flagdu.eq.0) then
          dudt=0.5d0*dti(pn)*(pdu_p(pn)+ndu_p(pn))
        else
          dudt=dti(pn)*ndu_p(pn)
        endif
#ifdef SF_EFD
! *** no cooling on feebdack particle flagc_p=-1 ***
        if(flagc_p(pn).eq.-1.and.Gsn_p(pn).gt.MININFESN*u_p(pn)) then
! once SNe II finish, cooling turned on. 
! *** only active particles since pv30.5
          asp(pn)=asp(pn)+dudt+gamrhogam1(pn)*Gsn_p(pn)*dti(pn)/dt_p(pn)
          p_p(pn)=rhogam(pn)*asp(pn)
          ram_p(pn)=0.0d0
          flagu(pn)=0
        else
#endif
! *** a_n+0.5*(da/dt)_n ***
          du0(pn)=asp(pn)+dudt &
#ifdef SF_EFD
! *** ram at step i
           -0.5d0*dti(pn)*gamrhogam1(pn)*adv(pn) &
! *** fraction of feedback in dti
           +gamrhogam1(pn)*Gsn_p(pn)*dti(pn)/dt_p(pn)
#else
           -0.5d0*dti(pn)*gamrhogam1(pn)*adv(pn)
#endif
          ram(pn) = 0.0d0
! *** check minimum entropy ***
          if(du0(pn).lt.ulim(pn)) then  
! *** note asp = du0
            asp(pn)=0.5d0*asp(pn)
            du0(pn)=asp(pn)
            flagu1(pn)=0
!            flagu1(pn)=-1
! *** check max entropy ***
          else if(du0(pn).gt.uhlim(pn)) then
            asp(pn)=uhlim(pn)
            flagu1(pn)=1
          else
            asp(pn) = du0(pn)
            flagu1(pn)=0
          endif
          hu(pn)=dlog10(asp(pn))
          lu(pn)=dlog10(asp(pn))
          fhx(pn)=-1.0d0
          fx(pn)=1.0d0
          palist(nlist)=pn
          flagu(pn)=1
          nlist=nlist+1
! no option for adiabatic feedback
#ifdef SF_EFD
        endif
#endif
      enddo
! * bisec iteration start from asp *
      ie=0
   70 if(nlist.eq.0) then
        goto 90
      endif
      ie = ie+1
      do i=0,nlist-1
        pn=palist(i)
! *** p is used in cool ***
        p_p(pn)=rhogam(pn)*asp(pn)      
      enddo

! *** set ram ***

! reallocate to clist length with nlist from palist
      allocate(clist(0:nlist))

      do i=0,nlist-1 
        clist(i)=palist(i)
      enddo

      call cool(nlist,clist,current_z,ngs,ram)

      deallocate(clist)

      if(ie.eq.1) then
        do i=0,nlist-1
! *** flagu: -1,1: changing u (-1:decrease, 1:increase)
!  2: bisection, 0: finish ***
          pn = palist(i)         
          if(ram(pn).gt.0.0d0) then
            flagu(pn) = -1
          else
            flagu(pn) = 1
          endif  
! *** check the lower temperature limit at 1st iteration ***
          if(flagu1(pn).eq.-1.and.du0(pn)-0.5d0*dti(pn)*ram(pn) &
           *gamrhogam1(pn).lt.asp(pn)) then
            asp(pn) = ulim(pn)
            flagu(pn) = 0
          endif
! *** check the highest temperature limit ***
          if(flagu1(pn).eq.1.and.du0(pn)-0.5d0*dti(pn)*ram(pn) &
           *gamrhogam1(pn).gt.asp(pn)) then
            asp(pn) = uhlim(pn)
            flagu(pn) = 0
#ifdef SF_EFD
!                                  123456789012345678901234
            write(6,'(a24,3I10)') 'Thlim1 myrank,id,flagfd=' &
             ,myrank,id_p(pn),flagfd_p(pn)
#endif
          endif
        enddo
      endif
! #pragma loop novrec hu,lu,u_p,p_p,rho_p,fx,flagu
      do i=0,nlist-1
! *** flagu: -1,1: changing u (-1:decrease, 1:increase)
!  2: bisection, 0: finish ***
        pn = palist(i)         
        if(flagu(pn).eq.1) then
! *** increasing u ***
          fhx(pn) = du0(pn)-0.5d0*dti(pn)*ram(pn)*gamrhogam1(pn)-asp(pn)
          if(dabs(fhx(pn))/asp(pn).lt.EPSC) then
            flagu(pn)=0
          endif
          if((ie.eq.1.and.flagu(pn).ne.0).and.fhx(pn).le.0.0d0) then
            write(6,*) ' initial fhx<0 T=',(as_p(pn)*myu_p(pn)/(GAM-1.0d0)) &
             *(rho_p(pn)**(GAM-1.0d0))*gkmm,(asp(pn)*myu_p(pn)/(GAM-1.0d0)) &
             *(rho_p(pn)**(GAM-1.0d0))*gkmm
            write(6,*) ' ram_p,ram,flagu=',ram_p(pn),ram(pn),flagu(pn)
          endif
        else if(flagu(pn).ne.0) then
! *** decreasing u ***
          fx(pn) = du0(pn)-0.5d0*dti(pn)*ram(pn)*gamrhogam1(pn)-asp(pn)
          if(dabs(fx(pn))/asp(pn).lt.EPSC) then
            flagu(pn)=0
          endif
          if((ie.eq.1.and.flagu(pn).ne.0).and.fx(pn).ge.0.0d0) then
            write(6,*) ' initial fx>0 T,fx=',(as_p(pn)*myu_p(pn)/(GAM-1.0d0)) &
             *(rho_p(pn)**(GAM-1.0d0))*gkmm,(asp(pn)*myu_p(pn)/(GAM-1.0d0)) &
             *(rho_p(pn)**(GAM-1.0d0))*gkmm,fx(pn),du0(pn) &
             ,du0(pn)-0.5d0*dt_p(pn)*ram(pn)*(gamrhogam1(pn)) &
             ,as_p(pn),asp(pn)
            write(6,*) ' ram_p,ram,flagu=',ram_p(pn),ram(pn),flagu(pn)
            stop
          endif
        endif
      enddo
      if(ie.gt.1) then
        do i=0,nlist-1
          pn = palist(i)         
! *** check the convergence: bisection ***
          if(flagu(pn).eq.2.and.dabs(hu(pn)-lu(pn)).lt.EPSC) then
            flagu(pn) = 0
          endif
          if(flagu(pn).eq.2.and.fx(pn)*fhx(pn).le.0.0d0) then
            lu(pn) = dlog10(asp(pn))
          else if(flagu(pn).eq.2) then
            hu(pn) = dlog10(asp(pn))
            fhx(pn) = fx(pn)
          else if(flagu(pn).ne.0.and.fx(pn)*fhx(pn).le.0.0d0) then
! *** start bisection ***
            flagu(pn) = 2
          endif
        enddo
      endif

! *** update fx, fhx, lu, hu ***
      do i=0,nlist-1
        pn = palist(i)         
        if(flagu(pn).eq.1) then
! *** increasing as ***
! *** check the highest temperature limit ***
          if(asp(pn).gt.uhlim(pn)) then

#ifdef SF_EFD
!                                  123456789012345678901234567890
            write(6,'(a30,5I10)') 'Thlim1 myrank,id,flagfd,ie,it=' &
             ,myrank,id_p(pn),flagfd_p(pn),ie,it
            write(6,*) 'flagdu,Gsn,ndu,ram,asp,uhlim=' &
             ,flagdu,Gsn_p(pn),ndu_p(pn),ram_p(pn),asp(pn),uhlim(pn)
            write(6,*) 'flagtc,flagu,dti,dt=' &
             ,flagtc(pn),flagu(pn),dti(pn),dt_p(pn)

            stop
#endif

            asp(pn) = uhlim(pn)
            flagu(pn) = 0
          else
            lu(pn) = hu(pn)
            hu(pn) = hu(pn)+DU_UPU
            asp(pn) = 10.0d0**(hu(pn))
            fx(pn) = fhx(pn)
          endif
        else if(flagu(pn).eq.-1) then
! *** decreasing a ***
! *** check the lower temperature limit ***
          if(asp(pn).lt.ulim(pn)) then
            asp(pn) = ulim(pn)
            flagu(pn) = 0
          else
            hu(pn) = lu(pn)
            lu(pn) = lu(pn)-DU_UPU
            asp(pn) = 10.0d0**(lu(pn))
            fhx(pn) = fx(pn)
          endif
        else if(flagu(pn).eq.2) then
          asp(pn) = 10.0d0**(0.5d0*(lu(pn)+hu(pn)))
        endif
      enddo

! * change the list of not finished particles *

      allocate(node(0:nlist))

      tnlist = nlist
      do i=0,tnlist-1
        node(i) = palist(i)
      enddo
      nlist = 0
      do i=0,tnlist-1
        pn=node(i)
        if(flagu(pn).ne.0) then
          palist(nlist)=node(i)
          nlist=nlist+1
        else
          if(asp(pn).lt.ulim(pn)) then
            asp(pn) = ulim(pn)
          else if(asp(pn).gt.uhlim(pn)) then
            asp(pn) = uhlim(pn)
          endif
          p_p(pn)=rhogam(pn)*asp(pn)
        endif
      enddo

      deallocate(node)

      if(ie.gt.mite.and.nlist.gt.0) then
        write(6,*) 'Error in esolve() : rank,flag,it=',myrank,flagdu,it
        write(6,*) '  No convergence ',nlist,' particles !'
        do i=0,nlist-1
          pn=palist(i)
          write(6,*) ' flagu,T(MYU),fx/as, ram,fx,fhx, Tl,Tu ='
          write(6,*) flagu(pn),(as_p(pn)*MYU/(GAM-1.0d0)) &
           *(rho_p(pn)**(GAM-1.0d0))*gkmm,fx(pn)/asp(pn)
          write(6,*) ram(pn),fx(pn),fhx(pn)
          write(6,*) (10.0**(lu(pn))/(GAM-1.0d0))*(rho_p(pn)**(GAM-1.0d0)) &
           *gkmm,(10.0**(hu(pn))/(GAM-1.0d0))*(rho_p(pn)**(GAM-1.0d0))*gkmm
          write(6,*) ' du, pdu, duram,nnb=',ndu_p(pn),pdu_p(pn) &
           ,(gamrhogam1(pn))*ram_p(pn),nnb_p(pn)
          write(6,*) ' alpv,u,u2,div,arot,h,omgh=',alpv_p(pn) &
           ,alpu_p(pn),d2u_p(pn),div_v_p(pn),arot_v_p(pn),h_p(pn),omgh_p(pn)
          write(6,*) ' id, original as,rho,nh,ram,myu,T,Y,Z=' &
           ,id_p(pn),as_p(pn),rho_p(pn) &
#ifdef METAL
           ,((m_p(pn)-((mzZ_p(pn)+mzHe_p(pn))/MUSM))/m_p(pn)) &
#else
           ,0.0d0 &
#endif
           ,rho_p(pn)*(DU/MP),ram_p(pn),tx(pn),(as_p(pn)*(tx(pn)/MYU) &
           /(GAM-1.0d0))*(rho_p(pn)**(GAM-1.0d0))*gkmm &
#ifdef METAL
           ,mzHe_p(pn)/(m_p(pn)*MUSM),mzZ_p(pn)/(m_p(pn)*MUSM) 
#else
           ,0.0d0,0.0d0
#endif
#ifdef SF_EFD
          write(6,*) ' Gsn,dtp,dtsf,TMd,flagfd,rfd,fc',Gsn_p(pn) &
          ,dt_p(pn),TMsf_dt,TM_dt,flagfd_p(pn),flagrfd_p(pn),flagc_p(pn)
#endif
        enddo
        stop
      endif
      goto 70
! *** one step of dti 
! *** update the list
   90 tnlist=nptc
    
      allocate(node(0:nptc))

      do i=0,tnlist-1
        node(i)=pn_nfp(i)
      enddo
      nptc=0
      do i=0,tnlist-1
        pn=node(i)
        if(flagtc(pn).eq.1) then
          pn_nfp(nptc)=pn
          nptc=nptc+1
! *** cooling rate at step i
          adv(pn)=ram(pn)
! *** update time
          av(pn)=av(pn)+dti(pn)
        endif
      enddo

      deallocate(node)

      if(nptc.gt.0) then
        goto 71
      endif

      deallocate(pn_nfp)
      deallocate(ulim)
      deallocate(uhlim)
      deallocate(tx)
      deallocate(adv)
      deallocate(av)
      deallocate(asp)
      deallocate(flagtc)
      deallocate(rhogam)
      deallocate(gamrhogam1)
      deallocate(dti)
      deallocate(du0)
      deallocate(ram)
      deallocate(hu)
      deallocate(lu)
      deallocate(fhx)
      deallocate(fx)
      deallocate(palist)
      deallocate(flagu)
      deallocate(flagu1)

end subroutine
#endif

