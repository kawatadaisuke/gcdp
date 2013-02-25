#include "gcdp.def"
! *****************************************
!     starfd.F for GCD+ ver. f03.0
!  27 Jan. 2013    written by D.KAWATA
! ******************************************
! *********************************
!  Definition star()
!    Process about star formation
! *********************************
#if defined(SF_EFD) || defined(SF_ZFD)
subroutine starfd(ng,ns,nag,nas,nagravg)
      use gcdp_const
      use gcdp_baryon
      use gcdp_metal
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(inout) :: ng,ns,nag,nas,nagravg
      integer i,j,pn,ngfd,ngfdt,dngfd,dngfdt,np,ierr
      integer newnas,newns
      integer na,nna,dns,dnsne,npsnia,dnsnet
      double precision fzms
! *** for check ***
      double precision desns,tdesns,dzsns(0:NYEL-1),tdzsns(0:NYEL-1)
! *** for time ***
      double precision ntime,ectime
! for work
! * for making star particle *
      integer,allocatable :: alistg(:),list(:),talist(:),nalist(:) &
       ,tivs(:),tivr(:)    
! * Number , Mass of SN *	  
      double precision,allocatable :: mfd(:),tdvs(:),tdvr(:)
! * for feedback  common *
! * Total *
      double precision,allocatable :: tmej_sn(:),tnsn_sn(:),tmzHe_sn(:) &
       ,tmzC_sn(:),tmzN_sn(:),tmzO_sn(:),tmzNe_sn(:),tmzMg_sn(:) &
       ,tmzSi_sn(:),tmzFe_sn(:),tmzZ_sn(:),nsw_sn(:)

! ***   Calculate basic Value   ***
! *** for check ***
      desns=0.0d0
      do i=0,NYEL-1
        dzsns(i)=0.0d0
      enddo
      dnsne=0
      np=ng+ns

! *** calculate yields within this time step ***
! *** set value to t*_sn() at TM_tot+TM_dt ***

      dngfd=0
      newns=0
      newnas=0
! *** should not initialise flagrfd_p, because it's setted in calc_dv_du ***
! allocation for feedback values
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

      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)
! *** from star to gas ***
! *** and make the candidate to gas particles (flagrfd.ne.0) ***
      allocate(mfd(0:np))
      allocate(list(0:ns))
      allocate(talist(0:ns))

      nna=0
      na=0
      newns=0
      do i=0,ns-1
        pn=list_ap(i+ng)
! *** flagrfd_p is 1, if star-> feedback gas at the next step ***
        if(flagrfd_p(pn).eq.1) then
! *** from star to gas ***
          list(dngfd)=pn
          dngfd=dngfd+1
          flagc_p(pn)=0
! since pv33.19: set flagc_p=-1 to set the particle no cooling
          if(flagfd_p(pn).le.SI_snii) then
            flagc_p(pn)=-1
          endif
          flagrfd_p(pn)=0
! keep flagfd_p positive, but updated in main.F
!          flagfd_p(pn)=-flagfd_p(pn)
! *** for time ***
          dt_p(pn)=TM_lt+TM_dt-lt_p(pn)
          dtg_p(pn)=dt_p(pn)
          ltg_p(pn)=lt_p(pn)
          flagt_p(pn)=2
          dtr_p(pn)=TM_dt
! *** check feedback energy rate ***
#ifdef SF_EFD
! *** number of supernovae ***
          mfd(pn)=tnsn_sn(pn)-tnsn_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          if(nsw_sn(pn).lt.0.0d0) then
            nsw_sn(pn)=0.0d0
          endif
! *** nsw_sn includes delta t ***
          Gsn_p(pn)=mfd(pn)*SNEU+SWEU*nsw_sn(pn)
#ifdef CONSTFDE
! *** erg/s unit 5.752e43
          Gsn_p(pn)=(1.739e-6)*TMsf_dt
#endif
! *** Gsn_p is energy gain per mass ***
          Gsn_p(pn)=Gsn_p(pn)/m_p(pn)
#endif
! not for dv?_p, since it is used correctly in updatev
!  not set p, cs, as
          ndu_p(pn)=0.0d0
          pdu_p(pn)=0.0d0
          myu_p(pn)=MYU
          alpv_p(pn)=V_ALPHAINI            
          alpu_p(pn)=0.0d0
          d2u_p(pn)=0.0d0
          arot_v_p(pn)=0.0d0
          hvsigdt_p(pn)=INF
#ifdef COOL
          ram_p(pn)=0.0d0
#endif
#ifdef METALDIFF
          vsig_p(pn)=0.0d0
#endif
        else
! *** active star particles ***
          talist(newns)=pn
          newns=newns+1
          if(i.lt.nas) then
            newnas=newnas+1
          endif
        endif
! dtyield only gives _sn values for age>tspiz, ie flagrfd=1
! but, here updates all. 
! *** update these values ***
        tnsn_p(pn)=tnsn_sn(pn)
        tmej_p(pn)=tmej_sn(pn)
! * Unit Solar Mass *
        tmzHe_p(pn)=tmzHe_sn(pn)
        tmzC_p(pn)=tmzC_sn(pn)
        tmzN_p(pn)=tmzN_sn(pn)
        tmzO_p(pn)=tmzO_sn(pn)
        tmzNe_p(pn)=tmzNe_sn(pn)
        tmzMg_p(pn)=tmzMg_sn(pn)
        tmzSi_p(pn)=tmzSi_sn(pn)
        tmzFe_p(pn)=tmzFe_sn(pn)
        tmzZ_p(pn)=tmzZ_sn(pn)
      enddo
! *** change the list ***
      if(dngfd.gt.0) then
        nna=0
        na=0

        allocate(alistg(0:nag))
        allocate(nalist(0:ng))

        do i=0,nag-1
          alistg(na)=list_ap(i)
          na=na+1
        enddo
        do i=nag,ng-1
          nalist(nna)=list_ap(i)
          nna=nna+1
        enddo
! *** re-arrange list_ap ***
! *** add new feeback gas to the active gas particle list ***
        do i=0,dngfd-1
          list_ap(i)=list(i)
        enddo
        nagravg=nagravg+dngfd
        nag=nag+dngfd
! original active gas particle
        do i=0,na-1
          list_ap(i+dngfd)=alistg(i)
        enddo
        do i=0,nna-1
          list_ap(i+nag)=nalist(i)
        enddo
        ng=ng+dngfd
        do i=0,newns-1
          list_ap(i+ng)=talist(i)
        enddo
        nas=newnas
        ns=newns

        deallocate(alistg)
        deallocate(nalist)

      endif

      deallocate(list)
      deallocate(talist)
      deallocate(mfd)

! *** set Feedback ***
! *** feedback from gas with flagfd_p<0 ***

      allocate(mfd(0:np))

      ngfd=0
      do i = 0,ng-1
        pn = list_ap(i)
        if(flagfd_p(pn).lt.0) then
! excluding new feedback gas particle flagfd>0
          ngfd=ngfd+1
! *** if age limit, flagrfd_p(pn)=-1 setted in dtyields ***
          if(flagrfd_p(pn).eq.-1) then
            flagfd_p(pn)=0
          endif
          flagrfd_p(pn)=1
! *** number of supernovae ***
          mfd(pn)=tnsn_sn(pn)-tnsn_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          tnsn_p(pn)=tnsn_p(pn)+mfd(pn)
#ifdef SF_EFD
          if(nsw_sn(pn).lt.0.0d0) then
            nsw_sn(pn)=0.0d0
          endif
          Gsn_p(pn)=mfd(pn)*SNEU+SWEU*nsw_sn(pn)
#ifdef CONSTFDE
! *** erg/s unit 5.752e43
          Gsn_p(pn)=(1.739e-6)*Tmsf_dt
#endif
! *** Gsn_p is energy gain per mass ***
          Gsn_p(pn)=Gsn_p(pn)/m_p(pn)
#endif
#ifdef SF_ZFD
! *** for Heavy Elements ***          
! *** He ***
          mfd(pn)=tmzHe_sn(pn)-tmzHe_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzHe_p(pn)=mzHe_p(pn)+mfd(pn)
          tmzHe_p(pn)=tmzHe_p(pn)+mfd(pn)
          if(mzHe_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzHe_p(pn)
            mzHe_p(pn)=0.0d0
          endif
          tmzHe_p(pn)=tmzHe_p(pn)+mfd(pn)
          dzsns(0)=dzsns(0)+mfd(pn)
! *** C ***
          mfd(pn)=tmzC_sn(pn)-tmzC_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzC_p(pn)=mzC_p(pn)+mfd(pn)
          if(mzC_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzC_p(pn)
            mzC_p(pn)=0.0d0
          endif
          dzsns(1)=dzsns(1)+mfd(pn)
          tmzC_p(pn)=tmzC_p(pn)+mfd(pn)
! *** N ***
          mfd(pn)=tmzN_sn(pn)-tmzN_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzN_p(pn)=mzN_p(pn)+mfd(pn)
          if(mzN_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzN_p(pn)
            mzN_p(pn)=0.0d0
          endif
          dzsns(2)=dzsns(2)+mfd(pn)
          tmzN_p(pn)=tmzN_p(pn)+mfd(pn)
! *** O ***
          mfd(pn)=tmzO_sn(pn)-tmzO_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzO_p(pn)=mzO_p(pn)+mfd(pn)
          if(mzO_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzO_p(pn)
            mzO_p(pn)=0.0d0
          endif
          dzsns(3)=dzsns(3)+mfd(pn)
          tmzO_p(pn)=tmzO_p(pn)+mfd(pn)
! *** Ne ***
          mfd(pn)=tmzNe_sn(pn)-tmzNe_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzNe_p(pn)=mzNe_p(pn)+mfd(pn)
          if(mzNe_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzNe_p(pn)
            mzNe_p(pn)=0.0d0
          endif
          dzsns(4)=dzsns(4)+mfd(pn)
          tmzNe_p(pn)=tmzNe_p(pn)+mfd(pn)
! *** Mg ***
          mfd(pn)=tmzMg_sn(pn)-tmzMg_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzMg_p(pn)=mzMg_p(pn)+mfd(pn)
          if(mzMg_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzMg_p(pn)
            mzMg_p(pn)=0.0d0
          endif
          dzsns(5)=dzsns(5)+mfd(pn)
          tmzMg_p(pn)=tmzMg_p(pn)+mfd(pn)
! *** Si ***
          mfd(pn)=tmzSi_sn(pn)-tmzSi_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzSi_p(pn)=mzSi_p(pn)+mfd(pn)
          if(mzSi_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzSi_p(pn)
            mzSi_p(pn)=0.0d0
          endif
          dzsns(6)=dzsns(6)+mfd(pn)
          tmzSi_p(pn)=tmzSi_p(pn)+mfd(pn)
! *** Fe ***
          mfd(pn)=tmzFe_sn(pn)-tmzFe_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          mzFe_p(pn)=mzFe_p(pn)+mfd(pn)
          if(mzFe_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzFe_p(pn)
            mzFe_p(pn)=0.0d0
          endif
          dzsns(7)=dzsns(7)+mfd(pn)
          tmzFe_p(pn)=tmzFe_p(pn)+mfd(pn)
! *** Z ***
          mfd(pn)=tmzZ_sn(pn)-tmzZ_p(pn)
          if(mfd(pn).lt.0.0d0) then
            mfd(pn)=0.0d0
          endif
          if(mzZ_p(pn)+mfd(pn).gt.0.95d0*MUSM*m_p(pn)) then
            write(6,*) ' Error: Z>0.95 myrank,id=',myrank,id_p(pn)
            write(6,*) ' flagfd,age,Z=',flagfd_p(pn) &
             ,(TM_tot+TM_dt-ts_p(pn))*TMUGYR,mzZ_p(pn)/(MUSM*m_p(pn))
            write(6,*) ' tmzsn,tmzp=',tmzZ_sn(pn),tmzZ_p(pn)
            stop
          endif
          mzZ_p(pn)=mzZ_p(pn)+mfd(pn)
          if(mzZ_p(pn).lt.0.0d0) then
            mfd(pn)=mfd(pn)-mzZ_p(pn)
            mzZ_p(pn)=0.0d0
          endif
          tmzZ_p(pn)=tmzZ_p(pn)+mfd(pn)
#endif
! for check
#ifdef SF_EFD
          desns = desns+Gsn_p(pn)*m_p(pn)
#else
! put a only mass data
          desns = desns+m_p(pn)
#endif
          dzsns(8)=dzsns(8)+mfd(pn)
! not used
          tmej_p(pn)=tmej_sn(pn)
        endif
      enddo

      deallocate(mfd)
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

      ngfdt=0
      allocate(tivs(0:1))
      allocate(tivr(0:1))

      tivs(0)=ngfd
      tivs(1)=dngfd
      tivr(0)=0
      tivr(1)=0
      call MPI_ALLREDUCE(tivs,tivr,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      ngfdt=tivr(0)
      dngfdt=tivr(1)

      deallocate(tivs)
      deallocate(tivr)

! *** flagrfd_p not 0, gas got feedback at the previous time step
!     and feedback gas in this time step 
!     they will be active particle and set dt=TM_dt ***
      if(ngfdt.gt.0) then
#ifdef MINDT_SFD

        allocate(talist(0:ng))
        allocate(nalist(0:ng))

        na=0
        nna=0
        do i=0,ng-1
! *** for the list ***
          pn=list_ap(i)
          if(flagrfd_p(pn).ne.0.and.i.ge.nag) then
! *** only change dt_p, but not dtg_p ***
            dt_p(pn)=TM_lt+TM_dt-lt_p(pn)
            dtr_p(pn)=TM_dt
            flagt_p(pn)=1
            talist(na)=pn
            na=na+1
          else if(i.ge.nag) then
            nalist(nna)=pn
            nna=nna+1
          endif
        enddo
        do i=0,na-1
          list_ap(nag+i)=talist(i)
        enddo
        nag=nag+na
        pn=nag
        do i=0,nna-1
          list_ap(pn)=nalist(i)
          pn=pn+1
        enddo

        deallocate(talist)
        deallocate(nalist)
#else
! *** require minum dt if feedback particles
        do i=0,ng-1
          pn=list_ap(i)
          if(flagrfd_p(pn).ne.0) then
!          if(flagfd_p(pn).ne.0.and.Gsn_p(pn).gt.0.0d0) then
            dtmnb_p(pn)=TM_dt
          endif
        enddo
#endif
! *** for check ***
        allocate(tdvs(0:NYEL))
        allocate(tdvr(0:NYEL))
        tdvs(0)=desns
        do i=1,NYEL
          tdvs(i)=dzsns(i-1)
        enddo
        do i=0,NYEL
          tdvr(i)=0.0d0
        enddo
        call MPI_ALLREDUCE(tdvs,tdvr,NYEL+1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
        SI_esns=SI_esns+tdvr(0)
        do i=1,NYEL
          SI_zsns(i-1)=SI_zsns(i-1)+tdvr(i)
        enddo

        deallocate(tdvs)
        deallocate(tdvr)

      endif


end subroutine
#endif
