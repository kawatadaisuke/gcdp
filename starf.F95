#include "gcdp.def"
! *****************************************
!     starf.F for GCD+ ver. f03.0
!  27 Jan. 2013    written by D.KAWATA
! ******************************************
! *********************************
!  Definition star()
!    Process about star formation
! *********************************

#ifdef SF
subroutine starf(ng,ns,nag,nas,nagravg)
      use gcdp_const
      use gcdp_baryon
      use gcdp_system
      use gcdp_metal

      implicit none
      include 'mpif.h'

      integer,intent(inout) :: ng,ns,nag,nas,nagravg
      integer i,j,pn,ifd,np,ierr
! * for criterion *
      double precision nhp,temp,tdyn
      double precision peff
!      double precision alpvir
! * for star formation *
! * for making star particle *
      integer na,nna,dns,npsnia,nagrav
      integer newns,newnas
      double precision zn
      double precision dtgc
! * for checking low metallicity gas *
      double precision fzms
! for work
      integer,allocatable :: slist(:),list(:),talist(:),nalist(:) &
       ,tivs(:),tivr(:)
! * for feedback *
      double precision,allocatable :: p(:)

      real ran1
      external ran1

! ***   Calculate basic Value   ***
! *** Initialization ***
      np=ng+ns

! *** gas->star ***
      allocate(p(0:np))
! these list is smaller than ng
      allocate(slist(0:ng))
      allocate(list(0:ng))
      allocate(talist(0:ng))
      allocate(nalist(0:ng))

      nagrav=0
      dns=0
      na=0
      nna=0
      do i=0,ng-1
        pn=list_ap(i)
        p(pn)=0.0d0
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).eq.0) then
#endif
! *** since pv34.6.7 to suppress accidental SF due to v field ~ 0 at low density
          nhp=XHSOL*rho_p(pn)*(DU/MP)
          if(nhp.gt.NSTHLIM) then
! since 34.6.8 
            peff=1.2d0*(ULJMFAC**(2.0d0/3.0d0))*G*(h_p(pn)**2) &
             *(rho_p(pn)**2)/GAM
            if(p_p(pn).lt.peff) then
! *** since pv34.6.7 to suppress accidental SF due to v field ~ 0 at low density
! calculate alpha_vir using Hopkins13 
!            alpvir=BETASF*(div_v_p(pn)**2+arot_v_p(pn)**2)/(G*rho_p(pn))
!            if(alpvir.lt.ALPHASF) then
              tdyn=dsqrt(3.0d0*M_PI/(16.0d0*G*rho_p(pn)))
              p(pn)=(1.0d0-exp(-CSEFF*TM_dt/tdyn))
              if(p(pn).ne.0.0d0.and.p(pn).gt.dble(ran1(idum))) then
                p(pn)=2.0d0
              endif
            endif
          endif
#if defined(SF_EFD) || defined(SF_ZFD)
        endif
#endif
        if(p(pn).gt.1.5d0) then
! * make star particle *
          slist(dns)=pn
          dns=dns+1
          flagc_p(pn) = 1
#if defined(SF_EFD) || defined(SF_ZFD)
! *** set flagfd_p
          ifd=int(ran1(idum)*dble(SI_nsp))+1
          flagfd_p(pn)=ifd
          flagrfd_p(pn)=0
#endif
          ts_p(pn)=TM_tot
! *** fixed softening ***
          h_p(pn)=SI_eps*(m_p(pn)**THIRD)
          omgh_p(pn)=1.0d0
          zetah_p(pn)=0.0d0
! *** keep the data of density and pressure, when star forms ***
          p_p(pn) = (GAM-1.0d0)*u_p(pn)*rho_p(pn)
          myu_p(pn)=MYU
#ifdef COOL
          ram_p(pn)=0.0d0
#endif
#ifdef METAL
! *** for metal ***
          mzHe0_p(pn)=mzHe_p(pn)
          mzC0_p(pn)=mzC_p(pn)
          mzN0_p(pn)=mzN_p(pn)
          mzO0_p(pn)=mzO_p(pn)
          mzNe0_p(pn)=mzNe_p(pn)
          mzMg0_p(pn)=mzMg_p(pn)
          mzSi0_p(pn)=mzSi_p(pn)
          mzFe0_p(pn)=mzFe_p(pn)
          mzZ0_p(pn)=mzZ_p(pn)
#endif
! *** for time ***
          dt_p(pn)=TM_lt+TM_dt-lt_p(pn)
! *** correction of velocity due to the change of dtg_p***
! *** overestimated delta t for initial kick ***
          dtgc=0.5d0*(TM_lt+TM_dt-(dtg_p(pn)+ltg_p(pn)))
! *** corrected vn?_p ***
          vnx_p(pn)=vnx_p(pn)+dtgc*ax_p(pn)
          vny_p(pn)=vny_p(pn)+dtgc*ay_p(pn)
          vnz_p(pn)=vnz_p(pn)+dtgc*az_p(pn)
          dtg_p(pn)=TM_lt+TM_dt-ltg_p(pn)
          ltg_p(pn)=lt_p(pn)
          flagt_p(pn)=1
          dtr_p(pn)=TM_dt
! *** keep the data of density and pressure, when star forms ***
        else if(i.lt.nagravg) then
          list(nagrav)=pn
          nagrav=nagrav+1
        else if(i.lt.nag) then
          talist(na)=pn
          na=na+1
        else 
          nalist(nna)=pn
          nna=nna+1
        endif
      enddo

      if(dns.gt.0) then
! * change list_ap() *
        ng=na+nna+nagrav
        nag=na+nagrav
        nagravg=nagrav
        ns=ns+dns
        nas=nas+dns
        j=0
        do i=0,nagrav-1
          list_ap(j)=list(i)
          j=j+1
        enddo
        do i=0,na-1
          list_ap(j)=talist(i)
          j=j+1
        enddo
        do i=0,nna-1
          list_ap(j)=nalist(i)
          j=j+1
        enddo
        do i=0,dns-1
          list_ap(j)=slist(i)
          j=j+1
        enddo
      endif

      deallocate(slist)
      deallocate(list)
      deallocate(talist)
      deallocate(nalist)

! *** for total number of star particles ***
      allocate(tivs(0:1))
      allocate(tivr(0:1))
      if(nprocs.gt.1) then
        tivs(0)=ns
        tivs(1)=dns
        call MPI_ALLREDUCE(tivs,tivr,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      else
        tivr(1)=dns
      endif
      SI_tnsf=SI_tnsf+tivr(1)
      deallocate(tivs)
      deallocate(tivr)

end subroutine
#endif
