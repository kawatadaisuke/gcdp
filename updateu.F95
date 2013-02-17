#include "gcdp.def"
! *****************************************************
!    updateu.F for GCD+ ver f03.0
!  27  Jan.  2013   written by D.Kawata
! ***************************************************** 
! *********************************************
!    Definition of Function for update data
!       updateu() & vertualp_cs()
! *********************************************

#ifdef GAS
! **********   Definition of updateu()   **********
subroutine updateu(np,ng)
      use gcdp_const
      use gcdp_baryon
#ifdef COOL
      use gcdp_system
#endif

      implicit none

      integer,intent(in) :: np,ng
      integer i,j,pn
#ifdef COOL
      integer nlist
      integer,allocatable :: plist(:),eslist(:)
      double precision zn
#endif
      character fileo*60

#ifdef COOL
      allocate(plist(0:ng))
      nlist = 0
      do i=0,ng-1
        pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).le.0) then
#endif
          plist(nlist)=pn
          nlist = nlist+1
#if defined(SF_EFD) || defined(SF_ZFD)
        endif
#endif
      enddo
      zn = 1.0d0/SI_a-1.0d0
! *** update p_p ***
      allocate(eslist(0:nlist))
      do i=0,nlist-1
        eslist(i)=plist(i)
      enddo
      call energy_solver(np,nlist,eslist,zn,0)
! *** update u_p ***
      do i=0,nlist-1
        pn=eslist(i)
        u_p(pn) = p_p(pn)/((GAM-1.0d0)*rho_p(pn))
      enddo
      deallocate(plist)
      deallocate(eslist)

#else
      do i=0,ng-1
        pn = list_ap(i)
! *** only active particles since pv30.5
#ifdef SF_EFD
          as_p(pn)=as_p(pn)+0.5d0*dt_p(pn)*(ndu_p(pn)+pdu_p(pn)) &
           +((GAM-1.0d0)/(rho_p(pn)**(GAM-1.0d0)))*Gsn_p(pn)
#else
          as_p(pn)=as_p(pn)+0.5d0*dt_p(pn)*(ndu_p(pn)+pdu_p(pn))
#endif
          u_p(pn)=(as_p(pn)/(GAM-1.0d0))*(rho_p(pn)**(GAM-1.0d0))
        enddo
#endif

end subroutine

! **********   Definition of virtualp_cs()   **********
subroutine virtualp_cs(np,ng)
      use gcdp_const
      use gcdp_baryon
#ifdef COOL
      use gcdp_system
#endif
      implicit none
      integer,intent(in) :: np,ng
      integer i,j,pn
#ifdef COOL
      integer nlist
      integer,allocatable :: palist(:),eslist(:)
      double precision zn
#endif

#ifdef COOL
      allocate(palist(0:ng))
      nlist = 0
      do i=0,ng-1
        pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).le.0) then
#endif
          palist(nlist)=pn
          nlist = nlist+1
#if defined(SF_EFD) || defined(SF_ZFD)
        endif
#endif
      enddo
      zn = 1.0d0/SI_a-1.0d0
! *** update p_p ***
      allocate(eslist(0:nlist))
      do i=0,nlist-1
        eslist(i)=palist(i)
      enddo
      call energy_solver(nlist,eslist,zn,1)
      deallocate(palist)
      deallocate(eslist)
#else
! *** no cooling ***
      do i=0,ng-1
        pn=list_ap(i)
! *** only for active particles since pv30.5 ***
        p_p(pn)=(rho_p(pn)**GAM)*(as_p(pn)+dt_p(pn)*ndu_p(pn) &
#ifdef SF_EFD
         +((GAM-1.0d0)/(rho_p(pn)**(GAM-1.0d0)))*Gsn_p(pn))
#else
         )
#endif
      enddo
#endif
! *** calculate cs ***
      do i=0,ng-1
        pn=list_ap(i)         
        cs_p(pn)=dsqrt(GAM*p_p(pn)/rho_p(pn))
      enddo
end subroutine
! endif for GAS
#endif

