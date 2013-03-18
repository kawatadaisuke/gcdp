#include "gcdp.def"
! ****************************************************
!    setsldmtr.F95 for GCD+ ver. f03.0
!  4  Jan., 2013    written by D.KAWATA
! ****************************************************

#ifdef DM
subroutine setsldmtr(nssl,nalist,ndm)
      use gcdp_const
      use gcdp_dm
      use gcdp_dmtree
#ifdef MULTI_LEV
      use gcdp_ldmtree
#endif
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: nssl,nalist(0:nssl),ndm
      integer i
! * Number of Notfinished particle, temp *      
      integer nlist,tnlist
! * Particle Number in this node *      
      integer pni,pn,nd
      double precision rc2
! * for calculate div(), rot() *
      double precision xij,yij,zij,tij,crit,crit2
! for work
     integer,allocatable :: list(:),node(:)

#ifdef SIM1D
      crit=0.5d0
#elif defined(SIM2D)
      crit=dsqrt(2.0d0)*0.5d0
#else
      crit=dsqrt(3.0d0)*0.5d0
#endif
      crit2=2.0d0*crit

! *** tree walk for all the proc ***
#ifdef MULTI_LEV
      if(np_ldmtr(0).le.1.or.proc_ldmtr(0).ne.myrank) then
#else
      if(np_dmtr(0).le.1.or.proc_dmtr(0).ne.myrank) then
#endif
        do i=0,nssl-1
          pni=nalist(i)
          h_dm(pni)=2.0d0*h_dm(pni)
        enddo
        return
      endif

      allocate(list(0:nssl))
      allocate(node(0:ndm))

      do i=0,nssl-1
        pni=nalist(i)
        list(i)=pni
        node(pni)=0
! *** set largest one ***
#ifdef MULTI_LEV
        h_dm(pni)=(crit2*l_ldmtr(nodess_ldmtr))**2
#else
        h_dm(pni)=(crit2*l_dmtr(nodess_dmtr))**2
#endif
      enddo 
      nlist=nssl
   77 if(nlist.eq.0) then

        deallocate(list)
        deallocate(node)

#ifdef MULTI_LEV
        goto 98
#else
        goto 99
#endif
      endif
      do i=0,nlist-1
        pni = list(i)
        nd = node(pni)         
        pn = pn_dmtr(nd)
        if(np_dmtr(nd).gt.1) then
          xij=x_dm(pni)-cx_dmtr(nd)
#ifdef BCX_PERIOD   
          tij=xij-SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
          tij=xij+SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
#endif
#ifdef SIM1D
          yij=0.0d0
          zij=0.0d0
#else
          yij=y_dm(pni)-cy_dmtr(nd)
#ifdef BCY_PERIOD   
          tij=yij-SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
          tij=yij+SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
#endif
#ifdef SIM2D
          zij=0.0d0
#else
          zij=z_dm(pni)-cz_dmtr(nd)
#ifdef BCZ_PERIOD   
          tij=zij-SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
          tij=zij+SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
#endif
#endif
#endif
          rc2=xij*xij+yij*yij+zij*zij
          if(rc2.lt.h_dm(pni) &
           .or.dsqrt(rc2)-crit*l_dmtr(nd).lt.dsqrt(h_dm(pni))) then
! *** possible to find a particle closer ***
            node(pni)=daughter_dmtr(nd)
          else
            node(pni)=next_dmtr(nd)
          endif
        else if(pn.ne.pni) then
! *** need to recalculate rij, because cx_dmtr is the center of the node. ***
          xij=x_dm(pni)-x_dm(pn)
#ifdef BCX_PERIOD   
          tij=xij-SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
          tij=xij+SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
#endif
#ifdef SIM1D
          yij=0.0d0
          zij=0.0d0
#else
          yij=y_dm(pni)-y_dm(pn)
#ifdef BCY_PERIOD   
          tij=yij-SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
          tij=yij+SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
#endif
#ifdef SIM2D
          zij=0.0d0
#else
          zij=z_dm(pni)-z_dm(pn)
#ifdef BCZ_PERIOD   
          tij=zij-SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
          tij=zij+SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
#endif
#endif
#endif
          rc2=xij*xij+yij*yij+zij*zij
          if(rc2.lt.h_dm(pni)) then
            h_dm(pni)=rc2
          endif
          node(pni)=next_dmtr(nd)
        else
          node(pni)=next_dmtr(nd)
        endif
      enddo
! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      do i=0,tnlist-1
        if(node(list(i)).ne.0) then
          list(nlist)=list(i)
          nlist=nlist+1
        endif
      enddo
!        write(6,*) ' tnlist,nlist=',nlist,tnlist
      goto 77
#ifdef MULTI_LEV
! *** for ldmtr ***
   98 do i=0,nssl-1
        pni=nalist(i)
        list(i)=pni
        node(pni)=0
      enddo 
      nlist=nssl
   78 if(nlist.eq.0) then
        goto 99
      endif
      do i=0,nlist-1
        pni = list(i)
        nd = node(pni)         
        pn = pn_ldmtr(nd)
        if(np_ldmtr(nd).gt.1) then
          xij=x_dm(pni)-cx_ldmtr(nd)
#ifdef BCX_PERIOD   
          tij=xij-SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
          tij=xij+SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
#endif
#ifdef SIM1D
          yij=0.0d0
          zij=0.0d0
#else
          yij=y_dm(pni)-cy_ldmtr(nd)
#ifdef BCY_PERIOD   
          tij=yij-SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
          tij=yij+SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
#endif
#ifdef SIM2D
          zij=0.0d0
#else
          zij=z_dm(pni)-cz_ldmtr(nd)
#ifdef BCZ_PERIOD   
          tij=zij-SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
          tij=zij+SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
#endif
#endif
#endif
        else if(pn.ne.pni) then
! *** need to recalculate rij, because cx_dmtr is the center of the node. ***
          xij=x_dm(pni)-x_dm(pn)
#ifdef BCX_PERIOD   
          tij=xij-SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
          tij=xij+SI_lbox(0)
          if(dabs(tij).lt.dabs(xij)) then
            xij=tij
          endif
#endif
#ifdef SIM1D
          yij=0.0d0
          zij=0.0d0
#else
          yij=y_dm(pni)-y_dm(pn)
#ifdef BCY_PERIOD   
          tij=yij-SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
          tij=yij+SI_lbox(1)
          if(dabs(tij).lt.dabs(yij)) then
            yij=tij
          endif
#endif
#ifdef SIM2D
          zij=0.0d0
#else
          zij=z_dm(pni)-z_dm(pn)
#ifdef BCZ_PERIOD   
          tij=zij-SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
          tij=zij+SI_lbox(2)
          if(dabs(tij).lt.dabs(zij)) then
            zij=tij
          endif
#endif
#endif
#endif
        endif
        rc2=xij*xij+yij*yij+zij*zij
        if(np_ldmtr(nd).eq.1) then
          if(rc2.lt.h_dm(pni)) then
            h_dm(pni)=rc2
          endif
          node(pni)=next_ldmtr(nd)
        else
          if(rc2.lt.h_dm(pni) &
           .or.dsqrt(rc2)-crit*l_ldmtr(nd).lt.dsqrt(h_dm(pni))) then
! *** possible to find a particle closer ***
            node(pni)=daughter_ldmtr(nd)
          else
            node(pni)=next_ldmtr(nd)
          endif
        endif
      enddo
! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      do i=0,tnlist-1
        if(node(list(i)).ne.0) then
          list(nlist)=list(i)
          nlist=nlist+1
        endif
      enddo
!        write(6,*) ' tnlist,nlist=',nlist,tnlist
      goto 78
#endif
! *** for ldmtr ***
! *** end itteration within the proc ***
   99 do i=0,nssl-1
        pni=nalist(i)
        if(h_dm(pni).gt.0.0d0) then
!          write(6,*) ' setsldmtr(): id,r_close=',pni,dsqrt(h_dm(pni))
          h_dm(pni)=ETAHDM*dsqrt(h_dm(pni))
        else
          write(6,*) ' Error in setsltr(): h is zero'
          write(6,*) ' nssl,id,rho=',nssl,pni,h_dm(pni)
          stop
        endif
      enddo

end subroutine
#endif
