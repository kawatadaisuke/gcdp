#include "gcdp.def"
! ****************************************************
!    setslgtr.F95 for GCD+ ver. f03.0
!  04  Jan. 2013    written by D.KAWATA
! ****************************************************
! flag: 0 for gas, otherwise for feedback

#ifdef GAS
subroutine setslgtr(nssl,nalist,npb,flag)
      use gcdp_const
      use gcdp_baryon
      use gcdp_btree
      use gcdp_gtree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: nssl,nalist(0:nssl),flag,npb
      integer i
! * Number of Notfinished particle, temp *      
      integer nlist,tnlist
! * Particle Number in this node *      
      integer pni,pn,nd
      double precision rc2,etasl
! * for calculate div(), rot() *
      double precision xij,yij,zij,tij,crit,crit2
! for work
      integer,allocatable :: list(:),node(:)
      double precision,allocatable :: tdvr(:)

#ifdef SIM1D
      crit=0.5d0
#elif defined(SIM2D)
      crit=dsqrt(2.0d0)*0.5d0
#else
      crit=dsqrt(3.0d0)*0.5d0
#endif
      crit2=2.0d0*crit

      if(flag.eq.0) then

        allocate(tdvr(0:nssl))

        etasl=ETAH
        do i=0,nssl-1
          pni=nalist(i)
          tdvr(pni)=h_p(pni)
        enddo
      else
        write(6,*) ' Error in setslgtr(): flag is not 0. myrank=',myrank
        write(6,*) '  but Feedback is off in gcdp.def.'
        stop
      endif

! *** tree walk for all the proc ***
      if(np_gtr(0).le.1.or.proc_gtr(0).ne.myrank) then
        do i=0,nssl-1
          pni=nalist(i)
!          if(tdvr(pni).lt.l_tr(0)) then
! *** use the box size of the node ***
!            tdvr(pni)=MGROOT*l_tr(0)
!          else
          tdvr(pni)=2.0d0*tdvr(pni)
!          endif
        enddo
        goto 999
      endif

      allocate(list(0:nssl))
      allocate(node(0:npb))

      do i=0,nssl-1
        pni=nalist(i)
        list(i)=pni
        node(pni)=nodess_gtr
! *** set 2nd largest one ***
        tdvr(pni)=(crit2*l_gtr(nodess_gtr))**2
      enddo 
      nlist=nssl
   77 if(nlist.eq.0) then
  
        deallocate(list)
        deallocate(node)

        goto 99
      endif

      do i=0,nlist-1
        pni=list(i)
        nd=node(pni)         
        pn=pn_gtr(nd)
        if(np_gtr(nd).gt.1) then
          xij=x_p(pni)-cx_gtr(nd)
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
          yij=y_p(pni)-cy_gtr(nd)
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
          zij=z_p(pni)-cz_gtr(nd)
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
          if(rc2.lt.tdvr(pni) &
           .or.dsqrt(rc2)-crit*l_gtr(nd).lt.dsqrt(tdvr(pni))) then
! *** possible to find a particle closer ***
! *** check if or not it's pseudo node
            if(proc_gtr(nd).ne.myrank) then
              node(pni)=next_gtr(nd)
            else
              node(pni)=daughter_gtr(nd)
            endif
          else
            node(pni)=next_gtr(nd)
          endif
        else if(pn.ne.pni) then
! *** need to recalculate rij, because cx_tr is the center of the node. ***
          xij=x_p(pni)-x_p(pn)
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
          yij=y_p(pni)-y_p(pn)
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
          zij=z_p(pni)-z_p(pn)
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
          if(rc2.gt.0.0d0.and.rc2.lt.tdvr(pni)) then
            tdvr(pni)=rc2
          endif
          node(pni)=next_gtr(nd)
        else
          node(pni)=next_gtr(nd)
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
! *** end itteration within the proc ***
   99 do i=0,nssl-1
        pni=nalist(i)
        if(tdvr(pni).gt.0.0d0) then
!          write(6,*) ' setslgtr(): id,r_close=',pni,dsqrt(tdvr(pni))
          tdvr(pni)=etasl*dsqrt(tdvr(pni))
        else
          write(6,*) ' Error in setslgtr(): h is zero'
          write(6,*) ' myrank,i,nssl,id,rho=',myrank,i,nssl,pni,tdvr(pni)
          write(6,*) ' l and n gtr0=',l_gtr(0),np_gtr(0)
          write(6,*) ' x=',x_p(pni),y_p(pni),z_p(pni),id_p(pni)
          write(6,*) ' vx=',vnx_p(pni),vny_p(pni),vnz_p(pni)
          if(i.lt.nssl-1) then 
            pni=pni+1
          else
            pni=0
          endif
          write(6,*) ' x=',x_p(pni),y_p(pni),z_p(pni),id_p(pni)
          write(6,*) ' vx=',vnx_p(pni),vny_p(pni),vnz_p(pni)
          stop
        endif
      enddo
  999 if(flag.eq.0) then
        do i=0,nssl-1
          pni=nalist(i)
          h_p(pni)=tdvr(pni)
        enddo
      endif

      deallocate(tdvr)

end subroutine
#endif
