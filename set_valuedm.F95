#include "gcdp.def"
! ****************************************************
!    set_value.F95 for gcd+ ver. f03.0
!  04  Jan. 2013    written by D.KAWATA
! ****************************************************

! calculate rho_dm,h_dm,omgh_dm,zetah_dm,minhv_dm
!
! if flag.ne.0, use vnx_dm for the veloicty at the step
! (vx_dm is not good, after ddecdm in some case.)

#ifdef DM
subroutine set_valuedm(ndm,nadm,flag)      
      use gcdp_const
      use gcdp_dm
      use gcdp_dmtree
#ifdef MULTI_LEV
      use gcdp_ldmtree
#endif
      use gcdp_system
      use gcdp_kernel

      implicit none
      include 'mpif.h'

      integer,intent(in) :: ndm,nadm,flag
#ifdef FIXEDDMEPS
      integer i,pni
      do i=0,ndm-1
        pni=list_adm(i)   
        nnb_dm(pni)=1
        div_v_dm(pni)=0.0d0
        omgh_dm(pni)=1.0d0
        zetah_dm(pni)=0.0d0
        h_dm(pni)=SI_eps*(m_dm(pni)**THIRD)
      enddo        
#else
      integer i,j,nit,mnit,mnitb,npj,srank,rrank,ip,nval,isend &
       ,irecv,nival,is,iproc,npj0,snval,rnval,snival,npjs,snvalmax
      integer flagnb,nssl
      integer ntw,maxntw,ierr
! * Number of Notfinished particle, temp *      
      integer nlist,tnlist
! * Particle Number in this node *      
      integer pni,pn,nd
      double precision r,s,rc,rh
! * Smoothing length *      
      double precision hsi,ihsi,wij
! * for calculate div(), rot() *
      double precision dwr,dwh,dphih
      double precision xij,yij,zij,tij
      double precision vxij,vyij,vzij,vr
      double precision crit,ipi
! *** information for particles need communication ***
      integer ncomp,ncompt,ndivcom,idivcom
      integer ncomptmax,npprocmax,ndivcompt
      integer ncompi
      integer pncompt(0:MNDM-1),proccompt(0:MNDM-1)
! *** for pseudo node 
      integer ndp,pnodess
! for test output
      character fileo*60
! *** flag for tree ***
      integer flagtr
! pn_nfp can be bigger than ng+ns in this proc, when receiving the particles
      integer pn_nfp(0:MNDM-1)
! for MPI
      integer,allocatable :: npproc(:),isproc(:),ieproc(:) &
       ,cproc(:),idisp(:),jjlen(:),ireqs(:),npjr(:)
! for work,allocatable
      integer,allocatable :: flagi(:)
      integer,allocatable :: pncomp(:),flagproc(:),flagcom(:)
! *** for pseudo node 
      integer,allocatable :: nextpnd(:),daupnd(:),ndpnd(:),idpnd(:)
      double precision,allocatable :: xpnd(:),ypnd(:),zpnd(:),lpnd(:)
      double precision,allocatable :: trbuf(:)
      integer,allocatable :: tibuf(:)
! *** for work ***
      integer,allocatable :: nd_nfp(:),list(:),node(:),tivr(:) &
        ,tivs(:),talist(:),nalist(:),sldmtrlist(:)
      double precision,allocatable :: tx(:),ty(:),tz(:),du0(:),tdvr(:) &
        ,tdvs(:),fhx(:),fx(:),hu(:),lu(:),av(:),adv(:)
#if defined(COSM) && defined(BOXSIM)
      double precision,allocatable :: int1oa(:),t0(:),t1(:)
#endif


!      write(fileo,'(a7,i3.3)') 'psetvdm',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ndm-1
!        pni=list_adm(i)
!        write(60,'(6(1pE13.5))') x_dm(pni),y_dm(pni),z_dm(pni) &
!        ,vx_dm(pni),vy_dm(pni),vz_dm(pni)
!      enddo
!      close(60)

! *** allocate work space for MPI ***
      allocate(idisp(0:nprocs))
      allocate(jjlen(0:nprocs))
      allocate(ireqs(0:nprocs))
      allocate(npproc(0:nprocs))
      allocate(isproc(0:nprocs))
      allocate(ieproc(0:nprocs))
      allocate(cproc(0:nprocs))
      allocate(npjr(0:nprocs))

! **** Basic Value ***
      ipi = 1.0d0/M_PI
#ifdef SIM1D
      crit=0.5d0
#elif defined(SIM2D)
      crit=dsqrt(2.0d0)*0.5d0
#else
      crit=dsqrt(3.0d0)*0.5d0
#endif
      nval = 11
      nival = 1
! *** initialization ***
      snval=1
      snival=1

! *** set tx,ty,tz for v* at the current time step ***

      allocate(tx(0:ndm))
      allocate(ty(0:ndm))
      allocate(tz(0:ndm))

      if(flag.eq.0) then

#if defined(COSM) && defined(BOXSIM)
! time at n+1/3
        allocate(t0(0:ndm))
        allocate(t1(0:ndm))
        allocate(int1oa(0:ndm))
! assign time range 
        do i=0,ndm-1
          t0(i)=tc_dm(i)
          t1(i)=tc_dm(i)+vdt_dm(i)
        enddo
! integrate 1/a from n to n+1/2 for each particle)
        call integrate1oa(ndm,t0,t1,int1oa)
#endif

        do i=0,ndm-1
! velocity at this step
#if defined(COSM) && defined(BOXSIM)
! store a*vpec only, add hubble flow for vxij
          tx(i)=(vnx_dm(i)+(asc_dm(i)**2)*dvx_dm(i)*int1oa(i))/SI_a
          ty(i)=(vny_dm(i)+(asc_dm(i)**2)*dvy_dm(i)*int1oa(i))/SI_a
          tz(i)=(vnz_dm(i)+(asc_dm(i)**2)*dvz_dm(i)*int1oa(i))/SI_a
#else
          tx(i)=vnx_dm(i)+vdt_dm(i)*dvx_dm(i)
          ty(i)=vny_dm(i)+vdt_dm(i)*dvy_dm(i)
          tz(i)=vnz_dm(i)+vdt_dm(i)*dvz_dm(i)
#endif
        enddo 

#if defined(COSM) && defined(BOXSIM)
        deallocate(t0)
        deallocate(t1)
        deallocate(int1oa)
#endif

      else
        do i=0,ndm-1
! physical velocity at this time step
#if defined(COSM) && defined(BOXSIM)
          tx(i)=SI_hub*x_dm(i)+vnx_dm(i)/SI_a
          ty(i)=SI_hub*y_dm(i)+vny_dm(i)/SI_a
          tz(i)=SI_hub*z_dm(i)+vnz_dm(i)/SI_a
#else
          tx(i)=vnx_dm(i)
          ty(i)=vny_dm(i)
          tz(i)=vnz_dm(i)
#endif
        enddo
      endif

! * Initialization *
      maxntw=ndm*3+nprocs*MAXNODESEND
! *** maximum iteration to switch to bi-section ***
      mnitb = 10
! *** maximum iteration ***
      mnit = 100
      nit = 0
! *** for communication ***
! *** only for active particles ***
      npj=nadm
! *** initialisation ***
      idivcom=0
      ndivcom=-1
      do i=0,nprocs-1
        npproc(i)=0
      enddo

! allocate work space for all dm paricles
      allocate(flagi(0:ndm))
      allocate(fhx(0:ndm))
      allocate(fx(0:ndm))
      allocate(hu(0:ndm))
      allocate(lu(0:ndm))
      allocate(du0(0:ndm))
      allocate(av(0:ndm))
      allocate(adv(0:ndm))
      allocate(nd_nfp(0:ndm))

      do i=0,npj-1
        pn_nfp(i) = list_adm(i)
      enddo
      do i=0,npj-1
        pn=pn_nfp(i)
! *** flag for iteration 1: ready for bi-section ***
        flagi(pn)=0
        fhx(pn)=INF
        fx(pn)=-INF
        hu(pn)=-INF
        lu(pn)=-INF
      enddo

! *** initialization ***
   71 do i=0,npj-1
        pni=pn_nfp(i)
        rho_dm(pni)=0.0d0
        omgh_dm(pni)=0.0d0
        zetah_dm(pni)=0.0d0
        div_v_dm(pni)=0.0d0
        nnb_dm(pni)=0
      enddo 

      flagtr=1
   72 allocate(tdvr(0:npj*nval))
      allocate(tivr(0:npj*nival))
      allocate(flagproc(0:npj))
      allocate(flagcom(0:npj))
      allocate(pncomp(0:npj))
      allocate(talist(0:npj))

      do i=0,npj-1
        pni=pn_nfp(i)
        tdvr(i)=x_dm(pni)
        tdvr(i+npj)=y_dm(pni)
        tdvr(i+npj*2)=z_dm(pni)
        tdvr(i+npj*3)=tx(pni)
        tdvr(i+npj*4)=ty(pni)
        tdvr(i+npj*5)=tz(pni)
        tdvr(i+npj*6)=h_dm(pni)
        talist(i)=pni
        flagproc(i)=myrank
        flagcom(i)=-1
        tivr(i)=0
! store original h
        du0(pni)=h_dm(pni)
      enddo
! rho_dm(8),zetah(9),div(10),omgh(11)
      do i=npj*7,npj*nval-1
        tdvr(i)=0.0d0
      enddo 
      ncomp=0
      ncompt=0
      iproc=0
      idivcom=0
      ndivcom=-1
      do i=0,nprocs-1
        npproc(i)=0
      enddo

      if(flagtr.eq.1) then
        goto 91
#ifdef MULTI_LEV
      else if(flagtr.eq.2) then
        goto 92
#endif
      else
        goto 94
      endif

   91 if(np_dmtr(0).eq.0) then
        goto 99
      endif
      if(iproc.ne.0.and.proc_dmtr(0).ne.myrank) then   
        goto 99
      endif

      allocate(list(0:npj))
      allocate(node(0:npj))

      do i=0,npj-1
        list(i)=i
! *** node has to start from 0
        node(i)=0
      enddo 
      nlist=npj
! **********    start tree walk   *********
      ntw=0
   77 if(nlist.eq.0) then

       deallocate(list)
       deallocate(node)

        goto 99
      endif
      ntw=ntw+1
      if(ntw.gt.maxntw) then
        write(6,*) ' Error set_valuedm(): ntw>',maxntw
        write(6,*) ' myrank,iproc,npj,nlist=',myrank,iproc,npj,nlist
        write(fileo,'(a7,i3.3)') 'errsvdm',myrank
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pni=list(i)
          write(60,'(4(1pE13.5),3I10)') tdvr(pni),tdvr(pni+npj) &
           ,tdvr(pni+npj*2),tdvr(pni+npj*6),pni,node(pni)
        enddo                  
        close(60)
        stop
      endif
      do i=0,nlist-1
        pni = list(i)
        nd = node(pni)         
        xij=tdvr(pni)-cx_dmtr(nd)
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
        yij=tdvr(pni+npj)-cy_dmtr(nd)
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
        zij=tdvr(pni+npj*2)-cz_dmtr(nd)
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
        rc=xij*xij+yij*yij+zij*zij
! *** only take into account gather condition ***
        rh=tdvr(pni+npj*6)+crit*l_dmtr(nd)
        rh=rh*rh
        if(np_dmtr(nd).eq.1) then
        if(rc.lt.rh) then
        if(proc_dmtr(nd).eq.myrank) then
          pn = pn_dmtr(nd)
          r=rc
          if(r.ne.0.0d0) then
            r=dsqrt(r)
          endif            
          hsi=tdvr(pni+npj*6)
          ihsi = 1.0d0/hsi
          s = r*ihsi
          if(s.le.1.0d0) then
! * satisfy the criterion for neigbour particles *
            is=int(s*dnktab)
            if(is.lt.0) then
              is=0
            else if(is.ge.NKTAB) then
              is=NKTAB-1
            endif
            wij=w_tb(is)+(w_tb(is+1)-w_tb(is))*(s-s_tb(is))*dnktab
! *** only 3D ***
            wij=wij*(ihsi**3)
! *** update density rho ***
            tdvr(pni+npj*7)=tdvr(pni+npj*7)+m_dm(pn)*wij
! *** dW/ds/r/h calculation of omega ***
            dwr=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
             *(s-s_tb(is))*dnktab
            dwr=dwr*(ihsi**5)
! *** calculating omgh_dm ***
! *** dW/dh ***
            dwh=-s*r*dwr-3.0d0*wij*ihsi
! *** for omgh_dm ***
            tdvr(pni+npj*8)=tdvr(pni+npj*8)+m_dm(pn)*dwh
#ifdef SGRAV
! *** calculate zetah_dm ***
! *** dphi/dh ***
            dphih=dphidh_tb(is)+(dphidh_tb(is+1)-dphidh_tb(is)) &
             *(s-s_tb(is))*dnktab
            dphih=dphih*(ihsi**2)
! *** for zetah_dm ***
            tdvr(pni+npj*10)=tdvr(pni+npj*10)+m_dm(pn)*dphih
#endif
            if(iproc.ne.0.or.talist(pni).ne.pn) then
! *** count neighbour particles 
              tivr(pni)=tivr(pni)+1
! *** velocity difference, only 3D ***
              vxij=tdvr(pni+npj*3)-tx(pn)
              vyij=tdvr(pni+npj*4)-ty(pn)
              vzij=tdvr(pni+npj*5)-tz(pn)
#if defined(COSM) && defined(BOXSIM)
              vxij=vxij+SI_hub*xij
              vyij=vyij+SI_hub*yij
              vzij=vzij+SI_hub*zij
#endif
! * div(v) *
              vr=(vxij*xij+vyij*yij+vzij*zij)
              tdvr(pni+npj*9)=tdvr(pni+npj*9)+m_dm(pn)*dwr*vr
            endif
          endif
        endif
        endif
        endif
! * update node *		
        if((np_dmtr(nd).eq.1.and.proc_dmtr(nd).eq.myrank).or.rc.ge.rh) then
          node(pni)=next_dmtr(nd)
        else
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_dmtr(nd).eq.-1) then
            if(proc_dmtr(nd).ne.myrank.and.flagproc(pni).ne.proc_dmtr(nd)) then
! proc_dmtr will be the same for different domain within the same proc

              flagproc(pni)=proc_dmtr(nd)

! *** store pn ***
              if(flagcom(pni).lt.0) then
                flagcom(pni)=proc_dmtr(nd)
                pncomp(ncomp)=pn_nfp(pni)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNDM) then
                pncompt(ncompt)=pn_nfp(pni)
                proccompt(ncompt)=proc_dmtr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_dmtr(nd))=npproc(proc_dmtr(nd))+1
            endif
            node(pni)=next_dmtr(nd)
          else
! *** information is in the other node ***
            node(pni)=daughter_dmtr(nd)
          endif
        endif          
      enddo
! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      if(iproc.eq.0.and.next_dmtr(0).ne.0) then
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          else if(node(list(i)).eq.0) then
            node(list(i))=next_dmtr(0)
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      else
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      endif
      goto 77
! *** end itteration within the proc ***

#ifdef MULTI_LEV
   92 if(np_ldmtr(0).eq.0) then
        goto 99
      endif
      if(iproc.ne.0.and.proc_ldmtr(0).ne.myrank) then   
        goto 99
      endif

      allocate(list(0:npj))
      allocate(node(0:npj))

      do i=0,npj-1
        list(i)=i
! *** node has to start from 0
        node(i)=0
      enddo 
      nlist=npj
! **********    start tree walk   *********
      ntw=0
   78 if(nlist.eq.0) then

        deallocate(list)
        deallocate(node)

        goto 99
      endif
      ntw=ntw+1
      if(ntw.gt.maxntw) then
        write(6,*) ' Error set_valuedm(): in ldmtr ntw>',maxntw
        write(6,*) ' myrank,iproc,npj,nlist=',myrank,iproc,npj,nlist
        write(fileo,'(a7,i3.3)') 'errsvdm',myrank
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pni=list(i)
          write(60,'(4(1pE13.5),3I10)') tdvr(pni),tdvr(pni+npj) &
           ,tdvr(pni+npj*2),tdvr(pni+npj*6),pni,node(pni)
        enddo                  
        close(60)
        call MPI_FINALIZE(ierr)
        stop
      endif
      do i=0,nlist-1
        pni = list(i)
        nd = node(pni)         
        xij=tdvr(pni)-cx_ldmtr(nd)
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
        yij=tdvr(pni+npj)-cy_ldmtr(nd)
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
        zij=tdvr(pni+npj*2)-cz_ldmtr(nd)
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
        rc=xij*xij+yij*yij+zij*zij
! *** only take into account gather condition ***
        rh=tdvr(pni+npj*6)+crit*l_ldmtr(nd)
        rh=rh*rh
        if(np_ldmtr(nd).eq.1) then
        if(rc.lt.rh) then
        if(proc_ldmtr(nd).eq.myrank) then
          pn = pn_ldmtr(nd)
          r=rc
          if(r.ne.0.0d0) then
            r=dsqrt(r)
          endif            
          hsi=tdvr(pni+npj*6)
          ihsi = 1.0d0/hsi
          s = r*ihsi
          if(s.le.1.0d0) then
! * satisfy the criterion for neigbour particles *
            is=int(s*dnktab)
            if(is.lt.0) then
              is=0
            else if(is.ge.NKTAB) then
              is=NKTAB-1
            endif
            wij=w_tb(is)+(w_tb(is+1)-w_tb(is))*(s-s_tb(is))*dnktab
! *** only 3D ***
            wij=wij*(ihsi**3)
! *** update density rho ***
            tdvr(pni+npj*7)=tdvr(pni+npj*7)+m_dm(pn)*wij
! *** dW/ds/r/h calculation of omega ***
            dwr=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
             *(s-s_tb(is))*dnktab
            dwr=dwr*(ihsi**5)
! *** calculating omgh_dm ***
! *** dW/dh ***
            dwh=-s*r*dwr-3.0d0*wij*ihsi
! *** for omgh_dm ***
            tdvr(pni+npj*8)=tdvr(pni+npj*8)+m_dm(pn)*dwh
#ifdef SGRAV
! *** calculate zetah_dm ***
! *** dphi/dh ***
            dphih=dphidh_tb(is)+(dphidh_tb(is+1)-dphidh_tb(is)) &
             *(s-s_tb(is))*dnktab
            dphih=dphih*(ihsi**2)
! *** for zetah_dm ***
            tdvr(pni+npj*10)=tdvr(pni+npj*10)+m_dm(pn)*dphih
#endif
            if(iproc.ne.0.or.talist(pni).ne.pn) then
! *** count neighbour particles 
              tivr(pni)=tivr(pni)+1
! *** velocity difference, only 3D ***
              vxij=tdvr(pni+npj*3)-tx(pn)
              vyij=tdvr(pni+npj*4)-ty(pn)
              vzij=tdvr(pni+npj*5)-tz(pn)
#if defined(COSM) && defined(BOXSIM)
              vxij=vxij+SI_hub*xij
              vyij=vyij+SI_hub*yij
              vzij=vzij+SI_hub*zij
#endif
! * div(v) *
              vr=(vxij*xij+vyij*yij+vzij*zij)
              tdvr(pni+npj*9)=tdvr(pni+npj*9)+m_dm(pn)*dwr*vr
            endif
          endif
        endif
        endif
        endif
! * update node *		
        if((np_ldmtr(nd).eq.1.and.proc_ldmtr(nd).eq.myrank).or.rc.ge.rh) then
          node(pni)=next_ldmtr(nd)
        else
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_ldmtr(nd).eq.-1) then
            if(proc_ldmtr(nd).ne.myrank &
             .and.flagproc(pni).ne.proc_ldmtr(nd)) then
! proc_dmtr will be the same for different domain within the same proc
              flagproc(pni)=proc_ldmtr(nd)
! *** store pn ***
              if(flagcom(pni).lt.0) then
                flagcom(pni)=proc_ldmtr(nd)
                pncomp(ncomp)=pn_nfp(pni)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNDM) then
                pncompt(ncompt)=pn_nfp(pni)
                proccompt(ncompt)=proc_ldmtr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_ldmtr(nd))=npproc(proc_ldmtr(nd))+1
            endif
            node(pni)=next_ldmtr(nd)
          else
! *** information is in the other node ***
            node(pni)=daughter_ldmtr(nd)
          endif
        endif          
      enddo
! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      if(iproc.eq.0.and.next_ldmtr(0).ne.0) then
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          else if(node(list(i)).eq.0) then
            node(list(i))=next_ldmtr(0)
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      else
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      endif
      goto 78
#endif
! *** end tree walk within the proc ***
   99 if(iproc.ne.0) then
! *** setting numbers for sending and receiving
! *** snval -> rnval: received nval
        rnval=snval
! *** new nval to send the results
        snval=4

        allocate(tivs(0:nprocs))

        do i=0,nprocs-1
! *** keep npjr: data received from each proc ***
          tivs(i)=npjr(i)
! *** ireqs(i)=jjlen(i)number of particles sent -> now receiving.
          npjr(i)=ireqs(i)
        enddo
! *** for double precision ***
        isend=0
        do i=0,nprocs-1
          idisp(i)=isend 
! *** tivs is numbr of particles received from each proc
          jjlen(i)=tivs(i)*snval
          isend=isend+jjlen(i)          
        enddo

        allocate(tdvs(0:npj*snval))

! *** setting the sending back data ***
        do i=0,npj-1
          tdvs(snval*i)=tdvr(i+npj*7)
          tdvs(snval*i+1)=tdvr(i+npj*8)
          tdvs(snval*i+2)=tdvr(i+npj*9)
          tdvs(snval*i+3)=tdvr(i+npj*10)
        enddo
! *** sending and reveiving the data ***
        irecv=0
! *** number of particles sent is ncomp
        npjs=npj
! *** ncompi: number of particles originally sent from iproc=0
        npj=ncompi

        deallocate(tdvr)
        allocate(tdvr(0:npj*nval))

        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tdvr(irecv+npj*7)=trbuf(snval*i)
            tdvr(irecv+npj*8)=trbuf(snval*i+1)
            tdvr(irecv+npj*9)=trbuf(snval*i+2)
            tdvr(irecv+npj*10)=trbuf(snval*i+3)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo    

        if(irecv.ne.npj) then
          write(6,*) ' Error in set_valuedm(): irecv,npj=',irecv,npj
          write(6,*) ' when sending back the results.'
          stop
        endif

        deallocate(tdvs)

! *** for integer values ***
        isend=0
        snival=1
        do i=0,nprocs-1
          idisp(i)=isend 
! *** tivs is numbr of particles received from each proc
          jjlen(i)=tivs(i)*snival
          isend=isend+jjlen(i)          
        enddo

        deallocate(tivs)
        allocate(tivs(0:npjs*snival))

! *** setting the sending back data, use npjs ***
        do i=0,npjs-1
          tivs(snival*i)=tivr(i)
        enddo
! *** sending and reveiving the data ***
        irecv=0
! *** ncompi: number of particles originally sent from iproc=0
        npj=ncompi

        deallocate(tivr)
        allocate(tivr(0:npj*nival))

        do ip=0,nprocs-1

          allocate(tibuf(0:npjr(ip)*snival))

          call MPI_SCATTERV(tivs,jjlen,idisp,MPI_INTEGER,tibuf &
           ,npjr(ip)*snival,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tivr(irecv)=tibuf(snival*i)
            irecv=irecv+1
          enddo

          deallocate(tibuf)

        enddo          

        deallocate(tivs)

      endif

      do i=0,npj-1
        pni = pn_nfp(i)
        rho_dm(pni)=rho_dm(pni)+tdvr(i+npj*7)
        nnb_dm(pni)=nnb_dm(pni)+tivr(i)
        div_v_dm(pni)=div_v_dm(pni)+tdvr(i+npj*9)
        omgh_dm(pni)=omgh_dm(pni)+tdvr(i+npj*8)
#ifdef SGRAV
        zetah_dm(pni)=zetah_dm(pni)+tdvr(i+npj*10)
#else
        zetah_dm(pni)=0.0d0
#endif
      enddo   
! *** end of neighbour search ***

      deallocate(tdvr)
      deallocate(tivr)

! *** restore the initail, i.e., all the list of pn_nfp
      if(iproc.ne.0.and.idivcom.gt.ndivcom) then
        do i=0,npj0-1
          pn_nfp(i)=nd_nfp(i)
        enddo   
        npj=npj0
      endif

      if(nprocs.le.1.or.(iproc.ne.0.and.idivcom.gt.ndivcom)) then
#ifdef MULTI_LEV
        flagtr=flagtr+1
        if(flagtr.eq.2) then

          deallocate(flagproc)
          deallocate(flagcom)
          deallocate(pncomp)
          deallocate(talist)

          goto 72
        else
#endif
          goto 94
#ifdef MULTI_LEV
        endif
#endif
      endif

      if(nprocs.gt.1.and.iproc.eq.0) then
        iproc=1
! *** check if need communication ***
        ncomptmax=0
        call MPI_ALLREDUCE(ncompt,ncomptmax,1,MPI_INTEGER &
         ,MPI_MAX,MPI_COMM_WORLD,ierr)
        if(ncomptmax.gt.0) then
! *** get how many particles each proc receives

          allocate(tivr(0:nprocs))

          do i=0,nprocs-1
            tivr(i)=0
          enddo
          call MPI_ALLREDUCE(npproc,tivr,nprocs,MPI_INTEGER &
           ,MPI_SUM,MPI_COMM_WORLD,ierr)
          npprocmax=0
          do i=0,nprocs-1
            if(tivr(i).gt.npprocmax) then
              npprocmax=tivr(i)
            endif
          enddo

          deallocate(tivr)

          if(npprocmax.gt.MNDM) then
! *** number of times receiving and calculating for the other nodes
            ndivcom=int((npprocmax+nprocs*nprocs)/MNDM)+1
          else
            ndivcom=1
          endif
          ndivcompt=0
          if(ncomptmax.gt.MNDM) then
            ndivcompt=int((ncomptmax+nprocs*nprocs)/MNDM)+1
            if(ndivcompt.gt.ndivcom) then
              ndivcom=ndivcompt
            endif
          endif
          idivcom=1
! *** store original pn_nfp in nd_nfp ***
          npj0=npj
          do i=0,npj0-1
            nd_nfp(i)=pn_nfp(i)
          enddo
        else
! *** no communication required ***
#ifdef MULTI_LEV
          flagtr=flagtr+1
          if(flagtr.eq.2) then

            deallocate(flagproc)
            deallocate(flagcom)
            deallocate(pncomp)
            deallocate(talist)

            goto 72
          else
#endif
            goto 94
#ifdef MULTI_LEV
          endif
#endif
        endif
      endif
! *** do communication ndivcom times
      if(idivcom.le.ndivcom) then
        if(ndivcom.eq.1) then
!        if(2.eq.1) then
! *** can use pncompt
! *** preparation for sending the data to the other procs ***
          isend=0          
! *** store particle list in the order of sending procs in list() ***
          do ip=0,nprocs-1
            idisp(ip)=isend           
            jjlen(ip)=0
            do i=0,ncompt-1
              if(proccompt(i).eq.ip) then
                pn_nfp(isend)=pncompt(i)
                jjlen(ip)=jjlen(ip)+1
                isend=isend+1
              endif
            enddo
! keep original jjlen 
            ireqs(ip)=jjlen(ip)
          enddo
          if(isend.ne.ncompt) then    
            write(6,*) ' Error in set_value(): isend.ne.ncomp'
            write(6,*) ' when counting N particles need communication'
            write(6,*) ' myrank,isend,ncomp=',myrank,isend,ncomp
            stop
          endif
        else
! *** set the range of number of particles for sending for each proc
          do i=0,nprocs-1
            call para_range(0,npproc(i)-1,ndivcom,idivcom-1,isproc(i),ieproc(i))
! initialise for counting particles for each proc
            cproc(i)=0
          enddo
! *** preparation for sending the data to the other procs ***
          isend=0          
          if(flagtr.eq.1) then
            if(proc_dmtr(0).eq.myrank) then
              pnodess=next_dmtr(0)
            else
! in case, if there is no local tree, but pseudo nodes
              pnodess=0
            endif
            ndp=num_dmtr-pnodess

            allocate(idpnd(0:num_dmtr))

#ifdef MULTI_LEV
          else if(flagtr.eq.2) then
            if(proc_ldmtr(0).eq.myrank) then
              pnodess=next_ldmtr(0)
            else
! in case, if there is no local tree, but pseudo nodes
              pnodess=0
            endif
            ndp=num_ldmtr-pnodess

            allocate(idpnd(0:num_ldmtr))

#endif      
          endif

          allocate(ndpnd(0:ndp))
          allocate(xpnd(0:ndp))
          allocate(ypnd(0:ndp))
          allocate(zpnd(0:ndp))
          allocate(lpnd(0:ndp))
          allocate(nextpnd(0:ndp))
          allocate(daupnd(0:ndp))

! *** store particle list in the order of sending procs in list() ***
          do ip=0,nprocs-1
            idisp(ip)=isend           
            jjlen(ip)=0
!            flagerr=0
            if(npproc(ip).gt.0.and.ip.ne.myrank) then
! *** get the coordinate and etc. for pseudo node for the proc
              ndp=0
              if(flagtr.eq.1) then
                do nd=pnodess,num_dmtr-1
                  if(proc_dmtr(nd).eq.ip) then
                    ndpnd(ndp)=nd
                    idpnd(nd)=ndp
                    xpnd(ndp)=cx_dmtr(nd)
                    ypnd(ndp)=cy_dmtr(nd)
                    zpnd(ndp)=cz_dmtr(nd)
                    lpnd(ndp)=l_dmtr(nd)
                    ndp=ndp+1
                  endif
                enddo
                do i=0,ndp-1
                  nd=ndpnd(i)                  
                  if(next_dmtr(nd).gt.pnodess &
                   .and.next_dmtr(nd).le.ndpnd(ndp-1)) then
                    nextpnd(i)=idpnd(next_dmtr(nd))      
                  else
                    nextpnd(i)=ndp
                  endif
                  if(daughter_dmtr(nd).ne.-1) then
                    daupnd(i)=idpnd(daughter_dmtr(nd))
                  else 
                    daupnd(i)=-1
                  endif
                enddo
#ifdef MULTI_LEV
              else if(flagtr.eq.2) then
                do nd=pnodess,num_ldmtr-1
                  if(proc_ldmtr(nd).eq.ip) then
                    ndpnd(ndp)=nd
                    idpnd(nd)=ndp
                    xpnd(ndp)=cx_ldmtr(nd)
                    ypnd(ndp)=cy_ldmtr(nd)
                    zpnd(ndp)=cz_ldmtr(nd)
                    lpnd(ndp)=l_ldmtr(nd)
                    ndp=ndp+1
                  endif
                enddo
                do i=0,ndp-1
                  nd=ndpnd(i)                  
                  if(next_ldmtr(nd).gt.pnodess &
                   .and.next_ldmtr(nd).le.ndpnd(ndp-1)) then
                    nextpnd(i)=idpnd(next_ldmtr(nd))      
                  else
                    nextpnd(i)=ndp
                  endif
                  if(daughter_ldmtr(nd).ne.-1) then
                    daupnd(i)=idpnd(daughter_ldmtr(nd))
                  else 
                    daupnd(i)=-1
                  endif
                enddo
#endif
              endif
! *** search particles need communication with ip
              do i=0,ncomp-1
                pn=pncomp(i)
                nd=0
   73           xij=x_dm(pn)-xpnd(nd)
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
                yij=y_dm(pn)-ypnd(nd)
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
                zij=z_dm(pn)-zpnd(nd)
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
                rc=xij*xij+yij*yij+zij*zij
! *** only take into account gather condition ***
                rh=h_dm(pn)+crit*lpnd(nd)
                rh=rh*rh
                if(rc.lt.rh) then
                  if(daupnd(nd).eq.-1) then
                    if(cproc(ip).ge.isproc(ip).and.cproc(ip) &
                     .le.ieproc(ip)) then
                      pn_nfp(isend)=pn
                      jjlen(ip)=jjlen(ip)+1
                      isend=isend+1
                    endif
                    cproc(ip)=cproc(ip)+1
                    goto 90
                  endif
                  nd=daupnd(nd)
                else
                  nd=nextpnd(nd)
                endif
                if(nd.ge.ndp) then
                  goto 90
                endif
                goto 73
   90         enddo

              if(cproc(ip).ne.npproc(ip)) then
                write(6,*) ' Error in set_valuedm():npproc,cproc,' &
                 ,'myrank,ip,ndp,nit,idiv,ndivcom,ncomp,flagtr=' &
                 ,npproc(ip),cproc(ip),myrank,ip,ndp,nit,idivcom,ndivcom,ncomp 
                write(fileo,'(a10,i3.3)') 'setvdmerrp',myrank,flagtr
                open(60,file=fileo,status='unknown')
                do i=0,ncomp-1
                  pn=pncomp(i)
                  write(60,'(4(1pE13.5))') x_dm(pn),y_dm(pn),z_dm(pn),h_dm(pn)
                enddo
                close(60)
                stop
              endif
            endif
! keep original jjlen 
            ireqs(ip)=jjlen(ip)
          enddo

          deallocate(idpnd)
          deallocate(ndpnd)
          deallocate(xpnd)
          deallocate(ypnd)
          deallocate(zpnd)
          deallocate(lpnd)
          deallocate(nextpnd)
          deallocate(daupnd)

        endif
! *** update idivcom
        idivcom=idivcom+1
! *** np for this communication
        ncompi=isend
! *** getting the total number of particles received at each proc ***
        npj=0
        do ip=0,nprocs-1
          irecv=0
          call MPI_SCATTER(jjlen,1,MPI_INTEGER &
           ,irecv,1,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
          npjr(ip)=irecv
          npj=npj+irecv
        enddo
        if(npj.gt.MNDM) then
          write(6,*) ' Error in set_valuedm(): npj receiving >MNDM'
          write(6,*) ' npj,MNB,ndiv,idiv=',npj,MNDM,ndivcom,idivcom
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
        snval=7
! ncompi = number of particles to be sent

        allocate(tdvs(0:ncompi*snval))

        do i=0,ncompi-1
          pni=pn_nfp(i)
          tdvs(snval*i)=x_dm(pni)
          tdvs(snval*i+1)=y_dm(pni)
          tdvs(snval*i+2)=z_dm(pni)
          tdvs(snval*i+3)=tx(pni)
          tdvs(snval*i+4)=ty(pni)
          tdvs(snval*i+5)=tz(pni)
          tdvs(snval*i+6)=h_dm(pni)
        enddo
! *** sending proc ***
        do i=0,nprocs-1
          idisp(i)=idisp(i)*snval
          jjlen(i)=ireqs(i)*snval
        enddo

! reallocate tdvr and tivr to receive npj particles
        allocate(tdvr(0:npj*nval))
        allocate(tivr(0:npj*nival))

! *** sending and reveiving the data ***
        irecv=0
        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! *** sethe data to tdvr ***
          do i=0,npjr(ip)-1
            tdvr(irecv)=trbuf(snval*i)
            tdvr(irecv+npj)=trbuf(snval*i+1)
            tdvr(irecv+npj*2)=trbuf(snval*i+2)
            tdvr(irecv+npj*3)=trbuf(snval*i+3)
            tdvr(irecv+npj*4)=trbuf(snval*i+4)
            tdvr(irecv+npj*5)=trbuf(snval*i+5)
            tdvr(irecv+npj*6)=trbuf(snval*i+6)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo       

        deallocate(tdvs)

        if(irecv.ne.npj) then
          write(6,*) ' Error in set_value(): irecv,npj=',irecv,npj
          write(6,*) ' after sending the data to the other pe'
          stop
        endif
! *** initialization ***
        do i=0,npj*nival-1
          tivr(i)=0
        enddo
        do i=npj*snval,npj*nval-1
          tdvr(i)=0.0d0
        enddo
        if(flagtr.eq.1) then
          goto 91
#ifdef MULTI_LEV
        else if(flagtr.eq.2) then
          goto 92
#endif
       else 
          if(myrank.eq.0) then
          write(*,*) ' Error in set_valuedm: flagtr=',flagtr
          endif
          call  MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif            
      endif

! *** set omgh_p and etc. ***
   94 do i=0,npj-1
        pni=pn_nfp(i)
        if(rho_dm(pni).gt.0.0d0) then
#ifdef SIM1D
          omgh_dm(pni)=1.0d0+(h_dm(pni)/rho_dm(pni))*omgh_dm(pni)
#elif defined(SIM2D)
          omgh_dm(pni)=1.0d0+0.5d0*(h_dm(pni)/rho_dm(pni))*omgh_dm(pni)
#else
          omgh_dm(pni)=1.0d0+THIRD*(h_dm(pni)/rho_dm(pni))*omgh_dm(pni)
#endif
#ifdef SGRAV
          zetah_dm(pni)=-THIRD*(h_dm(pni)/rho_dm(pni))*zetah_dm(pni)
#else
          zetah_dm(pni)=0.0d0
#endif
        endif
      enddo

! * Check smoothing length: tdvr(i+np*6)=h_dm(pni)
! * criterion h=etah (m/rho)^1/3

      allocate(list(0:npj))

      tnlist=npj
      do i=0,tnlist-1
        list(i)=pn_nfp(i)
      enddo
      nlist = 0
! *** nalist strores the particle whose density=0 ***
      nssl=0
      if(nit.lt.mnitb) then

        allocate(nalist(0:tnlist))

        do i = 0,tnlist-1
          pn = list(i)
          if(nnb_dm(pn).le.0) then
! *** density is zero. can happen for star ***
            nalist(nssl)=pn 
            nssl=nssl+1
            pn_nfp(nlist)=pn
            nlist = nlist+1
            flagi(pn)=0
          else
#ifdef SIM1D
! *** target density ***
            av(pn)=m_dm(pn)*(ETAHDM/h_dm(pn))
! *** difference between target and real density ***
            adv(pn)=av(pn)-rho_dm(pn)
! *** new h ***
            av(pn)=h_dm(pn)-adv(pn)/(-rho_dm(pn)*omgh_dm(pn)/h_dm(pn))
#elif defined(SIM2D)
! *** target density ***
            av(pn)=m_dm(pn)*((ETAHDM/h_dm(pn))**2)
! *** difference between target and real density ***
            adv(pn)=av(pn)-rho_dm(pn)
! *** new h ***
            av(pn)=h_dm(pn)-adv(pn)/(-2.0d0*rho_dm(pn)*omgh_dm(pn)/h_dm(pn))
#else
! *** target density ***
            av(pn)=m_dm(pn)*((ETAHDM/h_dm(pn))**3)
! *** difference between target and real density ***
            adv(pn)=av(pn)-rho_dm(pn)
! *** new h ***
            av(pn)=h_dm(pn)-adv(pn)/(-3.0d0*rho_dm(pn)*omgh_dm(pn)/h_dm(pn))
! tried below, but convergence is slower. 8/7 2008
!           av(pn)=h_dm(pn)-adv(pn)/(-3.0d0*av(pn)*omgh_dm(pn)/h_dm(pn))
             
#endif
            if(dabs(h_dm(pn)-av(pn))/h_dm(pn).gt.ERRH) then
! tried below, but too strict condtion 8/7 2008
!              if(dabs(h_dm(pn)-av(pn))/du0(pn).gt.ERRH) then
#ifdef SGRAV
! *** check for minimum h ***
              if(flagi(pn).eq.10) then
                omgh_dm(pn)=1.0d0
                zetah_dm(pn)=0.0d0
              else
#endif
! *** store the data for bi-section
              if(adv(pn).gt.0.0d0.and.adv(pn).lt.fhx(pn)) then
                hu(pn)=h_dm(pn)
                fhx(pn)=adv(pn)
              else if(adv(pn).lt.0.0d0.and.adv(pn).gt.fx(pn)) then
                lu(pn)=h_dm(pn)
                fx(pn)=adv(pn)
              endif
! *** update ***
              if(av(pn).gt.DHFLIM*h_dm(pn)) then
                h_dm(pn)=DHFLIM*h_dm(pn)
              else if(av(pn).gt.0.0d0) then
                 h_dm(pn)=av(pn)
              else
                h_dm(pn)=0.5d0*h_dm(pn)
              endif
              pn_nfp(nlist)=pn
              nlist = nlist+1
#ifdef SGRAV
            endif
#ifdef DMEPSLIM_G
            else if(h_dm(pn).le.SI_eps*(m_dm(pn)**THIRD) &
             .and.flagi(pn).ne.10) then
              pn_nfp(nlist)=pn
              nlist = nlist+1
              flagi(pn)=10
              h_dm(pn)=SI_eps*(m_dm(pn)**THIRD)
            else if(flagi(pn).eq.10) then
              omgh_dm(pn)=1.0d0
              zetah_dm(pn)=0.0d0
#endif
#endif
            endif
          endif
        enddo
        if(nssl.gt.0) then
! *** set h for particles whose density is zero, using local _dmtr ***

          allocate(sldmtrlist(0:nssl))

          do i=0,nssl-1
            sldmtrlist(i)=nalist(i)
          enddo

          call setsldmtr(nssl,sldmtrlist,ndm)

          deallocate(sldmtrlist)

        endif

!        deallocate(nalist)

      else
        do i = 0,tnlist-1
          pn = list(i)
#ifdef SIM1D
! *** target density ***
          av(pn)=m_dm(pn)*(ETAHDM/h_dm(pn))
! *** difference between target and real density ***
          adv(pn)=av(pn)-rho_dm(pn)
#elif defined(SIM2D)
! *** target density ***
          av(pn)=m_dm(pn)*((ETAHDM/h_dm(pn))**2)
! *** difference between target and real density ***
          adv(pn)=av(pn)-rho_dm(pn)
#else
! *** target density ***
          av(pn)=m_dm(pn)*((ETAHDM/h_dm(pn))**3)
! *** difference between target and real density ***
          adv(pn)=av(pn)-rho_dm(pn)
#endif
        enddo
        if(nit.eq.mnitb) then
          if(nssl.gt.0) then
            write(6,*) ' Error in set_value(): density is zero.after iteration.'
            write(6,*) 'flag,myrank,nssl=',flag,myrank,nssl
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
            stop
          endif
          do i=0,tnlist-1
            pn=list(i)
            if(hu(pn).gt.0.0d0.and.lu(pn).gt.0.0d0) then
              flagi(pn)=0
            else if(lu(pn).lt.0.0d0) then
              flagi(pn)=-1
            else if(hu(pn).lt.0.0d0) then
              h_dm(pn)=lu(pn)
              flagi(pn)=1
              if(lu(pn).lt.0.0d0) then
                write(6,*) ' Error in set_valuedm(): rank,id=',myrank,id_dm(pn)
                write(6,*) ' no guess for initial h'
                write(6,*) ' h,rho,flagi=',h_dm(pn),rho_dm(pn),flagi(pn)
                write(6,*) ' lu,hu,av,adv=',lu(pn),hu(pn),av(pn),adv(pn)
                stop
              endif
            else
              write(6,*) ' Error in set_value(): rank,id=',myrank,id_dm(pn)
              write(6,*) ' no guess for initial h'
              write(6,*) ' h,rho=',h_dm(pn),rho_dm(pn)
              call MPI_ABORT(MPI_COMM_WORLD,ierr)
              stop
            endif
          enddo
        endif
! *** bi-section ***
        do i=0,tnlist-1
          pn=list(i)
          if(flagi(pn).eq.0) then
            if(adv(pn)*fhx(pn).lt.0.0d0) then
              lu(pn)=h_dm(pn)
            else
              hu(pn)=h_dm(pn)
              fhx(pn)=adv(pn)
            endif
            av(pn)=0.5d0*(lu(pn)+hu(pn))
          else if(flagi(pn).eq.1) then
! *** need to find f(hu)>0 ***
            if(adv(pn).gt.0.0d0.and.nnb_dm(pn).gt.0) then
              fhx(pn)=adv(pn)
              hu(pn)=h_dm(pn)
              av(pn)=0.5d0*(lu(pn)+hu(pn))
              flagi(pn)=0
            else if(nnb_dm(pn).le.0) then
              av(pn)=1.01d0*h_dm(pn)
            else
! *** believe that decreasing h would increase f(h) ***
              av(pn)=0.8d0*h_dm(pn)
            endif
          else if(flagi(pn).eq.-1) then
! **** Note: hu < lu ***
            if(adv(pn).lt.0.0d0.and.nnb_dm(pn).gt.0) then
              lu(pn)=h_dm(pn)
              if(hu(pn).gt.0.0d0) then
                av(pn)=0.5d0*(lu(pn)+hu(pn))
                flagi(pn)=0
              else
                fx(pn)=adv(pn)
! *** now turned to flagi(pn) and decreasing h to find hu
                flagi(pn)=1
              endif
            else
              av(pn)=1.2d0*h_dm(pn)
            endif
          endif
          if(dabs(h_dm(pn)-av(pn))/h_dm(pn).gt.ERRH &
! tried below, but too strict condtion 8/7 2008
!            if(dabs(h_dm(pn)-av(pn))/du0(pn).gt.ERRH
           .or.flagi(pn).ne.0) then
#ifdef SGRAV
            if(flagi(pn).eq.10) then
              omgh_dm(pn)=1.0d0 
              zetah_dm(pn)=0.0d0
            else
#endif
! *** update ***
            h_dm(pn)=av(pn)
            pn_nfp(nlist)=pn
            nlist = nlist+1
#ifdef SGRAV
            endif
#ifdef DMEPSLIM_G
          else if(h_dm(pn).le.SI_eps*(m_dm(pn)**THIRD) &
           .and.flagi(pn).ne.10) then
            pn_nfp(nlist)=pn
            nlist = nlist+1
            flagi(pn)=10
            h_dm(pn)=SI_eps*(m_dm(pn)**THIRD)
          else if(flagi(pn).eq.10) then
            omgh_dm(pn)=1.0d0 
            zetah_dm(pn)=0.0d0
#endif
#endif
          endif
        enddo
      endif


      flagnb=0
      call MPI_ALLREDUCE(nlist,flagnb,1,MPI_INTEGER &
       ,MPI_SUM,MPI_COMM_WORLD,ierr)

      deallocate(list)
      if(nit.lt.mnitb) then
        deallocate(nalist)
      endif

      if(flagnb.gt.0) then
        if(nit.gt.mnit) then
          if(myrank.eq.0) then
            write(6,*) myrank &
             ,' Warning: iteration does not work in set_valuedm',nlist,flagnb,flag
          endif
          write(fileo,'(a7,i3.3)') 'errordm',myrank
          open(60,file=fileo,status='unknown')
          do i=0,nlist-1
            pn=pn_nfp(i)
            write(60,'(3I6,3(1pE13.5))') myrank,pn,nnb_dm(pn),h_dm(pn) &
             ,rho_dm(pn),ETAHDM*((m_dm(pn)/rho_dm(pn))**THIRD)
            write(60,'(I3,7(1pE13.5))') flagi(pn),av(pn),adv(pn) &
             ,hu(pn),lu(pn),fhx(pn),fx(pn),du0(pn)
          enddo
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
! *** since gas h not changed, neighbour search is needed only for 
! *** stars or gas which changed h.
        npj=nlist
        nit = nit+1
! *** initialze proc iteration
        iproc=0
! *** initialisation ***
        idivcom=0
        ndivcom=-1
        do i=0,nprocs-1
          npproc(i)=0
        enddo

        deallocate(flagproc)
        deallocate(flagcom)
        deallocate(pncomp)
        deallocate(talist)

        goto 71
      endif
! * define div(v) & |rot(v)| *
  999 do i=0,nadm-1 
        pn=list_adm(i)
        div_v_dm(pn)=-div_v_dm(pn)/rho_dm(pn)
      enddo

#ifdef DMPSLIM_G
! *** check
      do i=0,np-1      
        pn=list_adm(i)
        if(h_dm(pn).lt.SI_eps*(m_dm(pn)**THIRD)*0.9999d0) then
          write(6,*) ' Error in set_value(): too small h,nit=',nit
          write(6,*) 'myrank,h,h0,av,flag=',myrank,h_dm(pn),du0(pn) &
           ,av(pn),SI_eps*(m_dm(pn)**THIRD),flagi(pn)
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
      enddo
#endif

! deallocate
      deallocate(idisp)
      deallocate(jjlen)
      deallocate(ireqs)
      deallocate(npproc)
      deallocate(isproc)
      deallocate(ieproc)
      deallocate(cproc)
      deallocate(npjr)

      deallocate(flagi)
      deallocate(fhx)
      deallocate(fx)
      deallocate(hu)
      deallocate(lu)
      deallocate(du0)
      deallocate(av)
      deallocate(adv)
      deallocate(nd_nfp)
      deallocate(flagproc)
      deallocate(flagcom)
      deallocate(pncomp)
      deallocate(talist)
      deallocate(tx)
      deallocate(ty)
      deallocate(tz)


#endif
end subroutine
#endif

