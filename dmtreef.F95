#include "gcdp.def"
! ***************************************
!     dmtreef.F for GCD+ ver. f03.0
!   3 Feb. 2013    written by D. Kawata
! ***************************************
! *** Definition of treeforce() ***
! * calculate acceraration *

#if defined(DM) && defined (SGRAV)
subroutine dmtreeforce(ndm,nadm,np,flags)
      use gcdp_const
      use gcdp_system
      use gcdp_kernel
      use gcdp_dm
      use gcdp_dmtree
#ifdef MULTI_LEV
      use gcdp_ldmtree
#endif
#if defined(GAS) && defined(STAR)
      use gcdp_baryon
      use gcdp_btree
#endif

      implicit none
      include 'mpif.h'

      integer,intent(in) :: ndm,nadm,np,flags
      integer i,level,ip,is
      integer nval,isend,npj,srank,rrank,irecv
! * Number of Notfinished particle, temp *      
      integer nlist,tnlist
! * Node searching for each Particle *      
      integer nd,pn,pnj,ntf,flagp
      double precision r2,ir2,rij,ir1,ir3
! * for calculate gradients *
      double precision dwi,dwj
! * for test far away *
      double precision theta2
      double precision l_r2
      double precision xij,yij,zij
      double precision dphiir,dphijr,si,sj,hsi,hsj
#ifdef QPTREEF
      double precision dp_ij,d2phii,d2phij,d3phii,d3phij,d2p_ij,d3p_ij
      double precision x2ij,y2ij,z2ij,x3ij,y3ij,z3ij
      double precision ir5,ir7,ir4
      double precision dpijr
#endif
! *** potential and correction adaptive softening ***
      double precision pot,dhc
! *** information particles need communication ***
      integer iproc,snval,rnval,snival,npjs
      integer ncomp,ncompt,ndivcom,idivcom
      integer ncomptmax,npprocmax,ndivcompt
      integer ncompi
! *** for pseudo node 
      integer ndp,pnodess
! *** flag for tree ***
      integer flagtr
! *** for test file *** 
      character fileo*60
! for work for communication fixed size array
      integer pncompt(0:MNDM-1),proccompt(0:MNDM-1)
      integer ierr
! pn_nfp can be bigger than ng+ns in this proc, when receiving the particles
      integer pn_nfp(0:MNDM-1)
! for work allocatable
      integer,allocatable :: npltf(:),flago(:)
      integer,allocatable :: pncomp(:),flagproc(:),flagcom(:)
      integer,allocatable :: npproc(:),isproc(:),ieproc(:),cproc(:)
      integer,allocatable :: nppnd(:),nextpnd(:),daupnd(:) &
       ,ndpnd(:),idpnd(:)
      integer,allocatable :: idisp(:),jjlen(:),ireqs(:),npjr(:)
      integer,allocatable :: list(:),node(:)
      integer,allocatable :: tivs(:),tivr(:)
      double precision,allocatable :: xijp(:),yijp(:),zijp(:) &
       ,r2p(:),fap(:)
      double precision,allocatable :: xpnd(:),ypnd(:),zpnd(:) &
       ,lpnd(:),deltapnd(:),mpnd(:),cxpnd(:),cypnd(:),czpnd(:) &
       ,trbuf(:)
      double precision,allocatable :: tdvs(:),tdvr(:)

! ***** Initialization *****
      theta2 = THETA*THETA
      nval=11

!allocate work space for MPI
      allocate(idisp(0:nprocs))
      allocate(jjlen(0:nprocs))
      allocate(ireqs(0:nprocs))
      allocate(npproc(0:nprocs))
      allocate(isproc(0:nprocs))
      allocate(ieproc(0:nprocs))
      allocate(cproc(0:nprocs))
      allocate(npjr(0:nprocs))

      allocate(fap(0:ndm))

      do i=0,nadm-1
        pn=list_adm(i)
        fap(pn)=dsqrt(pdvx_dm(pn)**2+pdvy_dm(pn)**2+pdvz_dm(pn)**2)
      enddo

      allocate(tdvr(0:nval*nadm))
      allocate(flagproc(0:nadm))
      allocate(flagcom(0:nadm))
      allocate(pncomp(0:nadm))

#ifdef GAS
      flagtr=0
#else
      flagtr=1
#endif
   71 do i=0,nadm-1
        pn=list_adm(i)
        pn_nfp(i)=pn
        tdvr(i)=x_dm(pn)
        tdvr(i+nadm)=y_dm(pn)
        tdvr(i+nadm*2)=z_dm(pn)
        tdvr(i+nadm*3)=h_dm(pn)
        tdvr(i+nadm*4)=omgh_dm(pn)
        tdvr(i+nadm*5)=zetah_dm(pn)
        tdvr(i+nadm*9)=dt_dm(pn)
        tdvr(i+nadm*10)=fap(pn)
        flagproc(i)=myrank
        flagcom(i)=-1
      enddo
! *** initialization for dv?_dm ***
      do i=nadm*6,(nval-2)*nadm-1 
        tdvr(i)=0.0d0
      enddo

      npj=nadm
      ncomp=0
      ncompt=0
      iproc=0
      idivcom=0
      ndivcom=-1
      do i=0,nprocs-1
        npproc(i)=0
      enddo

#ifdef GAS    
      if(flagtr.eq.0) then
        goto 70
      else if(flagtr.eq.1) then
#else
      if(flagtr.eq.1) then
#endif
        goto 91
#ifdef MULTI_LEV
      else if(flagtr.eq.2) then
        goto 92
#endif
      else
        goto 94
      endif

#ifdef GAS
   70 if(np_tr(0).eq.0) then
       goto 93
      endif        
      if(iproc.ne.0.and.proc_tr(0).ne.myrank) then
        goto 93
      endif
! *** initialization ***

      allocate(list(0:npj))
      allocate(node(0:npj))

      do i=0,npj-1
        list(i)=i
        node(i)=0 
      enddo
      nlist=npj
      level=0

! *** start tree walk within a proc ***
   77 if(nlist.le.0) then

        deallocate(list)
        deallocate(node)

        goto 93
      endif
      level = level+1
      if(level.gt.MAXNODEDM) then
        write(6,*) ' Error in dmtreef(): failure in walking tree'
        write(6,*) 'myrank,iproc.idivcom,ndivcom=' &
         ,myrank,iproc,idivcom,ndivcom
        write(fileo,'(a7,i3.3)') 'errdtfb',myrank
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pn=list(i)
          write(60,'(6(1pE13.5),2I10)') tdvr(pn),tdvr(pn+npj),tdvr(pn+npj*2)
        enddo
        stop
      endif

      allocate(npltf(0:nlist))
      allocate(xijp(0:nlist))
      allocate(yijp(0:nlist))
      allocate(zijp(0:nlist))
      allocate(r2p(0:nlist))
      allocate(flago(0:nlist))

      ntf=0
      do i=0,nlist-1
        pn = list(i)         
        nd = node(pn)
        flago(pn)=0
        if(np_tr(nd).gt.1) then
! *** check for particles within a node ***
          xij=tdvr(pn)-cx_tr(nd)
          yij=tdvr(pn+npj)-cy_tr(nd)		
          zij=tdvr(pn+npj*2)-cz_tr(nd)
          if(dabs(xij).lt.0.6d0*l_tr(nd)) then
          if(dabs(yij).lt.0.6d0*l_tr(nd)) then
          if(dabs(zij).lt.0.6d0*l_tr(nd)) then
! *** go to daughter node ***
            flago(pn)=1
          endif
          endif
          endif
        endif
        if(flago(pn).eq.0) then
          xij=tdvr(pn)-cmx_tr(nd)
          yij=tdvr(pn+npj)-cmy_tr(nd)		
          zij=tdvr(pn+npj*2)-cmz_tr(nd)
          r2=xij*xij+yij*yij+zij*zij
          if(r2.gt.0.0d0) then
! *** register the particles need distance check ***
            npltf(ntf)=pn
! *** store distance from the center of the mass ***
            xijp(ntf)=xij
            yijp(ntf)=yij
            zijp(ntf)=zij
            r2p(ntf)=r2
            ntf=ntf+1
          else if(np_tr(nd).eq.1) then
            node(pn)=next_tr(nd)		
          else
            flago(pn)=1
          endif
        endif
      enddo
      if(flags.eq.0) then
        do i=0,ntf-1
          pn=npltf(i)
          nd=node(pn)
          if(np_tr(nd).eq.1) then
            flago(pn)=2
          else
#ifdef QPTREEF
            l_r2=(G*mass_tr(nd))*((l_tr(nd)**4)/(r2p(i)**3))
#else
            l_r2=(G*mass_tr(nd))*((l_tr(nd)**2)/(r2p(i)**2))
#endif
            if(l_r2.le.ALPHAT*tdvr(pn+npj*10)) then
              flago(pn)=3
            else
! *** go to daughter ***
              flago(pn)=1
            endif
          endif
        enddo
      else
! *** for step 0 or continue ***
        do i=0,ntf-1
          pn=npltf(i)
          nd=node(pn)
          if(np_tr(nd).eq.1) then
            flago(pn)=2
          else
            l_r2=(l_tr(nd)+delta_tr(nd))*(l_tr(nd)+delta_tr(nd))/r2p(i)
            if(l_r2.le.theta2) then
              flago(pn)=3
            else
! *** go to daughter ***
              flago(pn)=1
            endif
          endif
        enddo
      endif

! *** tree force calculation
      do i=0,ntf-1
        pn=npltf(i)
        nd=node(pn)
        if(flago(pn).eq.2) then
          if(iproc.eq.0.or.nd.gt.nodese_tr) then
            pnj = pn_tr(nd)
            rij=dsqrt(r2p(i))
            ir2=1.0d0/r2p(i)
            ir3=ir2/rij
! *** dphi/dr(hi) and dwi ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
                *(si-s_tb(is))/ds_tb
              dphiir=dphiir/(hsi**3)
! dw/dr/r
              dwi=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dwi=dwi/(hsi**5)
            else
              dphiir=ir3
              dwi=0.0d0
            endif
! *** dphi/dr(hj) ***
            hsj=h_p(pnj)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dphijr=dphijr/(hsj**3)
! dw/dr/r
              dwj=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dwj=dwj/(hsj**5)
            else
              dphijr=ir3
              dwj=0.0d0
            endif
! *** check for timestep ***
            if(si.lt.1.0d0.or.sj.lt.1.0d0) then
              if(tdvr(pn+npj*9).lt.dtmnb_p(pnj)) then
                dtmnb_p(pnj)=tdvr(pn+npj*9)
              endif
            endif
! *** force calculation ***
! *** dvx_p,dvy_p,dvz_p ***
            pot=m_p(pnj)*0.5d0*G*(dphiir+dphijr)
            dhc=m_p(pnj)*0.5d0*G*((tdvr(pn+npj*5)/tdvr(pn+npj*4))*dwi &
             +(zetah_p(pnj)/omgh_p(pnj))*dwj)
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-xijp(i)*(pot+dhc)
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-yijp(i)*(pot+dhc)
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-zijp(i)*(pot+dhc)
          endif
! * update node *
          node(pn)=next_tr(nd)
        else if(flago(pn).eq.3) then
          if(iproc.eq.0.or.nd.gt.nodese_tr) then
            ir2=1.0d0/r2p(i)
            r2=r2p(i)
            rij=dsqrt(r2)       
            ir1=1.0d0/rij
            ir3=ir2*ir1
! *** use softening of particle i
! *** dphi/dr(hi) ***
            hsi=tdvr(pn+npj*3)
            rij=dsqrt(r2)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dphiir=dphiir/(hsi**3)
#ifdef QPTREEF
! d2phi/dr2(hi)
              d2phii=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is)) &
               *(si-s_tb(is))/ds_tb
              d2phii=d2phii/(hsi**3)
! d3phi/dr3(hi)
              d3phii=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is)) &
               *(si-s_tb(is))/ds_tb
              d3phii=d3phii/(hsi**4)
#endif
            else
              dphiir=ir3
#ifdef QPTREEF

! d2phi/dr2(hi)
              d2phii=-2.0d0*ir3
! d3phi/dr3(hi)
              d3phii=6.0d0*ir2*ir2
#endif
            endif
! *** dphi/dr(hj) ***
            hsj=hm_tr(nd)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dphijr=dphijr/(hsj**3)
#ifdef QPTREEF
! d2phi/dr2(hj)
              d2phij=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              d2phij=d2phij/(hsj**3)
! d3phi/dr3(hj)
              d3phij=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              d3phij=d3phij/(hsj**4)
#endif
            else
              dphijr=ir3
#ifdef QPTREEF
! d2pj/dr2(hj)
              d2phij=-2.0d0*ir3
! d3phj/dr3(hj)
              d3phij=6.0d0*ir2*ir2
#endif
            endif
#ifdef QPTREEF
! *** force calculation ***
! *** dvx_p,dvy_p,dvz_p ***
            xij=xijp(i)
            yij=yijp(i)
            zij=zijp(i)
            dpijr=0.5d0*(dphiir+dphijr)
            pot=mass_tr(nd)*G*dpijr
            dp_ij=dpijr*rij
            d2p_ij=0.5d0*(d2phii+d2phij)
            d3p_ij=0.5d0*(d3phii+d3phij)
            ir4=ir2**2
            ir5=ir3*ir2
            ir7=ir5*ir2
            x2ij=xij**2
            y2ij=yij**2
            z2ij=zij**2
            x3ij=xij*x2ij
            y3ij=yij*y2ij
            z3ij=zij*z2ij
! *** ax ***
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-(xij*pot &
! dipole
             -mx_tr(nd)*((ir1-x2ij*ir3)*dp_ij+x2ij*ir2*d2p_ij) &
             -my_tr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij) &
             -mz_tr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_tr(nd)*((-3.0d0*xij*ir3+3.0d0*x3ij*ir5)*dp_ij &
             +(3.0d0*xij*ir2-3.0d0*x3ij*ir4)*d2p_ij+x3ij*ir3*d3p_ij) &
!
             +0.5d0*myy_tr(nd)*((-xij*ir3+3.0d0*xij*y2ij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*xij*y2ij*ir4)*d2p_ij+xij*y2ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_tr(nd)*((-xij*ir3+3.0d0*xij*z2ij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*xij*z2ij*ir4)*d2p_ij+xij*z2ij*ir3*d3p_ij) &
!
             +mxy_tr(nd)*((-yij*ir3+3.0d0*x2ij*yij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*x2ij*yij*ir4)*d2p_ij+x2ij*yij*ir3*d3p_ij) &
!
             +mzx_tr(nd)*((-zij*ir3+3.0d0*x2ij*zij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*x2ij*zij*ir4)*d2p_ij+x2ij*zij*ir3*d3p_ij) &
!
             +myz_tr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
              +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij))
! *** ay ***
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-(yij*pot &
! dipole
             -mx_tr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij) &
             -my_tr(nd)*((ir1-y2ij*ir3)*dp_ij+y2ij*ir2*d2p_ij) &
             -mz_tr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_tr(nd)*((-yij*ir3+3.0d0*yij*x2ij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*yij*x2ij*ir4)*d2p_ij+yij*x2ij*ir3*d3p_ij) &
!
             +0.5d0*myy_tr(nd)*((-3.0d0*yij*ir3+3.0d0*y3ij*ir5)*dp_ij &
             +(3.0d0*yij*ir2-3.0d0*y3ij*ir4)*d2p_ij+y3ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_tr(nd)*((-yij*ir3+3.0d0*yij*z2ij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*yij*z2ij*ir4)*d2p_ij+yij*z2ij*ir3*d3p_ij) &
!
             +mxy_tr(nd)*((-xij*ir3+3.0d0*y2ij*xij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*y2ij*xij*ir4)*d2p_ij+y2ij*xij*ir3*d3p_ij) &
!
             +mzx_tr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij)) &
!
             +myz_tr(nd)*((-zij*ir3+3.0d0*y2ij*zij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*y2ij*zij*ir4)*d2p_ij+y2ij*zij*ir3*d3p_ij)
! *** az ***
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-(zij*pot &
! dipole
             -mx_tr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij) &
             -my_tr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij) &
             -mz_tr(nd)*((ir1-z2ij*ir3)*dp_ij+z2ij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_tr(nd)*((-zij*ir3+3.0d0*zij*x2ij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*zij*x2ij*ir4)*d2p_ij+zij*x2ij*ir3*d3p_ij) &
!
             +0.5d0*myy_tr(nd)*((-zij*ir3+3.0d0*zij*y2ij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*zij*y2ij*ir4)*d2p_ij+zij*y2ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_tr(nd)*((-3.0d0*zij*ir3+3.0d0*z3ij*ir5)*dp_ij &
             +(3.0d0*zij*ir2-3.0d0*z3ij*ir4)*d2p_ij+z3ij*ir3*d3p_ij) &
!
             +mxy_tr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij)) &
!
             +mzx_tr(nd)*((-xij*ir3+3.0d0*z2ij*xij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*z2ij*xij*ir4)*d2p_ij+z2ij*xij*ir3*d3p_ij) &
!
             +myz_tr(nd)*((-yij*ir3+3.0d0*z2ij*yij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*z2ij*yij*ir4)*d2p_ij+z2ij*yij*ir3*d3p_ij)
#else
! *** force calculation ***
! *** dvx_p,dvy_p,dvz_p ***
            pot=mass_tr(nd)*0.5d0*G*(dphiir+dphijr)
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-(xijp(i)*pot)
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-(yijp(i)*pot)
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-(zijp(i)*pot)
#endif
          endif
! * update node *			
          node(pn)=next_tr(nd)
        endif
      enddo
      do i=0,nlist-1
        pn=list(i)
        if(flago(pn).eq.1) then
          nd=node(pn)
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_tr(nd).eq.-1) then
            if(proc_tr(nd).ne.myrank.and.flagproc(pn).ne.proc_tr(nd)) then
! *** proc_tr will be the same for different domain within the same proc
              flagproc(pn)=proc_tr(nd)
! *** store pn ***
              if(flagcom(pn).lt.0) then
                flagcom(pn)=proc_tr(nd)
                pncomp(ncomp)=pn_nfp(pn)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNDM) then
                pncompt(ncompt)=pn_nfp(pn)
                proccompt(ncompt)=proc_tr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_tr(nd))=npproc(proc_tr(nd))+1
            endif
            node(pn)=next_tr(nd)
          else
            node(pn) = daughter_tr(nd)
          endif
        endif
      enddo

      deallocate(npltf)
      deallocate(xijp)
      deallocate(yijp)
      deallocate(zijp)
      deallocate(r2p)
      deallocate(flago)

! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      if(iproc.eq.0.and.next_tr(0).ne.0) then
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          else if(node(list(i)).eq.0) then
            node(list(i))=next_tr(0)
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
#else
   70 goto 93
#endif

! ***   Calculate force contributed from DM particle ***
! * Initialization *
   91 if(np_dmtr(0).eq.0) then
        goto 93
      endif
      if(iproc.ne.0.and.proc_dmtr(0).ne.myrank) then
        goto 93
      endif

      allocate(list(0:npj))
      allocate(node(0:npj))

      do i=0,npj-1
        list(i)=i
        node(i)=0
      enddo
      nlist=npj
      level=0
  777 if(nlist.le.0) then

        deallocate(list)
        deallocate(node)

        goto 93
      endif
      level = level+1
      if(level.gt.MAXNODEDM) then
        write(6,*) ' Error in dmtreef(): failure in walking dmtree'
        write(6,*) 'myrank,iproc.idivcom,ndivcom=',myrank,iproc,idivcom,ndivcom
        write(fileo,'(a8,i3.3)') 'errdtfdt',myrank
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pn=list(i)
          write(60,'(6(1pE13.5),2I10)') tdvr(pn),tdvr(pn+npj),tdvr(pn+npj*2)
        enddo
        stop
      endif

      allocate(npltf(0:nlist))
      allocate(xijp(0:nlist))
      allocate(yijp(0:nlist))
      allocate(zijp(0:nlist))
      allocate(r2p(0:nlist))
      allocate(flago(0:nlist))

      ntf=0
      do i=0,nlist-1
        pn = list(i)         
        nd = node(pn)
        flago(pn)=0
        if(np_dmtr(nd).gt.1) then
          xij=tdvr(pn)-cx_dmtr(nd)
          yij=tdvr(pn+npj)-cy_dmtr(nd)		
          zij=tdvr(pn+npj*2)-cz_dmtr(nd)
          if(dabs(xij).lt.0.6d0*l_dmtr(nd)) then
          if(dabs(yij).lt.0.6d0*l_dmtr(nd)) then
          if(dabs(zij).lt.0.6d0*l_dmtr(nd)) then
            flago(pn)=1
          endif
          endif
          endif
        endif
        if(flago(pn).eq.0) then
          xij=tdvr(pn)-cmx_dmtr(nd)
          yij=tdvr(pn+npj)-cmy_dmtr(nd)		
          zij=tdvr(pn+npj*2)-cmz_dmtr(nd)
          r2=xij*xij+yij*yij+zij*zij
          if(r2.gt.0.0d0) then
! *** register the particles need distance calculation ***
            npltf(ntf)=pn
! *** store distance from the center of the mass ***
            xijp(ntf)=xij
            yijp(ntf)=yij
            zijp(ntf)=zij
            r2p(ntf)=r2
            ntf=ntf+1
          else if(np_dmtr(nd).eq.1) then
            node(pn)=next_dmtr(nd)		
          else
            flago(pn)=1
          endif
        endif
      enddo
      if(flags.eq.0) then
        do i=0,ntf-1
          pn=npltf(i)
          nd=node(pn)
          if(np_dmtr(nd).eq.1) then
            flago(pn)=2
          else
#ifdef QPTREEF
            l_r2=(G*mass_dmtr(nd))*((l_dmtr(nd)**4)/(r2p(i)**3))
#else
            l_r2=(G*mass_dmtr(nd))*((l_dmtr(nd)**2)/(r2p(i)**2))
#endif
            if(l_r2.le.ALPHAT*tdvr(pn+npj*10)) then
              flago(pn)=3
            else
! *** go to daughter ***
              flago(pn)=1
            endif
          endif
        enddo
      else
! *** for step 0 or continue ***
        do i=0,ntf-1
          pn=npltf(i)
          nd=node(pn)
          if(np_dmtr(nd).eq.1) then
            flago(pn)=2
          else
            l_r2=(l_dmtr(nd)+delta_dmtr(nd))*(l_dmtr(nd)+delta_dmtr(nd))/r2p(i)
            if(l_r2.le.theta2) then
              flago(pn)=3
            else
! *** go to daughter ***
              flago(pn)=1
            endif
          endif
        enddo
      endif
! *** tree force calculation
      do i=0,ntf-1
        pn=npltf(i)
        nd=node(pn)
        if(flago(pn).eq.2) then
          if(iproc.eq.0.or.nd.gt.nodese_dmtr) then
            pnj = pn_dmtr(nd)
            rij=dsqrt(r2p(i))
            ir2=1.0d0/r2p(i)
            ir3=ir2/rij
! *** dphi/dr(hi) and dwi ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dphiir=dphiir/(hsi**3)
! dw/dr/r
              dwi=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dwi=dwi/(hsi**5)
            else
              dphiir=ir3
              dwi=0.0d0
            endif
! *** dphi/dr(hj) ***
            hsj=h_dm(pnj)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dphijr=dphijr/(hsj**3)
! dw/dr/r
              dwj=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dwj=dwj/(hsj**5)
            else
              dphijr=ir3
              dwj=0.0d0
            endif
! *** check for timestep ***
            if(si.lt.1.0d0.or.sj.lt.1.0d0) then
              if(tdvr(pn+npj*9).lt.dtmnb_dm(pnj)) then
                dtmnb_dm(pnj)=tdvr(pn+npj*9)
              endif
            endif
! *** force calculation ***
! *** dvx_p,dvy_p,dvz_p ***
            pot=m_dm(pnj)*0.5d0*G*(dphiir+dphijr)
            dhc=m_dm(pnj)*0.5d0*G*((tdvr(pn+npj*5)/tdvr(pn+npj*4))*dwi &
             +(zetah_dm(pnj)/omgh_dm(pnj))*dwj)
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-xijp(i)*(pot+dhc)
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-yijp(i)*(pot+dhc)
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-zijp(i)*(pot+dhc)
          endif
! * update node *		  
          node(pn)=next_dmtr(nd)
        else if(flago(pn).eq.3) then
          if(iproc.eq.0.or.nd.gt.nodese_dmtr) then
            ir2=1.0d0/r2p(i)
            r2=r2p(i)
            rij=dsqrt(r2)       
            ir1=1.0d0/rij
            ir3=ir2*ir1
! *** use softening of particle i
! *** dphi/dr(hi) ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dphiir=dphiir/(hsi**3)
#ifdef QPTREEF
! d2phi/dr2(hi)
              d2phii=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is)) &
               *(si-s_tb(is))/ds_tb
              d2phii=d2phii/(hsi**3)
! d3phi/dr3(hi)
              d3phii=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is)) &
               *(si-s_tb(is))/ds_tb
              d3phii=d3phii/(hsi**4)
#endif
            else
              dphiir=ir3
#ifdef QPTREEF
! d2phi/dr2(hi)
              d2phii=-2.0d0*ir3
! d3phi/dr3(hi)
              d3phii=6.0d0*ir2*ir2
#endif
            endif
! *** dphi/dr(hj) ***
            hsj=hm_dmtr(nd)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dphijr=dphijr/(hsj**3)
#ifdef QPTREEF
! d2phi/dr2(hj)
              d2phij=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              d2phij=d2phij/(hsj**3)
! d3phi/dr3(hj)
              d3phij=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              d3phij=d3phij/(hsj**4)
#endif
            else
              dphijr=ir3
#ifdef QPTREEF
! d2phj/dr2(hj)
              d2phij=-2.0d0*ir3
! d3phj/dr3(hj)
              d3phij=6.0d0*ir2*ir2
#endif
            endif
! *** force calculation ***
! *** dvx_p,dvy_p,dvz_p ***
#ifdef QPTREEF
            xij=xijp(i)
            yij=yijp(i)
            zij=zijp(i)
            dpijr=0.5d0*(dphiir+dphijr)
            pot=mass_dmtr(nd)*G*dpijr
            dp_ij=dpijr*rij
            d2p_ij=0.5d0*(d2phii+d2phij)
            d3p_ij=0.5d0*(d3phii+d3phij)
            ir4=ir2**2
            ir5=ir3*ir2
            ir7=ir5*ir2
            x2ij=xij**2
            y2ij=yij**2
            z2ij=zij**2
            x3ij=xij*x2ij
            y3ij=yij*y2ij
            z3ij=zij*z2ij
! *** ax ***
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-(xij*pot &
! dipole
             -mx_dmtr(nd)*((ir1-x2ij*ir3)*dp_ij+x2ij*ir2*d2p_ij) &
             -my_dmtr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij) &
             -mz_dmtr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_dmtr(nd)*((-3.0d0*xij*ir3+3.0d0*x3ij*ir5)*dp_ij &
             +(3.0d0*xij*ir2-3.0d0*x3ij*ir4)*d2p_ij+x3ij*ir3*d3p_ij) &
!
             +0.5d0*myy_dmtr(nd)*((-xij*ir3+3.0d0*xij*y2ij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*xij*y2ij*ir4)*d2p_ij+xij*y2ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_dmtr(nd)*((-xij*ir3+3.0d0*xij*z2ij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*xij*z2ij*ir4)*d2p_ij+xij*z2ij*ir3*d3p_ij) &
!
             +mxy_dmtr(nd)*((-yij*ir3+3.0d0*x2ij*yij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*x2ij*yij*ir4)*d2p_ij+x2ij*yij*ir3*d3p_ij) &
!
             +mzx_dmtr(nd)*((-zij*ir3+3.0d0*x2ij*zij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*x2ij*zij*ir4)*d2p_ij+x2ij*zij*ir3*d3p_ij) &
!
             +myz_dmtr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij))
! *** ay ***
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-(yij*pot &
! dipole
             -mx_dmtr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij) &
             -my_dmtr(nd)*((ir1-y2ij*ir3)*dp_ij+y2ij*ir2*d2p_ij) &
             -mz_dmtr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_dmtr(nd)*((-yij*ir3+3.0d0*yij*x2ij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*yij*x2ij*ir4)*d2p_ij+yij*x2ij*ir3*d3p_ij) &
!
             +0.5d0*myy_dmtr(nd)*((-3.0d0*yij*ir3+3.0d0*y3ij*ir5)*dp_ij &
             +(3.0d0*yij*ir2-3.0d0*y3ij*ir4)*d2p_ij+y3ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_dmtr(nd)*((-yij*ir3+3.0d0*yij*z2ij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*yij*z2ij*ir4)*d2p_ij+yij*z2ij*ir3*d3p_ij) &
!
             +mxy_dmtr(nd)*((-xij*ir3+3.0d0*y2ij*xij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*y2ij*xij*ir4)*d2p_ij+y2ij*xij*ir3*d3p_ij) &
!
             +mzx_dmtr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij)) &
!
             +myz_dmtr(nd)*((-zij*ir3+3.0d0*y2ij*zij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*y2ij*zij*ir4)*d2p_ij+y2ij*zij*ir3*d3p_ij)
! *** az ***
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-(zij*pot &
! dipole
             -mx_dmtr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij) &
             -my_dmtr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij) &
             -mz_dmtr(nd)*((ir1-z2ij*ir3)*dp_ij+z2ij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_dmtr(nd)*((-zij*ir3+3.0d0*zij*x2ij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*zij*x2ij*ir4)*d2p_ij+zij*x2ij*ir3*d3p_ij) &
!
             +0.5d0*myy_dmtr(nd)*((-zij*ir3+3.0d0*zij*y2ij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*zij*y2ij*ir4)*d2p_ij+zij*y2ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_dmtr(nd)*((-3.0d0*zij*ir3+3.0d0*z3ij*ir5)*dp_ij &
             +(3.0d0*zij*ir2-3.0d0*z3ij*ir4)*d2p_ij+z3ij*ir3*d3p_ij) &
!
             +mxy_dmtr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij)) &
!
             +mzx_dmtr(nd)*((-xij*ir3+3.0d0*z2ij*xij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*z2ij*xij*ir4)*d2p_ij+z2ij*xij*ir3*d3p_ij) &
!
             +myz_dmtr(nd)*((-yij*ir3+3.0d0*z2ij*yij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*z2ij*yij*ir4)*d2p_ij+z2ij*yij*ir3*d3p_ij)
#else
            pot=mass_dmtr(nd)*0.5d0*G*(dphiir+dphijr) 
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-(xijp(i)*pot)
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-(yijp(i)*pot)
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-(zijp(i)*pot)
#endif
          endif
! * update node *		  
          node(pn)=next_dmtr(nd)
        endif
      enddo
      do i=0,nlist-1
        pn = list(i)
        if(flago(pn).eq.1) then
          nd = node(pn)		
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_dmtr(nd).eq.-1) then
            if(proc_dmtr(nd).ne.myrank.and.flagproc(pn).ne.proc_dmtr(nd)) then
! *** proc_tr will be the same for different domain within the same proc
              flagproc(pn)=proc_dmtr(nd)
! *** store pn ***
              if(flagcom(pn).lt.0) then
                flagcom(pn)=proc_dmtr(nd)
                pncomp(ncomp)=pn_nfp(pn)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNDM) then
                pncompt(ncompt)=pn_nfp(pn)
                proccompt(ncompt)=proc_dmtr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_dmtr(nd))=npproc(proc_dmtr(nd))+1
            endif
            node(pn)=next_dmtr(nd)
          else
            node(pn) = daughter_dmtr(nd)
          endif
        endif
      enddo       

      deallocate(npltf)
      deallocate(xijp)
      deallocate(yijp)
      deallocate(zijp)
      deallocate(r2p)
      deallocate(flago)
 
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
      goto 777        

! ***   Calculate force contributed from low resolution DM particle ***
#ifdef MULTI_LEV
   92 if(np_ldmtr(0).eq.0) then
        goto 93
      endif
      if(iproc.ne.0.and.proc_ldmtr(0).ne.myrank) then
        goto 93
      endif

      allocate(list(0:npj))
      allocate(node(0:npj))

      do i=0,npj-1
        list(i)=i
        node(i)=0 
      enddo
      nlist=npj
      level=0
 7777 if(nlist.le.0) then

        deallocate(list)
        deallocate(node)

        goto 93
      endif
      level = level+1
      if(level.gt.MAXNODELDM) then
        write(6,*) ' Error in dmtreef(): failure in walking ldmtree'
        write(6,*) 'myrank,iproc.idivcom,ndivcom=',myrank,iproc,idivcom,ndivcom
        write(fileo,'(a9,i3.3)') 'errdtfldt',myrank
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pn=list(i)
          write(60,'(6(1pE13.5),2I10)') tdvr(pn),tdvr(pn+npj),tdvr(pn+npj*2)
        enddo
        stop
      endif

      allocate(npltf(0:nlist))
      allocate(xijp(0:nlist))
      allocate(yijp(0:nlist))
      allocate(zijp(0:nlist))
      allocate(r2p(0:nlist))
      allocate(flago(0:nlist))

      ntf=0
      do i=0,nlist-1
        pn = list(i)         
        nd = node(pn)
        flago(pn)=0
        if(np_ldmtr(nd).gt.1) then
          xij=tdvr(pn)-cx_ldmtr(nd)
          yij=tdvr(pn+npj)-cy_ldmtr(nd)		
          zij=tdvr(pn+npj*2)-cz_ldmtr(nd)
          if(dabs(xij).lt.0.6d0*l_ldmtr(nd)) then
          if(dabs(yij).lt.0.6d0*l_ldmtr(nd)) then
          if(dabs(zij).lt.0.6d0*l_ldmtr(nd)) then
            flago(pn)=1
          endif
          endif
          endif
        endif
        if(flago(pn).eq.0) then
          xij=tdvr(pn)-cmx_ldmtr(nd)
          yij=tdvr(pn+npj)-cmy_ldmtr(nd)		
          zij=tdvr(pn+npj*2)-cmz_ldmtr(nd)
          r2=xij*xij+yij*yij+zij*zij
          if(r2.gt.0.0d0) then
! *** register the particles need distance calculation ***
            npltf(ntf)=pn
! *** store distance from the center of the mass ***
            xijp(ntf)=xij
            yijp(ntf)=yij
            zijp(ntf)=zij
            r2p(ntf)=r2
            ntf=ntf+1
          else if(np_ldmtr(nd).eq.1) then
            node(pn) = next_ldmtr(nd)		
          else
            flago(pn)=1
          endif
        endif
      enddo
      if(flags.eq.0) then
        do i=0,ntf-1
          pn=npltf(i)
          nd=node(pn)
          if(np_ldmtr(nd).eq.1) then
            flago(pn)=2
          else
#ifdef QPTREEF
            l_r2=(G*mass_ldmtr(nd))*((l_ldmtr(nd)**4)/(r2p(i)**3))
#else
            l_r2=(G*mass_ldmtr(nd))*((l_ldmtr(nd)**2)/(r2p(i)**2))
#endif
            if(l_r2.le.ALPHAT*tdvr(pn+npj*10)) then
              flago(pn)=3
            else
! *** go to daughter ***
              flago(pn)=1
            endif
          endif
        enddo
      else
! *** for step 0 or continue ***
        do i=0,ntf-1
          pn=npltf(i)
          nd=node(pn)
          if(np_ldmtr(nd).eq.1) then
            flago(pn)=2
          else
            l_r2=(l_ldmtr(nd)+delta_ldmtr(nd)) &
             *(l_ldmtr(nd)+delta_ldmtr(nd))/r2p(i)
            if(l_r2.le.theta2) then
              flago(pn)=3
            else
! *** go to daughter ***
              flago(pn)=1
            endif
          endif
        enddo
      endif

! *** tree force calculation
      do i=0,ntf-1
        pn=npltf(i)
        nd=node(pn)
        if(flago(pn).eq.2) then
          if(iproc.eq.0.or.nd.gt.nodese_ldmtr) then
            pnj = pn_ldmtr(nd)
            rij=dsqrt(r2p(i))
            ir2=1.0d0/r2p(i)
            ir3=ir2/rij
! *** dphi/dr(hi) and dwi ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dphiir=dphiir/(hsi**3)
! dw/dr/r
              dwi=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dwi=dwi/(hsi**5)
            else
              dphiir=ir3
              dwi=0.0d0
            endif
! *** dphi/dr(hj) ***
            hsj=h_dm(pnj)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dphijr=dphijr/(hsj**3)
! dw/dr/r
              dwj=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dwj=dwj/(hsj**5)
            else
              dphijr=ir3
              dwj=0.0d0
            endif
! *** check for timestep ***
            if(si.lt.1.0d0.or.sj.lt.1.0d0) then
              if(tdvr(pn+npj*9).lt.dtmnb_dm(pnj)) then
                dtmnb_dm(pnj)=tdvr(pn+npj*9)
              endif
            endif
! *** force calculation ***
! *** dvx_p,dvy_p,dvz_p ***
            pot=m_dm(pnj)*0.5d0*G*(dphiir+dphijr)
            dhc=m_dm(pnj)*0.5d0*G*((tdvr(pn+npj*5)/tdvr(pn+npj*4))*dwi &
             +(zetah_dm(pnj)/omgh_dm(pnj))*dwj)
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-xijp(i)*(pot+dhc)
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-yijp(i)*(pot+dhc)
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-zijp(i)*(pot+dhc)
          endif
! * update node *		  
          node(pn)=next_ldmtr(nd)
        else if(flago(pn).eq.3) then
          if(iproc.eq.0.or.nd.gt.nodese_ldmtr) then
            ir2=1.0d0/r2p(i)
            r2=r2p(i)
            rij=dsqrt(r2)
            ir1=1.0d0/rij
            ir3=ir2*ir1
! *** use softening of particle i
! *** dphi/dr(hi) ***
             hsi=tdvr(pn+npj*3)
             si=rij/hsi
             if(si.lt.1.0d0) then
              is=int(si/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))/ds_tb
              dphiir=dphiir/(hsi**3)
#ifdef QPTREEF
! d2phi/dr2(hi)
              d2phii=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is)) &
               *(si-s_tb(is))/ds_tb
              d2phii=d2phii/(hsi**3)
! d3phi/dr3(hi)
              d3phii=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is)) &
               *(si-s_tb(is))/ds_tb
              d3phii=d3phii/(hsi**4)
#endif
            else
              dphiir=ir3
#ifdef QPTREEF
! d2phi/dr2(hi)
              d2phii=-2.0d0*ir3
! d3phi/dr3(hi)
              d3phii=6.0d0*ir2*ir2
#endif
            endif
! *** dphi/dr(hj) ***
            hsj=hm_ldmtr(nd)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj/ds_tb)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              dphijr=dphijr/(hsj**3)
#ifdef QPTREEF
! d2phi/dr2(hj)
              d2phij=d2phidr2_tb(is)+(d2phidr2_tb(is+1)-d2phidr2_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              d2phij=d2phij/(hsj**3)
! d3phi/dr3(hj)
              d3phij=d3phidr3_tb(is)+(d3phidr3_tb(is+1)-d3phidr3_tb(is)) &
               *(sj-s_tb(is))/ds_tb
              d3phij=d3phij/(hsj**4)
#endif
            else
              dphijr=ir3
#ifdef QPTREEF
! d2phj/dr2(hj)
              d2phij=-2.0d0*ir3
! d3phj/dr3(hj)
              d3phij=6.0d0*ir2*ir2
#endif
            endif
! *** force calculation ***
! *** dvx_p,dvy_p,dvz_p ***
#ifdef QPTREEF
            xij=xijp(i)
            yij=yijp(i)
            zij=zijp(i)
            dpijr=0.5d0*(dphiir+dphijr)
            pot=mass_ldmtr(nd)*G*dpijr
            dp_ij=dpijr*rij
            d2p_ij=0.5d0*(d2phii+d2phij)
            d3p_ij=0.5d0*(d3phii+d3phij)
            ir4=ir2**2
            ir5=ir3*ir2
            ir7=ir5*ir2
            x2ij=xij**2
            y2ij=yij**2
            z2ij=zij**2
            x3ij=xij*x2ij
            y3ij=yij*y2ij
            z3ij=zij*z2ij
! *** ax ***
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-(xij*pot &
! dipole
             -mx_ldmtr(nd)*((ir1-x2ij*ir3)*dp_ij+x2ij*ir2*d2p_ij) &
             -my_ldmtr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij) &
             -mz_ldmtr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_ldmtr(nd)*((-3.0d0*xij*ir3+3.0d0*x3ij*ir5)*dp_ij &
             +(3.0d0*xij*ir2-3.0d0*x3ij*ir4)*d2p_ij+x3ij*ir3*d3p_ij) &
!
             +0.5d0*myy_ldmtr(nd)*((-xij*ir3+3.0d0*xij*y2ij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*xij*y2ij*ir4)*d2p_ij+xij*y2ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_ldmtr(nd)*((-xij*ir3+3.0d0*xij*z2ij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*xij*z2ij*ir4)*d2p_ij+xij*z2ij*ir3*d3p_ij) &
!
             +mxy_ldmtr(nd)*((-yij*ir3+3.0d0*x2ij*yij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*x2ij*yij*ir4)*d2p_ij+x2ij*yij*ir3*d3p_ij) &
!
             +mzx_ldmtr(nd)*((-zij*ir3+3.0d0*x2ij*zij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*x2ij*zij*ir4)*d2p_ij+x2ij*zij*ir3*d3p_ij) &
!
             +myz_ldmtr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij)) 
! *** ay ***
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-(yij*pot &
! dipole
             -mx_ldmtr(nd)*(-xij*yij*ir3*dp_ij+xij*yij*ir2*d2p_ij) &
             -my_ldmtr(nd)*((ir1-y2ij*ir3)*dp_ij+y2ij*ir2*d2p_ij) &
             -mz_ldmtr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij) &
!n quadrupole
             +0.5d0*mxx_ldmtr(nd)*((-yij*ir3+3.0d0*yij*x2ij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*yij*x2ij*ir4)*d2p_ij+yij*x2ij*ir3*d3p_ij) &
!
             +0.5d0*myy_ldmtr(nd)*((-3.0d0*yij*ir3+3.0d0*y3ij*ir5)*dp_ij &
             +(3.0d0*yij*ir2-3.0d0*y3ij*ir4)*d2p_ij+y3ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_ldmtr(nd)*((-yij*ir3+3.0d0*yij*z2ij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*yij*z2ij*ir4)*d2p_ij+yij*z2ij*ir3*d3p_ij) &
!
             +mxy_ldmtr(nd)*((-xij*ir3+3.0d0*y2ij*xij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*y2ij*xij*ir4)*d2p_ij+y2ij*xij*ir3*d3p_ij) &
!
             +mzx_ldmtr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij)) &
!
             +myz_ldmtr(nd)*((-zij*ir3+3.0d0*y2ij*zij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*y2ij*zij*ir4)*d2p_ij+y2ij*zij*ir3*d3p_ij)
! *** az ***
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-(zij*pot &
! dipole
             -mx_ldmtr(nd)*(-xij*zij*ir3*dp_ij+xij*zij*ir2*d2p_ij) &
             -my_ldmtr(nd)*(-yij*zij*ir3*dp_ij+yij*zij*ir2*d2p_ij) &
             -mz_ldmtr(nd)*((ir1-z2ij*ir3)*dp_ij+z2ij*ir2*d2p_ij) &
! quadrupole
             +0.5d0*mxx_ldmtr(nd)*((-zij*ir3+3.0d0*zij*x2ij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*zij*x2ij*ir4)*d2p_ij+zij*x2ij*ir3*d3p_ij) &
!
             +0.5d0*myy_ldmtr(nd)*((-zij*ir3+3.0d0*zij*y2ij*ir5)*dp_ij &
             +(zij*ir2-3.0d0*zij*y2ij*ir4)*d2p_ij+zij*y2ij*ir3*d3p_ij) &
!
             +0.5d0*mzz_ldmtr(nd)*((-3.0d0*zij*ir3+3.0d0*z3ij*ir5)*dp_ij &
             +(3.0d0*zij*ir2-3.0d0*z3ij*ir4)*d2p_ij+z3ij*ir3*d3p_ij) &
!
             +mxy_ldmtr(nd)*(3.0d0*xij*yij*zij*ir5*dp_ij &
             +3.0d0*xij*yij*zij*ir4*d2p_ij+xij*yij*zij*ir3*d3p_ij)) &
!
             +mzx_ldmtr(nd)*((-xij*ir3+3.0d0*z2ij*xij*ir5)*dp_ij &
             +(xij*ir2-3.0d0*z2ij*xij*ir4)*d2p_ij+z2ij*xij*ir3*d3p_ij) &
!
             +myz_ldmtr(nd)*((-yij*ir3+3.0d0*z2ij*yij*ir5)*dp_ij &
             +(yij*ir2-3.0d0*z2ij*yij*ir4)*d2p_ij+z2ij*yij*ir3*d3p_ij)
#else
            pot=mass_ldmtr(nd)*0.5d0*G*(dphiir+dphijr)
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-(xijp(i)*pot)
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-(yijp(i)*pot)
            tdvr(pn+npj*8)=tdvr(pn+npj*8)-(zijp(i)*pot)
#endif
          endif
! * update node *		  
          node(pn)=next_ldmtr(nd)
        endif
      enddo
      do i=0,nlist-1
        pn=list(i)
        if(flago(pn).eq.1) then
          nd=node(pn)
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_ldmtr(nd).eq.-1) then
            if(proc_ldmtr(nd).ne.myrank.and.flagproc(pn).ne.proc_ldmtr(nd)) then
! *** proc_tr will be the same for different domain within the same proc
              flagproc(pn)=proc_ldmtr(nd)
! *** store pn ***
              if(flagcom(pn).lt.0) then
                flagcom(pn)=proc_ldmtr(nd)
                pncomp(ncomp)=pn_nfp(pn)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNDM) then
                pncompt(ncompt)=pn_nfp(pn)
                proccompt(ncompt)=proc_ldmtr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_ldmtr(nd))=npproc(proc_ldmtr(nd))+1
            endif
            node(pn)=next_ldmtr(nd)
          else
            node(pn) = daughter_ldmtr(nd)
          endif
        endif
      enddo      

      deallocate(npltf)
      deallocate(xijp)
      deallocate(yijp)
      deallocate(zijp)
      deallocate(r2p)
      deallocate(flago)
  
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
      goto 7777
#endif

   93 if(iproc.ne.0) then
! *** sending back the results ***
! *** new snval to send the data back
        snval=3

        allocate(tivs(0:nprocs))

        do i=0,nprocs-1
! *** keep npjr: data received from each proc -> now sending***
          tivs(i)=npjr(i)
! *** ireqs=jjlen number of particles sent -> now receiving.
          npjr(i)=ireqs(i)
        enddo
! *** send and receive double precision values ***
        isend=0
        do i=0,nprocs-1
          idisp(i)=isend 
! *** tivs is numbr of particles received from each proc -> now sending
          jjlen(i)=tivs(i)*snval
          isend=isend+jjlen(i)
        enddo

        deallocate(tivs)
        allocate(tdvs(0:npj*snval))

! *** setting the sending back data ***
        do i=0,npj-1
! *** ax, ay, az
          tdvs(snval*i)=tdvr(i+npj*6)
          tdvs(snval*i+1)=tdvr(i+npj*7)
          tdvs(snval*i+2)=tdvr(i+npj*8)
        enddo
! *** sending and reveiving the data ***
        irecv=0
! *** number of particles sent is ncomp 
        npjs=npj
        npj=ncompi

        deallocate(tdvr)
        allocate(tdvr(0:npj*nval))

        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION &
           ,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tdvr(irecv+npj*6)=trbuf(snval*i)
            tdvr(irecv+npj*7)=trbuf(snval*i+1)
            tdvr(irecv+npj*8)=trbuf(snval*i+2)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo

        deallocate(tdvs)    

      endif
! *** update the values ***
      do i=0,npj-1
        pn=pn_nfp(i)
        dvx_dm(pn)=dvx_dm(pn)+tdvr(i+npj*6)
        dvy_dm(pn)=dvy_dm(pn)+tdvr(i+npj*7)
        dvz_dm(pn)=dvz_dm(pn)+tdvr(i+npj*8)
      enddo

      deallocate(tdvr)

! *** end of neighbour search ***
      if(nprocs.le.1.or.(iproc.ne.0.and.idivcom.gt.ndivcom)) then
        flagtr=flagtr+1
        if(flagtr.eq.1) then
          goto 71
#ifdef MULTI_LEV
        else if(flagtr.eq.2) then
          goto 71
#endif
        else
          goto 94
        endif
      endif

      if(nprocs.gt.1.and.iproc.eq.0) then
        iproc=1
! *** check if need communication ***
        ncomptmax=0
        call MPI_ALLREDUCE(ncompt,ncomptmax,1,MPI_INTEGER &
         ,MPI_MAX,MPI_COMM_WORLD,ierr)
        if(ncomptmax.gt.0) then

          allocate(tivr(0:nprocs))

! *** get how many particles each proc receives
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
        else
! no communication is required    
          flagtr=flagtr+1
          if(flagtr.eq.1) then
            goto 71
#ifdef MULTI_LEV
          else if(flagtr.eq.2) then
            goto 71
#endif
          else
            goto 94
          endif
        endif
      endif
! *** do communication ndivcom times
      if(idivcom.le.ndivcom) then
        if(ndivcom.eq.1) then
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
        else
! *** set the range of number of particles for sending for each proc
          do i=0,nprocs-1
            call para_range(0,npproc(i)-1,ndivcom,idivcom-1,isproc(i),ieproc(i))
            cproc(i)=0
          enddo
! *** preparation for sending the data to the other procs ***
          isend=0          
! *** store particle list in the order of sending procs in list() ***
          do ip=0,nprocs-1
            idisp(ip)=isend           
            jjlen(ip)=0
            if(npproc(ip).gt.0.and.ip.ne.myrank) then
! *** get the coordinate and etc. for pseudo node for the proc

              ndp=0
              if(flagtr.eq.0) then
#ifdef GAS
                if(proc_tr(0).eq.myrank) then
                  pnodess=next_tr(0)
                else
                  pnodess=0
                endif

                ndp=num_tr-pnodess-1
                allocate(idpnd(0:num_tr))
                allocate(nppnd(0:ndp))
                allocate(ndpnd(0:ndp))
                allocate(xpnd(0:ndp))
                allocate(ypnd(0:ndp))
                allocate(zpnd(0:ndp))
                allocate(lpnd(0:ndp))
                allocate(mpnd(0:ndp))
                allocate(cxpnd(0:ndp))
                allocate(cypnd(0:ndp))
                allocate(czpnd(0:ndp))
                allocate(deltapnd(0:ndp))
                allocate(nextpnd(0:ndp))
                allocate(daupnd(0:ndp))

                ndp=0
                do nd=pnodess,num_tr-1
                  if(proc_tr(nd).eq.ip) then
                    ndpnd(ndp)=nd
                    idpnd(nd)=ndp
                    nppnd(ndp)=np_tr(nd)
                    xpnd(ndp)=cmx_tr(nd)
                    ypnd(ndp)=cmy_tr(nd)
                    zpnd(ndp)=cmz_tr(nd)
                    lpnd(ndp)=l_tr(nd)
                    mpnd(ndp)=mass_tr(nd)
                    cxpnd(ndp)=cx_tr(nd)
                    cypnd(ndp)=cy_tr(nd)
                    czpnd(ndp)=cz_tr(nd)
                    deltapnd(ndp)=delta_tr(nd)
                    ndp=ndp+1
                  endif
                enddo
                do i=0,ndp-1
                  nd=ndpnd(i)                  
                  if(next_tr(nd).gt.pnodess &
                   .and.next_tr(nd).le.ndpnd(ndp-1)) then
                    nextpnd(i)=idpnd(next_tr(nd))      
                  else
                    nextpnd(i)=ndp
                  endif
                  if(daughter_tr(nd).ne.-1) then
                    daupnd(i)=idpnd(daughter_tr(nd))
                  else 
                    daupnd(i)=-1
                  endif
                enddo
#else
                if(myrank.eq.0) then
                  write(6,*) ' Error in dmtreef: flagtr=1, but GAS is off'
                  call MPI_FINALIZE()
                  stop
                endif
#endif
#ifdef DM
              else if(flagtr.eq.1) then
                if(proc_dmtr(0).eq.myrank) then
                  pnodess=next_dmtr(0)
                else
                  pnodess=0
                endif

                ndp=num_dmtr-pnodess-1
                allocate(idpnd(0:num_dmtr))
                allocate(nppnd(0:ndp))
                allocate(ndpnd(0:ndp))
                allocate(xpnd(0:ndp))
                allocate(ypnd(0:ndp))
                allocate(zpnd(0:ndp))
                allocate(lpnd(0:ndp))
                allocate(mpnd(0:ndp))
                allocate(cxpnd(0:ndp))
                allocate(cypnd(0:ndp))
                allocate(czpnd(0:ndp))
                allocate(deltapnd(0:ndp))
                allocate(nextpnd(0:ndp))
                allocate(daupnd(0:ndp))

                ndp=0
                do nd=pnodess,num_dmtr-1
                  if(proc_dmtr(nd).eq.ip) then
                    ndpnd(ndp)=nd
                    idpnd(nd)=ndp
                    nppnd(ndp)=np_dmtr(nd)
                    xpnd(ndp)=cmx_dmtr(nd)
                    ypnd(ndp)=cmy_dmtr(nd)
                    zpnd(ndp)=cmz_dmtr(nd)
                    lpnd(ndp)=l_dmtr(nd)
                    mpnd(ndp)=mass_dmtr(nd)
                    cxpnd(ndp)=cx_dmtr(nd)
                    cypnd(ndp)=cy_dmtr(nd)
                    czpnd(ndp)=cz_dmtr(nd)
                    deltapnd(ndp)=delta_dmtr(nd)
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
                if(proc_ldmtr(0).eq.myrank) then
                  pnodess=next_ldmtr(0)
                else
                  pnodess=0
                endif

                ndp=num_ldmtr-pnodess-1
                allocate(idpnd(0:num_ldmtr))
                allocate(nppnd(0:ndp))
                allocate(ndpnd(0:ndp))
                allocate(xpnd(0:ndp))
                allocate(ypnd(0:ndp))
                allocate(zpnd(0:ndp))
                allocate(lpnd(0:ndp))
                allocate(mpnd(0:ndp))
                allocate(cxpnd(0:ndp))
                allocate(cypnd(0:ndp))
                allocate(czpnd(0:ndp))
                allocate(deltapnd(0:ndp))
                allocate(nextpnd(0:ndp))
                allocate(daupnd(0:ndp))

                ndp=0
                do nd=pnodess,num_ldmtr-1
                  if(proc_ldmtr(nd).eq.ip) then
                    ndpnd(ndp)=nd
                    idpnd(nd)=ndp
                    nppnd(ndp)=np_ldmtr(nd)
                    xpnd(ndp)=cmx_ldmtr(nd)
                    ypnd(ndp)=cmy_ldmtr(nd)
                    zpnd(ndp)=cmz_ldmtr(nd)
                    lpnd(ndp)=l_ldmtr(nd)
                    mpnd(ndp)=mass_ldmtr(nd)
                    cxpnd(ndp)=cx_ldmtr(nd)
                    cypnd(ndp)=cy_ldmtr(nd)
                    czpnd(ndp)=cz_ldmtr(nd)
                    deltapnd(ndp)=delta_ldmtr(nd)
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
#endif
              endif
! *** search particles need communication with ip
              do i=0,ncomp-1
                pn=pncomp(i)
                nd=0
   72           flagp=0
! *** check for particles within a node ***
                if(nppnd(nd).gt.1) then
                  xij=x_dm(pn)-cxpnd(nd)
                  yij=y_dm(pn)-cypnd(nd)
                  zij=z_dm(pn)-czpnd(nd)
                  if(dabs(xij).lt.0.6d0*lpnd(nd)) then
                  if(dabs(yij).lt.0.6d0*lpnd(nd)) then
                  if(dabs(zij).lt.0.6d0*lpnd(nd)) then
! *** go to daughter node ***
                    flagp=1
                  endif
                  endif
                  endif
                  if(flagp.eq.0) then
                    xij=x_dm(pn)-xpnd(nd)
                    yij=y_dm(pn)-ypnd(nd)
                    zij=z_dm(pn)-zpnd(nd)
                    r2=xij*xij+yij*yij+zij*zij
                    if(flags.eq.0) then
                      if(r2.gt.0.0d0) then
#ifdef QPTREEF
                        l_r2=(G*mpnd(nd))*((lpnd(nd)**4)/(r2**3))
#else
                        l_r2=(G*mpnd(nd))*((lpnd(nd)**2)/(r2**2))
#endif
                      else
                        l_r2=0.0d0
                      endif
                      if(l_r2.gt.ALPHAT*fap(pn)) then
                        flagp=1
                      endif
                    else
                      if(r2.gt.0.0d0) then
                        l_r2=(lpnd(nd)+deltapnd(nd))*(lpnd(nd)+deltapnd(nd))/r2
                      else
                        l_r2 = 0.0d0
                      endif
                      if(l_r2.gt.theta2) then
                        flagp=1
                      endif
                    endif
                  endif
                endif
                if(flagp.eq.1) then
                  if(daupnd(nd).eq.-1) then
                    if(cproc(ip).ge.isproc(ip).and.cproc(ip).le.ieproc(ip)) then
                      pn_nfp(isend)=pn
                      jjlen(ip)=jjlen(ip)+1
                      isend=isend+1
                    endif
                    cproc(ip)=cproc(ip)+1
                    goto 95
                  endif
                  nd=daupnd(nd)
                else
                  nd=nextpnd(nd)
                endif
                if(nd.ge.ndp) then
                  goto 95
                endif
                goto 72
   95         enddo

              if(cproc(ip).ne.npproc(ip)) then
                write(6,*) ' Error in dmtreef():npproc,cproc' &
                 ,',myrank,idiv,ndiv,ip,flagtr=' &
                 ,npproc(ip),cproc(ip),myrank,idivcom,ndivcom,ip,flagtr
                stop
              endif

              deallocate(idpnd)
              deallocate(nppnd)
              deallocate(ndpnd)
              deallocate(xpnd)
              deallocate(ypnd)
              deallocate(zpnd)
              deallocate(lpnd)
              deallocate(mpnd)
              deallocate(cxpnd)
              deallocate(cypnd)
              deallocate(czpnd)
              deallocate(deltapnd)
              deallocate(nextpnd)
              deallocate(daupnd)

            endif
! keep original jjlen 
            ireqs(ip)=jjlen(ip)
          enddo
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
! *** update npj ***
          npj=npj+irecv
        enddo
        if(npj.gt.MNDM) then
          write(6,*) ' Error in dmtreef(): npj > MNDM'
          write(6,*) ' npj,MNDM=',npj,MNDM
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
        snval=8

        allocate(tdvs(0:ncompi*snval))

        do i = 0,ncompi-1
          pn=pn_nfp(i)
! *** pn_nfp is not used in the other proc, so keep this ***
          tdvs(snval*i)=x_dm(pn)
          tdvs(snval*i+1)=y_dm(pn)
          tdvs(snval*i+2)=z_dm(pn)             
          tdvs(snval*i+3)=h_dm(pn)
          tdvs(snval*i+4)=omgh_dm(pn)
          tdvs(snval*i+5)=zetah_dm(pn)
          tdvs(snval*i+6)=dt_dm(pn)        
          tdvs(snval*i+7)=fap(pn)
        enddo
! *** reset sending parameters ***
        do i=0,nprocs-1
          idisp(i)=idisp(i)*snval
          jjlen(i)=ireqs(i)*snval
        enddo

        allocate(tdvr(0:npj*nval))

! *** sending and reveiving the data ***
        irecv=0
        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tdvr(irecv)=trbuf(snval*i)
            tdvr(irecv+npj)=trbuf(snval*i+1)
            tdvr(irecv+npj*2)=trbuf(snval*i+2)
            tdvr(irecv+npj*3)=trbuf(snval*i+3)
            tdvr(irecv+npj*4)=trbuf(snval*i+4)
            tdvr(irecv+npj*5)=trbuf(snval*i+5)
            tdvr(irecv+npj*9)=trbuf(snval*i+6)
            tdvr(irecv+npj*10)=trbuf(snval*i+7)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo       

        deallocate(tdvs)

! *** check ***
        if(irecv.ne.npj) then
          write(6,*) ' Error in dmtreef(): irecv,npj=',irecv,npj
          write(6,*) ' after sending the data to the other pe'
          stop
        endif
! *** initialization ***
! *** setting 0 for ax,ay,az ***
        do i=npj*6,npj*(nval-2)-1
          tdvr(i)=0.0d0
        enddo
        if(flagtr.eq.0) then
          goto 70
#ifdef DM
        else if(flagtr.eq.1) then
          goto 91
#ifdef MULTI_LEV
        else if(flagtr.eq.2) then
          goto 92
#endif
#endif
        else 
          if(myrank.eq.0) then
          write(*,*) ' Error in treef(): flagtr=',flagtr
          endif
          call  MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif            
      endif


   94 deallocate(idisp)
      deallocate(jjlen)
      deallocate(ireqs)
      deallocate(npproc)
      deallocate(isproc)
      deallocate(ieproc)
      deallocate(cproc)
      deallocate(npjr)

      deallocate(fap)
      deallocate(flagproc)
      deallocate(flagcom)
      deallocate(pncomp)

end subroutine
#endif      

