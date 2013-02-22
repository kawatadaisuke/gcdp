#include "gcdp.def"
! ****************************************************
!    calc_du_dv.F95 for gcd+ ver. f03.0
!  23 Jan., 2013    written by D.KAWATA
! **************************************************** 
! ************************************
!    Definition of calc_dv_du()
!  This program set dv/dt & du/dt
! ************************************

subroutine calc_dv_du(np,ng,ndm,nag,nadm,nas,nagravg,nagt &
 ,nadmt,nast,nagravgt,flagv,flagu,flags)
      use gcdp_const
      use gcdp_system
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#ifdef GAS
      use gcdp_gtree
      use gcdp_kernel
#endif
#ifdef METALDIFF
      use gcdp_metal
#endif
#endif
#ifdef DM
      use gcdp_dm
#endif
#ifdef MEXT
      use mext
#endif
      implicit none
      include 'mpif.h'

      integer,intent(in) :: np,ng,ndm,nag,nadm,nas,nagravg,nagt,nadmt,nast &
       ,nagravgt,flagv,flagu,flags
      integer i,j,is
      integer nd
      integer npj,srank,rrank,ip,nval,nival,isend,irecv,isend0
      integer nsn,irecv1
! number of particles in a list
      integer nlist,tnlist,nbp
! list for active and their neighbour paritcles
      integer,allocatable :: listb(:),nblistb(:)
! particle id
      integer pn,pni
! * distance, smoothing length
      double precision r,rc,rh,si,sj
      double precision hsi,ihsi,hsj,ihsj
! for calculate SPH weights
      double precision mdwi,mdwj,mwi
      double precision xij,yij,zij,tij
      double precision vxij,vyij,vzij
      double precision,allocatable :: mdwij(:),mdwic(:),mdwjc(:),rijp(:) &
       ,xijp(:),yijp(:),zijp(:),ihsip(:),ihsjp(:)
! * for artificial viscocity *
      double precision,allocatable :: vsigav(:),v_x_r(:),dvdiss(:) &
       ,dudissv(:),dudissu(:)
      double precision vsigu,alpij,alpuij,vsigdt,hvsigdt
! mean density
      double precision,allocatable :: rhoij(:)
! for Kinetic Equation *
      double precision,allocatable :: mpqdwj(:)
! u for j particles using updated p and cs, peff
      double precision,allocatable :: upj(:),peffp(:)
      double precision plim
! * etc *	
      double precision crit,tmp,zn
! *** for time ***
      double precision ntime,ectime
! *** information for particles need communication ***
      integer iproc,snval,rnval,snival,npjs
      integer ncomp,ncompt,ndivcom,idivcom
      integer ncomptmax,npprocmax,ndivcompt
      integer ncompi
      integer pncompt(0:MNB-1),proccompt(0:MNB-1)
      integer,allocatable :: pncomp(:),flagproc(:),flagcom(:)
      integer,allocatable :: npproc(:),isproc(:),ieproc(:),cproc(:)
! *** for pseudo node 
      integer ndp,pnodess
      integer,allocatable :: nextpnd(:),daupnd(:),ndpnd(:),idpnd(:)
      double precision,allocatable :: xpnd(:),ypnd(:),zpnd(:),lpnd(:),hmpnd(:)
      double precision,allocatable :: trbuf(:)
      integer,allocatable :: tibuf(:)
! *** for work ***
      integer ierr
      integer,allocatable :: idisp(:),jjlen(:),ireqs(:)
! pn_nfp can be bigger than ng+ns in this proc, when receiving the particles
      integer pn_nfp(0:MNB-1)
      integer,allocatable :: list(:),node(:),tivr(:),tivs(:),npjr(:),talist(:)
      double precision,allocatable :: tx(:),ty(:),tz(:),tdvr(:),tdvs(:)
#ifdef COOL   
      double precision,allocatable :: ram(:)
#endif
! *** for test file *** 
      character fileo*60
#ifdef GLASS
      real ran1
      external ran1
#endif
#ifdef ARTFDM
      double precision beta1,beta2,facaf,vvn2
      double precision,allocatable :: txf,tyf,tzf
#endif
#ifdef WTIMESUB
! *** to track time spent ***
      double precision nstime,estime,stime(4)
#endif
#ifdef METALDIFF
      integer nvalzd,czdi
      double precision,allocatable :: kij(:),av(:),tdzd(:)
#endif
#ifdef MEXT
      integer ir
      double precision rp,mr,pot
#endif
#ifdef B85TEST
      double precision rp,mpert,rpert,pot
      mpert=123.45679*0.05d0
      rpert=2.0d0*0.1d0
#endif

!      write(fileo,'(a6,i3.3)') 'pcalvu',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ng-1
!        pni=list_ap(i)
!        write(60,'(4(1pE13.5))') x_p(pni),y_p(pni),z_p(pni),ndu_p(pni)
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

! *** Basic Value ***
#ifdef SIM1D
      crit=0.5d0
#elif defined(SIM2D)
      crit=dsqrt(2.0d0)*0.5d0
#else
      crit=dsqrt(3.0d0)*0.5d0
#endif
      if(flagv.eq.0) then
        nval=22
      else
        nval=16
      endif
#if defined(SF_EFD) || defined(SF_ZFD)
      if(flagv.eq.0) then
        nival=2
      else
        nival=1
      endif
#else
      nival=0
#endif
#ifdef METALDIFF
      nvalzd=11
#endif
! *** initialization ***
      snval=1
      snival=1

      if(flagv.eq.0) then
! *** reset dt minimum of neighbours ***
#ifdef GAS
        do i=0,np-1
          dtmnb_p(i)=INF
        enddo
#endif
#ifdef DM 
        do i=0,ndm-1
          dtmnb_dm(i)=INF
        enddo
#endif
      endif
! *** Calculate Acceleration ***
      if(flagv.eq.0) then
#if defined(GEPSLIM_VH) || defined(GEPSLIM_G)
#ifdef GAS
! *** set min limited h temporary ***
        allocate(tx(0:ng))
        allocate(ty(0:ng))
        allocate(tz(0:ng))
        do i=0,ng-1
          pn=list_ap(i)
          tx(i)=h_p(pn)
          ty(i)=omgh_p(pn)
          tz(i)=zetah_p(pn)
          if(h_p(pn).lt.SI_eps*(m_p(pn)**THIRD)) then
            h_p(pn)=SI_eps*(m_p(pn)**THIRD)
            omgh_p(pn)=1.0d0
            zetah_p(pn)=0.0d0
          endif
        enddo
#endif
#endif
! *** initialization ***
#ifdef DM
        do i=0,nadm-1
          pn = list_adm(i)
          pdvx_dm(pn)=dvx_dm(pn)
          pdvy_dm(pn)=dvy_dm(pn)
          pdvz_dm(pn)=dvz_dm(pn)
          dvx_dm(pn)=0.0d0
          dvy_dm(pn)=0.0d0
          dvz_dm(pn)=0.0d0
        enddo
#endif
! *** add artificial force following Widrow (2008)
#ifdef ARTFDM
        if(TM_tot.lt.TG_AF) then
          facaf=3.0d0*(TM_tot/TG_AF)**2-2.0d0*(TM_tot/TG_AF)**3
          beta1=B1MAX_AF*facaf
          beta2=B2MAX_AF*facaf
        else if(TM_tot.lt.TMC_AF) then
          beta1=B1MAX_AF
          beta2=B2MAX_AF
        else if(TM_tot.lt.TME_AF) then
          facaf=1.0d0-(3.0d0*((TM_tot-TMC_AF)/TD_AF)**2 &
           -2.0d0*((TM_tot-TMC_AF)/TG_AF)**3)
          beta1=B1MAX_AF
          beta2=B2MAX_AF
        else
          beta1=0.0d0
          beta2=0.0d0
        endif
        if(TM_tot.lt.TME_AF) then
          do i=0,nadm-1
            pn=list_adm(i)
            if(id_dm(pn).ge.SI_artfdm(0).and.id_dm(pn).le.SI_artfdm(1)) then
! *** velocity at tn
              txf=vnx_dm(pn)+pdvx_dm(pn)*dt_dm(pn)
              tyf=vny_dm(pn)+pdvy_dm(pn)*dt_dm(pn)
              tzf=vnz_dm(pn)+pdvz_dm(pn)*dt_dm(pn)
              vvn2=(txf**2+tyf**2+tzf**2)
              if(vvn2.gt.0.0d0) then
                dvx_dm(pn)=dvx_dm(pn) &
                 +((beta1-beta2)*(tyf**2)+beta1*(tzf**2))*txf/vvn2
                dvy_dm(pn)=dvy_dm(pn) &
                 +((beta2-beta1)*(txf**2)+beta2*(tzf**2))*tyf/vvn2
                dvz_dm(pn)=dvz_dm(pn) &
                 -(beta1*(txf**2)+beta2*(tyf**2))*tzf/vvn2
              endif
            endif
          enddo
        endif
#endif

#ifdef GAS
        do i=0,nagravg-1
          pn=list_ap(i)
          pax_p(pn)=ax_p(pn)
          pay_p(pn)=ay_p(pn)
          paz_p(pn)=az_p(pn)
          ax_p(pn)=0.0d0
          ay_p(pn)=0.0d0
          az_p(pn)=0.0d0
        enddo
        do i=0,nag-1
          pn=list_ap(i)
          dvx_p(pn)=0.0d0
          dvy_p(pn)=0.0d0
          dvz_p(pn)=0.0d0
        enddo
#endif
#ifdef STAR
        if(nas.gt.0) then
          do i=ng,nas+ng-1
            pn=list_ap(i)           
            pax_p(pn)=ax_p(pn)
            pay_p(pn)=ay_p(pn)
            paz_p(pn)=az_p(pn)
            ax_p(pn)=0.0d0
            ay_p(pn)=0.0d0
            az_p(pn)=0.0d0
            dvx_p(pn)=0.0d0
            dvy_p(pn)=0.0d0
            dvz_p(pn)=0.0d0
          enddo
        endif
#endif
#ifdef DM
        if(nadmt.gt.0) then
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ntime=MPI_WTIME()
#endif
#ifdef SGRAV
          call dmtreeforce(ndm,nadm,np,flags)
#endif

! *** Set dv in DM Particle ***
#ifdef MEXT
          do i=0,nadm-1
            pn=list_adm(i)
            rp=dsqrt(x_dm(pn)**2+y_dm(pn)**2+z_dm(pn)**2)
            if(rp.gt.10.0d0**SI_lro) then
              pot=mextr(SI_nmext-1)/(rp**3)
            else if(rp.gt.0.0d0) then
              ir=int((dlog10(rp)-SI_lri)/SI_dlr)
              if(ir.lt.0) then
                ir=0
              else if(ir.ge.SI_nmext-1) then
                ir=SI_nmext-2
              endif
              pot=(mextr(ir)+(rp-rmext(ir)) &
               *(mextr(ir+1)-mextr(ir))/(rmext(ir+1)-rmext(ir)))/(rp**3)
            else
              pot=0.0d0
            endif
            dvx_dm(pn)=dvx_dm(pn)-x_dm(pn)*pot
            dvy_dm(pn)=dvy_dm(pn)-y_dm(pn)*pot
            dvz_dm(pn)=dvz_dm(pn)-z_dm(pn)*pot
          enddo
#endif
#ifdef B85TEST
          do i=0,nadm-1
            pn=list_adm(i)
            rp=dsqrt(x_dm(pn)**2+y_dm(pn)**2+z_dm(pn)**2)
            if(rp.lt.rpert) then
              pot=mpert/(rpert**3)
            else
              pot=mpert/(rp**3)
            endif
            dvx_dm(pn)=dvx_dm(pn)-x_dm(pn)*pot
            dvy_dm(pn)=dvy_dm(pn)-y_dm(pn)*pot
            dvz_dm(pn)=dvz_dm(pn)-z_dm(pn)*pot
          enddo
#endif
#ifdef COSM
          do i=0,nadm-1
            pn = list_adm(i)
            dvx_dm(pn)=G*dvx_dm(pn)+SI_lam0*x_dm(pn)
            dvy_dm(pn)=G*dvy_dm(pn)+SI_lam0*y_dm(pn)
            dvz_dm(pn)=G*dvz_dm(pn)+SI_lam0*z_dm(pn)
          enddo
#endif
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ectime=MPI_WTIME()
          WTM_dmf=WTM_dmf+ectime-ntime
#endif
        endif
#endif

! *** TreeForce for Gas and Stars ***
#if defined(GAS) || defined(STAR)
        if(nagravgt.gt.0 .or. nast.gt.0) then
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ntime=MPI_WTIME()
#endif
#ifdef SGRAV
          call treeforce(ng,nagravg,np-ng,nas,ndm,flags)
#endif
#ifdef MEXT
          do i=0,nagravg-1
            pn=list_ap(i)
            rp=dsqrt(x_p(pn)**2+y_p(pn)**2+z_p(pn)**2)
            if(rp.gt.10.0d0**SI_lro) then
              pot=mextr(SI_nmext-1)/(rp**3)
            else if(rp.gt.0.0d0) then
              ir=int((dlog10(rp)-SI_lri)/SI_dlr)
              if(ir.lt.0) then
                ir=0
              else if(ir.ge.SI_nmext-1) then
                ir=SI_nmext-2
              endif
              pot=(mextr(ir)+(rp-rmext(ir)) &
               *(mextr(ir+1)-mextr(ir))/(rmext(ir+1)-rmext(ir)))/(rp**3)
            else
              pot=0.0d0
            endif
            ax_p(pn)=ax_p(pn)-x_p(pn)*pot
            ay_p(pn)=ay_p(pn)-y_p(pn)*pot
            az_p(pn)=az_p(pn)-z_p(pn)*pot
          enddo
#endif
#ifdef B85TEST
          do i=0,nagravg-1
            pn=list_ap(i)
            rp=dsqrt(x_p(pn)**2+y_p(pn)**2+z_p(pn)**2)
            if(rp.lt.rpert) then
              pot=mpert/(rpert**3)
            else
              pot=mpert/(rp**3)
            endif
            ax_p(pn)=ax_p(pn)-x_p(pn)*pot
            ay_p(pn)=ay_p(pn)-y_p(pn)*pot
            az_p(pn)=az_p(pn)-z_p(pn)*pot
          enddo
#endif
          do i=0,nagravg-1
            pn=list_ap(i)           
#ifdef COSM
            ax_p(pn)=G*ax_p(pn)+SI_lam0*x_p(pn)
            ay_p(pn)=G*ay_p(pn)+SI_lam0*y_p(pn)
            az_p(pn)=G*az_p(pn)+SI_lam0*z_p(pn)
#else
            ax_p(pn)=G*ax_p(pn)
            ay_p(pn)=G*ay_p(pn)
            az_p(pn)=G*az_p(pn)
#endif
          enddo
#ifdef STAR
! *** Set dv in Star Particle ***
          if(nas.gt.0) then
#ifdef MEXT
            do i=ng,ng+nas-1
              pn=list_ap(i)
              rp=dsqrt(x_p(pn)**2+y_p(pn)**2+z_p(pn)**2)
              if(rp.gt.10.0d0**SI_lro) then
                pot=mextr(SI_nmext-1)/(rp**3)
              else if(rp.gt.0.0d0) then
                ir=int((dlog10(rp)-SI_lri)/SI_dlr)
                if(ir.lt.0) then
                  ir=0
                else if(ir.ge.SI_nmext-1) then
                  ir=SI_nmext-2
                endif
                pot=(mextr(ir)+(rp-rmext(ir))*(mextr(ir+1)-mextr(ir)) &
                 /(rmext(ir+1)-rmext(ir)))/(rp**3)
              else
                pot=0.0d0
              endif
              ax_p(pn)=ax_p(pn)-x_p(pn)*pot
              ay_p(pn)=ay_p(pn)-y_p(pn)*pot
              az_p(pn)=az_p(pn)-z_p(pn)*pot
            enddo
#endif
            do i=ng,nas+ng-1
              pn=list_ap(i)           
#ifdef COSM
              ax_p(pn)=G*ax_p(pn)+SI_lam0*x_p(pn)
              ay_p(pn)=G*ay_p(pn)+SI_lam0*y_p(pn)
              az_p(pn)=G*az_p(pn)+SI_lam0*z_p(pn)
#else
              ax_p(pn)=G*ax_p(pn)
              ay_p(pn)=G*ay_p(pn)
              az_p(pn)=G*az_p(pn)
#endif
              dvx_p(pn)=ax_p(pn)
              dvy_p(pn)=ay_p(pn)
              dvz_p(pn)=az_p(pn)
            enddo
          endif
#endif
#ifdef WTIME
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          ectime=MPI_WTIME()
          WTM_bf=WTM_bf+ectime-ntime
#endif
        endif
! end for defined(GAS) || defined(STAR)
#endif


#if defined(GEPSLIM_VH) || defined(GEPSLIM_G)
#ifdef GAS
! *** recover original h ***
        do i=0,ng-1
          pn=list_ap(i)
          h_p(pn)=tx(i)
          omgh_p(pn)=ty(i)
          zetah_p(pn)=tz(i)
        enddo

        deallocate(tx)
        deallocate(ty)
        deallocate(tz)

#endif
#endif

      endif      

#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ntime=MPI_WTIME()
#endif
#ifdef WTIMESUB
      do i=1,4
        stime(i)=0.0d0
      enddo
#endif
! ***** Set du & dv in SPH Particle *****
#ifdef GAS
      if(flagu.ne.1) then
        do i=0,nag-1
          pdu_p(list_ap(i))=ndu_p(list_ap(i))
        enddo          
        do i=0,nag-1
          ndu_p(list_ap(i))=0.0d0
        enddo
      endif
      do i=0,nag-1
        d2u_p(list_ap(i))=0.0d0
      enddo
      if(flagv.eq.0) then      
        do i=0,nag-1
          pn=list_ap(i)         
          dvx_p(pn)=0.0d0
          dvy_p(pn)=0.0d0
          dvz_p(pn)=0.0d0
          hvsigdt_p(pn)=INF
        enddo

! allocate for working space for all the baryon particles
        allocate(peffp(0:np))

! *** set effective pressure ***
        do i=0,ng-1
          pn=list_ap(i)
          peffp(pn)=p_p(pn)
#ifdef UTHLIM_JM
! from eq. (1) of Hopkins et al. (2011)
          plim=1.2d0*(ULJMFAC**(2.0d0/3.0d0))*G*(h_p(pn)**2)*(rho_p(pn)**2)
          if(peffp(pn).lt.plim) then
            peffp(pn)=plim
          endif
#endif
        enddo
      endif
      npj=nag
! alocate for npj particles
      allocate(tdvr(0:npj*nval))
      allocate(flagproc(0:npj))
      allocate(flagcom(0:npj))
      allocate(pncomp(0:npj))
      allocate(talist(0:npj))
#if defined(SF_EFD) || defined(SF_ZFD)
      allocate(tivr(0:npj*nival))
#endif
#ifdef METALDIFF
      allocate(tdzd(0:npj*nvalzd))
      allocate(av(0:np))
#endif

! ****  Set Data  ****
      do i = 0,npj-1
        pni=list_ap(i)
        pn_nfp(i)=pni
        tdvr(i)=x_p(pni)
        tdvr(i+npj)=y_p(pni)
        tdvr(i+npj*2)=z_p(pni)             
        tdvr(i+npj*3)=vvnx_p(pni)
        tdvr(i+npj*4)=vvny_p(pni)
        tdvr(i+npj*5)=vvnz_p(pni)
        tdvr(i+npj*6)=h_p(pni)        
        tdvr(i+npj*7)=cs_p(pni)
        tdvr(i+npj*8)=rho_p(pni)        
        tdvr(i+npj*9)=div_v_p(pni)
        tdvr(i+npj*10)=arot_v_p(pni)
        tdvr(i+npj*11)=omgh_p(pni)
        tdvr(i+npj*12)=alpv_p(pni)        
        tdvr(i+npj*13)=alpu_p(pni)
        tdvr(i+npj*14)=ndu_p(pni)        
        tdvr(i+npj*15)=d2u_p(pni)
#if defined(SF_EFD) || defined(SF_ZFD)
! *** need this also for flagu=2 case
        tivr(i)=flagfd_p(pni)
#endif
        flagproc(i)=myrank
        flagcom(i)=-1
! *** keep the id(in proc) for particles ***
        talist(i)=pni
      enddo
      if(flagv.eq.0) then
        do i = 0,npj-1
          pni = pn_nfp(i)
          tdvr(i+npj*16)=dvx_p(pni)
          tdvr(i+npj*17)=dvy_p(pni)
          tdvr(i+npj*18)=dvz_p(pni)
          tdvr(i+npj*19)=hvsigdt_p(pni)
          tdvr(i+npj*20)=dt_p(pni)
          tdvr(i+npj*21)=peffp(pni)
#if defined(SF_EFD) || defined(SF_ZFD)
! *** need this also for flagu=2 case
          tivr(i+npj)=flagrfd_p(pni)
#endif
#ifdef METALDIFF
! *** for metal diffusion ***
          tdzd(i)=vsig_p(pni)
! *** initialization
          zdA_p(pni)=0.0d0
          zdHeB_p(pni)=0.0d0
          zdCB_p(pni)=0.0d0
          zdNB_p(pni)=0.0d0
          zdOB_p(pni)=0.0d0
          zdNeB_p(pni)=0.0d0
          zdMgB_p(pni)=0.0d0
          zdSiB_p(pni)=0.0d0
          zdFeB_p(pni)=0.0d0
          zdZB_p(pni)=0.0d0
#endif
        enddo
#ifdef METALDIFF
        do i=npj,npj*nvalzd-1
          tdzd(i)=0.0d0
        enddo
! *** set metal diffusion coefficient ***        
        do i=0,ng-1
          pn=list_ap(i)
#ifdef CSZDIFF
          av(pn)=ZDFCFAC*cs_p(pn)*rho_p(pn)*h_p(pn)
#else
          av(pn)=ZDFCFAC*vsig_p(pn)*h_p(pn)*rho_p(pn)
#endif
        enddo
#endif      
      endif

      ncomp=0
      ncompt=0
      iproc=0
      idivcom=0
      ndivcom=-1
      do i=0,nprocs-1
        npproc(i)=0
      enddo
      
! *** tree walk for all the proc ***
   70 if(np_gtr(0).eq.0) then
        goto 90
      endif
      if(iproc.ne.0.and.proc_gtr(0).ne.myrank) then
        goto 90
      endif
! * Initialization *

      allocate(list(0:npj))
      allocate(node(0:npj))

      do i=0,npj-1
        list(i)=i
        node(i)=0
      enddo
      nlist=npj

! * start tree walk *
   77 if(nlist.eq.0) then

        deallocate(list)
        deallocate(node)

        goto 90
      endif

! allocate temporary arrays should be smaller than npj
      allocate(xijp(0:nlist))
      allocate(yijp(0:nlist))
      allocate(zijp(0:nlist))
      allocate(rijp(0:nlist))
      allocate(listb(0:nlist))
      allocate(nblistb(0:nlist))

! *** nabp: number of particles who has neighbour in this loop
      nbp=0   
      do i=0,nlist-1
        pni=list(i)
        nd=node(pni)        
        xij=tdvr(pni)-cx_gtr(nd)
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
        yij=tdvr(pni+npj)-cy_gtr(nd)
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
        zij=tdvr(pni+npj*2)-cz_gtr(nd)
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
! *** from gather ***
        rh=tdvr(pni+npj*6)+crit*l_gtr(nd)
        if(rh.lt.hm_gtr(nd)) then
! *** from scatter ***
          rh=hm_gtr(nd)
        endif
        rh=rh**2
        if(np_gtr(nd).eq.1) then
        if(rc.lt.rh) then
        if(proc_gtr(nd).eq.myrank) then
        if(rc.gt.0.0d0) then
          pn=pn_gtr(nd)
          if(iproc.ne.0.or.talist(pni).ne.pn) then
#if defined(SF_EFD) || defined(SF_ZFD)
! *** set flag for the particles affected by feedback ***
            if(flagv.eq.0.and.tivr(pni).ne.0) then
              tivr(pni+npj)=1
              flagrfd_p(pn)=1
            endif
#endif

#if defined(SF_EFD) || defined(SF_ZFD)
            if(tivr(pni).le.0) then
#endif
              xijp(nbp)=xij  
              yijp(nbp)=yij  
              zijp(nbp)=zij   
              rijp(nbp)=dsqrt(rc)
              listb(nbp)=pni
              nblistb(nbp)=pn
              nbp=nbp+1
#if defined(SF_EFD) || defined(SF_ZFD)
            else if(flagu.eq.0.and.tivr(pni).gt.0) then
              r=dsqrt(rc)
              if(r.lt.tdvr(pni+npj*6)) then
! *** for gather ***
                hsi=tdvr(pni+npj*6)
                ihsi=1.0d0/hsi
                si=r*ihsi
                is=int(si/ds_tb)
                if(is.lt.0) then
                  is=0
                else if(is.ge.NKTAB) then
                  is=NKTAB-1
                endif
                mwi=w_tb(is)+(w_tb(is+1)-w_tb(is))*(si-s_tb(is))/ds_tb
! *** only 3D ***
                mwi=m_p(pn)*mwi*(ihsi**3)
! *** store u for new feebdack particles ***
                tdvr(pni+npj*14)=tdvr(pni+npj*14)+u_p(pn)*mwi
              endif
            endif
#endif
          endif
        endif
        endif
        endif
        endif
! * update node *
        if((np_gtr(nd).eq.1.and.proc_gtr(nd).eq.myrank).or.rc.ge.rh) then
          node(pni)=next_gtr(nd)
        else
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_gtr(nd).eq.-1) then
            if(proc_gtr(nd).ne.myrank.and.flagproc(pni).ne.proc_gtr(nd)) then
! *** proc_gtr will be the same for different domain within the same proc
              flagproc(pni)=proc_gtr(nd)
! *** store pn ***
              if(flagcom(pni).lt.0) then
                flagcom(pni)=proc_gtr(nd)
                pncomp(ncomp)=pn_nfp(pni)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNB) then
                pncompt(ncompt)=pn_nfp(pni)
                proccompt(ncompt)=proc_gtr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_gtr(nd))=npproc(proc_gtr(nd))+1
            endif
            node(pni)=next_gtr(nd)
          else
            node(pni)=daughter_gtr(nd)
          endif
        endif
      enddo

! for temporary array
      allocate(mdwij(0:nbp))
      allocate(mdwic(0:nbp))
      allocate(mdwjc(0:nbp))
      allocate(v_x_r(0:nbp))
      allocate(ihsip(0:nbp))
      allocate(ihsjp(0:nbp))

      do i=0,nbp-1
        pni=listb(i)
        pn=nblistb(i)
! *** for gather ***
        hsi=tdvr(pni+npj*6)
        ihsi=1.0d0/hsi
        si=rijp(i)*ihsi
! *** for scatter ***
        hsj=h_p(pn)
        ihsj=1.0d0/hsj
        sj=rijp(i)*ihsj
        if(si.lt.1.0d0) then
          is=int(si/ds_tb)
          if(is.lt.0) then
            is=0
          else if(is.ge.NKTAB) then
            is=NKTAB-1
          endif
! * calculate m_j dw_i/dr/r *
          mdwi=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
           *(si-s_tb(is))/ds_tb
! * calculate m_j dw_i/dr *
          mdwic(i)=dwdsc_tb(is)+(dwdsc_tb(is+1)-dwdsc_tb(is)) &
           *(si-s_tb(is))/ds_tb
          mdwic(i)=mdwic(i)/rijp(i)
#ifdef SIM1D
          mdwi=m_p(pn)*mdwi*(ihsi**3)
          mdwic(i)=m_p(pn)*mdwic(i)*(ihsi**2)
#elif defined(SIM2D)
          mdwi=m_p(pn)*mdwi*(ihsi**4)
          mdwic(i)=m_p(pn)*mdwic(i)*(ihsi**3)
#else
          mdwi=m_p(pn)*mdwi*(ihsi**5)
          mdwic(i)=m_p(pn)*mdwic(i)*(ihsi**4)
#endif
        else 
          mdwi=0.0d0
          mdwic(i)=0.0d0
        endif
! * calculate m_j dw_j/dr/r *
        if(sj.lt.1.0d0) then
          is=int(sj/ds_tb)
          if(is.lt.0) then
            is=0
          else if(is.ge.NKTAB) then
            is=NKTAB-1
          endif
          mdwj=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is))*(sj-s_tb(is))/ds_tb
          mdwjc(i)=dwdsc_tb(is)+(dwdsc_tb(is+1)-dwdsc_tb(is)) &
           *(sj-s_tb(is))/ds_tb
          mdwjc(i)=mdwjc(i)/rijp(i)
#ifdef SIM1D
          mdwj=m_p(pn)*mdwj*(ihsj**3)
          mdwjc(i)=m_p(pn)*mdwjc(i)*(ihsj**2)
#elif defined(SIM2D)
          mdwj=m_p(pn)*mdwj*(ihsj**4)
          mdwjc(i)=m_p(pn)*mdwjc(i)*(ihsj**3)
#else
          mdwj=m_p(pn)*mdwj*(ihsj**5)
          mdwjc(i)=m_p(pn)*mdwjc(i)*(ihsj**4)
#endif
        else 
          mdwj=0.0d0
          mdwjc(i)=0.0d0
        endif
        vxij=tdvr(pni+npj*3)-vvnx_p(pn)
#ifdef SIM1D
        vyij=0.0d0
        vzij=0.0d0
#else
        vyij=tdvr(pni+npj*4)-vvny_p(pn)
#ifdef SIM2D
        vzij=0.0d0
#else
        vzij=tdvr(pni+npj*5)-vvnz_p(pn)
#endif
#endif
! *** store the data only for particles who have neighbour
! *** v*xij/rij
        v_x_r(i)=(xijp(i)*vxij+yijp(i)*vyij+zijp(i)*vzij)/rijp(i)
! RP07
!                mdwij=0.5d0*(mdwi/tdvr(pni+npj*14)+mdwj/omgh_p(pn))
! no omgh for w since 28/09/09
! Price 08  0.5 mj (|dw_i/dr|/r +|dw_j/dr|/r )
        mdwij(i)=0.5d0*(mdwi+mdwj)
        ihsip(i)=ihsi
        ihsjp(i)=ihsj
      enddo

! *** gas particles who has nenighbours
      allocate(tz(0:nbp))
      allocate(upj(0:nbp))
      allocate(rhoij(0:nbp))
      allocate(vsigav(0:nbp))
      allocate(dvdiss(0:nbp))
      allocate(dudissv(0:nbp))
      allocate(dudissu(0:nbp))
      do i=0,nbp-1
! *  u=p/((gam-1)*rho), p=cs_p^2*rhop/gam -> u=cs^2/(gam*(gam-1))
        tz(i)=(tdvr(listb(i)+npj*7)**2)/(GAM*(GAM-1.0d0))
      enddo
      do i=0,nbp-1
! *** for j particles
        upj(i)=p_p(nblistb(i))/(rho_p(nblistb(i))*(GAM-1.0d0))
      enddo
      do i=0,nbp-1
        rhoij(i)=(tdvr(listb(i)+npj*8)+rho_p(nblistb(i)))*0.5d0
      enddo
! *** calculate d2u_p ***
      do i=0,nbp-1
        tdvr(listb(i)+npj*15)=tdvr(listb(i)+npj*15) &
         +dabs(mdwij(i))*(tz(i)-upj(i))/(rho_p(nblistb(i)))
      enddo
! * calculate alpha, rho, vsig for artificial viscosity 
! * and conductivity 
      do i=0,nbp-1
        if(v_x_r(i).lt.0.0d0) then
! *** artificial viscosity parameter ***
#ifdef NOAVSWITCH
          alpij=V_ALPHAINI
#else
          alpij=(tdvr(listb(i)+npj*12)+alpv_p(nblistb(i)))*0.25d0 &
           *(dabs(tdvr(listb(i)+npj*9))/(dabs(tdvr(listb(i)+npj*9)) &
           +SHAREV*tdvr(listb(i)+npj*10) &
           +0.0002d0*tdvr(listb(i)+npj*7)*ihsip(i)) &
           +dabs(div_v_p(nblistb(i))) &
           /(dabs(div_v_p(nblistb(i)))+SHAREV*arot_v_p(nblistb(i)) &
           +0.0002d0*cs_p(nblistb(i))*ihsjp(i)))
#endif
! *** artificial conductivity parameter ***
          alpuij=0.5d0*(tdvr(listb(i)+npj*13)+alpu_p(nblistb(i)))
! *** signal velocity ***
          vsigav(i)=(tdvr(listb(i)+npj*7)+cs_p(nblistb(i)) &
           -V_BETA*v_x_r(i))*0.5d0
          dvdiss(i)=(alpij*vsigav(i)*v_x_r(i)/rhoij(i))*mdwij(i)
          dudissv(i)=0.5d0*vsigav(i)*alpij*(v_x_r(i)**2)
#ifdef NOAC
          dudissu(i)=0.0d0
#else
          dudissu(i)=alpuij*vsigav(i)*(tz(i)-upj(i))
#endif
        else
! A.C. for positive v_x_r
! Price & Monaghan (2005)
          vsigav(i)=(tdvr(listb(i)+npj*7)+cs_p(nblistb(i)) &
           -V_BETA*v_x_r(i))*0.5d0
          if(vsigav(i).lt.0.0d0) then
            vsigav(i)=0.0d0
          endif
#ifdef NOAC
          alpuij=0.0d0
          dudissu(i)=0.0d0
#else
          alpuij=0.5d0*(tdvr(listb(i)+npj*13)+alpu_p(nblistb(i)))
          dudissu(i)=alpuij*vsigav(i)*(tz(i)-upj(i))
#endif
! set positive vsigav value even when v_x_r>0, not on since pv34.2
!          vsigav(i)=(tdvr(listb(i)+npj*7)+cs_p(nblistb(i))
!     &      +v_x_r(i))*0.5d0
! since pv33.6
          dvdiss(i)=0.0d0
          dudissv(i)=0.0d0
        endif
      enddo

      deallocate(tz)
      deallocate(upj)
      deallocate(v_x_r)
      deallocate(ihsip)
      deallocate(ihsjp)

      if(flagu.ne.1) then
!DEC$ IVDEP
        do i=0,nbp-1
! * Sum for du/dt (artificial viscosity and thermal conductivity ) *
          tdvr(listb(i)+npj*14)=tdvr(listb(i)+npj*14) &
           -(rijp(i)*mdwij(i)/rhoij(i))*(dudissv(i)-dudissu(i))
        enddo
      endif

      deallocate(rhoij)
      deallocate(dudissv)
      deallocate(dudissu)

! * Sum dv/dt for active particle (pressure gradient) *
      if(flagv.eq.0) then

        allocate(tx(0:nbp))
        allocate(mpqdwj(0:nbp))

        do i=0,nbp-1
! *** check for timestep ***
          if(tdvr(listb(i)+npj*20).lt.dtmnb_p(nblistb(i))) then
            dtmnb_p(nblistb(i))=tdvr(listb(i)+npj*20)
          endif
        enddo
        do i=0,nbp-1
! *  p/rho^2
          tx(i)=tdvr(listb(i)+npj*21)/(tdvr(listb(i)+npj*8)**2)
        enddo
        do i=0,nbp-1
! *** dvdiss pation term ***
! *** use core dw/dr only for pressure gradient since pv32.14
          mpqdwj(i)=(tx(i)*mdwic(i)/tdvr(listb(i)+npj*11) &
           +peffp(nblistb(i))*mdwjc(i)/(omgh_p(nblistb(i)) &
           *(rho_p(nblistb(i))**2))-dvdiss(i))
        enddo
        do i=0,nbp-1
          tdvr(listb(i)+npj*16)=tdvr(listb(i)+npj*16)-mpqdwj(i)*xijp(i)
        enddo
        do i=0,nbp-1
          tdvr(listb(i)+npj*17)=tdvr(listb(i)+npj*17)-mpqdwj(i)*yijp(i)
        enddo
        do i=0,nbp-1
          tdvr(listb(i)+npj*18)=tdvr(listb(i)+npj*18)-mpqdwj(i)*zijp(i)
        enddo
        do i=0,nbp-1
! calculate hvsigdt for time step constraints
! Price thesis
          if(vsigav(i).gt.0.0d0) then
            vsigdt=vsigav(i)
            if(h_p(nblistb(i)).gt.tdvr(listb(i)+npj*6)) then
              hvsigdt=0.5d0*(tdvr(listb(i)+npj*6)/vsigdt)
            else
              hvsigdt=0.5d0*(h_p(nblistb(i))/vsigdt)
            endif
            if(hvsigdt.lt.tdvr(listb(i)+npj*19)) then
              tdvr(listb(i)+npj*19)=hvsigdt
            endif
          endif
        enddo

! deallocate 
        deallocate(tx)
        deallocate(mpqdwj)

      endif

      deallocate(mdwic)
      deallocate(mdwjc)
      deallocate(dvdiss)
      deallocate(vsigav)

#ifdef METALDIFF
      if(flagv.eq.0) then

        allocate(kij(0:nbp))

        do i=0,nbp-1
#ifdef CSZDIFF
! cs rho h
          czdi=ZDFCFAC*tdvr(listb(i)+npj*8)*tdvr(listb(i)+npj*6) &
            *tdvr(listb(i)+npj*7)
#else
! calculate metal diffusion parameters: Greif et al. (2009)
! metal diffusion coefficient
! 2 rho vsig h -> f rho vsig h (since pv32.8 30/04/2010)
          czdi=ZDFCFAC*tdvr(listb(i)+npj*8)*tdvr(listb(i)+npj*6) &
           *tdzd(listb(i))
#endif
          kij(i)=mdwij(i)*4.0d0*czdi*av(nblistb(i)) &
           /((czdi+av(nblistb(i)))*tdvr(listb(i)+npj*8) &
           *rho_p(nblistb(i)))
        enddo
        do i=0,nbp-1
! *** for zdA_p
          tdzd(listb(i)+npj)=tdzd(listb(i)+npj)+kij(i)
        enddo
        do i=0,nbp-1
! *** zdHeB_p ... zdZB_p
          tdzd(listb(i)+npj*2)=tdzd(listb(i)+npj*2)+kij(i)*mzHe_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*3)=tdzd(listb(i)+npj*3)+kij(i)*mzC_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*4)=tdzd(listb(i)+npj*4)+kij(i)*mzN_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*5)=tdzd(listb(i)+npj*5)+kij(i)*mzO_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*6)=tdzd(listb(i)+npj*6)+kij(i)*mzNe_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*7)=tdzd(listb(i)+npj*7)+kij(i)*mzMg_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*8)=tdzd(listb(i)+npj*8)+kij(i)*mzSi_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*9)=tdzd(listb(i)+npj*9)+kij(i)*mzFe_p(nblistb(i))
        enddo
        do i=0,nbp-1
          tdzd(listb(i)+npj*10)=tdzd(listb(i)+npj*10)+kij(i)*mzZ_p(nblistb(i))
        enddo

        deallocate(kij)
      endif
#endif
     
      deallocate(mdwij)
      deallocate(xijp)
      deallocate(yijp)
      deallocate(zijp)
      deallocate(rijp)
      deallocate(listb)
      deallocate(nblistb)

! * update not-finished particle list *
      tnlist=nlist
      nlist=0
      if(iproc.eq.0.and.next_gtr(0).ne.0) then
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          else if(node(list(i)).eq.0) then
            node(list(i))=next_gtr(0)
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
! *** update variables ***
   90 if(iproc.ne.0) then
! *** sending back the results ***
! *** new snval to send the data back
        if(flagv.eq.0) then
          snval=6
        else
          snval=2
        endif

        allocate(tivs(0:nprocs))          

        do i=0,nprocs-1
! *** keep npjr: data received from each proc ***
          tivs(i)=npjr(i)
! *** ireqs=jjlen number of particles sent -> now receiving.
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
! *** setting the sending back data ***

        allocate(tdvs(0:npj*snval))

        do i=0,npj-1
! *** ndu 
          tdvs(snval*i)=tdvr(i+npj*14)
! *** d2u 
          tdvs(snval*i+1)=tdvr(i+npj*15)
        enddo
        if(flagv.eq.0) then
          do i=0,npj-1
            tdvs(snval*i+2)=tdvr(i+npj*16)
            tdvs(snval*i+3)=tdvr(i+npj*17)
            tdvs(snval*i+4)=tdvr(i+npj*18)
!  hvsigdt_p
            tdvs(snval*i+5)=tdvr(i+npj*19)
          enddo
        endif
! *** sending and reveiving the data ***
        irecv=0
        irecv1=0
! *** number of particles sent is ncomp 
        npjs=npj
! *** ncompi: number of particles originally sent from iproc=0
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
            tdvr(irecv+npj*14)=trbuf(snval*i)
            tdvr(irecv+npj*15)=trbuf(snval*i+1)
            irecv=irecv+1
          enddo
          if(flagv.eq.0) then
            do i=0,npjr(ip)-1
              tdvr(irecv1+npj*16)=trbuf(snval*i+2)
              tdvr(irecv1+npj*17)=trbuf(snval*i+3)
              tdvr(irecv1+npj*18)=trbuf(snval*i+4)
              tdvr(irecv1+npj*19)=trbuf(snval*i+5)
              irecv1=irecv1+1
            enddo
          endif

          deallocate(trbuf)

        enddo

        deallocate(tdvs)

#ifdef METALDIFF
        isend=0
        snval=10

        allocate(tdvs(0:npj*snval))

! *** setting the number of particles to send
        do i=0,nprocs-1
          idisp(i)=isend 
          jjlen(i)=tivs(i)*snval
          isend=isend+jjlen(i)          
        enddo
! *** setting the sending back data ***
        do i=0,npjs-1
          tdvs(snval*i)=tdzd(i+npjs)
          tdvs(snval*i+1)=tdzd(i+npjs*2)
          tdvs(snval*i+2)=tdzd(i+npjs*3)
          tdvs(snval*i+3)=tdzd(i+npjs*4)
          tdvs(snval*i+4)=tdzd(i+npjs*5)
          tdvs(snval*i+5)=tdzd(i+npjs*6)
          tdvs(snval*i+6)=tdzd(i+npjs*7)
          tdvs(snval*i+7)=tdzd(i+npjs*8)
          tdvs(snval*i+8)=tdzd(i+npjs*9)
          tdvs(snval*i+9)=tdzd(i+npjs*10)
        enddo
! *** sending and reveiving the data ***
        irecv=0
! *** ncompi: number of particles originally sent from iproc=0
        npj=ncompi

        deallocate(tdzd)
        allocate(tdzd(0:npj*nvalzd))

        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION &
           ,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tdzd(irecv+npj)=trbuf(snval*i)
            tdzd(irecv+npj*2)=trbuf(snval*i+1)
            tdzd(irecv+npj*3)=trbuf(snval*i+2)
            tdzd(irecv+npj*4)=trbuf(snval*i+3)
            tdzd(irecv+npj*5)=trbuf(snval*i+4)
            tdzd(irecv+npj*6)=trbuf(snval*i+5)
            tdzd(irecv+npj*7)=trbuf(snval*i+6)
            tdzd(irecv+npj*8)=trbuf(snval*i+7)
            tdzd(irecv+npj*9)=trbuf(snval*i+8)
            tdzd(irecv+npj*10)=trbuf(snval*i+9)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo

        deallocate(tdvs)

#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** for integer values ***
        if(flagv.eq.0) then
          isend=0
          snival=1
          do i=0,nprocs-1
            idisp(i)=isend 
! *** tivs is numbr of particles received from each proc
            jjlen(i)=tivs(i)*snival
            isend=isend+jjlen(i)          
          enddo
! *** setting the sending back data, use npjs ***

          deallocate(tivs)
          allocate(tivs(0:npjs*snival))

          do i=0,npjs-1
! *** flagrfd_p
            tivs(snival*i)=tivr(i+npjs)
          enddo
! *** sending and reveiving the data ***
          irecv=0
! *** ncompi: number of particles originally sent from iproc=0
          npj=ncompi

          deallocate(tivr)
          allocate(tivr(0:npj*nival))

          do ip=0,nprocs-1

            allocate(tibuf(0:npjr(ip)*snival))

            call MPI_SCATTERV(tivs,jjlen,idisp,MPI_INTEGER &
             ,tibuf,npjr(ip)*snival,MPI_INTEGER,ip &
             ,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
            do i=0,npjr(ip)-1
              tivr(irecv+npj)=tibuf(snival*i)
              irecv=irecv+1
            enddo

            deallocate(tibuf)

          enddo          
        endif
#endif

          deallocate(tivs)

      endif

      if(flagu.ne.1) then
        do i = 0,npj-1
          pni=pn_nfp(i)
! *** ndu_p is used for u_p for flagfd_p>0
          ndu_p(pni)=ndu_p(pni)+tdvr(i+npj*14)
        enddo
      endif
      do i = 0,npj-1
        pni=pn_nfp(i)
        d2u_p(pni)=d2u_p(pni)+tdvr(i+npj*15)
      enddo
      if(flagv.eq.0) then      
! *** only active particles ***
        do i=0,npj-1
          pni=pn_nfp(i)
          dvx_p(pni)=dvx_p(pni)+tdvr(i+npj*16)
          dvy_p(pni)=dvy_p(pni)+tdvr(i+npj*17)
          dvz_p(pni)=dvz_p(pni)+tdvr(i+npj*18)
          if(tdvr(i+npj*19).lt.hvsigdt_p(pni)) then
            hvsigdt_p(pni)=tdvr(i+npj*19)
          endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** original flagrfd_p may be already 1 ***
          if(tivr(i+npj).lt.flagrfd_p(pni)) then
            flagrfd_p(pni)=tivr(i+npj)
          endif
#endif
#ifdef METALDIFF
          zdA_p(pni)=zdA_p(pni)+tdzd(i+npj)
          zdHeB_p(pni)=zdHeB_p(pni)+tdzd(i+npj*2)
          zdCB_p(pni)=zdCB_p(pni)+tdzd(i+npj*3)
          zdNB_p(pni)=zdNB_p(pni)+tdzd(i+npj*4)
          zdOB_p(pni)=zdOB_p(pni)+tdzd(i+npj*5)
          zdNeB_p(pni)=zdNeB_p(pni)+tdzd(i+npj*6)
          zdMgB_p(pni)=zdMgB_p(pni)+tdzd(i+npj*7)
          zdSiB_p(pni)=zdSiB_p(pni)+tdzd(i+npj*8)
          zdFeB_p(pni)=zdFeB_p(pni)+tdzd(i+npj*9)
          zdZB_p(pni)=zdZB_p(pni)+tdzd(i+npj*10)
#endif
        enddo
      endif

! *** end of neighbour search ***
      if(nprocs.le.1.or.(iproc.ne.0.and.idivcom.gt.ndivcom)) then
        goto 91
      endif

      if(nprocs.gt.1.and.iproc.eq.0) then
        iproc=1
! *** check if need communication ***
        ncomptmax=0
        call MPI_ALLREDUCE(ncompt,ncomptmax,1,MPI_INTEGER &
         ,MPI_MAX,MPI_COMM_WORLD,ierr)
        if(ncomptmax.gt.0) then
! *** get how many particles each proc receives

          deallocate(tivr)
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

          if(npprocmax.gt.MNB) then
! *** number of times receiving and calculating for the other nodes
            ndivcom=int((npprocmax+nprocs*nprocs)/MNB)+1
          else
            ndivcom=1
          endif
          ndivcompt=0
          if(ncomptmax.gt.MNB) then
            ndivcompt=int((ncomptmax+nprocs*nprocs)/MNB)+1
            if(ndivcompt.gt.ndivcom) then
              ndivcom=ndivcompt
            endif
          endif
          idivcom=1
        else
! *** no communication required ***
          goto 91
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
          if(isend.ne.ncompt) then    
            write(6,*) ' Error in calc_dv_du(): isend.ne.ncomp'
            write(6,*) ' when counting N particles need communication'
            write(6,*) ' myrank,isend,ncomp=',myrank,isend,ncomp
            stop
          endif
        else
! *** set the range of number of particles for sending for each proc
          do i=0,nprocs-1
            call para_range(0,npproc(i)-1,ndivcom,idivcom-1 &
             ,isproc(i),ieproc(i))
            cproc(i)=0
          enddo
! *** preparation for sending the data to the other procs ***
          isend=0          
! gtree id for starting pseudo node
          if(proc_gtr(0).eq.myrank) then
            pnodess=next_gtr(0)
          else
! in case, if there is no local tree, but pseudo nodes
            pnodess=0
          endif
          ndp=num_gtr-pnodess-1
          if(allocated(ndpnd)) then
            deallocate(ndpnd)
            deallocate(idpnd)
            deallocate(xpnd)
            deallocate(ypnd)
            deallocate(zpnd)
            deallocate(lpnd)
            deallocate(hmpnd)
            deallocate(nextpnd)
            deallocate(daupnd)
          endif
          allocate(idpnd(0:num_gtr))
          allocate(ndpnd(0:ndp))
          allocate(xpnd(0:ndp))
          allocate(ypnd(0:ndp))
          allocate(zpnd(0:ndp))
          allocate(lpnd(0:ndp))
          allocate(hmpnd(0:ndp))
          allocate(nextpnd(0:ndp))
          allocate(daupnd(0:ndp))

! *** store particle list in the order of sending procs in list() ***
          do ip=0,nprocs-1
            idisp(ip)=isend           
            jjlen(ip)=0
            if(npproc(ip).gt.0.and.ip.ne.myrank) then
! *** get the coordinate and etc. for pseudo node for the proc
              ndp=0
              do nd=pnodess,num_gtr-1
                if(proc_gtr(nd).eq.ip) then
                  ndpnd(ndp)=nd
                  idpnd(nd)=ndp
                  xpnd(ndp)=cx_gtr(nd)
                  ypnd(ndp)=cy_gtr(nd)
                  zpnd(ndp)=cz_gtr(nd)
                  lpnd(ndp)=l_gtr(nd)
                  hmpnd(ndp)=hm_gtr(nd)
                  ndp=ndp+1
                endif
              enddo
              do i=0,ndp-1
                nd=ndpnd(i)                  
                if(next_gtr(nd).gt.pnodess &
                 .and.next_gtr(nd).le.ndpnd(ndp-1)) then
                  nextpnd(i)=idpnd(next_gtr(nd))      
                else
                  nextpnd(i)=ndp
                endif
                if(daughter_gtr(nd).ne.-1) then
                  daupnd(i)=idpnd(daughter_gtr(nd))
                else 
                  daupnd(i)=-1
                endif
              enddo
! *** search particles need communication with ip
              do i=0,ncomp-1
                pn=pncomp(i)
                nd=0
   73           xij=x_p(pn)-xpnd(nd)
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
                yij=y_p(pn)-ypnd(nd)
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
                zij=z_p(pn)-zpnd(nd)
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
                rh=h_p(pn)+crit*lpnd(nd)
! *** hm_gtr is not only maximum one but including crit*l_gtr
                if(rh.lt.hmpnd(nd)) then
! *** from gather ***
                  rh=hmpnd(nd)
                endif
                rh=rh*rh
                if(rc.lt.rh) then
                  if(daupnd(nd).eq.-1) then
                    if(cproc(ip).ge.isproc(ip).and.cproc(ip).le.ieproc(ip)) then
                      pn_nfp(isend)=pn
                      jjlen(ip)=jjlen(ip)+1
                      isend=isend+1
                    endif
                    cproc(ip)=cproc(ip)+1
                    goto 92
                  endif
                  nd=daupnd(nd)
                else
                  nd=nextpnd(nd)
                endif
                if(nd.ge.ndp) then
                  goto 92
                endif
                goto 73
   92         enddo

              if(cproc(ip).ne.npproc(ip)) then
                write(6,*) ' Error in calc_dvdu():npproc,cproc' &
                 ,',myrank,idiv,ip,ndp,ncompt,ndivcom,ncomp=' &
                 ,npproc(ip),cproc(ip),myrank,idivcom,ip,ndp &
                 ,ncompt,ndivcom,ncomp
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
          deallocate(hmpnd)
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
! *** update npj ***
          npj=npj+irecv
        enddo
        if(npj.gt.MNB) then
          write(6,*) ' Error in calc_dv_du(): npj > MNB'
          write(6,*) ' npj,MNB=',npj,MNB
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
        if(flagv.eq.0) then
#ifdef METALDIFF
          snval=17
#else
          snval=16
#endif
        else
          snval=14
        endif

        if(allocated(tdvs)) then
          deallocate(tdvs)
        endif
        allocate(tdvs(0:ncompi*snval))

        do i = 0,ncompi-1
          pni=pn_nfp(i)
! *** pn_nfp is not used in the other proc, so keep this ***
          tdvs(snval*i)=x_p(pni)
          tdvs(snval*i+1)=y_p(pni)
          tdvs(snval*i+2)=z_p(pni)             
          tdvs(snval*i+3)=vvnx_p(pni)
          tdvs(snval*i+4)=vvny_p(pni)
          tdvs(snval*i+5)=vvnz_p(pni)
          tdvs(snval*i+6)=h_p(pni)        
          tdvs(snval*i+7)=cs_p(pni)
          tdvs(snval*i+8)=rho_p(pni)        
          tdvs(snval*i+9)=div_v_p(pni)
          tdvs(snval*i+10)=arot_v_p(pni)
          tdvs(snval*i+11)=omgh_p(pni)
          tdvs(snval*i+12)=alpv_p(pni)        
          tdvs(snval*i+13)=alpu_p(pni)
        enddo
        if(flagv.eq.0) then
          do i=0,ncompi-1
            pni=pn_nfp(i)
            tdvs(snval*i+14)=dt_p(pni)
            tdvs(snval*i+15)=peffp(pni)
#ifdef METALDIFF
            tdvs(snval*i+16)=vsig_p(pni)
#endif
          enddo
        endif
! *** reset sending parameters ***
        do i=0,nprocs-1
          idisp(i)=idisp(i)*snval
          jjlen(i)=ireqs(i)*snval
        enddo

! reallocate tdvr and tivr to receive npj particles
        deallocate(tdvr)
        allocate(tdvr(0:npj*nval))
#if defined(SF_EFD) || defined(SF_ZFD)
        if(allocated(tivr)) then
          deallocate(tivr)
        endif
        allocate(tivr(0:npj*nival))
#endif

! *** sending and reveiving the data ***
        irecv=0
        irecv1=0
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
            tdvr(irecv+npj*6)=trbuf(snval*i+6)
            tdvr(irecv+npj*7)=trbuf(snval*i+7)
            tdvr(irecv+npj*8)=trbuf(snval*i+8)
            tdvr(irecv+npj*9)=trbuf(snval*i+9)
            tdvr(irecv+npj*10)=trbuf(snval*i+10)
            tdvr(irecv+npj*11)=trbuf(snval*i+11)
            tdvr(irecv+npj*12)=trbuf(snval*i+12)
            tdvr(irecv+npj*13)=trbuf(snval*i+13)
            irecv=irecv+1
          enddo

          if(flagv.eq.0) then
            do i=0,npjr(ip)-1
! *** dt_p
              tdvr(irecv1+npj*20)=trbuf(snval*i+14)
! *** peffp
              tdvr(irecv1+npj*21)=trbuf(snval*i+15)
#ifdef METALDIFF             
              tdzd(irecv1)=trbuf(snval*i+16)
#endif
              irecv1=irecv1+1
            enddo
          endif

          deallocate(trbuf)

        enddo       
        deallocate(tdvs)

        if(irecv.ne.npj) then
          write(6,*) ' Error in calc_dv_du(): irecv,npj=',irecv,npj
          write(6,*) ' after sending the data to the other pe'
          stop
        endif
! *** initialization ***
! *** setting 0 for ndu_p and d2u_p ***
        do i=npj*14,npj*16-1
          tdvr(i)=0.0d0
        enddo
        if(flagv.eq.0) then
! *** setting 0 for dvx, dvy, dvz ***
          do i=npj*16,npj*19-1
            tdvr(i)=0.0d0
          enddo
! *** setting hvsigdt ***
          do i=npj*19,npj*20-1
            tdvr(i)=INF
          enddo
#ifdef METALDIFF
          do i=npj,npj*nvalzd-1
            tdzd(i)=0.0d0
          enddo
#endif
        endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** for integer values ***
        isend=0
        snival=1
! *** reset sending parameters ***
        do i=0,nprocs-1
          idisp(i)=(idisp(i)/snval)*snival
          jjlen(i)=ireqs(i)*snival
        enddo
! *** setting the sending data ***
        allocate(tivs(0:ncompi*snival))
        do i = 0,ncompi-1
          pni=pn_nfp(i)
          tivs(snival*i)=flagfd_p(pni)
        enddo
! *** sending and reveiving the data ***
        irecv=0
! *** number of particles sent is ncomp
        do ip=0,nprocs-1
          if(allocated(tibuf)) then
            deallocate(tibuf)
          endif
          allocate(tibuf(0:npjr(ip)*snival))
          call MPI_SCATTERV(tivs,jjlen,idisp,MPI_INTEGER &
           ,tibuf,npjr(ip)*snival,MPI_INTEGER,ip &
           ,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
! *** flagfd_p
            tivr(irecv)=tibuf(snival*i)
            irecv=irecv+1
          enddo
          deallocate(tibuf)
          deallocate(tivs)
        enddo
        if(flagv.eq.0) then
! *** set flagrfd=0 ***
          do i=0,npj-1             
            tivr(i+npj)=0
          enddo
        endif
#endif
        goto 70
      endif

! *** Definition dv/dt ***
   91 if(flagv.eq.0) then

        deallocate(peffp)

! *** contribution from hydro ****
        do i=0,nag-1
          pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
          if(flagfd_p(pn).ne.0) then
            dvx_p(pn)=0.0d0
            dvy_p(pn)=0.0d0
            dvz_p(pn)=0.0d0
          endif
#endif
#if defined(FIXDVYZ) || defined(SIM1D)
          dvy_p(pn)=0.0d0
          dvz_p(pn)=0.0d0
#endif
#ifdef SIM2D
          dvz_p(pn)=0.0d0
#endif
        enddo
      endif
! *** factor 2 for d2u_p ***
      do i=0,nag-1
        pn = list_ap(i)            
        d2u_p(pn)=2.0d0*d2u_p(pn)
      enddo

      allocate(listb(0:nag))

! *** Definition of du/dt for n step ***
      if(flagu.ne.1) then
! since ver.31.15
        do i=0,nag-1
          pn = list_ap(i)            
#if defined(SF_EFD) || defined(SF_ZFD)
          if(flagfd_p(pn).le.0) then
#endif
          ndu_p(pn)=((GAM-1.d0)/(rho_p(pn)**(GAM-1.0d0)))*ndu_p(pn)
#if defined(SF_EFD) || defined(SF_ZFD)
          endif
#endif
        enddo
#ifdef COOL
        if(flagu.eq.0) then
          nlist=0
          do i=0,nag-1
            pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
            if(flagfd_p(pn).le.0) then
#endif
              listb(nlist)=pn
              nlist=nlist+1
#if defined(SF_EFD) || defined(SF_ZFD)
            else
              ram_p(pn)=0.0d0
            endif
#endif
          enddo

          zn = 1.0d0/SI_a-1.0d0
! reallocate to nlist length list from nblist(0:nag-1)
          allocate(list(0:nlist))
          do i=0,nlist-1
            list(i)=listb(i)
          enddo

          allocate(ram(0:np))
! this is done within cool
!          do i=0,np-1
!            ram(i)=0.0d0
!          enddo

          call cool(nlist,list,zn,np,ram)
! receive ram(0:nlist-1)
          do i=0,nlist-1
            pn=list(i)
            ram_p(pn)=ram(pn)
          enddo
          deallocate(list)
          deallocate(ram)
        endif
#endif
      endif
! *** set u,p,cs,as for new feedback particles ***
#if defined(SF_EFD) || defined(SF_ZFD)
      if(flagu.eq.0) then
        nlist=0
        do i=0,nag-1
          pn=list_ap(i)
          if(flagfd_p(pn).gt.0) then
            u_p(pn)=ndu_p(pn)/rho_p(pn)
            listb(nlist)=pn
            nlist=nlist+1
            p_p(pn)=(GAM-1.0d0)*rho_p(pn)*u_p(pn)
            cs_p(pn)=dsqrt(GAM*p_p(pn)/rho_p(pn))
            as_p(pn)=p_p(pn)/(rho_p(pn)**GAM)
            ndu_p(pn)=0.0d0
            pdu_p(pn)=0.0d0
            d2u_p(pn)=0.0d0
            hvsigdt_p(pn)=INF
#ifdef METALDIFF
            zdA_p(pn)=0.0d0
            zdHeB_p(pn)=0.0d0
            zdCB_p(pn)=0.0d0
            zdNB_p(pn)=0.0d0
            zdOB_p(pn)=0.0d0
            zdNeB_p(pn)=0.0d0
            zdMgB_p(pn)=0.0d0
            zdSiB_p(pn)=0.0d0
            zdFeB_p(pn)=0.0d0
            zdZB_p(pn)=0.0d0
#endif
          endif
        enddo
#ifdef COOL
! reallocate to nlist length list from nblist(0:nag-1)
        allocate(list(0:nlist))
        allocate(ram(0:np))
        do i=0,nlist-1
          list(i)=listb(i)
        enddo
! this is done within cool
!        do i=0,np-1
!          ram(i)=0.0d0
!        enddo

! *** set myu
        call cool(nlist,list,zn,np,ram)
        call cool(nlist,list,zn,np,ram)
        do i=0,nlist-1
          pn=list(i)
          ram_p(pn)=ram(pn)          
        enddo
        deallocate(list) 
        deallocate(ram)
#endif
      endif
#endif
      deallocate(listb)

! end of GAS
#endif      

#ifdef GLASS
#ifdef GAS
        do i=0,nag-1
          pn=list_ap(i)
          dvx_p(pn)=dvx_p(pn)+FGLASS*(dble(ran1(idum))-0.5d0) &
!           *dvx_p(pn)
           *(vnx_p(pn)/dt_p(pn))
          dvy_p(pn)=dvy_p(pn)+FGLASS*(dble(ran1(idum))-0.5d0) &
!           *dvy_p(pn)
           *(vny_p(pn)/dt_p(pn))
          dvz_p(pn)=dvz_p(pn)+FGLASS*(dble(ran1(idum))-0.5d0) &
!           *dvz_p(pn)
           *(vnz_p(pn)/dt_p(pn))
        enddo
#endif
#endif
#ifdef WTIME
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ectime=MPI_WTIME()
        WTM_calcvu=WTM_calcvu+ectime-ntime
#endif
#ifdef WTIMESUB
        estime=MPI_WTIME()
        stime(4)=stime(4)+(estime-nstime)
        if(myrank.eq.0) then
!                                    1234567890123456
          write(6,'(a16,5(1pE13.5),2I3,I10)') 'calc_dv_du time=' &
          ,ectime-ntime,stime(1),stime(2),stime(3),stime(4),flagv,flagu
        endif
#endif
      deallocate(idisp)
      deallocate(jjlen)
      deallocate(ireqs)
      deallocate(npproc)
      deallocate(isproc)
      deallocate(ieproc)
      deallocate(cproc)
      deallocate(npjr)

      deallocate(tdvr)
      deallocate(flagproc)
      deallocate(flagcom)
      deallocate(pncomp)
      deallocate(talist)
#if defined(SF_EFD) || defined(SF_ZFD)
      deallocate(tivr)
#endif
#ifdef METALDIFF
      deallocate(tdzd)
      deallocate(av)
#endif

!      write(fileo,'(a5,i3.3)') 'calvu',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ng-1
!        pni=list_ap(i)
!        write(60,'(21(1pE13.5),3I10)') x_p(pni),y_p(pni),z_p(pni) &
!         ,ndu_p(pni),dvx_p(pni),dvy_p(pni),dvz_p(pni),div_v_p(pni) &
!         ,arot_v_p(pni),d2u_p(pni),p_p(pni),u_p(pni),ax_p(pni),ay_p(pni) &
!         ,az_p(pni),rho_p(pni),omgh_p(pni),hvsigdt_p(pni),alpv_p(pni) &
!         ,ram_p(pni),myu_p(pni),id_p(pni),nnb_p(pni),flagfd_p(pni)
!      enddo
!      close(60)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      stop

end subroutine


