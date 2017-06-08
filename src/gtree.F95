#include "gcdp.def"
! *******************************************
!    gtree.F  for GCD+  ver. f03.0
!   30  Jan. 2013  written by D. Kawata
! *******************************************

! **************************************************
!    Definition of function about octtree
!  This program build octtree and
!    compute the center of mass and total mass and 
! **************************************************

#ifdef GAS
subroutine gtreebuild(ng,ns)
      use gcdp_const
      use gcdp_baryon
      use gcdp_gtree
      use gcdp_system

      implicit none
      include 'mpif.h'
    
      integer,intent(in) :: ng,ns
      integer i,j,pn,nd,pnd,level,npn,np,ierr
! * for finish tree (tr) *
      integer numgtr,maxntr
! * for not-finished particle (nfp) *
! * number of not-finished particle *
      integer numnfp
! * for subcell (sc) *
! * start subcell in tree *
      integer stsc
! * number of subcells *
      integer numsc
! * diameter of subcell *
      double precision l_sc,l_sc2,crit
! *** for calculating hm_gtr ***
      integer d,num,npare
      double precision hmnd
! * for work *
      integer numtmp
      integer ix,iy,iz
      integer,allocatable :: clist(:),slist(:),s_sc(:)
! *** for Paeno-Hilbert ordering ***
      integer,allocatable :: ixp(:),iyp(:),izp(:)
! *** to find root node for domain ***
      integer numtr0,nval,nvali,ip,nrnd,numtri
      integer numtrs
      integer noders,nodere
      character fileo*60
      double precision dx(0:1)
      data dx/1.0d0,-1.0d0/
! * for work for constructing tree
      integer,allocatable :: pn_nfp(:),nd_nfp(:),c_nfp(:),label_nfp(:),snfp(:)
      integer,allocatable :: list(:),talist(:),nalist(:),node(:)
      integer,allocatable :: np_sc(:),pare_tr(:),nd_sc(:),pare_sc(:),c_sc(:) &
       ,flag_pc(:)
      integer,allocatable :: tivr(:),npjr(:),nd0(:)
      double precision,allocatable :: tdvr(:)

! tree is always 3D
      crit=0.25d0*dsqrt(3.0d0)
      level=0

! *** Make root ***
      np=ng+ns

! an expected maximum number of tree
      maxntr=ng*3+nprocs*MAXNODESEND

! allocate memory for tree
      call allocate_gtree(maxntr)

! allocate list for keeping the top node id for each level
      allocate(list(0:maxntr-1))
      allocate(flag_pc(0:maxntr-1))
      allocate(pare_tr(0:maxntr-1))

! allocate the arrays for _nfp
      allocate(pn_nfp(0:np-1))
      allocate(nd_nfp(0:np-1))
      allocate(c_nfp(0:np-1))
      allocate(label_nfp(0:np-1))
      allocate(snfp(0:np-1))

      do i=0,ng+ns-1
        pn_nfp(i)=i
      enddo
! the last 0 for node 0
! gmakeroot use pn_nfp
      call gmakeroot(ng+ns,pn_nfp)
! not set np_gtr(0) in gmakeroot, but set = ng after gmakeroot,
! because gmakeroot use ng+ns
      np_gtr(0)=ng
      next_gtr(0)=0
      pare_tr(0)=-1
      daughter_gtr(0)=1
      numgtr=1
      if(np_gtr(0).eq.0) then
        numgtr=0
        goto 94
      endif
      numnfp = 0
      numtr0=numgtr
      do i=0,ng-1
        pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).le.0) then
#endif
          pn_nfp(numnfp)=pn
          nd_nfp(numnfp)=numtr0
          label_nfp(numnfp)=0
          snfp(numnfp)=0
          c_nfp(numnfp)=0
          numnfp=numnfp+1
#if defined(SF_EFD) || defined(SF_ZFD)
        endif
#endif
      enddo
! note that np_gtr(0) can be 0 even if numnfp=0, if there is only transition 
! (star->gas) particle, flagfd_p>0
      if(numnfp.eq.0) then
        np_gtr(0)=0
        numgtr=0
        goto 94
      else
        np_gtr(0)=numnfp
        pn_gtr(0)=pn_nfp(0)
      endif
! set root for numgtr=1
! np,l,cx,cy,cz, pn is setted in gmakeroot
      l_gtr(numgtr)=l_gtr(0)
      cx_gtr(numgtr)=cx_gtr(0)
      cy_gtr(numgtr)=cy_gtr(0)
      cz_gtr(numgtr)=cz_gtr(0)
      np_gtr(numgtr)=numnfp
      pn_gtr(numgtr)=pn_nfp(0)

      next_gtr(numgtr)=0
      pare_tr(numgtr)=0
      numgtr=numgtr+1
! *** subcell starts from 0
      numsc=8
      stsc=numgtr
      l_sc=l_gtr(numtr0)*0.5d0

      allocate(np_sc(0:numsc-1))
      allocate(pare_sc(0:numsc-1))
      allocate(c_sc(0:numsc-1))

      do i=0,7
        np_sc(i)=0
        pare_sc(i)=numtr0
        c_sc(i)=i
      enddo
      flag_pc(numtr0)=0
! * in subcell *
      if(numnfp.gt.1) then
        daughter_gtr(numtr0)=0
      else
        daughter_gtr(numtr0)=-1
        daughter_gtr(numtr0+1)=-1
      endif
      level = 0
      list(0)=numtr0
!      nd_sc(0) = 0
! * start iteration in level *
   77 if(numnfp.le.1) then
        goto 99
      endif
! * initialization for the case np_gtr < 8 *

      allocate(nd_sc(0:numgtr-1))

      do i=list(level),numgtr-1
        nd_sc(i)=0
      enddo
      level=level+1	  
      list(level) = numgtr
      l_sc2=l_sc
      l_sc = 0.5d0*l_sc
! * find out which subcell it is in *

      allocate(ixp(0:numnfp-1))
      allocate(iyp(0:numnfp-1))
      allocate(izp(0:numnfp-1))

      do i=0,numnfp-1
        pn=pn_nfp(i)
        nd=nd_nfp(i)
        if(x_p(pn)-cx_gtr(nd).ge.0.0d0) then
          ixp(i)=0 
        else
          ixp(i)=1
        endif 
        if(y_p(pn)-cy_gtr(nd).ge.0.0d0) then
          iyp(i)=0 
        else
          iyp(i)=1
        endif 
        if(z_p(pn)-cz_gtr(nd).ge.0.0d0) then
          izp(i)=0 
        else
          izp(i)=1
        endif 
      enddo

      call phcurven(0,numnfp-1,snfp(0:numnfp-1),c_nfp(0:numnfp-1) &
        ,ixp,iyp,izp,level)

      deallocate(ixp)
      deallocate(iyp)
      deallocate(izp)

      allocate(nalist(0:numnfp-1))
      allocate(talist(0:numnfp-1))
      npn = 0
      do i=0,numnfp-1
        if(np_gtr(nd_nfp(i)).gt.8) then
! nd_nfp: node in subcell
          nd_nfp(i)=daughter_gtr(nd_nfp(i))+c_nfp(i)
        else
! nd_nfp: node in _gtr
          nalist(npn)=nd_nfp(i)
          talist(npn)=i
          npn = npn+1
        endif
      enddo                

      do i=0,npn-1
! nd_nfp: node in subcell
        nd_nfp(talist(i))=daughter_gtr(nalist(i))+nd_sc(nalist(i))
        nd_sc(nalist(i))=nd_sc(nalist(i))+1
      enddo

      deallocate(nd_sc)

      do i=0,npn-1
! nd_nfp: node in subcell < numsc
        c_sc(nd_nfp(talist(i)))=c_nfp(talist(i))  
      enddo

      deallocate(nalist)
      deallocate(talist)

      allocate(s_sc(0:numsc-1))

! * update info of subcell *
      do i=0,numnfp-1
        np_sc(nd_nfp(i))=np_sc(nd_nfp(i))+1
        s_sc(nd_nfp(i))=snfp(i)
      enddo
      if(numgtr.gt.maxntr-numnfp) then
        write(6,*) ' Error in gbuildtree():Node is overflow!'
        write(6,*) ' This level is ',level,'numnfp=',numnfp
        write(6,*) ' rank=',myrank,' ng=',ng
!        do i=0,numnfp
        do i=0,10
          pn=pn_gtr(i)
          write(6,*) id_p(pn),x_p(pn),vnx_p(pn),dvx_p(pn) &
           ,ax_p(pn),ndu_p(pn),pdu_p(pn),u_p(pn),rho_p(pn) &
           ,h_p(pn),nnb_p(pn) &
#if defined(SF_EFD) || defined(SF_ZFD) 
           ,m_p(pn),flagfd_p(pn),alpv_p(pn),alpu_p(pn)
#else
           ,m_p(pn),alpv_p(pn),alpu_p(pn)
#endif
        enddo
        stop
      endif

! *** subcell is connected to tree ***
      allocate(talist(0:numsc-1))
      allocate(nalist(0:numsc-1))
      allocate(node(0:numsc-1))
      allocate(clist(0:numsc-1))
      allocate(slist(0:numsc-1))

      npn = 0
      do i=0,numsc-1
        if(np_sc(i).ge.1) then
          nalist(npn)=i
          talist(npn)=numgtr
          node(npn)=pare_sc(i)
          clist(npn)=c_sc(i)
          slist(npn)=s_sc(i)
          numgtr=numgtr+1
          npn=npn+1
        endif
      enddo

      deallocate(s_sc)

      do i=0,npn-1
! node=pare_sc: tree node id of parent node.
        if(flag_pc(node(i)).eq.0) then
! there is a case no particle in nalist(i)-1, but same parent
!          pare_sc(nalist(i)).ne.pare_sc(nalist(i)-1)) then
          daughter_gtr(node(i))=talist(i)
          flag_pc(node(i)) = 1
        endif
      enddo

! check
      do i=0,numsc-1
        if(np_gtr(pare_sc(i)).gt.1.and.flag_pc(pare_sc(i)).eq.0) then
          write(6,*) ' Error in gtree(): no daughter found for np>1'
          write(6,*) ' myrank,node,np,i=',myrank,pare_sc(i),np_gtr(pare_sc(i)) &
           ,i,np_sc(i)
          stop
        endif
      enddo

      allocate(nd_sc(0:numsc-1))

      do i=0,npn-1
! nd_sc: for subsell < numsc
        nd_sc(nalist(i)) = talist(i)
        np_gtr(talist(i))=np_sc(nalist(i))
        l_gtr(talist(i))=l_sc2
        pare_tr(talist(i))=node(i)
        call phixyzp(slist(i),clist(i),ix,iy,iz,level)
        cx_gtr(talist(i))=dx(ix)*l_sc+cx_gtr(node(i))
        cy_gtr(talist(i))=dx(iy)*l_sc+cy_gtr(node(i))
        cz_gtr(talist(i))=dx(iz)*l_sc+cz_gtr(node(i))
      enddo

      deallocate(talist)
      deallocate(nalist)
      deallocate(node)
      deallocate(clist)
      deallocate(slist)

! *** Set label not-finished particle ***
      do i=0,numnfp-1
! d_nfp: position in _sc -> _gtr
        nd_nfp(i)=nd_sc(nd_nfp(i))
        pn_gtr(nd_nfp(i))=pn_nfp(i)
! *  this node is leaf *		
        if(np_gtr(nd_nfp(i)).eq.1) then
          label_nfp(i)=1
        endif
      enddo

      deallocate(nd_sc)

! *** rebuild not finished particle list ***
      numtmp = 0
      do i=0,numnfp-1
        if(label_nfp(i).eq.0) then
          pn_nfp(numtmp)=pn_nfp(i)
          nd_nfp(numtmp)=nd_nfp(i)
! *** update c and s
          snfp(numtmp)=snfp(i)
          c_nfp(numtmp)=c_nfp(i)
          label_nfp(numtmp)=0
          numtmp=numtmp+1
        endif
      enddo
      numnfp=numtmp

! *** rebuild subcell ***
      numtmp=numgtr-stsc

      deallocate(np_sc)
      deallocate(pare_sc)
      deallocate(c_sc)
      if(numtmp.gt.0) then
        allocate(np_sc(0:numtmp*8-1))
        allocate(pare_sc(0:numtmp*8-1))
        allocate(c_sc(0:numtmp*8-1))
      endif

      numsc = 0
      do i=0,numtmp-1
        if(np_gtr(stsc).ge.2) then
          daughter_gtr(stsc)=numsc
          flag_pc(stsc)=0

          np_sc(numsc)=0
          pare_sc(numsc)=stsc
          c_sc(numsc) = 0

          np_sc(numsc+1)=0
          pare_sc(numsc+1)=stsc
          c_sc(numsc+1) = 1

          np_sc(numsc+2)=0
          pare_sc(numsc+2)=stsc
          c_sc(numsc+2) = 2

          np_sc(numsc+3)=0
          pare_sc(numsc+3)=stsc
          c_sc(numsc+3) = 3

          np_sc(numsc+4)=0
          pare_sc(numsc+4)=stsc
          c_sc(numsc+4) = 4

          np_sc(numsc+5)=0
          pare_sc(numsc+5)=stsc
          c_sc(numsc+5) = 5

          np_sc(numsc+6)=0
          pare_sc(numsc+6)=stsc
          c_sc(numsc+6) = 6

          np_sc(numsc+7)=0
          pare_sc(numsc+7)=stsc
          c_sc(numsc+7) = 7

          numsc=numsc+8
        else if(np_gtr(stsc).eq.1) then
          daughter_gtr(stsc)=-1
        endif
        stsc=stsc+1
      enddo
      stsc=numgtr
      goto 77


   99 if(allocated(np_sc)) then
        deallocate(np_sc)
        deallocate(pare_sc)
        deallocate(c_sc)
      endif
      deallocate(flag_pc)
      deallocate(pn_nfp)
      deallocate(nd_nfp)
      deallocate(c_nfp)
      deallocate(label_nfp)
      deallocate(snfp)

! *** set next node ***      
      next_gtr(numtr0)=0
      if(numgtr.gt.1) then
        allocate(pare_sc(0:numgtr-2))
        do i=numtr0+1,numgtr-2 
          pare_sc(i)=pare_tr(i+1)
        enddo		  
        do i=numtr0+1,numgtr-2
          if(pare_sc(i).eq.pare_tr(i)) then
            next_gtr(i)=i+1
          else
            next_gtr(i) = next_gtr(pare_tr(i))
          endif
        enddo
        next_gtr(numgtr-1) = next_gtr(pare_tr(numgtr-1))
        deallocate(pare_sc)
      endif

! *** set hm_gtr ***
      do i=numtr0,numgtr-1
        if(np_gtr(i).eq.1) then
          pn=pn_gtr(i)
! *** if np=1, cx should be the position of the particle ***
          cx_gtr(i)=x_p(pn)
          cy_gtr(i)=y_p(pn)
          cz_gtr(i)=z_p(pn)
          l_gtr(i)=0.0d0
! *** no need for adding cx-cmx term, because cx_gtr set to x_p
          hm_gtr(i)=h_p(pn)
        else
          hm_gtr(i) = 0.0d0
        endif
      enddo

      num = numgtr-1
      do j=level,0,-1

        allocate(pare_sc(0:num-list(j)))
        allocate(nd_sc(0:num-list(j)))

        npare = 0
        do i=num,list(j),-1
          if(pare_tr(i).ne.pare_tr(i-1)) then
            pare_sc(npare) = pare_tr(i)
            nd_sc(npare) = i
            npare = npare+1
          endif
        enddo
        do d=0,7
          do i=0,npare-1
            nd=nd_sc(i)+d
            pnd=pare_sc(i)
            if(pnd.eq.pare_tr(nd).and.nd.le.num) then
              if(np_gtr(nd).eq.1) then
                hmnd=hm_gtr(nd)+crit*dsqrt((cx_gtr(nd)-cx_gtr(pnd))**2 &
                 +(cy_gtr(nd)-cy_gtr(pnd))**2+(cz_gtr(nd)-cz_gtr(pnd))**2)
              else
                hmnd=hm_gtr(nd)+crit*l_gtr(pnd)
              endif
              if(hmnd.gt.hm_gtr(pnd)) then
                hm_gtr(pnd)= hmnd
              endif
            endif
          enddo
        enddo


        deallocate(pare_sc)
        deallocate(nd_sc)

        num=list(j)-1
      enddo

! *** send the data to the other nodes ***
! *** get number of subdomains from each proc ***
   91 do i=0,numgtr-1
! *** set procid ***
        proc_gtr(i)=myrank
      enddo

! *** adjust cx,cy,cz,l_tr ***
      num = numgtr-1
      do j=level-1,0,-1
        do i=num,list(j),-1
          if(daughter_gtr(i).gt.0.and.(np_gtr(i) &
            .eq.np_gtr(daughter_gtr(i)))) then
            nd=daughter_gtr(i)
            l_gtr(i)=l_gtr(nd)
            cx_gtr(i)=cx_gtr(nd)
            cy_gtr(i)=cy_gtr(nd)
            cz_gtr(i)=cz_gtr(nd)   
            hm_gtr(i)=hm_gtr(nd)
            daughter_gtr(i)=daughter_gtr(nd)
          endif
        enddo
        num=list(j)-1
      enddo
      hm_gtr(0)=hm_gtr(1)
      daughter_gtr(0)=daughter_gtr(1)
      np_gtr(0)=np_gtr(1)
      pn_gtr(0)=pn_gtr(1)
      l_gtr(0)=l_gtr(1)
      cx_gtr(0)=cx_gtr(1)
      cy_gtr(0)=cy_gtr(1)
      cz_gtr(0)=cz_gtr(1)

! *** nodess_gtr should be starting from 1, if numgtr.eq.0, it will be corrected later
   94 nodess_gtr=1
      if(nprocs.gt.1) then  
        numtr0=numgtr
        numtrs=0
! *** list(0)=1=numtr0 ***
        i=0
! list(level) is the first node for the last level 
! if level=1, this do loop will be skipped, but nodess_gtr=1
        do i=0,level-2
!          if(nodess_gtr.eq.1.and.np_gtr(list(i)).ne.np_gtr(list(i+1))) then
!            nodess_gtr=list(i)
!          endif
          if(list(i+2)-nodess_gtr+1.gt.MAXNODESEND) then
! *** sending up to level i+1
            goto 92
          endif
        enddo
        if(i.eq.level-1.and.numgtr-nodess_gtr+1.gt.MAXNODESEND) then
          i=level-1
! *** sending up to level-1 ***
          goto 92
        endif
! *** sending all the node ***
        nodese_gtr=numgtr-1
        goto 93
   92   nodese_gtr=list(i+1)-1

!        write(6,*) ' gtree: myrank,nodess,nodese,np=' &
!          ,myrank,nodess_gtr,nodese_gtr,np

   93   if(numgtr.eq.0) then
          numtrs=0
          nodess_gtr=0
          nodese_gtr=0
        else if(nodess_gtr.gt.nodese_gtr) then
          write(6,*) ' Error in gbuildtree(): nodess > nodese'
          write(6,*) ' nodess,nodese=',nodess_gtr,nodese_gtr
          write(6,*) ' rank=',myrank,' ng,numgtr=',ng,numgtr
          stop
        else
          numtrs=nodese_gtr-nodess_gtr+1
        endif

        allocate(npjr(0:nprocs-1))

        do ip=0,nprocs-1
          if(ip.eq.myrank) then
            npjr(ip)=numtrs
          endif
          call MPI_BCAST(npjr(ip),1,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
!        write(6,*) ' npjr in gtree=',npjr(ip),myrank
        enddo

        allocate(nd0(0:maxntr-1))

! *** get number of particles ***
        nvali=4
        nval=5
        do ip=0,nprocs-1
          if(ip.eq.myrank) then
            allocate(tivr(0:numtrs*nvali-1))
            allocate(tdvr(0:numtrs*nval-1))
            do i=0,numtrs-1
              nd=nodess_gtr+i
              tivr(i)=np_gtr(nd)
              tivr(i+numtrs)=next_gtr(nd)
              tivr(i+numtrs*2)=daughter_gtr(nd)
              if(tivr(i+numtrs*2).gt.nodese_gtr) then
                tivr(i+numtrs*2)=-1
              endif
              tivr(i+numtrs*3)=nd
! *** double precision values for neighbour search ***
              tdvr(i)=l_gtr(nd)
              tdvr(i+numtrs)=cx_gtr(nd)
              tdvr(i+numtrs*2)=cy_gtr(nd)
              tdvr(i+numtrs*3)=cz_gtr(nd)
              tdvr(i+numtrs*4)=hm_gtr(nd)
            enddo
          else
            allocate(tivr(0:npjr(ip)*nvali-1))
            allocate(tdvr(0:npjr(ip)*nval-1))
          endif
          call MPI_BCAST(tivr,npjr(ip)*nvali,MPI_INTEGER,ip &
           ,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(tdvr,npjr(ip)*nval,MPI_DOUBLE_PRECISION &
           ,ip,MPI_COMM_WORLD,ierr)
! even if numgtr=0 np_gtr(0)=0, there may be a new feedback particle and ng>=0
          if(ip.ne.myrank.and.ng.gt.0) then
! *** add tree ***
            nrnd=npjr(ip)
            numtri=numgtr
            do i=0,npjr(ip)-1
! even if numgtr=0, it should work.
              np_gtr(numgtr)=tivr(i)
              pn_gtr(numgtr)=0
! *** set the other parameter for node ***
              proc_gtr(numgtr)=ip
              next_gtr(numgtr)=tivr(i+nrnd)
              daughter_gtr(numgtr)=tivr(i+nrnd*2)
! *** link between original node id and id in this proc ***
              nd0(tivr(i+nrnd*3))=numgtr
              pare_tr(numgtr)=0
! *** adding nodes from the other procs
              l_gtr(numgtr)=tdvr(i)
              cx_gtr(numgtr)=tdvr(i+nrnd)
              cy_gtr(numgtr)=tdvr(i+nrnd*2)
              cz_gtr(numgtr)=tdvr(i+nrnd*3)
              hm_gtr(numgtr)=tdvr(i+nrnd*4)
              numgtr=numgtr+1
            enddo
            if(nrnd.gt.0) then
              noders=tivr(nrnd*3)
              nodere=tivr(nrnd*4-1)
            else
              noders=numgtr
              nodere=numgtr
            endif

            do i=numtri,numgtr-1
! *** set next_gtr ***
              if(next_gtr(i).ge.noders.and.next_gtr(i).le.nodere) then
                next_gtr(i)=nd0(next_gtr(i))
              else
                next_gtr(i)=numgtr
              endif
! *** set daughter_gtr ***     
              if(daughter_gtr(i).ne.-1) then
                daughter_gtr(i)=nd0(daughter_gtr(i))
              endif
            enddo
          endif

          deallocate(tivr)
          deallocate(tdvr)
      
        enddo

        deallocate(npjr)
        deallocate(nd0)

! *** check there is any imported nodes if yes, change next_gtr(0) ***
        if(numtr0.ne.numgtr) then
! numtr0: number of the local tree
          if(numtr0.ne.0) then
            next_gtr(0)=numtr0
          endif
! *** set end of tree ***
          do i=numtri,numgtr-1
! *** set next_tr ***
            if(next_gtr(i).eq.numgtr) then
              next_gtr(i)=-1
            endif
          enddo
        endif
      endif

! *** set common 
      if(numgtr.gt.maxntr) then
        write(6,*) &
         ' Error in gbuildtree():Node is overflow after combined!'
        write(6,*) ' rank=',myrank,' ng=',ng,'numgtr=',numgtr
        stop
      endif

!      write(fileo,'(a5,i3.3,a1,i3.3)') 'gtree',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,numgtr-1
!        if(np_gtr(i).eq.1.and.proc_gtr(i).eq.myrank) then
!        write(60,'(6I10,5(1pE13.5))') i,np_gtr(i),pn_gtr(i) &
!         ,next_gtr(i),daughter_gtr(i),proc_gtr(i),l_gtr(i) &
!        ,cx_gtr(i),cy_gtr(i),cz_gtr(i),hm_gtr(i)
!        endif
!      enddo
!      close(60)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      stop


      num_gtr=numgtr

      deallocate(list)
      deallocate(pare_tr)

end subroutine
	  
! *** Definition of gmakeroot() ***
subroutine gmakeroot(ng,pn_nfp)
      use gcdp_const
      use gcdp_baryon
      use gcdp_gtree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: ng
      integer,intent(in) :: pn_nfp(0:ng-1)
      integer i,pn,ierr
! * max coordinate *	  
      double precision max_x,max_y,max_z
      double precision min_x,min_y,min_z
! * max,temp length *	  
      double precision maxl,tl,tdvr(0:2),tdvs(0:2)

! *** Define root node ***
      max_x=-INF
      max_y=-INF
      max_z=-INF
      min_x=INF
      min_y=INF
      min_z=INF
      do i=0,ng-1
        pn=pn_nfp(i)
        if(x_p(pn).lt.min_x) then
          min_x = x_p(pn)
        endif
        if(y_p(pn).lt.min_y) then
          min_y = y_p(pn)
        endif
        if(z_p(pn).lt.min_z) then
          min_z = z_p(pn)
        endif
        if(x_p(pn).gt.max_x) then
          max_x = x_p(pn)
        endif
        if(y_p(pn).gt.max_y) then
          max_y = y_p(pn)
        endif	
        if(z_p(pn).gt.max_z) then
          max_z = z_p(pn)
        endif
      enddo
! *** get the maximum and minimum for all the particles ***
! not get largest box size since pv33.2
! *** maximum ***
      tdvs(0)=max_x
      tdvs(1)=max_y
      tdvs(2)=max_z
      call MPI_ALLREDUCE(tdvs,tdvr,3,MPI_DOUBLE_PRECISION &
        ,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_x=tdvr(0)
      max_y=tdvr(1)
      max_z=tdvr(2)
! *** minimum ***
      tdvs(0)=min_x
      tdvs(1)=min_y
      tdvs(2)=min_z
      call MPI_ALLREDUCE(tdvs,tdvr,3,MPI_DOUBLE_PRECISION &
        ,MPI_MIN,MPI_COMM_WORLD,ierr)
      min_x=tdvr(0)
      min_y=tdvr(1)
      min_z=tdvr(2)
! *** check x range ***
      tl=max_x-min_x
      if(tl.lt.0.0d0) then
        tl = -tl
      endif
      maxl=tl
! *** check y range ***
      tl=max_y-min_y
      if(tl.lt.0.0d0) then
        tl = -tl
      endif
      if(tl.gt.maxl) then
        maxl = tl
      endif
! *** check z range ***
      tl=max_z - min_z
      if(tl.lt.0.0d0) then
        tl = -tl
      endif
      if(tl.gt.maxl) then
        maxl = tl
      endif
! *** Set root node ***
      l_gtr(0)=MGROOT*maxl
      cx_gtr(0)=(max_x+min_x)*0.5d0
      cy_gtr(0)=(max_y+min_y)*0.5d0
      cz_gtr(0)=(max_z+min_z)*0.5d0
      pn_gtr(0)=list_ap(0)

end subroutine
#endif
