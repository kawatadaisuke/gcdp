#include "gcdp.def"
! ***************************************
!     ldmtree.F95  for GCD+ ver f03.0
! copy from dmtree.F and changing 0,np-1 -> SI_ndm1,np-1
!  3 Feb. 2013   written by D. Kawata
!  *************************************** 

#if defined(DM) && defined(MULTI_LEV)
subroutine ldmtreebuild(np)      
      use gcdp_const
      use gcdp_dm
      use gcdp_ldmtree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np

      integer pn,nd 
      integer i,ierr
      integer level,npn
! * for finished tree (dmtr) *
      integer numdmtr,maxntr
! * for not-finished particle (nfp) *
! * number of not-finished particle *      
      integer numnfp
! * for subcell (sc) *
! * start subcell in tree *      
      integer stsc
! * number of subcells *      
      integer numsc
! * diameter of subcell *      
      double precision l_sc,l_sc2
! * for work *
      integer numtmp
! *** for calculating hm_ldmtr ***
      integer d,num,npare,pnd,j
      double precision xsd,ysd,zsd
! *** for Paeno-Hilbert ordering ***
      integer ix,iy,iz
! *** to find root node for domain ***
      integer numtr0,nval,nvali,ip,nrnd,numtri
      integer noders,nodere,numtrs
      character fileo*60
      double precision dx(0:1)
      data dx/1.0d0,-1.0d0/
! for allocatable array
      integer,allocatable :: clist(:),slist(:),s_sc(:),lvlist(:),paretrcm(:)
! *** for Paeno-Hilbert ordering ***
      integer,allocatable :: ixp(:),iyp(:),izp(:)
      integer,allocatable :: pn_nfp(:),nd_nfp(:),c_nfp(:),label_nfp(:),snfp(:)
      integer,allocatable :: list(:),talist(:),nalist(:),node(:)
      integer,allocatable :: np_sc(:),pare_tr(:),nd_sc(:),pare_sc(:),c_sc(:) &
       ,flag_pc(:)
      integer,allocatable :: tivr(:),npjr(:),nd0(:)
      double precision,allocatable :: tdvr(:)

! an expected maximum number of tree
      maxntr=np*3+nprocs*MAXNODESEND

! allocate memory for tree
      call allocate_ldmtree(maxntr)

! *** Make root ***
! allocate list for keeping the top node id for each level
      allocate(list(0:maxntr))
      allocate(flag_pc(0:maxntr))
      allocate(pare_tr(0:maxntr))
! allocate the arrays for _nfp
      allocate(pn_nfp(0:np))
      allocate(nd_nfp(0:np))
      allocate(c_nfp(0:np))
      allocate(label_nfp(0:np))
      allocate(snfp(0:np))

      do i=SI_ndm1,np-1
        pn_nfp(i-SI_ndm1)=i
      enddo
      level=0
      call ldmmakeroot(np)
      next_ldmtr(0)=0
      pare_tr(0)=-1
      daughter_ldmtr(0)=1
      numdmtr=1
      if(np-SI_ndm1.eq.0) then
        np_ldmtr(0)=0
        numdmtr=0
        cmx_ldmtr(1)=0.0d0
        cmy_ldmtr(1)=0.0d0
        cmz_ldmtr(1)=0.0d0
        hm_ldmtr(1)=0.0d0
#ifdef QPTREEF
        mx_ldmtr(1)=0.0d0
        my_ldmtr(1)=0.0d0
        mz_ldmtr(1)=0.0d0
        mxx_ldmtr(1)=0.0d0
        myy_ldmtr(1)=0.0d0
        mzz_ldmtr(1)=0.0d0
        mxy_ldmtr(1)=0.0d0
        myz_ldmtr(1)=0.0d0
        mzx_ldmtr(1)=0.0d0
#endif
        goto 94
      endif

      numnfp=0
      numtr0=numdmtr
      do i=SI_ndm1,np-1
        pn_nfp(numnfp)=i
        nd_nfp(numnfp)=numtr0
        label_nfp(numnfp)=0
        snfp(numnfp)=0
        c_nfp(numnfp)=0
        numnfp=numnfp+1
      enddo
! *** set root for numtr=1 ***
      l_ldmtr(numdmtr)=l_ldmtr(0)
      cx_ldmtr(numdmtr)=cx_ldmtr(0)
      cy_ldmtr(numdmtr)=cy_ldmtr(0)
      cz_ldmtr(numdmtr)=cz_ldmtr(0)
      np_ldmtr(numdmtr)=np-SI_ndm1
      pn_ldmtr(numdmtr)=pn_nfp(0)
      next_ldmtr(numdmtr)=0
      pare_tr(numdmtr)=0
      numdmtr=numdmtr+1
! *** subcell starts from 0
      numsc=8
      stsc=numdmtr
      l_sc=l_ldmtr(numtr0)*0.5d0

      allocate(np_sc(0:numsc))
      allocate(pare_sc(0:numsc))
      allocate(c_sc(0:numsc))

      do i=0,7
        np_sc(i)=0
        pare_sc(i)=numtr0
        c_sc(i)=i
      enddo
      flag_pc(numtr0)=0
! * in subcell *      
      if(numnfp.gt.1) then
        daughter_ldmtr(numtr0)=0 
      else
        daughter_ldmtr(numtr0)=0 
        daughter_ldmtr(numtr0+1)=-1
      endif
      level=0
      list(0)=numtr0
! * start iteration in level *
   77 if(numnfp.le.1) then
        goto 99
      endif
! * initialization for the case np_ldmtr < 8 *

      allocate(nd_sc(0:numdmtr))

      do i=list(level),numdmtr-1
        nd_sc(i)=0
      enddo
      level=level+1	  
      list(level)=numdmtr
      l_sc2=l_sc
      l_sc=0.5d0*l_sc2
! * find out which subcell it is in *

      allocate(ixp(0:numnfp))
      allocate(iyp(0:numnfp))
      allocate(izp(0:numnfp))

      do i=0,numnfp-1
        pn=pn_nfp(i)
        nd=nd_nfp(i)
        if(x_dm(pn)-cx_ldmtr(nd).ge.0.0d0) then
          ixp(i)=0 
        else
          ixp(i)=1
        endif 
        if(y_dm(pn)-cy_ldmtr(nd).ge.0.0d0) then
          iyp(i)=0 
        else
          iyp(i)=1
        endif 
        if(z_dm(pn)-cz_ldmtr(nd).ge.0.0d0) then
          izp(i)=0 
        else
          izp(i)=1
        endif 
        call phcurve(snfp(i),c_nfp(i),ixp(i),iyp(i),izp(i),level)
      enddo

      deallocate(ixp)
      deallocate(iyp)
      deallocate(izp)

      allocate(nalist(0:numnfp))
      allocate(talist(0:numnfp))

      npn = 0
      do i=0,numnfp-1
        if(np_ldmtr(nd_nfp(i)).gt.8) then
          nd_nfp(i)=daughter_ldmtr(nd_nfp(i))+c_nfp(i)
        else
          nalist(npn)=nd_nfp(i)
          talist(npn)=i
          npn = npn+1
        endif
      enddo   
      do i=0,npn-1
        nd_nfp(talist(i)) = daughter_ldmtr(nalist(i))+nd_sc(nalist(i))
        nd_sc(nalist(i)) = nd_sc(nalist(i))+1
      enddo

      deallocate(nd_sc)

      do i=0,npn-1
        c_sc(nd_nfp(talist(i))) = c_nfp(talist(i))
      enddo

      deallocate(nalist)
      deallocate(talist)

      allocate(s_sc(0:numsc))

! * update info of subcell *
      do i=0,numnfp-1
        np_sc(nd_nfp(i))=np_sc(nd_nfp(i))+1
        s_sc(nd_nfp(i))=snfp(i)
      enddo
      
      if(numdmtr.gt.maxntr-numnfp) then
        write(6,*) ' Error in ldmbuildtree() : Node is overflow!'
        write(6,*) '   This level is ',level,'numnfp=',numnfp
        stop
      endif

! *** subcell is conected to tree ***
      allocate(talist(0:numsc))
      allocate(nalist(0:numsc))
      allocate(node(0:numsc))
      allocate(clist(0:numsc))
      allocate(slist(0:numsc))

      npn = 0
      do i=0,numsc-1
        if(np_sc(i).ge.1) then
          nalist(npn)=i
          talist(npn)=numdmtr
          node(npn)=pare_sc(i)
          clist(npn)=c_sc(i)
          slist(npn)=s_sc(i)
          numdmtr=numdmtr+1
          npn=npn+1
        endif
      enddo

      deallocate(s_sc)

! *** setting daughter_ldmtr ***
! *** npn is the number of new nodes
      do i=0,npn-1
        if(flag_pc(node(i)).eq.0) then
! *** update daughter node name (now _ldmtr id) only by first one
          daughter_ldmtr(node(i))=talist(i)
          flag_pc(node(i))=1
        endif
      enddo

      allocate(nd_sc(0:numsc))

! *** set up the other values for _ldmtr ***
      do i=0,npn-1
! *** nalist is subcell id, talist is _ldmtr id
        nd_sc(nalist(i))=talist(i)
        np_ldmtr(talist(i))=np_sc(nalist(i))
        l_ldmtr(talist(i))=l_sc2
        pare_tr(talist(i))=node(i)
        call phixyzp(slist(i),clist(i),ix,iy,iz,level)
        cx_ldmtr(talist(i))=dx(ix)*l_sc+cx_ldmtr(node(i))
        cy_ldmtr(talist(i))=dx(iy)*l_sc+cy_ldmtr(node(i))
        cz_ldmtr(talist(i))=dx(iz)*l_sc+cz_ldmtr(node(i))
      enddo

      deallocate(talist)
      deallocate(nalist)
      deallocate(node)
      deallocate(clist)
      deallocate(slist)

! *** Set label not-finished particle ***
      do i=0,numnfp-1
        nd_nfp(i)=nd_sc(nd_nfp(i))
        pn_ldmtr(nd_nfp(i))=pn_nfp(i)
! *  this node is leaf *        
        if(np_ldmtr(nd_nfp(i)).eq.1) then
          label_nfp(i)=1
        endif
      enddo        

      deallocate(nd_sc)

! *** rebuild not finished particle list ***
      numtmp = 0
      do i=0,numnfp-1
        if(label_nfp(i).eq.0) then
! *** need further refinement
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

! *** stsc is starting node at this level
! *** numtmp is the number of new node
      numtmp = numdmtr-stsc

      deallocate(np_sc)
      deallocate(pare_sc)
      deallocate(c_sc)
      if(numtmp.gt.0) then
        allocate(np_sc(0:numtmp*8))
        allocate(pare_sc(0:numtmp*8))
        allocate(c_sc(0:numtmp*8))
      endif

      numsc = 0
      do i=0,numtmp-1
        if(np_ldmtr(stsc).ge.2) then
          daughter_ldmtr(stsc)=numsc
          flag_pc(stsc)=0

          np_sc(numsc)=0
          pare_sc(numsc)=stsc
          c_sc(numsc)=0

          np_sc(numsc+1)=0
          pare_sc(numsc+1)=stsc
          c_sc(numsc+1)=1

          np_sc(numsc+2)=0
          pare_sc(numsc+2)=stsc
          c_sc(numsc+2) = 2

          np_sc(numsc+3)=0
          pare_sc(numsc+3)=stsc
          c_sc(numsc+3)=3

          np_sc(numsc+4)=0
          pare_sc(numsc+4)=stsc
          c_sc(numsc+4)=4

          np_sc(numsc+5)=0
          pare_sc(numsc+5)=stsc
          c_sc(numsc+5)=5

          np_sc(numsc+6)=0
          pare_sc(numsc+6)=stsc
          c_sc(numsc+6)=6

          np_sc(numsc+7)=0
          pare_sc(numsc+7)=stsc
          c_sc(numsc+7) = 7

          numsc=numsc+8
        else if(np_ldmtr(stsc).eq.1) then
          daughter_ldmtr(stsc)=-1
        endif          
        stsc=stsc+1
      enddo
! *** store starting _ldmtr id 
      stsc = numdmtr
      if(numsc.gt.maxntr) then
        write(6,*) 'Error in ldmbuildtree():Subcell is overflow!'
        write(6,*) ' This level is ',level
        stop
      endif
      goto 77
! *** set next node ***
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

      next_ldmtr(numtr0)=0
      if(numdmtr.gt.1) then

        allocate(pare_sc(0:numdmtr-2))

        do i=numtr0,numdmtr-2
          pare_sc(i)=pare_tr(i+1)
        enddo
        do i=numtr0+1,numdmtr-2
          if(pare_sc(i).eq.pare_tr(i)) then
            next_ldmtr(i)=i+1
          else
            next_ldmtr(i)=next_ldmtr(pare_tr(i))
          endif
        enddo
        next_ldmtr(numdmtr-1)=next_ldmtr(pare_tr(numdmtr-1))

        deallocate(pare_sc)

      endif

! *** compute mass ***
      allocate(lvlist(0:level))
      allocate(paretrcm(0:numdmtr))
      do i=0,level
        lvlist(i)=list(i)
      enddo
      do i=0,numdmtr-1
        paretrcm(i)=pare_tr(i)
      enddo

      call ldmcompute_mass(numtr0,numdmtr,level,lvlist,paretrcm)

      deallocate(lvlist)
      deallocate(paretrcm)

! *** send the data to the other nodes ***
! *** get number of subdomains from each proc ***
   91 do i=0,numdmtr-1
! *** set procid ***
        proc_ldmtr(i)=myrank
      enddo 

! *** adjust cx,cy,cz,l_tr ***
      num = numdmtr-1
      do j=level-1,0,-1
        do i=num,list(j),-1
          if(daughter_ldmtr(i).gt.0.and.(np_ldmtr(i) &
           .eq.np_ldmtr(daughter_ldmtr(i)))) then
            nd=daughter_ldmtr(i)
            l_ldmtr(i)=l_ldmtr(nd)
            cx_ldmtr(i)=cx_ldmtr(nd)
            cy_ldmtr(i)=cy_ldmtr(nd)
            cz_ldmtr(i)=cz_ldmtr(nd)
            daughter_ldmtr(i)=daughter_ldmtr(nd)
          endif
        enddo
        num=list(j)-1
      enddo

! set mass, cx and etc. for node 0
      daughter_ldmtr(0)=daughter_ldmtr(1)
      mass_ldmtr(0)=mass_ldmtr(1)
      cmx_ldmtr(0)=cmx_ldmtr(1)
      cmy_ldmtr(0)=cmy_ldmtr(1)
      cmz_ldmtr(0)=cmz_ldmtr(1)
      hm_ldmtr(0)=hm_ldmtr(1)
      l_ldmtr(0)=l_ldmtr(1)
      cx_ldmtr(0)=cx_ldmtr(1)
      cy_ldmtr(0)=cy_ldmtr(1)
      cz_ldmtr(0)=cz_ldmtr(1)
#ifdef QPTREEF
      mx_ldmtr(0)=mx_ldmtr(1)
      my_ldmtr(0)=my_ldmtr(1)
      mz_ldmtr(0)=mz_ldmtr(1)
      mxx_ldmtr(0)=mxx_ldmtr(1)
      myy_ldmtr(0)=myy_ldmtr(1)
      mzz_ldmtr(0)=mzz_ldmtr(1)
      mxy_ldmtr(0)=mxy_ldmtr(1)
      myz_ldmtr(0)=myz_ldmtr(1)
      mzx_ldmtr(0)=mzx_ldmtr(1)
#endif

   94 nodess_ldmtr=1
      nodese_ldmtr=1
      if(nprocs.gt.1) then
        numtr0=numdmtr
        numtrs=0
! *** list(0)=1=numtr0 ***
        i=0
        do i=0,level-2
!          if(nodess_ldmtr.eq.0.and.np_ldmtr(list(i)) &
!           .ne.np_ldmtr(list(i+1))) then
!            nodess_ldmtr=list(i)
!          endif
          if(list(i+2)-nodess_ldmtr+1.gt.MAXNODESEND) then
            goto 92
          endif
        enddo
        if(i.eq.level-1.and.numdmtr-nodess_ldmtr+1.gt.MAXNODESEND) then
          i=level-1
! *** sending up to level-1 ***
          goto 92
        endif
! *** sending all the node ***
        nodese_ldmtr=numdmtr-1
        goto 93
   92   nodese_ldmtr=list(i+1)-1

   93   if(numdmtr.eq.0) then
          numtrs=0
          nodess_ldmtr=0
          nodese_ldmtr=0
        else if(nodess_ldmtr.gt.nodese_ldmtr) then
          write(6,*) ' Error in ldmbuildtree(): nodess > nodese'
          write(6,*) ' nodess,nodese=',nodess_ldmtr,nodese_ldmtr
          write(6,*) ' rank=',myrank,' np,numgtr=',np,numdmtr
          stop
        else
          numtrs=nodese_ldmtr-nodess_ldmtr+1
        endif

        allocate(npjr(0:nprocs))

        do ip=0,nprocs-1
          if(ip.eq.myrank) then
            npjr(ip)=numtrs
          endif
          call MPI_BCAST(npjr(ip),1,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
        enddo

        allocate(nd0(0:maxntr))

        nvali=4
#ifdef QPTREEF
        nval=18
#else
        nval=9
#endif
        do ip=0,nprocs-1
          if(ip.eq.myrank) then

            allocate(tivr(0:numtrs*nvali))
            allocate(tdvr(0:numtrs*nval))

            do i=0,numtrs-1
              nd=nodess_ldmtr+i
              tivr(i)=np_ldmtr(nd)
              tivr(i+numtrs)=next_ldmtr(nd)
              tivr(i+numtrs*2)=daughter_ldmtr(nd)
              if(tivr(i+numtrs*2).gt.nodese_ldmtr) then
                tivr(i+numtrs*2)=-1
              endif
              tivr(i+numtrs*3)=nd
! *** double precision values for neighbour search ***
              tdvr(i)=l_ldmtr(nd)
              tdvr(i+numtrs)=cx_ldmtr(nd)
              tdvr(i+numtrs*2)=cy_ldmtr(nd)
              tdvr(i+numtrs*3)=cz_ldmtr(nd)
              tdvr(i+numtrs*4)=hm_ldmtr(nd)
              tdvr(i+numtrs*5)=mass_ldmtr(nd)
              tdvr(i+numtrs*6)=cmx_ldmtr(nd)
              tdvr(i+numtrs*7)=cmy_ldmtr(nd)
              tdvr(i+numtrs*8)=cmz_ldmtr(nd)
#ifdef QPTREEF
              tdvr(i+numtrs*9)=mx_ldmtr(nd)
              tdvr(i+numtrs*10)=my_ldmtr(nd)
              tdvr(i+numtrs*11)=mz_ldmtr(nd)
              tdvr(i+numtrs*12)=mxx_ldmtr(nd)
              tdvr(i+numtrs*13)=myy_ldmtr(nd)
              tdvr(i+numtrs*14)=mzz_ldmtr(nd)
              tdvr(i+numtrs*15)=mxy_ldmtr(nd)
              tdvr(i+numtrs*16)=myz_ldmtr(nd)
              tdvr(i+numtrs*17)=mzx_ldmtr(nd)
#endif
            enddo
          else

            allocate(tivr(0:npjr(ip)*nvali))
            allocate(tdvr(0:npjr(ip)*nval))

          endif
          call MPI_BCAST(tivr,npjr(ip)*nvali,MPI_INTEGER,ip &
           ,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(tdvr,npjr(ip)*nval,MPI_DOUBLE_PRECISION,ip &
           ,MPI_COMM_WORLD,ierr)

          if(ip.ne.myrank) then
! *** add tree ***
            nrnd=npjr(ip)
            numtri=numdmtr
            do i=0,nrnd-1
              np_ldmtr(numdmtr)=tivr(i)
! *** set the other parameter for node ***
              proc_ldmtr(numdmtr)=ip
              next_ldmtr(numdmtr)=tivr(i+nrnd)
              daughter_ldmtr(numdmtr)=tivr(i+nrnd*2)
! *** link between original node id and id in this proc ***
              nd0(tivr(i+nrnd*3))=numdmtr
              pare_tr(numdmtr)=0
! *** adding nodes from the other procs
              l_ldmtr(numdmtr)=tdvr(i)
              cx_ldmtr(numdmtr)=tdvr(i+nrnd)
              cy_ldmtr(numdmtr)=tdvr(i+nrnd*2)
              cz_ldmtr(numdmtr)=tdvr(i+nrnd*3)
              hm_ldmtr(numdmtr)=tdvr(i+nrnd*4)
              mass_ldmtr(numdmtr)=tdvr(i+nrnd*5)
              cmx_ldmtr(numdmtr)=tdvr(i+nrnd*6)
              cmy_ldmtr(numdmtr)=tdvr(i+nrnd*7)
              cmz_ldmtr(numdmtr)=tdvr(i+nrnd*8)
#ifdef QPTREEF
              mx_ldmtr(numdmtr)=tdvr(i+nrnd*9)
              my_ldmtr(numdmtr)=tdvr(i+nrnd*10)
              mz_ldmtr(numdmtr)=tdvr(i+nrnd*11)
              mxx_ldmtr(numdmtr)=tdvr(i+nrnd*12)
              myy_ldmtr(numdmtr)=tdvr(i+nrnd*13)
              mzz_ldmtr(numdmtr)=tdvr(i+nrnd*14)
              mxy_ldmtr(numdmtr)=tdvr(i+nrnd*15)
              myz_ldmtr(numdmtr)=tdvr(i+nrnd*16)
              mzx_ldmtr(numdmtr)=tdvr(i+nrnd*17)
#endif
              numdmtr=numdmtr+1
            enddo
            if(nrnd.ne.0) then
              noders=tivr(nrnd*3)
              nodere=tivr(nrnd*4-1)
            else
              noders=numdmtr
              nodere=numdmtr
            endif

            do i=numtri,numdmtr-1
! *** set next_ldmtr ***
              if(next_ldmtr(i).ge.noders.and.next_ldmtr(i).le.nodere) then
                next_ldmtr(i)=nd0(next_ldmtr(i))
              else
                next_ldmtr(i)=numdmtr
              endif
! *** set daughter_ldmtr ***     
              if(daughter_ldmtr(i).ne.-1) then
                daughter_ldmtr(i)=nd0(daughter_ldmtr(i))
              endif
            enddo
          endif

          deallocate(tivr)
          deallocate(tdvr)

        enddo

        deallocate(npjr)
        deallocate(nd0)

! *** check there is any imported nodes if yes, change next_gtr(0) ***
        if(numtr0.ne.numdmtr) then
          if(numtr0.ne.0) then
            next_ldmtr(0)=numtr0
          endif
! *** set end of tree ***
          do i=numtri,numdmtr-1
! *** set next_ldmtr ***
            if(next_ldmtr(i).eq.numdmtr) then
              next_ldmtr(i)=-1
            endif
          enddo
        endif
      endif

! *** compute delta ***
      do i=0,numdmtr-1
        delta_ldmtr(i)=dsqrt((cmx_ldmtr(i)-cx_ldmtr(i))**2 &
         +(cmy_ldmtr(i)-cy_ldmtr(i))**2+(cmz_ldmtr(i)-cz_ldmtr(i))**2)
      enddo		

      num_ldmtr=numdmtr

      deallocate(list)
      deallocate(pare_tr)	

end subroutine
#endif

#if defined(DM) && defined(MULTI_LEV)	  
! *** Definition of makeroot() ***
subroutine ldmmakeroot(np)
      use gcdp_const
      use gcdp_dm
      use gcdp_ldmtree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np
      integer i,ierr
! * max coordinate *      
      double precision max_x,max_y,max_z
! * min coordinate <0 *      
      double precision min_x,min_y,min_z
! * max,temp length *      
      double precision maxl,tl
! for work
      double precision tdvr(0:2),tdvs(0:2)
      
! *** Define root node ***
      max_x=-INF
      max_y=-INF
      max_z=-INF
      min_x=INF
      min_y=INF
      min_z = INF
      do i=SI_ndm1,np-1
        if(x_dm(i).lt.min_x) then
          min_x = x_dm(i)
        endif          
        if(y_dm(i).lt.min_y) then
          min_y = y_dm(i)
        endif
        if(z_dm(i).lt.min_z) then
          min_z = z_dm(i)
        endif
        if(x_dm(i).gt.max_x) then
          max_x = x_dm(i)
        endif
        if(y_dm(i).gt.max_y) then
          max_y = y_dm(i)
        endif
        if(z_dm(i).gt.max_z) then
          max_z = z_dm(i)
        endif
      enddo
! *** get the maximum and minimum for all the particles ***
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
      tl=max_z-min_z
      if(tl.lt.0.0d0) then
        tl = -tl
      endif
      if(tl.gt.maxl) then
        maxl=tl
      endif

! *** Set root node for node "flag" ***
      l_ldmtr(0)=MGROOT*maxl
      cx_ldmtr(0)=(max_x+min_x)*0.5d0
      cy_ldmtr(0)=(max_y+min_y)*0.5d0
      cz_ldmtr(0)=(max_z+min_z)*0.5d0
      np_ldmtr(0)=np-SI_ndm1
      pn_ldmtr(0)=SI_ndm1

end subroutine
#endif


#if defined(DM) && defined(MULTI_LEV)
! *** Definition of compute_mass() ***
! * This function compute the center of mass and total mass *
! * tree, center of mass, total mass, node No. *	       
subroutine ldmcompute_mass(numtr0,numdmtr,level,list,pare_tr)
      use gcdp_const
      use gcdp_dm
      use gcdp_ldmtree

      implicit none

      integer,intent(in) :: numtr0,numdmtr,level
      integer,intent(in) :: list(0:level),pare_tr(0:numdmtr-1)
      integer i
      integer j,d,npare,num,nd,pnd
      double precision xsd,ysd,zsd
! for work
      integer,allocatable :: pare_sc(:),nd_sc(:)

! * compute total mass in node *
      do i=numtr0,numdmtr-1
        pnd=pn_ldmtr(i)
        if(np_ldmtr(i).eq.1) then
          mass_ldmtr(i)=m_dm(pnd)
          cmx_ldmtr(i)=m_dm(pnd)*x_dm(pnd)
          cmy_ldmtr(i)=m_dm(pnd)*y_dm(pnd)
          cmz_ldmtr(i)=m_dm(pnd)*z_dm(pnd)
! *** if np=1, cx should be the position of the particle ***
          cx_ldmtr(i)=x_dm(pnd)
          cy_ldmtr(i)=y_dm(pnd)
          cz_ldmtr(i)=z_dm(pnd)
          l_ldmtr(i)=0.0d0
          hm_ldmtr(i)=m_dm(pnd)*h_dm(pnd)
        else
          mass_ldmtr(i)=0.0d0
          cmx_ldmtr(i)=0.0d0
          cmy_ldmtr(i)=0.0d0
          cmz_ldmtr(i)=0.0d0
          hm_ldmtr(i) = 0.0d0
        endif
      enddo

      num = numdmtr-1	
      do j=level,1,-1
        npare = 0

        if(num.ge.list(j)) then
          allocate(pare_sc(0:num-list(j)))
          allocate(nd_sc(0:num-list(j)))    
        endif

        do i=num,list(j),-1
          if(pare_tr(i).ne.pare_tr(i-1)) then
            pare_sc(npare)=pare_tr(i)
            nd_sc(npare)=i
            npare=npare+1
          endif
        enddo
        do d=0,7
          do i=0,npare-1
            pnd=pare_sc(i)             
            nd=nd_sc(i)+d
            if(pnd.eq.pare_tr(nd).and.nd.le.num) then
              mass_ldmtr(pnd)=mass_ldmtr(pnd)+mass_ldmtr(nd)
              cmx_ldmtr(pnd)=cmx_ldmtr(pnd)+cmx_ldmtr(nd)
              cmy_ldmtr(pnd)=cmy_ldmtr(pnd)+cmy_ldmtr(nd)
              cmz_ldmtr(pnd)=cmz_ldmtr(pnd)+cmz_ldmtr(nd)
              hm_ldmtr(pnd)=hm_ldmtr(pnd)+hm_ldmtr(nd)
            endif
          enddo            
        enddo

        if(num.ge.list(j)) then
          deallocate(pare_sc)
          deallocate(nd_sc)
        endif

        num = list(j)-1        
      enddo         

! * compute center of mass *
      do i=numtr0,numdmtr-1
        cmx_ldmtr(i)=cmx_ldmtr(i)/mass_ldmtr(i)
        cmy_ldmtr(i)=cmy_ldmtr(i)/mass_ldmtr(i)
        cmz_ldmtr(i)=cmz_ldmtr(i)/mass_ldmtr(i)
       hm_ldmtr(i)=hm_ldmtr(i)/mass_ldmtr(i)
      enddo
#ifdef QPTREEF
! * compute Multipole Momentum *
      do i=numtr0,numdmtr-1
        mx_ldmtr(i)=0.0d0
        my_ldmtr(i)=0.0d0
        mz_ldmtr(i)=0.0d0
        mxx_ldmtr(i)=0.0d0
        myy_ldmtr(i)=0.0d0
        mzz_ldmtr(i)=0.0d0
        mxy_ldmtr(i)=0.0d0
        myz_ldmtr(i)=0.0d0
        mzx_ldmtr(i) = 0.0d0
      enddo
      num=numdmtr-1
      do j=level,1,-1
        npare = 0

        if(num.ge.list(j)) then
          allocate(pare_sc(0:num-list(j)))
          allocate(nd_sc(0:num-list(j)))    
        endif

        do i=num,list(j),-1
          if(pare_tr(i).ne.pare_tr(i-1)) then
            pare_sc(npare)=pare_tr(i)
            nd_sc(npare)=i
            npare=npare+1
          endif
        enddo
        do d=0,7
          do i=0,npare-1
            nd=nd_sc(i)+d
            pnd=pare_sc(i)
            if(pnd.eq.pare_tr(nd).and.nd.le.num) then
              xsd=cmx_ldmtr(nd)-cmx_ldmtr(pnd)
              ysd=cmy_ldmtr(nd)-cmy_ldmtr(pnd)
              zsd=cmz_ldmtr(nd)-cmz_ldmtr(pnd)
              mx_ldmtr(pnd)=mx_ldmtr(pnd)+(mx_ldmtr(nd)-xsd*mass_ldmtr(nd))
              my_ldmtr(pnd)=my_ldmtr(pnd)+(my_ldmtr(nd)-ysd*mass_ldmtr(nd))
              mz_ldmtr(pnd)=mz_ldmtr(pnd)+(mz_ldmtr(nd)-zsd*mass_ldmtr(nd))
              mxx_ldmtr(pnd)=mxx_ldmtr(pnd)+(mxx_ldmtr(nd) &
               -2.0d0*xsd*mx_ldmtr(nd)+xsd*xsd*mass_ldmtr(nd))
              myy_ldmtr(pnd)=myy_ldmtr(pnd)+(myy_ldmtr(nd) &
               -2.0d0*ysd*my_ldmtr(nd)+ysd*ysd*mass_ldmtr(nd))
              mzz_ldmtr(pnd)=mzz_ldmtr(pnd)+(mzz_ldmtr(nd) &
               -2.0d0*zsd*mz_ldmtr(nd)+zsd*zsd*mass_ldmtr(nd))
              mxy_ldmtr(pnd)=mxy_ldmtr(pnd)+(mxy_ldmtr(nd) &
               -xsd*my_ldmtr(nd)-ysd*mx_ldmtr(nd)+xsd*ysd*mass_ldmtr(nd))
              myz_ldmtr(pnd)=myz_ldmtr(pnd)+(myz_ldmtr(nd) &
               -ysd*mz_ldmtr(nd)-zsd*my_ldmtr(nd)+ysd*zsd*mass_ldmtr(nd))
              mzx_ldmtr(pnd)=mzx_ldmtr(pnd)+(mzx_ldmtr(nd) &
               -zsd*mx_ldmtr(nd)-xsd*mz_ldmtr(nd)+zsd*xsd*mass_ldmtr(nd))
            endif
          enddo
        enddo

        if(num.ge.list(j)) then
          deallocate(pare_sc)
          deallocate(nd_sc)
        endif

        num = list(j)-1
      enddo
#endif

end subroutine
#endif


