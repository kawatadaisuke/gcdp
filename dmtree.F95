#include "gcdp.def"
! ***************************************
!     dmtree.F95  for GCD+ ver f03.0
!  22 Mar. 2013   written by D. Kawata
!  *************************************** 

#ifdef DM
subroutine dmtreebuild(np)      
      use gcdp_const
      use gcdp_dm
      use gcdp_dmtree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np
      integer pn,nd 
      integer i,ierr
      integer level,npn
! * for finished tree (dmtr) *
      integer numdmtr
! * for not-finished particle (nfp) *
! * number of not-finished particle *      
      integer numnfp,maxntr
! * for subcell (sc) *
! * start subcell in tree *      
      integer stsc
! * number of subcells *      
      integer numsc
! * diameter of subcell *      
      double precision l_sc,l_sc2
! * for work *
      integer numtmp
! *** for calculating hm_dmtr ***
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
      call allocate_dmtree(maxntr)

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

      do i=0,np-1
        pn_nfp(i)=i
      enddo
      level=0
      call dmmakeroot(np,0.0d0)
      next_dmtr(0)=0
      pare_tr(0)=-1
      daughter_dmtr(0)=1
      numdmtr=1
      if(np.eq.0) then
        np_dmtr(0)=0
        numdmtr=0
        cmx_dmtr(1)=0.0d0
        cmy_dmtr(1)=0.0d0
        cmz_dmtr(1)=0.0d0
        hm_dmtr(1)=0.0d0
#ifdef QPTREEF
        mx_dmtr(1)=0.0d0
        my_dmtr(1)=0.0d0
        mz_dmtr(1)=0.0d0
        mxx_dmtr(1)=0.0d0
        myy_dmtr(1)=0.0d0
        mzz_dmtr(1)=0.0d0
        mxy_dmtr(1)=0.0d0
        myz_dmtr(1)=0.0d0
        mzx_dmtr(1)=0.0d0
#endif
        goto 94
      endif

      numnfp=0
      numtr0=numdmtr
      do i=0,np-1
        pn_nfp(numnfp)=i
        nd_nfp(numnfp)=numtr0
        label_nfp(numnfp)=0
        snfp(numnfp)=0
        c_nfp(numnfp)=0
        numnfp=numnfp+1
      enddo
      if(numnfp.eq.0) then
        goto 94
      endif
! *** set root for numtr=1 ***
      l_dmtr(numdmtr)=l_dmtr(0)
      cx_dmtr(numdmtr)=cx_dmtr(0)
      cy_dmtr(numdmtr)=cy_dmtr(0)
      cz_dmtr(numdmtr)=cz_dmtr(0)
      np_dmtr(numdmtr)=numnfp
      pn_dmtr(numdmtr)=pn_nfp(0)
      next_dmtr(numdmtr)=0
      pare_tr(numdmtr)=0
      numdmtr=numdmtr+1
! *** subcell starts from 0
      numsc=8
      stsc=numdmtr
      l_sc=l_dmtr(numtr0)*0.5d0

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
! *** daughter_tr keeps subcell id.
        daughter_dmtr(numtr0)=0 
      else
        daughter_dmtr(numtr0)=0 
        daughter_dmtr(numtr0+1)=-1
      endif
      level=0
! *** need to set to get hm_dmtr(0) and etc. ***
      list(0)=numtr0
! * start iteration in level *
   77 if(numnfp.le.1) then
        goto 99
      endif
! * initialization for the case np_dmtr < 8 *

      allocate(nd_sc(0:numdmtr))

      do i=list(level),numdmtr-1
        nd_sc(i)=0
      enddo
      level=level+1	  
! beginning of this level
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
        if(x_dm(pn)-cx_dmtr(nd).ge.0.0d0) then
          ixp(i)=0 
        else
          ixp(i)=1
        endif 
        if(y_dm(pn)-cy_dmtr(nd).ge.0.0d0) then
          iyp(i)=0 
        else
          iyp(i)=1
        endif 
        if(z_dm(pn)-cz_dmtr(nd).ge.0.0d0) then
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
        if(np_dmtr(nd_nfp(i)).gt.8) then
! *** nd_nfp is subcell id (starting from 0)
!     because daughter_tr keeps subcell id
          nd_nfp(i)=daughter_dmtr(nd_nfp(i))+c_nfp(i)
        else
! *** nd_nfp is still tree id 
          nalist(npn)=nd_nfp(i)
          talist(npn)=i
          npn = npn+1
        endif
      enddo   
      do i=0,npn-1
! *** nd_sc counting number of subcell for nalist(i)
! *** nd_nfp is now subcell id
        nd_nfp(talist(i)) = daughter_dmtr(nalist(i))+nd_sc(nalist(i))
        nd_sc(nalist(i)) = nd_sc(nalist(i))+1
      enddo

      deallocate(nd_sc)

      do i=0,npn-1
! *** set c for the case np_dmtr is less than 0
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
        write(6,*) ' Error in dmbuildtree() : Node is overflow!'
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

! *** setting daughter_dmtr ***
! *** npn is the number of new nodes
      do i=0,npn-1
! *** flag for parent node, initially 0 
        if(flag_pc(node(i)).eq.0) then
! *** update daughter node name (now _dmtr id) only by first one
          daughter_dmtr(node(i))=talist(i)
          flag_pc(node(i))=1
        endif
      enddo

      allocate(nd_sc(0:numsc))

! *** set up the other values for _dmtr ***
      do i=0,npn-1
! *** nalist is subcell id, talist is _dmtr id
        nd_sc(nalist(i))=talist(i)
        np_dmtr(talist(i))=np_sc(nalist(i))
        l_dmtr(talist(i))=l_sc2
! *** node is parent _tr id
        pare_tr(talist(i))=node(i)
        call phixyzp(slist(i),clist(i),ix,iy,iz,level)
        cx_dmtr(talist(i))=dx(ix)*l_sc+cx_dmtr(node(i))
        cy_dmtr(talist(i))=dx(iy)*l_sc+cy_dmtr(node(i))
        cz_dmtr(talist(i))=dx(iz)*l_sc+cz_dmtr(node(i))
      enddo

      deallocate(talist)
      deallocate(nalist)
      deallocate(node)
      deallocate(clist)
      deallocate(slist)

! *** Set label not-finished particle ***
      do i=0,numnfp-1
! *** nd_nfp from subcell id to _dmtr id ***
        nd_nfp(i)=nd_sc(nd_nfp(i))
        pn_dmtr(nd_nfp(i))=pn_nfp(i)
! *  this node is leaf *        
        if(np_dmtr(nd_nfp(i)).eq.1) then
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
        if(np_dmtr(stsc).ge.2) then
          daughter_dmtr(stsc)=numsc
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
        else if(np_dmtr(stsc).eq.1) then
          daughter_dmtr(stsc)=-1
        endif          
        stsc=stsc+1
      enddo
! *** store starting _dmtr id 
      stsc = numdmtr
      if(numsc.gt.maxntr) then
        write(6,*) 'Error in buildtree():Subcell is overflow!'
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

! *** set next node for node numtr0=1 (always go back to node 0, then look at node 0s next ***
      next_dmtr(numtr0)=0
      if(numdmtr.gt.1) then

        allocate(pare_sc(0:numdmtr-2))

        do i=numtr0,numdmtr-2
          pare_sc(i)=pare_tr(i+1)
        enddo
        do i=numtr0+1,numdmtr-2
! *** set next_tr, compare the parent of the next tr ***
          if(pare_sc(i).eq.pare_tr(i)) then
            next_dmtr(i)=i+1
          else
            next_dmtr(i)=next_dmtr(pare_tr(i))
          endif
        enddo
        next_dmtr(numdmtr-1)=next_dmtr(pare_tr(numdmtr-1))

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

      call dmcompute_mass(numtr0,numdmtr,level,lvlist,paretrcm)

      deallocate(lvlist)
      deallocate(paretrcm)

! *** send the data to the other nodes ***
! *** get number of subdomains from each proc ***
   91 do i=0,numdmtr-1
! *** set procid ***
        proc_dmtr(i)=myrank
      enddo 
! *** adjust cx,cy,cz,l_tr ***
      num = numdmtr-1
      do j=level-1,0,-1
        do i=num,list(j),-1
          if(daughter_dmtr(i).gt.0.and.(np_dmtr(i) &
           .eq.np_dmtr(daughter_dmtr(i)))) then
            nd=daughter_dmtr(i)
            l_dmtr(i)=l_dmtr(nd)
            cx_dmtr(i)=cx_dmtr(nd)
            cy_dmtr(i)=cy_dmtr(nd)
            cz_dmtr(i)=cz_dmtr(nd)
            daughter_dmtr(i)=daughter_dmtr(nd)
          endif
        enddo
        num=list(j)-1
      enddo
! set mass, cx and etc. for node 0
      daughter_dmtr(0)=daughter_dmtr(1)
      mass_dmtr(0)=mass_dmtr(1)
      cmx_dmtr(0)=cmx_dmtr(1)
      cmy_dmtr(0)=cmy_dmtr(1)
      cmz_dmtr(0)=cmz_dmtr(1)
      hm_dmtr(0)=hm_dmtr(1)
      l_dmtr(0)=l_dmtr(1)
      cx_dmtr(0)=cx_dmtr(1)
      cy_dmtr(0)=cy_dmtr(1)
      cz_dmtr(0)=cz_dmtr(1)
#ifdef QPTREEF
      mx_dmtr(0)=mx_dmtr(1)
      my_dmtr(0)=my_dmtr(1)
      mz_dmtr(0)=mz_dmtr(1)
      mxx_dmtr(0)=mxx_dmtr(1)
      myy_dmtr(0)=myy_dmtr(1)
      mzz_dmtr(0)=mzz_dmtr(1)
      mxy_dmtr(0)=mxy_dmtr(1)
      myz_dmtr(0)=myz_dmtr(1)
      mzx_dmtr(0)=mzx_dmtr(1)
#endif

! *** nodess should be starting from 1, if numdmtr.eq.0, corrected later
   94 nodess_dmtr=1
      nodese_dmtr=1
      if(nprocs.gt.1) then
        numtr0=numdmtr
        numtrs=0
! *** list(0)=1=numtr0 ***
        i=0
        do i=0,level-2
          if(nodess_dmtr.eq.0.and.np_dmtr(list(i)).ne.np_dmtr(list(i+1))) then
! *** set starting node ***
            nodess_dmtr=list(i)
          endif
          if(list(i+2)-nodess_dmtr+1.gt.MAXNODESEND) then
! *** sending up to level i+1
            goto 92
          endif
        enddo
        if(i.eq.level-1.and.numdmtr-nodess_dmtr+1.gt.MAXNODESEND) then
          i=level-1
! *** sending up to level-1 ***
          goto 92
        endif
! *** sending all the node ***
        nodese_dmtr=numdmtr-1
        goto 93
   92   nodese_dmtr=list(i+1)-1

   93   if(numdmtr.eq.0) then
          numtrs=0
          nodess_dmtr=0
          nodese_dmtr=0
        else if(nodess_dmtr.gt.nodese_dmtr) then
          write(6,*) ' Error in dmbuildtree(): nodess > nodese'
          write(6,*) ' nodess,nodese=',nodess_dmtr,nodese_dmtr
          write(6,*) ' rank=',myrank,' np,numgtr=',np,numdmtr
          stop
        else
          numtrs=nodese_dmtr-nodess_dmtr+1
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
              nd=nodess_dmtr+i
              tivr(i)=np_dmtr(nd)
              tivr(i+numtrs)=next_dmtr(nd)
              tivr(i+numtrs*2)=daughter_dmtr(nd)
              if(tivr(i+numtrs*2).gt.nodese_dmtr) then
                tivr(i+numtrs*2)=-1
              endif
              tivr(i+numtrs*3)=nd
! *** double precision values for neighbour search ***
              tdvr(i)=l_dmtr(nd)
              tdvr(i+numtrs)=cx_dmtr(nd)
              tdvr(i+numtrs*2)=cy_dmtr(nd)
              tdvr(i+numtrs*3)=cz_dmtr(nd)
              tdvr(i+numtrs*4)=hm_dmtr(nd)
              tdvr(i+numtrs*5)=mass_dmtr(nd)
              tdvr(i+numtrs*6)=cmx_dmtr(nd)
              tdvr(i+numtrs*7)=cmy_dmtr(nd)
              tdvr(i+numtrs*8)=cmz_dmtr(nd)
#ifdef QPTREEF
              tdvr(i+numtrs*9)=mx_dmtr(nd)
              tdvr(i+numtrs*10)=my_dmtr(nd)
              tdvr(i+numtrs*11)=mz_dmtr(nd)
              tdvr(i+numtrs*12)=mxx_dmtr(nd)
              tdvr(i+numtrs*13)=myy_dmtr(nd)
              tdvr(i+numtrs*14)=mzz_dmtr(nd)
              tdvr(i+numtrs*15)=mxy_dmtr(nd)
              tdvr(i+numtrs*16)=myz_dmtr(nd)
              tdvr(i+numtrs*17)=mzx_dmtr(nd)
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

! *** should have all the pseudo node even if np=0, because DM may need it
          if(ip.ne.myrank) then
! *** add tree ***
            nrnd=npjr(ip)
            numtri=numdmtr
            do i=0,nrnd-1
              np_dmtr(numdmtr)=tivr(i)
! *** set the other parameter for node ***
              proc_dmtr(numdmtr)=ip
              next_dmtr(numdmtr)=tivr(i+nrnd)
              daughter_dmtr(numdmtr)=tivr(i+nrnd*2)
! *** link between original node id and id in this proc ***
              nd0(tivr(i+nrnd*3))=numdmtr
              pare_tr(numdmtr)=0
! *** adding nodes from the other procs
              l_dmtr(numdmtr)=tdvr(i)
              cx_dmtr(numdmtr)=tdvr(i+nrnd)
              cy_dmtr(numdmtr)=tdvr(i+nrnd*2)
              cz_dmtr(numdmtr)=tdvr(i+nrnd*3)
              hm_dmtr(numdmtr)=tdvr(i+nrnd*4)
              mass_dmtr(numdmtr)=tdvr(i+nrnd*5)
              cmx_dmtr(numdmtr)=tdvr(i+nrnd*6)
              cmy_dmtr(numdmtr)=tdvr(i+nrnd*7)
              cmz_dmtr(numdmtr)=tdvr(i+nrnd*8)
#ifdef QPTREEF
              mx_dmtr(numdmtr)=tdvr(i+nrnd*9)
              my_dmtr(numdmtr)=tdvr(i+nrnd*10)
              mz_dmtr(numdmtr)=tdvr(i+nrnd*11)
              mxx_dmtr(numdmtr)=tdvr(i+nrnd*12)
              myy_dmtr(numdmtr)=tdvr(i+nrnd*13)
              mzz_dmtr(numdmtr)=tdvr(i+nrnd*14)
              mxy_dmtr(numdmtr)=tdvr(i+nrnd*15)
              myz_dmtr(numdmtr)=tdvr(i+nrnd*16)
              mzx_dmtr(numdmtr)=tdvr(i+nrnd*17)
#endif
              numdmtr=numdmtr+1
            enddo
            if(nrnd.gt.0) then
              noders=tivr(nrnd*3)
              nodere=tivr(nrnd*4-1)
            else
              noders=numdmtr
              nodere=numdmtr
            endif

            do i=numtri,numdmtr-1
! *** set next_dmtr ***
              if(next_dmtr(i).ge.noders.and.next_dmtr(i).le.nodere) then
                next_dmtr(i)=nd0(next_dmtr(i))
              else
                next_dmtr(i)=numdmtr
              endif
! *** set daughter_dmtr ***     
              if(daughter_dmtr(i).ne.-1) then
                daughter_dmtr(i)=nd0(daughter_dmtr(i))
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
            next_dmtr(0)=numtr0
          endif
! *** set end of tree ***
          do i=numtri,numdmtr-1
! *** set next_dmtr ***
            if(next_dmtr(i).eq.numdmtr) then
              next_dmtr(i)=-1
            endif
          enddo
        endif
      endif

! *** compute delta ***
      do i=0,numdmtr-1
        delta_dmtr(i)=dsqrt((cmx_dmtr(i)-cx_dmtr(i))**2 &
         +(cmy_dmtr(i)-cy_dmtr(i))**2+(cmz_dmtr(i)-cz_dmtr(i))**2)
      enddo		

      num_dmtr=numdmtr

      deallocate(list)
      deallocate(pare_tr)	

end subroutine
#endif

#ifdef DM	  
! *** Definition of makeroot() ***
subroutine dmmakeroot(np,cz)
      use gcdp_const
      use gcdp_dm
      use gcdp_dmtree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np
      double precision,intent(in) :: cz
      integer i,ierr
! * max coordinate *      
      double precision max_x,max_y,max_z
#ifndef CENTRE0
! * min coordinate <0 *      
      double precision min_x,min_y,min_z
#endif
! * max,temp length *      
      double precision maxl,tl
! for work
      double precision tdvr(0:2),tdvs(0:2)
      
! *** Define root node ***
#ifdef CENTRE0
      max_x=-INF
      max_y=-INF
      max_z=-INF
      cz_dmtr(0)=cz
      do i=0,np-1
        if(dabs(x_dm(i)).gt.max_x) then
          max_x = dabs(x_dm(i))
        endif
        if(dabs(y_dm(i)).gt.max_y) then
          max_y = dabs(y_dm(i))
        endif
        if(dabs(z_dm(i)-cz).gt.max_z) then
          max_z = dabs(z_dm(i)-cz)
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
! *** check max range ***
      if(max_x.gt.max_y) then
        maxl=max_x
      else
        maxl=max_y
      endif
      if(max_z.gt.maxl) then
        maxl=max_z
      endif
      maxl=2.0d0*maxl
! *** Set root node for node "flag" ***
      l_dmtr(0)=MGROOT*maxl
      cx_dmtr(0)=0.0d0
      cy_dmtr(0)=0.0d0
      np_dmtr(0)=np
      pn_dmtr(0)=0
#else
      max_x=-INF
      max_y=-INF
      max_z=-INF
      min_x=INF
      min_y=INF
      min_z = INF
      do i=0,np-1
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
      l_dmtr(0)=MGROOT*maxl
      cx_dmtr(0)=(max_x+min_x)*0.5d0
      cy_dmtr(0)=(max_y+min_y)*0.5d0
      cz_dmtr(0)=(max_z+min_z)*0.5d0
      np_dmtr(0)=np
      pn_dmtr(0)=0
#endif

end subroutine
#endif

#ifdef DM
! *** Definition of compute_mass() ***
! * This function compute the center of mass and total mass *
! * tree, center of mass, total mass, node No. *	       
subroutine dmcompute_mass(numtr0,numdmtr,level,list,pare_tr)
      use gcdp_const
      use gcdp_dm
      use gcdp_dmtree

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
        pnd=pn_dmtr(i)
        if(np_dmtr(i).eq.1) then
          mass_dmtr(i)=m_dm(pnd)
          cmx_dmtr(i)=m_dm(pnd)*x_dm(pnd)
          cmy_dmtr(i)=m_dm(pnd)*y_dm(pnd)
          cmz_dmtr(i)=m_dm(pnd)*z_dm(pnd)
! *** if np=1, cx should be the position of the particle, used for set_valuedm ***
          cx_dmtr(i)=x_dm(pnd)
          cy_dmtr(i)=y_dm(pnd)
          cz_dmtr(i)=z_dm(pnd)
          l_dmtr(i)=0.0d0
          hm_dmtr(i)=m_dm(pnd)*h_dm(pnd)
        else
          mass_dmtr(i)=0.0d0
          cmx_dmtr(i)=0.0d0
          cmy_dmtr(i)=0.0d0
          cmz_dmtr(i)=0.0d0
          hm_dmtr(i) = 0.0d0
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
              mass_dmtr(pnd)=mass_dmtr(pnd)+mass_dmtr(nd)
              cmx_dmtr(pnd)=cmx_dmtr(pnd)+cmx_dmtr(nd)
              cmy_dmtr(pnd)=cmy_dmtr(pnd)+cmy_dmtr(nd)
              cmz_dmtr(pnd)=cmz_dmtr(pnd)+cmz_dmtr(nd)
              hm_dmtr(pnd)=hm_dmtr(pnd)+hm_dmtr(nd)
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
        cmx_dmtr(i)=cmx_dmtr(i)/mass_dmtr(i)
        cmy_dmtr(i)=cmy_dmtr(i)/mass_dmtr(i)
        cmz_dmtr(i)=cmz_dmtr(i)/mass_dmtr(i)
        hm_dmtr(i)=hm_dmtr(i)/mass_dmtr(i)
      enddo
#ifdef QPTREEF
! * compute Multipole Momentum *
      do i=numtr0,numdmtr-1
        mx_dmtr(i)=0.0d0
        my_dmtr(i)=0.0d0
        mz_dmtr(i)=0.0d0
        mxx_dmtr(i)=0.0d0
        myy_dmtr(i)=0.0d0
        mzz_dmtr(i)=0.0d0
        mxy_dmtr(i)=0.0d0
        myz_dmtr(i)=0.0d0
        mzx_dmtr(i) = 0.0d0
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
              xsd=cmx_dmtr(nd)-cmx_dmtr(pnd)
              ysd=cmy_dmtr(nd)-cmy_dmtr(pnd)
              zsd=cmz_dmtr(nd)-cmz_dmtr(pnd)
              mx_dmtr(pnd)=mx_dmtr(pnd)+(mx_dmtr(nd)-xsd*mass_dmtr(nd))
              my_dmtr(pnd)=my_dmtr(pnd)+(my_dmtr(nd)-ysd*mass_dmtr(nd))
              mz_dmtr(pnd)=mz_dmtr(pnd)+(mz_dmtr(nd)-zsd*mass_dmtr(nd))
              mxx_dmtr(pnd)=mxx_dmtr(pnd)+(mxx_dmtr(nd)-2.0d0*xsd*mx_dmtr(nd) &
               +xsd*xsd*mass_dmtr(nd))
              myy_dmtr(pnd)=myy_dmtr(pnd)+(myy_dmtr(nd)-2.0d0*ysd*my_dmtr(nd) &
               +ysd*ysd*mass_dmtr(nd))
              mzz_dmtr(pnd)=mzz_dmtr(pnd)+(mzz_dmtr(nd)-2.0d0*zsd*mz_dmtr(nd) &
               +zsd*zsd*mass_dmtr(nd))
              mxy_dmtr(pnd)=mxy_dmtr(pnd)+(mxy_dmtr(nd)-xsd*my_dmtr(nd) &
               -ysd*mx_dmtr(nd)+xsd*ysd*mass_dmtr(nd))
              myz_dmtr(pnd)=myz_dmtr(pnd)+(myz_dmtr(nd)-ysd*mz_dmtr(nd) &
               -zsd*my_dmtr(nd)+ysd*zsd*mass_dmtr(nd))
              mzx_dmtr(pnd)=mzx_dmtr(pnd)+(mzx_dmtr(nd)-zsd*mx_dmtr(nd) &
               -xsd*mz_dmtr(nd)+zsd*xsd*mass_dmtr(nd))
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
!      write(6,*) ' numdmtr = ',numdmtr      

end subroutine
#endif

