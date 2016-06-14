#include "gcdp.def"
! ***********************************************
!    tree.F for GCD+ ver. f03.0
!  31  Jan. 2013   written by D. Kawata
! ***********************************************
! ***********************************************
!    Definition of function about octtree
!  This program build octtree and
!    compute the center of mass and total mass and 
! ***********************************************

#ifdef SGRAV
#if defined(GAS) || defined(STAR)
subroutine treebuild(np)
      use gcdp_const
      use gcdp_baryon
      use gcdp_btree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np
      integer i,pn,nd,level,npn,ierr
! change test
! * for not-finished particle (nfp) *
! * number of not-finished particle *	  
      integer numnfp,maxntr
! * for subcell (sc) *
! * start subcell in tree *	  
      integer stsc
! * number of subcells *	  
      integer numsc
!
! * diameter of subcell *	  
      double precision l_sc,l_sc2
! * for work *
      integer numtmp
! *** for calculating hm_tr ***
      integer d,num,npare,pnd,j
      double precision xsd,ysd,zsd
! for work
      integer ix,iy,iz
      integer numtr0,nval,nvali,ip,nrnd,numtri
      integer nodess,numtrs
! *** to find root node for domain ***
      integer noders,nodere
      character fileo*60
      double precision dx(0:1)
      data dx/1.0d0,-1.0d0/
! for allocatable arrays
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
      call allocate_btree(maxntr)
!      write(6,*) ' myrank,maxntr=',myrank,maxntr

! *** Make root ***
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
      do i=0,np-1
        pn_nfp(i)=i
      enddo
      level=0


! the last 0 for node 0
! gmakeroot use pn_nfp
      call makeroot(np,pn_nfp,0.0d0)

!      write(6,*) 'myrank,np,l=',myrank,np,l_tr(0)

      next_tr(0)=0
      pare_tr(0)=-1
      daughter_tr(0)=1
      num_tr=1	
      if(np.eq.0) then
        np_tr(0)=0
        num_tr=0
        cmx_tr(1)=0.0d0
        cmy_tr(1)=0.0d0
        cmz_tr(1)=0.0d0
        hm_tr(1)=0.0d0
#ifdef QPTREEF
        mx_tr(1)=0.0d0
        my_tr(1)=0.0d0
        mz_tr(1)=0.0d0
        mxx_tr(1)=0.0d0
        myy_tr(1)=0.0d0
        mzz_tr(1)=0.0d0
        mxy_tr(1)=0.0d0
        myz_tr(1)=0.0d0
        mzx_tr(1)=0.0d0
#endif
        goto 94
      endif
      numnfp=0
      numtr0=num_tr
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
      l_tr(num_tr)=l_tr(0)
      cx_tr(num_tr)=cx_tr(0)
      cy_tr(num_tr)=cy_tr(0)
      cz_tr(num_tr)=cz_tr(0)
      np_tr(num_tr)=numnfp
      pn_tr(num_tr)=pn_nfp(0)
      next_tr(num_tr)=0
      pare_tr(num_tr)=0
      num_tr=num_tr+1
! *** subcell starts from 0
      numsc=8
      stsc=num_tr
      l_sc=l_tr(numtr0)*0.5d0	

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
! *** daughter_tr keeps subcell id.
        daughter_tr(numtr0)=0 
      else
        daughter_tr(numtr0)=0
        daughter_tr(numtr0+1)=-1
      endif
      level = 0
! *** need to set to get hm_tr(0) and etc. ***
      list(0)=numtr0
! *****   start iteration in level *****
   77 if(numnfp.le.1) then
        goto 99
      endif
! * initialization for the case np_tr < 8 *

      allocate(nd_sc(0:num_tr-1))

      do i=list(level),num_tr-1
        nd_sc(i)=0
      enddo
      level=level+1
! beginning of this level
      list(level) = num_tr
      l_sc2=l_sc
      l_sc=0.5d0*l_sc2
! * find out which subcell it is in *

      allocate(ixp(0:numnfp-1))
      allocate(iyp(0:numnfp-1))
      allocate(izp(0:numnfp-1))

      do i=0,numnfp-1
        pn=pn_nfp(i)
        nd=nd_nfp(i)
        if(x_p(pn)-cx_tr(nd).ge.0.0d0) then
          ixp(i)=0 
        else
          ixp(i)=1
        endif 
        if(y_p(pn)-cy_tr(nd).ge.0.0d0) then
          iyp(i)=0 
        else
          iyp(i)=1
        endif 
        if(z_p(pn)-cz_tr(nd).ge.0.0d0) then
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
        if(np_tr(nd_nfp(i)).gt.8) then
! *** nd_nfp is subcell id (starting from 0)
!     because daughter_tr keeps subcell id
          nd_nfp(i) = daughter_tr(nd_nfp(i))+c_nfp(i)
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
        nd_nfp(talist(i))=daughter_tr(nalist(i))+nd_sc(nalist(i))
        nd_sc(nalist(i))=nd_sc(nalist(i))+1
      enddo

      do i=0,npn-1
! *** set c for the case np_tr is less than 0
        c_sc(nd_nfp(talist(i)))=c_nfp(talist(i))  
      enddo

      deallocate(nalist)
      deallocate(talist)

      allocate(s_sc(0:numsc-1))

! * update number of particles for subcell *
      do i=0,numnfp-1
        np_sc(nd_nfp(i))=np_sc(nd_nfp(i))+1
        s_sc(nd_nfp(i))=snfp(i)
      enddo

      deallocate(nd_sc)

      if(num_tr.gt.maxntr-numnfp) then
        write(6,*) ' Error in buildtree():Node is overflow!'
        write(6,*) ' This level is ',level,'numnfp=',numnfp
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
          talist(npn)=num_tr
          node(npn)=pare_sc(i)
          clist(npn)=c_sc(i)
          slist(npn)=s_sc(i)
          num_tr=num_tr+1
          npn=npn+1
        endif
      enddo

      deallocate(s_sc)

! *** setting daughter_tr ***
! *** npn is the number of new nodes
      do i=0,npn-1
! *** flag for parent node, initially 0 
        if(flag_pc(node(i)).eq.0) then
! *** update daughter node name (now _tr id) only by first one
          daughter_tr(node(i))=talist(i)
          flag_pc(node(i))=1
        endif
      enddo

      allocate(nd_sc(0:numsc-1))

! *** set up the other values for _tr ***
      do i=0,npn-1
! *** nalist is subcell id, talist is _tr id
        nd_sc(nalist(i))=talist(i)
        np_tr(talist(i))=np_sc(nalist(i))
        l_tr(talist(i))=l_sc2
! *** node is parent _tr id
        pare_tr(talist(i))=node(i)
        call phixyzp(slist(i),clist(i),ix,iy,iz,level)
        cx_tr(talist(i))=dx(ix)*l_sc+cx_tr(node(i))
        cy_tr(talist(i))=dx(iy)*l_sc+cy_tr(node(i))
        cz_tr(talist(i))=dx(iz)*l_sc+cz_tr(node(i))
      enddo		

      deallocate(talist)
      deallocate(nalist)
      deallocate(node)
      deallocate(clist)
      deallocate(slist)

! *** Set label not-finished particle ***
      do i=0,numnfp-1
! *** nd_nfp from subcell id to _tr id ***
        nd_nfp(i)=nd_sc(nd_nfp(i))
        pn_tr(nd_nfp(i))=pn_nfp(i)
        if(np_tr(nd_nfp(i)).eq.1) then
! *  this node is leaf *
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
      numtmp = num_tr-stsc

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
! *** create subcell, numsc -> subcell id -> daughter_tr
!     for connection to parent node ***
        if(np_tr(stsc).ge.2) then
          daughter_tr(stsc)=numsc
          flag_pc(stsc) = 0
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
        else if(np_tr(stsc).eq.1) then
          daughter_tr(stsc)=-1
        endif
        stsc=stsc+1
      enddo
! *** store starting _tr id 
      stsc = num_tr
      if(numsc.gt.maxntr) then
        write(6,*) ' Error in buildtree() : '
        write(6,*) '  Subcell is overflow	!'
        write(6,*) '  This level is ',level
        stop
      endif
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

!      write(6,*) 'myrank,level,num=',myrank,level,num_tr

! *** set next node for node numtr0=1 (always go back to node 0, then look at node 0s next ***
      next_tr(numtr0)=0
      if(num_tr.gt.1) then

        allocate(pare_sc(0:num_tr-2))

        do i=numtr0,num_tr-2
          pare_sc(i)=pare_tr(i+1)
        enddo		
        do i=numtr0+1,num_tr-2
! *** set next_tr, compare the parent of the next tr ***
          if(pare_sc(i).eq.pare_tr(i)) then 
            next_tr(i)=i+1
          else
            next_tr(i) = next_tr(pare_tr(i))
          endif
        enddo
        next_tr(num_tr-1)=next_tr(pare_tr(num_tr-1))

        deallocate(pare_sc)

      endif

! *** compute mass and center of mass ***
      allocate(lvlist(0:level))
      allocate(paretrcm(0:num_tr-1))
      do i=0,level
        lvlist(i)=list(i)
      enddo
      do i=0,num_tr-1
        paretrcm(i)=pare_tr(i)
      enddo
      call compute_mass(numtr0,level,lvlist,paretrcm)
      deallocate(lvlist)
      deallocate(paretrcm)

! *** send the data to the other nodes ***
! *** get number of subdomains from each proc ***
   91 do i=0,num_tr-1
! *** set procid ***
        proc_tr(i)=myrank
      enddo 

! *** adjust cx,cy,cz,l_tr ***
      num = num_tr-1
      do j=level-1,0,-1
        do i=num,list(j),-1
          if(daughter_tr(i).gt.0.and.(np_tr(i).eq.np_tr(daughter_tr(i)))) then
            nd=daughter_tr(i)
            l_tr(i)=l_tr(nd)
            cx_tr(i)=cx_tr(nd)
            cy_tr(i)=cy_tr(nd)
            cz_tr(i)=cz_tr(nd)
            daughter_tr(i)=daughter_tr(nd)
          endif
        enddo
        num=list(j)-1
      enddo
! set mass, cx and etc. for node 0
      daughter_tr(0)=daughter_tr(1)
      mass_tr(0)=mass_tr(1)
      np_tr(0)=np_tr(1)
      pn_tr(0)=pn_tr(1)
      cmx_tr(0)=cmx_tr(1)
      cmy_tr(0)=cmy_tr(1)
      cmz_tr(0)=cmz_tr(1)
      hm_tr(0)=hm_tr(1)
      l_tr(0)=l_tr(1)
      cx_tr(0)=cx_tr(1)
      cy_tr(0)=cy_tr(1)
      cz_tr(0)=cz_tr(1)
#ifdef QPTREEF
      mx_tr(0)=mx_tr(1)
      my_tr(0)=my_tr(1)
      mz_tr(0)=mz_tr(1)
      mxx_tr(0)=mxx_tr(1)
      myy_tr(0)=myy_tr(1)
      mzz_tr(0)=mzz_tr(1)
      mxy_tr(0)=mxy_tr(1)
      myz_tr(0)=myz_tr(1)
      mzx_tr(0)=mzx_tr(1)
#endif

! *** nodess should be starting from 1, if numgtr.eq.0, corrected later
   94 nodese_tr=1      
      if(nprocs.gt.1) then
        numtr0=num_tr
        numtrs=0
        nodess=1
! *** list(0)=1=numtr0 ***
        i=0
        do i=0,level-2
!          if(nodess.eq.0.and.np_tr(list(i)).ne.np_tr(list(i+1))) then
!            nodess=list(i)
!          endif
          if(list(i+2)-nodess+1.gt.MAXNODESEND) then
! *** sending up to level i+1
            goto 92
          endif
        enddo
        if(i.eq.level-1.and.num_tr-nodess+1.gt.MAXNODESEND) then
          i=level-1
! *** sending up to level-1 ***
          goto 92
        endif
! *** sending all the node ***
        nodese_tr=num_tr-1
        goto 93
   92   nodese_tr=list(i+1)-1

!        write(6,*) ' myrank,nodess,nodese,np=',myrank,nodess,nodese_tr,np

   93   if(num_tr.eq.0) then
          numtrs=0
          nodess=0
          nodese_tr=0
        else if(nodess.gt.nodese_tr) then
          write(6,*) ' Error in buildtree(): nodess > nodese'
          write(6,*) ' nodess,nodese=',nodess,nodese_tr
          write(6,*) ' rank=',myrank,' np,numgtr=',np,num_tr
          stop
        else
          numtrs=nodese_tr-nodess+1
        endif

        allocate(npjr(0:nprocs-1))

        do ip=0,nprocs-1
          if(ip.eq.myrank) then
            npjr(ip)=numtrs
          endif
          call MPI_BCAST(npjr(ip),1,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
        enddo

        allocate(nd0(0:maxntr-1))

        nvali=4
#ifdef QPTREEF
        nval=18
#else
        nval=9
#endif
        do ip=0,nprocs-1
          if(ip.eq.myrank) then

            allocate(tivr(0:numtrs*nvali-1))
            allocate(tdvr(0:numtrs*nval-1))

            do i=0,numtrs-1
              nd=nodess+i
              tivr(i)=np_tr(nd)
              tivr(i+numtrs)=next_tr(nd)
              tivr(i+numtrs*2)=daughter_tr(nd)
              if(tivr(i+numtrs*2).gt.nodese_tr) then
                tivr(i+numtrs*2)=-1
              endif
              tivr(i+numtrs*3)=nd
! *** double precision values for neighbour search ***
              tdvr(i)=l_tr(nd)
              tdvr(i+numtrs)=cx_tr(nd)
              tdvr(i+numtrs*2)=cy_tr(nd)
              tdvr(i+numtrs*3)=cz_tr(nd)
              tdvr(i+numtrs*4)=hm_tr(nd)
              tdvr(i+numtrs*5)=mass_tr(nd)
              tdvr(i+numtrs*6)=cmx_tr(nd)
              tdvr(i+numtrs*7)=cmy_tr(nd)
              tdvr(i+numtrs*8)=cmz_tr(nd)
#ifdef QPTREEF
              tdvr(i+numtrs*9)=mx_tr(nd)
              tdvr(i+numtrs*10)=my_tr(nd)
              tdvr(i+numtrs*11)=mz_tr(nd)
              tdvr(i+numtrs*12)=mxx_tr(nd)
              tdvr(i+numtrs*13)=myy_tr(nd)
              tdvr(i+numtrs*14)=mzz_tr(nd)
              tdvr(i+numtrs*15)=mxy_tr(nd)
              tdvr(i+numtrs*16)=myz_tr(nd)
              tdvr(i+numtrs*17)=mzx_tr(nd)
#endif
            enddo
          else 

            allocate(tivr(0:npjr(ip)*nvali-1))
            allocate(tdvr(0:npjr(ip)*nval-1))

          endif
          call MPI_BCAST(tivr,npjr(ip)*nvali,MPI_INTEGER,ip &
           ,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(tdvr,npjr(ip)*nval,MPI_DOUBLE_PRECISION,ip &
           ,MPI_COMM_WORLD,ierr)

! *** should have all the pseudo node even if np=0, because DM may need it
          if(ip.ne.myrank) then
! *** add tree ***
            nrnd=npjr(ip)
            numtri=num_tr
            do i=0,nrnd-1
              np_tr(num_tr)=tivr(i)
              pn_tr(num_tr)=0
! *** set the other parameter for node ***
              proc_tr(num_tr)=ip
              next_tr(num_tr)=tivr(i+nrnd)
              daughter_tr(num_tr)=tivr(i+nrnd*2)
! *** link between original node id and id in this proc ***
              nd0(tivr(i+nrnd*3))=num_tr
              pare_tr(num_tr)=0
! *** adding nodes from the other procs
              l_tr(num_tr)=tdvr(i)
              cx_tr(num_tr)=tdvr(i+nrnd)
              cy_tr(num_tr)=tdvr(i+nrnd*2)
              cz_tr(num_tr)=tdvr(i+nrnd*3)
              hm_tr(num_tr)=tdvr(i+nrnd*4)
              mass_tr(num_tr)=tdvr(i+nrnd*5)
              cmx_tr(num_tr)=tdvr(i+nrnd*6)
              cmy_tr(num_tr)=tdvr(i+nrnd*7)
              cmz_tr(num_tr)=tdvr(i+nrnd*8)
#ifdef QPTREEF
              mx_tr(num_tr)=tdvr(i+nrnd*9)
              my_tr(num_tr)=tdvr(i+nrnd*10)
              mz_tr(num_tr)=tdvr(i+nrnd*11)
              mxx_tr(num_tr)=tdvr(i+nrnd*12)
              myy_tr(num_tr)=tdvr(i+nrnd*13)
              mzz_tr(num_tr)=tdvr(i+nrnd*14)
              mxy_tr(num_tr)=tdvr(i+nrnd*15)
              myz_tr(num_tr)=tdvr(i+nrnd*16)
              mzx_tr(num_tr)=tdvr(i+nrnd*17)
#endif
              num_tr=num_tr+1
            enddo
            if(nrnd.gt.0) then
              noders=tivr(nrnd*3)
              nodere=tivr(nrnd*4-1)
            else
              noders=num_tr
              nodere=num_tr
            endif

            do i=numtri,num_tr-1
! *** set next_tr ***
              if(next_tr(i).ge.noders.and.next_tr(i).le.nodere) then
                next_tr(i)=nd0(next_tr(i))
              else
                next_tr(i)=num_tr
              endif
! *** set daughter_tr ***     
              if(daughter_tr(i).ne.-1) then
                daughter_tr(i)=nd0(daughter_tr(i))
              endif
            enddo
          endif

          deallocate(tivr)
          deallocate(tdvr)

        enddo

        deallocate(npjr)
        deallocate(nd0)

! *** check there is any imported nodes if yes, change next_tr(0) ***
        if(numtr0.ne.num_tr) then
          if(numtr0.ne.0) then
            next_tr(0)=numtr0
          endif
! *** set end of tree ***
          do i=numtri,num_tr-1
! *** set next_tr ***
            if(next_tr(i).eq.num_tr) then
              next_tr(i)=-1
            endif
          enddo
        endif
      endif

! *** compute delta ***
      do i=0,num_tr-1
        delta_tr(i)=dsqrt((cmx_tr(i)-cx_tr(i))*(cmx_tr(i)-cx_tr(i)) &
         +(cmy_tr(i)-cy_tr(i))*(cmy_tr(i)-cy_tr(i)) &
         +(cmz_tr(i)-cz_tr(i))*(cmz_tr(i)-cz_tr(i)))
      enddo	

!      write(fileo,'(a4,i3.3,a1,i3.3)') 'tree',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,num_tr-1
!        write(60,'(7I10,9(1pE13.5))') i,np_tr(i),pn_tr(i) &
!        ,next_tr(i),daughter_tr(i),pare_tr(i),proc_tr(i),l_tr(i) &
!        ,cx_tr(i),cy_tr(i),cz_tr(i),hm_tr(i) &
!        ,mass_tr(i),cmx_tr(i),cmy_tr(i),cmz_tr(i) 
!        ,mx_tr(i),my_tr(i),mz_tr(i) &
!        ,mxx_tr(i),myy_tr(i),mzz_tr(i) &
!        ,mxy_tr(i),myz_tr(i),mzx_tr(i),delta_tr(i)
!      enddo
!      close(60)
!      stop

      deallocate(list)
      deallocate(pare_tr)	

end subroutine
#endif
#endif

#if defined(GAS) || defined(STAR)
! *** Definition of makeroot() ***
subroutine makeroot(np,pn_nfp,cz)
      use gcdp_const
      use gcdp_baryon
      use gcdp_btree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np
      integer,intent(in) :: pn_nfp(0:np-1)
      double precision,intent(in) :: cz

      integer i,pn,ierr
! * max coordinate *	  
      double precision max_x,max_y,max_z
#ifndef CENTRE0
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
! *** set z plane away from disk plane
      cz_tr(0)=cz
      do i=0,np-1
        pn=pn_nfp(i)
        if(dabs(x_p(pn)).gt.max_x) then
          max_x = dabs(x_p(pn))
        endif
        if(dabs(y_p(pn)).gt.max_y) then
          max_y = dabs(y_p(pn))
        endif
        if(dabs(z_p(pn)-cz).gt.max_z) then
          max_z = dabs(z_p(pn)-cz)
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
! *** check maxx range ***
      if(max_x.gt.max_y) then
        maxl=max_x
      else
        maxl=max_y
      endif
      if(max_z.gt.maxl) then
        maxl=max_z
      endif
      maxl=2.0d0*maxl
! *** Set the info for node "flag" ***
      l_tr(0)=MGROOT*maxl
      cx_tr(0)=0.0d0
      cy_tr(0)=0.0d0
      np_tr(0)=np
      pn_tr(0)=pn_nfp(0)
#else
      max_x=-INF
      max_y=-INF
      max_z=-INF
      min_x=INF
      min_y=INF
      min_z=INF
      do i=0,np-1
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
      l_tr(0)=MGROOT*maxl
      cx_tr(0)=(max_x+min_x)*0.5d0
      cy_tr(0)=(max_y+min_y)*0.5d0
      cz_tr(0)=(max_z+min_z)*0.5d0
      pn_tr(0)=pn_nfp(0)
      np_tr(0)=np
#endif

end subroutine
#endif

#ifdef SGRAV
#if defined(GAS) || defined(STAR)
! *** Definition of compute_mass() ***
! * This function compute the center of mass and total mass *
! * tree, center of mass, total mass, node No. *	 	  
subroutine compute_mass(numtr0,level,list,pare_tr)
      use gcdp_const
      use gcdp_baryon
      use gcdp_btree

      implicit none

      integer,intent(in) :: numtr0,level
      integer,intent(in) :: list(0:level),pare_tr(0:num_tr-1)
      integer i
      integer j,d,npare,num,nd,pnd
      double precision xsd,ysd,zsd
! for work
      integer,allocatable :: pare_sc(:),nd_sc(:)
      
! * compute total mass in node *
      do i=numtr0,num_tr-1
        if(np_tr(i).eq.1) then
          mass_tr(i)=m_p(pn_tr(i))
          cmx_tr(i)=m_p(pn_tr(i))*x_p(pn_tr(i))
          cmy_tr(i)=m_p(pn_tr(i))*y_p(pn_tr(i))
          cmz_tr(i)=m_p(pn_tr(i))*z_p(pn_tr(i))
          hm_tr(i)=m_p(pn_tr(i))*h_p(pn_tr(i))
        else 
          mass_tr(i)=0.0d0
          cmx_tr(i)=0.0d0
          cmy_tr(i)=0.0d0
          cmz_tr(i)=0.0d0
          hm_tr(i)=0.0d0
        endif
      enddo

      num=num_tr-1	
      do j=level,1,-1
! *** counting number of the parent nodes ***
        npare = 0

        if(num.ge.list(j)) then
          allocate(pare_sc(0:num-list(j)))
          allocate(nd_sc(0:num-list(j)))    
        endif

        do i=num,list(j),-1
          if(pare_tr(i).ne.pare_tr(i-1)) then
            pare_sc(npare) = pare_tr(i)
            nd_sc(npare) = i		  
            npare=npare+1
          endif
        enddo
        do d=0,7
! *** do loop for number of parent node.
          do i=0,npare-1
! *** node for potential child nodes ***
            nd=nd_sc(i)+d
! *** id for parent node ***
            pnd=pare_sc(i)  
            if(pnd.eq.pare_tr(nd).and.nd.le.num) then
              mass_tr(pnd)=mass_tr(pnd)+mass_tr(nd)
              cmx_tr(pnd)=cmx_tr(pnd)+cmx_tr(nd)
              cmy_tr(pnd)=cmy_tr(pnd)+cmy_tr(nd)
              cmz_tr(pnd)=cmz_tr(pnd)+cmz_tr(nd)
              hm_tr(pnd)=hm_tr(pnd)+hm_tr(nd)
            endif
          enddo
        enddo		  

        if(num.ge.list(j)) then
          deallocate(pare_sc)
          deallocate(nd_sc)
        endif

        num=list(j)-1
      enddo		
! * compute center of mass *
      do i=numtr0,num_tr-1
        cmx_tr(i)=cmx_tr(i)/mass_tr(i)
        cmy_tr(i)=cmy_tr(i)/mass_tr(i)
        cmz_tr(i)=cmz_tr(i)/mass_tr(i)
        hm_tr(i)=hm_tr(i)/mass_tr(i)      
      enddo
#ifdef QPTREEF
! * compute Multipole Momentum *
      do i=numtr0,num_tr-1
        mx_tr(i)=0.0d0
        my_tr(i)=0.0d0
        mz_tr(i)=0.0d0
        mxx_tr(i)=0.0d0
        myy_tr(i)=0.0d0
        mzz_tr(i)=0.0d0
        mxy_tr(i)=0.0d0
        myz_tr(i)=0.0d0
        mzx_tr(i)=0.0d0 
      enddo
      num = num_tr-1	
      do j=level,1,-1
        npare = 0

        if(num.ge.list(j)) then
          allocate(pare_sc(0:num-list(j)))
          allocate(nd_sc(0:num-list(j)))    
        endif

        do i=num,list(j),-1
          if(pare_tr(i).ne.pare_tr(i-1)) then
            pare_sc(npare) = pare_tr(i)
            nd_sc(npare) = i
            npare=npare+1
          endif
        enddo
        do d=0,7
          do i=0,npare-1
            pnd=pare_sc(i)
            nd=nd_sc(i)+d
            if(pnd.eq.pare_tr(nd).and.nd.le.num) then
              xsd=cmx_tr(nd)-cmx_tr(pnd)
              ysd=cmy_tr(nd)-cmy_tr(pnd)
              zsd=cmz_tr(nd)-cmz_tr(pnd)
              mx_tr(pnd)=mx_tr(pnd)+(mx_tr(nd)-xsd*mass_tr(nd))
              my_tr(pnd)=my_tr(pnd)+(my_tr(nd)-ysd*mass_tr(nd))
              mz_tr(pnd)=mz_tr(pnd)+(mz_tr(nd)-zsd*mass_tr(nd))
              mxx_tr(pnd)=mxx_tr(pnd)+(mxx_tr(nd)-2.0d0*xsd*mx_tr(nd) &
               +xsd*xsd*mass_tr(nd))
              myy_tr(pnd)=myy_tr(pnd)+(myy_tr(nd)-2.0d0*ysd*my_tr(nd) &
               +ysd*ysd*mass_tr(nd))
              mzz_tr(pnd)=mzz_tr(pnd)+(mzz_tr(nd)-2.0d0*zsd*mz_tr(nd) &
               +zsd*zsd*mass_tr(nd))
              mxy_tr(pnd)=mxy_tr(pnd)+(mxy_tr(nd)-xsd*my_tr(nd)-ysd*mx_tr(nd) &
               +xsd*ysd*mass_tr(nd))
              myz_tr(pnd)=myz_tr(pnd)+(myz_tr(nd)-ysd*mz_tr(nd)-zsd*my_tr(nd) &
               +ysd*zsd*mass_tr(nd))
              mzx_tr(pnd)=mzx_tr(pnd)+(mzx_tr(nd)-zsd*mx_tr(nd)-xsd*mz_tr(nd) &
               +zsd*xsd*mass_tr(nd))
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
#endif


