#include "gcdp.def"
! ***********************************************
!   ddecdm.F95 for GCD+ ver. f03.0
!  04 Feb. 2013    written by D.KAWATA
! ***********************************************

! ***********************************************
!  Domain decomposition for DM particles
!  using Paeno-Hilbert ordering 
! ***********************************************

#ifdef DM
subroutine ddecdm(ndmt,ndm,nadm,step)
      use gcdp_const
      use gcdp_dm
      use gcdp_dmtree
#ifdef MULTI_LEV
      use gcdp_ldmtree
#endif
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: ndmt,step
      integer,intent(inout) :: ndm,nadm
      integer i,j,k,level,maxlevel,level_mergin,ip,nin,ierr,pn
#ifdef FIXEDDMP
      integer nc
#endif
      integer idmr(0:1),idmrjr(0:1),idmrjs(0:1),it
! * diameter of subcell *	  
      double precision l_sc
! * center of the simulation box *	  
      double precision cx0,cy0,cz0
! *** for parallel process ***
      integer ndmj,nval
      integer crank,isend,jsta,jend,npnord
      integer ireqs(0:3),ireqr(0:3)
      double precision dx(0:1)
      data dx/1.0d0,-1.0d0/
! *** to find active particles ***
      integer na,nanp
      double precision lwt,upt
! *** for test ***
      character fileo*60
! for work
      integer,allocatable :: srank(:),rrank(:)
      integer,allocatable :: jstap(:),jendp(:)
      integer,allocatable :: npsp(:),nprp(:)
      integer,allocatable :: istatus(:),istatusr(:)
      integer,allocatable :: ord(:),pnord(:,:),ord1(:),nd_nfp(:)
      integer,allocatable :: list(:),talist(:),nalist(:)
      integer,allocatable :: tivs(:),tivr(:)
      integer*8,allocatable :: phkeyr(:),phkeys(:)
      integer*8,allocatable :: phkey(:),phkey1(:)
      double precision,allocatable :: tdvs(:),tdvr(:),cx0p(:),cy0p(:),cz0p(:)
! *** for Paeno-Hilbert ordering ***
      integer,allocatable :: ixp(:),iyp(:),izp(:),snfp(:),c_nfp(:)
! subcell
      integer numsc,pnumsc,nscnf
      integer,allocatable :: npsc(:),dpsc(:),scp(:),nppsc(:)
! for counting the time
      double precision ntime,ctime,c1time

! *** mergin of tree level to make a key for the particles ***
      level_mergin = 1

      allocate(ord(0:ndm))
      allocate(istatus(MPI_STATUS_SIZE))
      allocate(istatusr(MPI_STATUS_SIZE))

#ifdef MULTI_LEV

      allocate(pnord(0:ndm,0:1))

      do it=0,1
#else

      allocate(pnord(0:ndm,0:0))

      do it=0,0
#endif
        if(it.eq.0) then

          if(.not.allocated(np_dmtr)) then
            call allocate_dmtree(SI_ndm1)
          endif

          call dmmakeroot(SI_ndm1,0.0d0)
          idmr(0)=0
          idmr(1)=SI_ndm1-1
          cx0=cx_dmtr(0)
          cy0=cy_dmtr(0)
          cz0=cz_dmtr(0)
          l_sc=l_dmtr(0)*0.5d0	
#ifdef MULTI_LEV
        else
! *** for low-reso particles SI_ndm1-ndm ***

          if(.not.allocated(np_ldmtr)) then
            call allocate_ldmtree(ndm-SI_ndm1)
          endif

          call ldmmakeroot(ndm)
          idmr(0)=SI_ndm1
          idmr(1)=ndm-1
          cx0=cx_ldmtr(0)
          cy0=cy_ldmtr(0)
          cz0=cz_ldmtr(0)
          l_sc=l_ldmtr(0)*0.5d0	
#endif
        endif
        maxlevel = 0
! * initialization *
        allocate(phkey(idmr(0):idmr(1)+1))
        allocate(c_nfp(idmr(0):idmr(1)+1))
        allocate(snfp(idmr(0):idmr(1)+1))
        allocate(cx0p(idmr(0):idmr(1)+1))
        allocate(cy0p(idmr(0):idmr(1)+1))
        allocate(cz0p(idmr(0):idmr(1)+1))

! time
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        ntime=MPI_WTIME()

        do i=idmr(0),idmr(1)
          phkey(i)=0
          cx0p(i)=cx0
          cy0p(i)=cy0
          cz0p(i)=cz0
          c_nfp(i)=0
          snfp(i)=0
        enddo
        level = 0

! allocate daughter subcell
        numsc=NDTREE

        allocate(npsc(0:numsc))
        allocate(dpsc(0:numsc))

        do i=0,NDTREE-1
          npsc(i)=0
          dpsc(i)=0
        enddo

! record current subcell for particles
        
        allocate(scp(idmr(0):idmr(1)+1))

        do i=idmr(0),idmr(1)
          scp(i)=0
        enddo

! *****   start iteration in level *****
   77   level=level+1	  
        l_sc = 0.5d0*l_sc
! *** check for the level ***
        if(level.gt.ndmt.and.ndmt.gt.0) then
          write(6,*) ' Error in ddecdm(): too many tree level is required'
          write(6,*) ' myrank,ndm,it=',myrank,ndm,it
          write(6,*) ' level,ndmt=',level,ndmt
          stop
        endif
! * find out which subcell it is in *

        allocate(ixp(idmr(0):idmr(1)+1))
        allocate(iyp(idmr(0):idmr(1)+1))
        allocate(izp(idmr(0):idmr(1)+1))

        do i=idmr(0),idmr(1)
          if(x_dm(i)-cx0p(i).ge.0.0d0) then
            ixp(i)=0 
          else
             ixp(i)=1
          endif 
          if(y_dm(i)-cy0p(i).ge.0.0d0) then
            iyp(i)=0 
          else
            iyp(i)=1
          endif 
          if(z_dm(i)-cz0p(i).ge.0.0d0) then
            izp(i)=0 
          else
            izp(i)=1
          endif 
        enddo
! *** set c and state ***

        call phcurven(idmr(0),idmr(1),snfp,c_nfp &
         ,ixp,iyp,izp,level)

        do i=idmr(0),idmr(1)
! *** asign the ph key ***
! *** has to be phkey(i)*NDPH 
          phkey(i)=phkey(i)*NDPH+c_nfp(i)
! *** update cx for the cell ***
          cx0p(i)=dx(ixp(i))*l_sc+cx0p(i)
          cy0p(i)=dx(iyp(i))*l_sc+cy0p(i)
          cz0p(i)=dx(izp(i))*l_sc+cz0p(i)
        enddo

        deallocate(ixp)
        deallocate(iyp)
        deallocate(izp)

        do i=idmr(0),idmr(1)
! allocate subcell
          scp(i)=dpsc(scp(i))+c_nfp(i)
! counting number of particles in subcell
          npsc(scp(i))=npsc(scp(i))+1
        enddo

! make new subcell
        pnumsc=numsc
        numsc=0
        nscnf=0
! counting numsc for the next level
        do i=0,pnumsc-1
          if(npsc(i).gt.NDTREE) then
            nscnf=nscnf+1
          endif
          if(npsc(i).gt.0) then
            do j=0,NDTREE-1
              numsc=numsc+1
            enddo
          endif
        enddo

        allocate(nppsc(0:pnumsc))

        do i=0,pnumsc-1
          nppsc(i)=npsc(i)
        enddo

        deallocate(dpsc)
        deallocate(npsc)

        allocate(dpsc(0:pnumsc))
        allocate(npsc(0:numsc))

        numsc=0
        do i=0,pnumsc-1
          if(nppsc(i).gt.0) then
            dpsc(i)=numsc
            do j=0,NDTREE-1
              npsc(numsc)=0
              numsc=numsc+1
            enddo
          endif
        enddo

        deallocate(nppsc)

        if(maxlevel.eq.0) then
          if(nscnf.gt.0.and.numsc*NDTREE.lt.MAXNODE) then
! goto lower level
            goto 77
          endif

!          write(6,*) ' maxlev,myrank=',level,myrank
! *** get the max level among the processors ***
          call MPI_ALLREDUCE(level,maxlevel,1,MPI_INTEGER &
             ,MPI_MAX,MPI_COMM_WORLD,ierr)



        endif
        if(level.lt.maxlevel) then
          goto 77
        endif

        deallocate(c_nfp)
        deallocate(snfp)
        deallocate(cx0p)
        deallocate(cy0p)
        deallocate(cz0p)
        deallocate(dpsc)
        deallocate(npsc)
        deallocate(scp)

! *** ordering with indexxl ***
        npnord = idmr(1)-idmr(0)+1

        allocate(phkey1(1:npnord+1))
        allocate(ord1(1:npnord+1))


! time
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        ctime=MPI_WTIME()
!        if(myrank.eq.0) then
!          write(6,*) ' assign key=',ctime-ntime
!        endif   
!        ntime=ctime
  
        nin=0
        do i=idmr(0),idmr(1)
          nin=nin+1
          phkey1(nin)=phkey(i)
        enddo

        call indexxl(npnord,phkey1,ord1)

!        write(fileo,'(a5,i6.6)') 'orddm',myrank
!        open(60,file=fileo,status='unknown')
!        do i=1,npnord
!          pn=ord1(i)-1
!          write(60,'(3(1pE13.5),I10)') x_dm(pn),y_dm(pn),z_dm(pn) &
!           ,phkey1(pn)
!        enddo
!        close(60)

        deallocate(phkey1)

        do i=0,npnord-1
          pnord(i,it)=ord1(i+1)+idmr(0)-1
        enddo
        do i=0,npnord-1
          ord(pnord(i,it))=i
        enddo

! time
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        ctime=MPI_WTIME()
!        if(myrank.eq.0) then
!          write(6,*) ' sorting within the core=',ctime-ntime
!        endif   
!        ntime=ctime

        deallocate(ord1)
        allocate(phkeyr(idmr(0):idmr(1)+1))
        allocate(tivr(idmr(0):idmr(1)+1))
   
        do i=idmr(0),idmr(1)
          phkeyr(i)=phkey(i)
          tivr(i)=pnord(i-idmr(0),it)
        enddo
        idmrjr(0)=idmr(0)
        idmrjr(1)=idmr(1)

! *** define rank for sending and rank for recieving ***

        allocate(srank(0:0))
        allocate(rrank(0:0))

        srank(0) = myrank+1
        if(srank(0).ge.nprocs) then
          srank(0) = srank(0)-nprocs
        endif
        rrank(0) = myrank-1
        if(rrank(0).lt.0) then
          rrank(0)=rrank(0)+nprocs
        endif

! *** start ordering @ other PE *** 
        do j = 0,nprocs-1
          if(j.ne.0) then
            crank = myrank-j
            if(crank.lt.0) then
              crank = crank+nprocs
            endif
! *** new ordering using the sorted list ***
            nin=idmrjr(0)
            do k=0,npnord-1
              do i=nin,idmrjr(1)
                if(phkeyr(tivr(i)).gt.phkey(pnord(k,it))) then
                  goto 98
                else if(phkeyr(tivr(i)).eq.phkey(pnord(k,it)) &
                 .and.crank.ge.myrank) then 
                  goto 98
                endif
              enddo
   98         ord(pnord(k,it))=ord(pnord(k,it))+i-idmrjr(0)
              nin=i
            enddo
          endif

          allocate(tivs(idmrjr(0):idmrjr(1)+1))
          allocate(phkeys(idmrjr(0):idmrjr(1)+1))

! *** Message Passing ***
          do i=idmrjr(0),idmrjr(1)
            tivs(i) = tivr(i)
            phkeys(i) = phkeyr(i)
          enddo

          deallocate(tivr)
          deallocate(phkeyr)

          if(j.ne.nprocs-1) then

            idmrjs(0)=idmrjr(0)
            idmrjs(1)=idmrjr(1)

!          write(6,*) ' myrank,idmrjs=',myrank,idmrjs(0),idmrjs(1)
! *** range of the key and order ***
! *** cannot use MPI_GET_COUNT ***
            call MPI_ISEND(idmrjs,2,MPI_INTEGER,srank(0),1 &
             ,MPI_COMM_WORLD,ireqs(0),ierr)
            call MPI_IRECV(idmrjr,2,MPI_INTEGER,rrank(0),1 &
             ,MPI_COMM_WORLD,ireqr(0),ierr)
            call MPI_WAIT(ireqs(0),istatus,ierr)
            call MPI_WAIT(ireqr(0),istatusr,ierr)

            allocate(tivr(idmrjr(0):idmrjr(1)+1))
            allocate(phkeyr(idmrjr(0):idmrjr(1)+1))

! *** for order ***
            call MPI_ISEND(tivs(idmrjs(0)),idmrjs(1)-idmrjs(0)+1 &
             ,MPI_INTEGER,srank(0),2,MPI_COMM_WORLD,ireqs(1),ierr)
            call MPI_IRECV(tivr(idmrjr(0)),idmrjr(1)-idmrjr(0)+1 &
             ,MPI_INTEGER,rrank(0),2,MPI_COMM_WORLD,ireqr(1),ierr)
            call MPI_WAIT(ireqs(1),istatus,ierr)
            call MPI_WAIT(ireqr(1),istatusr,ierr)
! *** for phkey ***
            call MPI_ISEND(phkeys(idmrjs(0)),idmrjs(1)-idmrjs(0)+1 &
             ,MPI_INTEGER8,srank(0),3,MPI_COMM_WORLD,ireqs(2),ierr)
            call MPI_IRECV(phkeyr(idmrjr(0)),idmrjr(1)-idmrjr(0)+1 &
             ,MPI_INTEGER8,rrank(0),3,MPI_COMM_WORLD,ireqr(2),ierr)
            call MPI_WAIT(ireqs(2),istatus,ierr)
            call MPI_WAIT(ireqr(2),istatusr,ierr)
          endif

          deallocate(tivs)
          deallocate(phkeys)

        enddo

        deallocate(phkey)
        deallocate(srank)
        deallocate(rrank)

      enddo
! *** finish ordering for both high and low reso particles ***
! *** define srank and rrank ***

      allocate(srank(0:nprocs))
      allocate(rrank(0:nprocs))

! time
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      ctime=MPI_WTIME()
!      if(myrank.eq.0) then
!        write(6,*) ' global sorting =',ctime-ntime
!      endif   
!      ntime=ctime

      do j = 0,nprocs-1
        srank(j) = myrank+j
        if(srank(j).ge.nprocs) then
          srank(j) = srank(j)-nprocs
        endif
        rrank(j) = myrank-j
        if(rrank(j).lt.0) then
          rrank(j) = rrank(j)+nprocs
        endif
      enddo

! *** set order of low-resolution particles ***
      do i=0,ndm-SI_ndm1-1
        ord(pnord(i,1))=ord(pnord(i,1))+SI_ndm1t
      enddo
! *** sort talist for sending proc ***

      allocate(talist(0:ndm))
      allocate(jstap(0:nprocs))
      allocate(jendp(0:nprocs))
      allocate(npsp(0:nprocs))
      allocate(nprp(0:nprocs))

      na=0
      do ip=0,nprocs-1
        jstap(ip)=na
! for proc srank(ip)
        call para_range(0,SI_ndm1t-1,nprocs,srank(ip),jsta,jend)
        do i=0,SI_ndm1-1
          if(ord(pnord(i,0)).ge.jsta.and.ord(pnord(i,0)).le.jend) then
            talist(na)=pnord(i,0)
            na=na+1
          endif
        enddo
! *** set the asigned rank for low-reso particles ***
        call para_range(SI_ndm1t,ndmt-1,nprocs,srank(ip),jsta,jend)
        do i=0,ndm-SI_ndm1-1
          if(ord(pnord(i,1)).ge.jsta.and.ord(pnord(i,1)).le.jend) then
            talist(na)=pnord(i,1)
            na=na+1
          endif
        enddo
        jendp(ip)=na-1
        npsp(ip) = jendp(ip)-jstap(ip)+1
        nprp(ip) = 0
      enddo

      if(ndm.ne.na) then
        write(6,*) ' Error in ddecdm(): npnord not equal ndm'
        write(6,*) '  myrank,ndm,npnord=',myrank,ndm,na
        write(6,*) ' ndm1t,ndmt-ndm1t,ndmt=',SI_ndm1t,ndmt-SI_ndm1t,ndmt
        write(fileo,'(a9,i6.6)') 'errddecdm',myrank
        open(60,file=fileo,status='unknown')
        do i=0,SI_ndm1-1
          write(60,'(3(1pE13.5),2I10)') x_dm(i),y_dm(i),z_dm(i) &
           ,ord(pnord(i,0)),pnord(i,0)
        enddo
        close(60)
        write(fileo,'(a10,i6.6)') 'errddecdm1',myrank
        open(60,file=fileo,status='unknown')
        do i=SI_ndm1,ndm-1
          write(60,'(3(1pE13.5),2I10)') x_dm(i),y_dm(i),z_dm(i) &
           ,ord(pnord(i-SI_ndm1,1)),pnord(i-SI_ndm1,1)
        enddo
        close(60)        
        stop
      endif

! *** set pnord, jstap, jendp in order of asigned rank ***
      npnord=na 
      do i=0,ndm-1
        pnord(i,0)=talist(i)
      enddo 

      deallocate(talist)

! time
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      ctime=MPI_WTIME()
!      if(myrank.eq.0) then
!        write(6,*) ' assign the core to be sent =',ctime-ntime
!      endif   
!      ntime=ctime

! *** start data transfer ***
! get number of particles -> nprp and np
      ndm=0
      do j = 0,nprocs-1
        call MPI_ISEND(npsp(j),1,MPI_INTEGER,srank(j) &
         ,1,MPI_COMM_WORLD,ireqs(0),ierr)
        call MPI_IRECV(nprp(j),1,MPI_INTEGER,rrank(j) &
         ,1,MPI_COMM_WORLD,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)
        ndm=ndm+nprp(j)
      enddo

! *** ord ***

      allocate(tivs(0:npnord))

      do i=0,npnord-1
        tivs(i)=ord(pnord(i,0))
      enddo

      deallocate(ord)
      allocate(ord(0:ndm))

      ndm = 0
      do i=jstap(0),jendp(0)
        ord(ndm)=tivs(i)
        ndm=ndm+1
      enddo
      do j=1,nprocs-1
        call MPI_ISEND(tivs(jstap(j)),npsp(j),MPI_INTEGER,srank(j) &
         ,1,MPI_COMM_WORLD,ireqs(0),ierr)
        call MPI_IRECV(ord(ndm),nprp(j),MPI_INTEGER,rrank(j) &
         ,1,MPI_COMM_WORLD,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)
        ndm=ndm+nprp(j)
      enddo

      deallocate(tivs)

! *** sort ***
      call para_range(0,SI_ndm1t-1,nprocs,myrank,jsta,jend)
      ndmj=0
      do i=0,ndm-1
        if(ord(i).le.jend) then
          ord(i)=ord(i)-jsta
          ndmj=ndmj+1
        else if(ord(i).lt.jsta) then
          write(6,*) ' Error in ddecdm(): wrong ord myrank=',myrank
          write(6,*) ' ord,jsta,jend=',ord(i),jsta,jend 
        endif
      enddo
      call para_range(SI_ndm1t,ndmt-1,nprocs,myrank,jsta,jend)
      do i=0,ndm-1
        if(ord(i).ge.jsta) then
          ord(i)=ord(i)-jsta+ndmj
        else if(ord(i).gt.jend) then
          write(6,*) ' Error in ddecdm() for low-reso particles:'
          write(6,*) '   wrong ord myrank=',myrank
          write(6,*) ' ord,jsta,jend=',ord(i),jsta,jend 
        endif
      enddo

! *** integer values ***  
      nval=2

      allocate(tivr(0:ndm*nval))

      ndm=0
      do ip=0,nprocs-1
        ndmj=0

        allocate(tivs(0:npsp(ip)*nval))

        do i=jstap(ip),jendp(ip)
          tivs(ndmj)=flagt_dm(pnord(i,0))
          tivs(ndmj+npsp(ip))=id_dm(pnord(i,0))
          ndmj=ndmj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tivs,npsp(ip)*nval,MPI_INTEGER &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(0),ierr)
          call MPI_IRECV(tivr(ndm*nval),nprp(ip)*nval,MPI_INTEGER &
           ,rrank(ip),1,MPI_COMM_WORLD,ireqr(0),ierr)
          call MPI_WAIT(ireqs(0),istatus,ierr)
          call MPI_WAIT(ireqr(0),istatusr,ierr)
        else
          nprp(ip)=npsp(ip)
          do j=0,nval-1
            do i=0,nprp(ip)-1
              tivr(i+nprp(ip)*j)=tivs(i+npsp(ip)*j)
            enddo
          enddo
        endif

        deallocate(tivs)

        ndm=ndm+nprp(ip)
      enddo

! reallocate DM variables
      call reallocate_dm_int(ndm)

      ndmj=0
      ndm=0
      do ip=0,nprocs-1
        do i=ndm*nval,ndm*nval+nprp(ip)-1
          flagt_dm(ord(ndmj))=tivr(i)
          id_dm(ord(ndmj))=tivr(i+nprp(ip))
          ndmj=ndmj+1
        enddo
        ndm=ndm+nprp(ip)
      enddo

      deallocate(tivr)

! *** x?_dm,vn?_dm
      nval=6

      allocate(tdvr(0:ndm*nval))

      ndm=0
      do ip=0,nprocs-1
        ndmj=0

        allocate(tdvs(0:npsp(ip)*nval))

        do i=jstap(ip),jendp(ip)
          tdvs(ndmj)=xc_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip))=yc_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*2)=zc_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*3)=vnx_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*4)=vny_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*5)=vnz_dm(pnord(i,0))
          ndmj=ndmj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(ndm*nval),nprp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
        else
          nprp(ip)=npsp(ip)
          do j=0,nval-1         
            do i=0,nprp(ip)-1
              tdvr(i+nprp(ip)*j)=tdvs(i+npsp(ip)*j)
            enddo
          enddo
        endif

        deallocate(tdvs)

        ndm=ndm+nprp(ip)
      enddo

! reallocate DM variables
      call reallocate_dm_d1(ndm)

      ndmj=0
      ndm=0
      do ip=0,nprocs-1
        do i=ndm*nval,ndm*nval+nprp(ip)-1
          xc_dm(ord(ndmj))=tdvr(i)
          yc_dm(ord(ndmj))=tdvr(i+nprp(ip))
          zc_dm(ord(ndmj))=tdvr(i+nprp(ip)*2)
          vnx_dm(ord(ndmj))=tdvr(i+nprp(ip)*3)
          vny_dm(ord(ndmj))=tdvr(i+nprp(ip)*4)
          vnz_dm(ord(ndmj))=tdvr(i+nprp(ip)*5)
          ndmj=ndmj+1
        enddo
        ndm=ndm+nprp(ip)
      enddo

      deallocate(tdvr)

      if(ndm.ne.ndmj) then
        write(6,*) ' Error in ddecdm ndm,ndmj=',ndm,ndmj
        stop
      endif

! *** m_dm,h_dm,dv?_dm,dt_dm,lt_dm
      nval=7

      allocate(tdvr(0:ndm*nval))

      ndm=0
      do ip=0,nprocs-1
        ndmj=0

        allocate(tdvs(0:npsp(ip)*nval))

        do i=jstap(ip),jendp(ip)
          tdvs(ndmj)=m_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip))=h_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*2)=dvx_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*3)=dvy_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*4)=dvz_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*5)=dt_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*6)=lt_dm(pnord(i,0))
          ndmj=ndmj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(ndm*nval),nprp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
        else
          nprp(ip)=npsp(ip)
          do j=0,nval-1         
            do i=0,nprp(ip)-1
              tdvr(i+nprp(ip)*j)=tdvs(i+npsp(ip)*j)
            enddo
          enddo
        endif

        deallocate(tdvs)

        ndm=ndm+nprp(ip)
      enddo

! reallocate DM variables
      call reallocate_dm_d2(ndm)

      ndmj=0
      ndm=0
      do ip=0,nprocs-1
        do i=ndm*nval,ndm*nval+nprp(ip)-1
          m_dm(ord(ndmj))=tdvr(i)
          h_dm(ord(ndmj))=tdvr(i+nprp(ip))
          dvx_dm(ord(ndmj))=tdvr(i+nprp(ip)*2)
          dvy_dm(ord(ndmj))=tdvr(i+nprp(ip)*3)
          dvz_dm(ord(ndmj))=tdvr(i+nprp(ip)*4)
          dt_dm(ord(ndmj))=tdvr(i+nprp(ip)*5)
          lt_dm(ord(ndmj))=tdvr(i+nprp(ip)*6)
          ndmj=ndmj+1
        enddo
        ndm=ndm+nprp(ip)
      enddo

      deallocate(tdvr)

      if(ndm.ne.ndmj) then
        write(6,*) ' Error in ddecdm ndm,ndmj=',ndm,ndmj,'myrank=',myrank
        stop
      endif

! *** rho_dm,zetah_dm,omgh_dm,div_v_dm,dtr_dm
! *** because only active particles will update these values.
#if defined(COSM) && defined(BOXSIM)
      nval=8
#else
      nval=5
#endif

      allocate(tdvr(0:ndm*nval))

      ndm=0
      do ip=0,nprocs-1
        ndmj=0

        allocate(tdvs(0:npsp(ip)*nval)) 

        do i=jstap(ip),jendp(ip)
          tdvs(ndmj)=rho_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip))=zetah_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*2)=omgh_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*3)=div_v_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*4)=dtr_dm(pnord(i,0))
#if defined(COSM) && defined(BOXSIM)
          tdvs(ndmj+npsp(ip)*5)=hub_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*6)=asc_dm(pnord(i,0))
          tdvs(ndmj+npsp(ip)*7)=tc_dm(pnord(i,0))
#endif
          ndmj=ndmj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(ndm*nval),nprp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
        else
          nprp(ip)=npsp(ip)
          do j=0,nval-1         
            do i=0,nprp(ip)-1
              tdvr(i+nprp(ip)*j)=tdvs(i+npsp(ip)*j)
            enddo
          enddo
        endif

        deallocate(tdvs)

        ndm=ndm+nprp(ip)
      enddo

! reallocate DM variables
      call reallocate_dm_d3(ndm)

      ndmj=0
      ndm=0
      do ip=0,nprocs-1
        do i=ndm*nval,ndm*nval+nprp(ip)-1
          rho_dm(ord(ndmj))=tdvr(i)
          zetah_dm(ord(ndmj))=tdvr(i+nprp(ip))
          omgh_dm(ord(ndmj))=tdvr(i+nprp(ip)*2)
          div_v_dm(ord(ndmj))=tdvr(i+nprp(ip)*3)
          dtr_dm(ord(ndmj))=tdvr(i+nprp(ip)*4)
#if defined(COSM) && defined(BOXSIM)
          hub_dm(ord(ndmj))=tdvr(i+nprp(ip)*5)
          asc_dm(ord(ndmj))=tdvr(i+nprp(ip)*6)
          tc_dm(ord(ndmj))=tdvr(i+nprp(ip)*7)
#endif
          ndmj=ndmj+1
        enddo
        ndm=ndm+nprp(ip)
      enddo

      deallocate(tdvr)
      deallocate(istatus)
      deallocate(istatusr)

      if(ndm.ne.ndmj) then
        write(6,*) ' Error in ddecdm ndm,ndmj=',ndm,ndmj,'myrank=',myrank
        stop
      endif

! *** change the order to set high-reso particle first ***

      allocate(talist(0:ndm))
      allocate(nalist(0:ndm))
      allocate(tivr(0:ndm))
      allocate(tdvr(0:ndm))

      na = 0
      nanp = 0
      do i=0,ndm-1
        if(id_dm(i).lt.SI_ndm1t) then
          talist(na)=i
          na=na+1
        else
          nalist(nanp)=i
          nanp=nanp+1
        endif
      enddo
      SI_ndm1=na
!      write(6,*) 'rank',myrank,'ndm,ndm1=',ndm,SI_ndm1
      if(ndm.ne.na+nanp) then
        write(6,*) ' Error in ddecdm(): ndm is inconsistent at rank',myrank 
        write(6,*) '  ndm,n_high,low-reso=',ndm,na,nanp
        stop
      endif
      nanp=0
      do i=na,ndm-1
        talist(i)=nalist(nanp)
        nanp=nanp+1
      enddo
! *** id_dm ***
      do i=0,ndm-1
        tivr(i)=id_dm(i)
      enddo
      do i=0,ndm-1
        id_dm(i)=tivr(talist(i))
      enddo
! *** flagt_dm ***
      do i=0,ndm-1
        tivr(i)=flagt_dm(i)
      enddo
      do i=0,ndm-1
        flagt_dm(i)=tivr(talist(i))
      enddo
! *** xc_dm ***
      do i=0,ndm-1
        tdvr(i)=xc_dm(i)
      enddo
      do i=0,ndm-1
        xc_dm(i)=tdvr(talist(i))
      enddo
! *** yc_dm ***
      do i=0,ndm-1
        tdvr(i)=yc_dm(i)
      enddo
      do i=0,ndm-1
        yc_dm(i)=tdvr(talist(i))
      enddo
! *** zc_dm ***
      do i=0,ndm-1
        tdvr(i)=zc_dm(i)
      enddo
      do i=0,ndm-1
        zc_dm(i)=tdvr(talist(i))
      enddo
! *** vnx_dm ***
      do i=0,ndm-1
        tdvr(i)=vnx_dm(i)
      enddo
      do i=0,ndm-1
        vnx_dm(i)=tdvr(talist(i))
      enddo
! *** vny_dm ***
      do i=0,ndm-1
        tdvr(i)=vny_dm(i)
      enddo
      do i=0,ndm-1
        vny_dm(i)=tdvr(talist(i))
      enddo
! *** vnz_dm ***
      do i=0,ndm-1
        tdvr(i)=vnz_dm(i)
      enddo
      do i=0,ndm-1
        vnz_dm(i)=tdvr(talist(i))
      enddo
! *** vnz_dm ***
      do i=0,ndm-1
        tdvr(i)=m_dm(i)
      enddo
      do i=0,ndm-1
        m_dm(i)=tdvr(talist(i))
      enddo
! *** dvx_dm ***
      do i=0,ndm-1
        tdvr(i)=dvx_dm(i)
      enddo
      do i=0,ndm-1
        dvx_dm(i)=tdvr(talist(i))
      enddo
! *** dvy_dm ***
      do i=0,ndm-1
        tdvr(i)=dvy_dm(i)
      enddo
      do i=0,ndm-1
        dvy_dm(i)=tdvr(talist(i))
      enddo
! *** dvz_dm ***
      do i=0,ndm-1
        tdvr(i)=dvz_dm(i)
      enddo
      do i=0,ndm-1
        dvz_dm(i)=tdvr(talist(i))
      enddo
! *** h_dm ***
      do i=0,ndm-1
        tdvr(i)=h_dm(i)
      enddo
      do i=0,ndm-1
        h_dm(i)=tdvr(talist(i))
      enddo
! *** lt_dm ***
      do i=0,ndm-1
        tdvr(i)=lt_dm(i)
      enddo
      do i=0,ndm-1
        lt_dm(i)=tdvr(talist(i))
      enddo
! *** dt_dm ***
      do i=0,ndm-1
        tdvr(i)=dt_dm(i)
      enddo
      do i=0,ndm-1
        dt_dm(i)=tdvr(talist(i))
      enddo
! *** rho_dm ***
      do i=0,ndm-1
        tdvr(i)=rho_dm(i)
      enddo
      do i=0,ndm-1
        rho_dm(i)=tdvr(talist(i))
      enddo
! *** zetah_dm ***
      do i=0,ndm-1
        tdvr(i)=zetah_dm(i)
      enddo
      do i=0,ndm-1
        zetah_dm(i)=tdvr(talist(i))
      enddo
! *** omgh_dm ***
      do i=0,ndm-1
        tdvr(i)=omgh_dm(i)
      enddo
      do i=0,ndm-1
        omgh_dm(i)=tdvr(talist(i))
      enddo
! *** div_v_dm ***
      do i=0,ndm-1
        tdvr(i)=div_v_dm(i)
      enddo
      do i=0,ndm-1
        div_v_dm(i)=tdvr(talist(i))
      enddo
! *** dtr_dm ***
      do i=0,ndm-1
        tdvr(i)=dtr_dm(i)
      enddo
      do i=0,ndm-1
        dtr_dm(i)=tdvr(talist(i))
      enddo

      deallocate(talist)
      deallocate(nalist)
      deallocate(tivr)
      deallocate(tdvr)

      deallocate(srank)
      deallocate(rrank)
      deallocate(jstap)
      deallocate(jendp)
      deallocate(npsp)
      deallocate(nprp)


! *** after care ***
! *** make active particle list ***
      if(step.ge.0) then

        allocate(nalist(0:ndm))

! *** find active particles ***
        lwt =TM_lt+MGTL*TM_dt
        upt=TM_lt+MGTU*TM_dt
        nadm=0
        nanp=0
        do i=0,ndm-1
          if(lwt.lt.lt_dm(i)+dt_dm(i).and.lt_dm(i)+dt_dm(i).lt.upt) then
            list_adm(nadm)=i
            nadm=nadm+1
          else
            nalist(nanp)=i
            nanp=nanp+1
          endif
        enddo
! *** list_adm for non-active particles ***
        do i=0,nanp-1
          list_adm(i+nadm)=nalist(i)
        enddo

        deallocate(nalist)

! time
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      ctime=MPI_WTIME()
!      if(myrank.eq.0) then
!        write(6,*) ' data transfer =',ctime-ntime
!      endif   
!      ntime=ctime

! *** setup virtual timestep ***
        do i=0,ndm-1
#ifdef FIXEDDMP
        if(id_dm(i).lt.SI_fiddm(0).or.id_dm(i).gt.SI_fiddm(1)) then
#endif
          vdt_dm(i)=TM_lt+TM_dt-lt_dm(i)
#ifdef FIXEDDMP
          else
            vdt_dm(i)=0.0d0
          endif
#endif
        enddo
      else
! *** step = 0
        do i=0,ndm-1
          list_adm(i)=i
          x_dm(i)=xc_dm(i)
          y_dm(i)=yc_dm(i)
          z_dm(i)=zc_dm(i)
        enddo
        nadm = ndm

#ifdef FIXEDDMP
! *** change nag and nagravg

        allocate(nalist(0:ndm))

        nadm=0
        nc=0
        do i=0,ndm-1
          if(id_dm(i).ge.SI_fiddm(0).and.id_dm(i).le.SI_fiddm(1)) then
            nalist(nc)=i
            nc=nc+1
          else
            list(nadm)=i
            nadm=nadm+1
          endif
        enddo
! *** change list_ap putting fixedp at the end
        do i=0,nadm-1
          list_adm(i)=list(i)
        enddo
        do i=0,nc-1
          list_adm(nadm+i)=nalist(i)
        enddo

        deallocate(nalist)

#endif

      endif

#ifdef FIXEDDMP
      do i=0,ndm-1
        if(id_dm(i).ge.SI_fiddm(0).and.id_dm(i).le.SI_fiddm(1)) then
          x_dm(i)=xc_dm(i)
          y_dm(i)=yc_dm(i)
          z_dm(i)=zc_dm(i)
! *** set dvx and etc.
          dvx_dm(i)=0.0d0
          dvy_dm(i)=0.0d0
          dvz_dm(i)=0.0d0
        endif
      enddo
#endif

!      write(fileo,'(a6,i6.6)') 'ddecdm',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ndm-1
!        write(60,'(3(1pE13.5),2I10)') x_dm(i),y_dm(i),z_dm(i) &
!         ,ord(pnord(i,0)),pnord(i,0)
!      enddo
!      close(60)

end subroutine
#endif



