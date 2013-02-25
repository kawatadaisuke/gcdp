#include "gcdp.def"
! ***********************************************
!    ddecb.F95  for GCD+  ver. f03.0
!  03  Feb. 2013    written by D.KAWATA
! ***********************************************

! ***********************************************
!  Domain decomposition for baryon particles
!  using Paeno-Hilbert ordering 
! ***********************************************

#if defined(GAS) || defined(STAR)
subroutine ddecb(npt,np,ng,ns,nag,nas,nagravg,step,flagcont)
      use gcdp_const
      use gcdp_baryon
      use gcdp_btree
      use gcdp_system
#ifdef METAL
      use gcdp_metal
#endif
#ifdef IMSL
      use svrgp_init
      use umach_init
#endif

      implicit none
      include 'mpif.h'

      integer,intent(in) :: npt,step,flagcont
      integer,intent(inout) :: np,ng,ns,nag,nas,nagravg
      double precision kmmp
      integer i,j,k,level,maxlevel,level_mergin,level_limit,ip,nin,ierr
#ifdef FIXEDP  
      integer nc
#endif    
! *** diameter of subcell ***	  
      double precision l_sc
!  center of the simulation box 	  
      double precision cx0,cy0,cz0
! ** for parallel process ***
      integer npj,crank,isend
      integer nvali,nval,jsta,jend
      integer ireqs(0:3),ireqr(0:3)
      double precision dx(0:1)
      data dx/1.0d0,-1.0d0/
! *** for domain decomposition ver. 2 ***
! *** to find active particles ***
      integer pn,na,nanp
      double precision lwt,upt
! *** for check ***
      integer nptc
! *** for test ***
      character fileo*60
! for work
      integer,allocatable :: istatus(:),istatusr(:)
      integer,allocatable :: srank(:),rrank(:)
      integer,allocatable :: jstap(:),jendp(:)
      integer,allocatable :: npsp(:),nprp(:)
      integer,allocatable :: idisp(:),jjlen(:)
      integer,allocatable :: ord(:),ord1(:),pn_nfp(:),nd_nfp(:)
      integer,allocatable :: list(:),talist(:),nalist(:)
      integer,allocatable :: tivs(:),tivr(:)
      double precision,allocatable :: tdvs(:),tdvr(:)
! ** for Paeno-Hilbert ordering **
      integer,allocatable :: ixp(:),iyp(:),izp(:),snfp(:),c_nfp(:)
      integer*8,allocatable :: phkey(:)
      integer*8,allocatable :: phkeyr(:),phkeys(:)
#ifdef NUMREC
      integer*8,allocatable :: phkey1(:)
#endif
#ifdef IMSL
      double precision,allocatable :: phkey1d(:),sortkey1d(:)
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! * for feedback  common *
! * Total *
      double precision,allocatable :: tmej_sn(:),tnsn_sn(:),tmzHe_sn(:) &
       ,tmzC_sn(:),tmzN_sn(:),tmzO_sn(:),tmzNe_sn(:),tmzMg_sn(:) &
       ,tmzSi_sn(:),tmzFe_sn(:),tmzZ_sn(:),nsw_sn(:)
#endif

! *** mergin of tree level to make a key for the particles ***
      level_mergin = 1
! *** phkey < 2**64 ***
      level_limit=20-level_mergin

!      if(step.gt.0) then
!        write(fileo,'(a5,i3.3,i3.3)') 'pdecb',myrank,step
!      else
!        write(fileo,'(a5,i3.3,i3.3)') 'pdecb',myrank,0
!      endif
!       write(fileo,'(a5,i3.3)') 'pdecb',myrank
!      open(60,file=fileo,status='unknown')
!      write(60,*) ng,nag,np
!      do i=0,np-1
!        write(60,'(13(1pE13.5),2I10)') xc_p(i),yc_p(i),zc_p(i) &
!        ,vnx_p(i),vny_p(i),vnz_p(i),dvx_p(i),dvy_p(i),dvz_p(i) &
!        ,rho_p(i),vdt_p(i),dt_p(i),lt_p(i),id_p(i),i
!      enddo
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      close(60)
!      write(6,*) step,' ng,nag,ns,nas,nagravg,myrank=',ng,nag,ns,nas,nagravg,myrank

! ** Make root in tree.f **	
! *** since ver.33. need to list pn_nfp

      allocate(pn_nfp(0:np))      

      do i=0,np-1
        pn_nfp(i)=i
      enddo
      call makeroot(np,pn_nfp,0.0d0)
      cx0=cx_tr(0)
      cy0=cy_tr(0)
      cz0=cz_tr(0)
      maxlevel=0

      deallocate(pn_nfp)
      allocate(phkey(0:np))
      allocate(c_nfp(0:np))
      allocate(snfp(0:np))

!  initialization 
      do i=0,np-1
        phkey(i)=0
! *** use cx_tr for the position of the cell for the particles 
        cx_tr(i)=cx0
        cy_tr(i)=cy0
        cz_tr(i)=cz0
        c_nfp(i)=0
        snfp(i)=0
      enddo
      l_sc=l_tr(0)*0.5d0	
      level = 0
! ****   start iteration in level ****
   77 level=level+1	  
      l_sc = 0.5d0*l_sc
!  find out which subcell it is in 

      allocate(ixp(0:np))
      allocate(iyp(0:np))
      allocate(izp(0:np))

      do i=0,np-1
        if(x_p(i)-cx_tr(i).ge.0.0d0) then
          ixp(i)=0 
        else
          ixp(i)=1
        endif 
        if(y_p(i)-cy_tr(i).ge.0.0d0) then
          iyp(i)=0 
        else
          iyp(i)=1
        endif 
        if(z_p(i)-cz_tr(i).ge.0.0d0) then
          izp(i)=0 
        else
          izp(i)=1
        endif 
      enddo
      do i=0,np-1
! *** set c and state ***
        call phcurve(snfp(i),c_nfp(i),ixp(i),iyp(i),izp(i),level)
! *** asign the ph key ***
! *** has to be phkey(i)*NDPH 
        phkey(i)=phkey(i)*NDPH+c_nfp(i)
! *** update cx for the cell ***
        cx_tr(i)=dx(ixp(i))*l_sc+cx_tr(i)
        cy_tr(i)=dx(iyp(i))*l_sc+cy_tr(i)
        cz_tr(i)=dx(izp(i))*l_sc+cz_tr(i)
      enddo

      deallocate(ixp)
      deallocate(iyp)
      deallocate(izp)

! *** check the number of particles whose key overlap>8 ***
      if(maxlevel.eq.0) then
        if(level.ge.level_limit) then
          write(6,*) ' ddecb(): rank',myrank,' hit level limit'
          goto 90
        endif

        allocate(nd_nfp(0:np))

        do i=0,np-1
          nd_nfp(i)=0
          do j=i+1,np-1
            if(phkey(i).eq.phkey(j)) then
              nd_nfp(i)=nd_nfp(i)+1
            endif
          enddo
          if(nd_nfp(i).gt.int(NDPH)) then

            deallocate(nd_nfp)

            goto 77
          endif
        enddo

        deallocate(nd_nfp)

! *** get the max level among the processors ***
!        write(6,*) myrank,' np,level=',np,level
   90   call MPI_ALLREDUCE(level,maxlevel,1,MPI_INTEGER &
         ,MPI_MAX,MPI_COMM_WORLD,ierr)
        maxlevel = maxlevel+level_mergin
!        write(6,*) myrank,' level,maxlevel=',level,maxlevel
        goto 77
      else if(level+1.le.maxlevel) then
        goto 77
      endif

      deallocate(c_nfp)
      deallocate(snfp)

#ifdef NUMREC

      allocate(phkey1(1:np))
      allocate(ord1(1:np))

      do i=0,np-1
        phkey1(i+1)=phkey(i)
      enddo
      call indexxl(np,phkey1,ord1)

      deallocate(phkey1)

#endif
#ifdef IMSL

      allocate(phkey1d(1:np))
      allocate(sortkey1d(1:np))
      allocate(ord1(1:np))

      do i=1,np
        phkey1d(i)=dble(phkey(i-1))
        ord1(i)=i
      enddo
      call DSVRGP(np,phkey1d,sortkey1d,ord1)

      deallocate(phkey1d)
      deallocate(sortkey1d)

#endif

      allocate(list(0:np))
      allocate(ord(0:np))

      do i=0,np-1
        list(i)=ord1(i+1)-1
      enddo
      do i=0,np-1
        ord(list(i))=i
      enddo

!      write(fileo,'(a4,i3.3)') 'ord1',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np-1
!        write(60,'(3(1pE13.5),4I20)') x_p(i),y_p(i),z_p(i) &
!         ,id_p(i),phkey(i),ord(i),list(i)
!      enddo
!      close(60)

      deallocate(ord1)
      allocate(phkeyr(0:np))
      allocate(tivr(0:np))

      do i=0,np-1
        phkeyr(i)=phkey(i)
        tivr(i)=list(i)
      enddo
      npj=np
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
! ** start ordering @ other PE ** 

      allocate(istatus(MPI_STATUS_SIZE))
      allocate(istatusr(MPI_STATUS_SIZE))

      do j = 0,nprocs-1
        if(j.ne.0) then
          crank = myrank-j
          if(crank.lt.0) then
            crank = crank+nprocs
          endif
! *** new ordering using the sorted list ***
          nin=0
          do k=0,np-1
            do i=nin,npj-1
              if(phkeyr(tivr(i)).gt.phkey(list(k))) then
                goto 98
              else if(phkeyr(tivr(i)).eq.phkey(list(k)) &
               .and.crank.ge.myrank) then 
                goto 98
              endif
            enddo
   98       ord(list(k))=ord(list(k))+i
            nin=i
          enddo
        endif
! ** Message Passing **
        isend = npj

        allocate(tivs(0:npj))                
        allocate(phkeys(0:npj))                

        do i=0,isend-1
          tivs(i) = tivr(i)
          phkeys(i) = phkeyr(i)
        enddo
  
        deallocate(tivr)
        deallocate(phkeyr)

!        write(6,*) j,' myrank,isend,npj=',myrank,isend,npj

! get npj
        call MPI_ISEND(isend,1,MPI_INTEGER,srank(0),1,MPI_COMM_WORLD &
         ,ireqs(0),ierr)
        call MPI_IRECV(npj,1,MPI_INTEGER,rrank(0),1,MPI_COMM_WORLD &
         ,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)

!        write(6,*) ' myrank,isend,npj=',myrank,isend,npj

        if(j.ne.nprocs-1) then

          allocate(tivr(0:npj))
          allocate(phkeyr(0:npj))

! ** for order **
          call MPI_ISEND(tivs,isend,MPI_INTEGER,srank(0),2 &
           ,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tivr,npj,MPI_INTEGER,rrank(0),2 &
           ,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
! ** for phkey **
          call MPI_ISEND(phkeys,isend,MPI_INTEGER8,srank(0),3 &
           ,MPI_COMM_WORLD,ireqs(2),ierr)
          call MPI_IRECV(phkeyr,npj,MPI_INTEGER8,rrank(0),3 &
           ,MPI_COMM_WORLD,ireqr(2),ierr)
          call MPI_WAIT(ireqs(2),istatus,ierr)
          call MPI_WAIT(ireqr(2),istatusr,ierr)

!          write(fileo,'(a6,i3.3,i3.3)') 'phkeyr',myrank,rrank(0)
!          open(60,file=fileo,status='unknown')
!          do i=0,npj-1
!            write(60,'(2I20)') phkeyr(i),tivr(i)
!           enddo
!           close(60)

        endif

        deallocate(tivs)
        deallocate(phkeys)

      enddo

!      write(fileo,'(a3,i3.3)') 'ord',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np-1
!        write(60,'(3(1pE13.5),3I20)') x_p(i),y_p(i),z_p(i) &
!         ,id_p(i),phkey(i),ord(i)
!      enddo
!      close(60)
!      stop

      deallocate(phkey)

      deallocate(srank)
      deallocate(rrank)

! ** send and recieve data ***
! ** define srank and rrank **

      allocate(srank(0:nprocs))
      allocate(rrank(0:nprocs))

      do j=0,nprocs-1
        srank(j) = myrank+j
        if(srank(j).ge.nprocs) then
          srank(j) = srank(j)-nprocs
        endif
        rrank(j) = myrank-j
        if(rrank(j).lt.0) then
          rrank(j) = rrank(j)+nprocs
        endif
      enddo
! ** initialization **

      allocate(idisp(0:nprocs))
      allocate(jjlen(0:nprocs))

      do i=0,nprocs-1
        call para_range(0,npt-1,nprocs,i,jsta,jend)
        jjlen(i)=jend-jsta+1
        idisp(i)=jsta
!        write(*,*) ' @ ',myrank,' jsta,jend=',jsta,jend        
      enddo
      npj = np

      allocate(jstap(0:nprocs))
      allocate(jendp(0:nprocs))
      allocate(npsp(0:nprocs))
      allocate(nprp(0:nprocs))

! ** get start and end id, and number of particles ***
      do j = 0,nprocs-1
! ** define array and number sended to srank **
        jstap(j) = 0
        do i = 0,npj-1
          if(ord(list(i)).lt.idisp(srank(j))) then
            jstap(j) = i+1
          else if(ord(list(i)).gt.idisp(srank(j))+jjlen(srank(j))-1) then
            jendp(j) = i-1
            goto 999
          endif
        enddo
        jendp(j) = npj-1
! *** particle number to send ***
  999   npsp(j) = jendp(j)-jstap(j)+1
        nprp(j) = 0
      enddo
! *** start data transfer ***
! get number of particles -> nprp and np
      np=0
      do j = 0,nprocs-1
        call MPI_ISEND(npsp(j),1,MPI_INTEGER,srank(j) &
         ,1,MPI_COMM_WORLD,ireqs(0),ierr)
        call MPI_IRECV(nprp(j),1,MPI_INTEGER,rrank(j) &
         ,1,MPI_COMM_WORLD,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)
        np=np+nprp(j)
!        write(6,*) ' myrank,j,npsp,nprp=',myrank,j,npsp(j),nprp(j)
      enddo

! *** ord ***

! size of the original particles
      allocate(tivs(0:npj))

      do i=0,npj-1
        tivs(i)=ord(list(i))
      enddo

      deallocate(ord)
! size with the new set of particles
      allocate(ord(0:np))

      np = 0
      do i=jstap(0),jendp(0)
        ord(np)=tivs(i)
        np=np+1
      enddo
      do j = 1,nprocs-1
        call MPI_ISEND(tivs(jstap(j)),npsp(j),MPI_INTEGER,srank(j) &
         ,1,MPI_COMM_WORLD,ireqs(0),ierr)
        call MPI_IRECV(ord(np),nprp(j),MPI_INTEGER,rrank(j) &
         ,1,MPI_COMM_WORLD,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)
        np=np+nprp(j)
      enddo

      deallocate(tivs)

! *** sort ***
      do i=0,np-1
        ord(i)=ord(i)-idisp(srank(0))
      enddo

! *** for integer values ***
#if defined(SF_EFD) || defined(SF_ZFD)
      nval=4
#else
      nval=3
#endif

      allocate(tivr(0:np*nval))

      np=0
      do ip = 0,nprocs-1
        npj=0

        allocate(tivs(0:npsp(ip)*nval))

        do i=jstap(ip),jendp(ip)
          tivs(npj)=flagc_p(list(i))
          tivs(npj+npsp(ip))=flagt_p(list(i))
          tivs(npj+npsp(ip)*2)=id_p(list(i))
#if defined(SF_EFD) || defined(SF_ZFD)
          tivs(npj+npsp(ip)*3)=flagfd_p(list(i))
#endif
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tivs,npsp(ip)*nval,MPI_INTEGER &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tivr(np*nval),nprp(ip)*nval &
           ,MPI_INTEGER,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
        else
          nprp(ip)=npsp(ip)
          do j=0,nval-1         
            do i=0,nprp(ip)-1
              tivr(i+nprp(ip)*j)=tivs(i+npsp(ip)*j)
            enddo
          enddo
        endif

        deallocate(tivs)

        np=np+nprp(ip)
      enddo
      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          flagc_p(ord(npj))=tivr(i)
          flagt_p(ord(npj))=tivr(i+nprp(ip))
          id_p(ord(npj))=tivr(i+nprp(ip)*2)
#if defined(SF_EFD) || defined(SF_ZFD)
          flagfd_p(ord(npj))=tivr(i+nprp(ip)*3)
#endif
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tivr)
  
! *** ?c_p,vn?_p,m_p,u_p(7),h_p,Gsn_p
! *** dv?_p(14),ndu_p,pdu_p,alpv_p,alpu_p ***
      nval=19

      allocate(tdvr(0:np*nval))

      np=0
      do ip = 0,nprocs-1
        npj=0

        allocate(tdvs(0:npsp(ip)*nval))

        do i=jstap(ip),jendp(ip)
          tdvs(npj)=xc_p(list(i))
          tdvs(npj+npsp(ip))=yc_p(list(i))
          tdvs(npj+npsp(ip)*2)=zc_p(list(i))
          tdvs(npj+npsp(ip)*3)=vnx_p(list(i))
          tdvs(npj+npsp(ip)*4)=vny_p(list(i))
          tdvs(npj+npsp(ip)*5)=vnz_p(list(i))
          tdvs(npj+npsp(ip)*6)=m_p(list(i))
          tdvs(npj+npsp(ip)*7)=u_p(list(i))
          tdvs(npj+npsp(ip)*8)=h_p(list(i))
          tdvs(npj+npsp(ip)*9)=dvx_p(list(i))
          tdvs(npj+npsp(ip)*10)=dvy_p(list(i))
          tdvs(npj+npsp(ip)*11)=dvz_p(list(i))
          tdvs(npj+npsp(ip)*12)=ndu_p(list(i))
          tdvs(npj+npsp(ip)*13)=pdu_p(list(i))
          tdvs(npj+npsp(ip)*14)=div_v_p(list(i))
          tdvs(npj+npsp(ip)*15)=lt_p(list(i))
          tdvs(npj+npsp(ip)*16)=dt_p(list(i))
          tdvs(npj+npsp(ip)*17)=alpv_p(list(i))
          tdvs(npj+npsp(ip)*18)=alpu_p(list(i))
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(np*nval),nprp(ip)*nval &
           ,MPI_DOUBLE_PRECISION,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
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

        np=np+nprp(ip)
      enddo
      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          xc_p(ord(npj))=tdvr(i)
          yc_p(ord(npj))=tdvr(i+nprp(ip))
          zc_p(ord(npj))=tdvr(i+nprp(ip)*2)
          vnx_p(ord(npj))=tdvr(i+nprp(ip)*3)
          vny_p(ord(npj))=tdvr(i+nprp(ip)*4)
          vnz_p(ord(npj))=tdvr(i+nprp(ip)*5)
          m_p(ord(npj))=tdvr(i+nprp(ip)*6)
          u_p(ord(npj))=tdvr(i+nprp(ip)*7)
          h_p(ord(npj))=tdvr(i+nprp(ip)*8)
          dvx_p(ord(npj))=tdvr(i+nprp(ip)*9)
          dvy_p(ord(npj))=tdvr(i+nprp(ip)*10)
          dvz_p(ord(npj))=tdvr(i+nprp(ip)*11)
          ndu_p(ord(npj))=tdvr(i+nprp(ip)*12)
          pdu_p(ord(npj))=tdvr(i+nprp(ip)*13)
          div_v_p(ord(npj))=tdvr(i+nprp(ip)*14)
          lt_p(ord(npj))=tdvr(i+nprp(ip)*15)
          dt_p(ord(npj))=tdvr(i+nprp(ip)*16)
          alpv_p(ord(npj))=tdvr(i+nprp(ip)*17)
          alpu_p(ord(npj))=tdvr(i+nprp(ip)*18)
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tdvr)

      if(np.ne.npj) then
        write(6,*) ' Error in ddecb np,npj=',np,npj
        stop
      endif

! *** since pv30.5 individual time step applied for thermal energy update
! *** which requires to update the followings ***
#ifdef COOL
      nval=13
#else
      nval=12
#endif

      allocate(tdvr(0:np*nval))

      np=0
      do ip = 0,nprocs-1

        allocate(tdvs(0:npsp(ip)*nval))

        npj=0
        do i=jstap(ip),jendp(ip)
          tdvs(npj)=rho_p(list(i))
          tdvs(npj+npsp(ip))=omgh_p(list(i))
          tdvs(npj+npsp(ip)*2)=zetah_p(list(i))
          tdvs(npj+npsp(ip)*3)=arot_v_p(list(i))
          tdvs(npj+npsp(ip)*4)=d2u_p(list(i))
          tdvs(npj+npsp(ip)*5)=dtr_p(list(i))
          tdvs(npj+npsp(ip)*6)=myu_p(list(i))
          tdvs(npj+npsp(ip)*7)=dtg_p(list(i))
          tdvs(npj+npsp(ip)*8)=ltg_p(list(i))
          tdvs(npj+npsp(ip)*9)=ax_p(list(i))
          tdvs(npj+npsp(ip)*10)=ay_p(list(i))
          tdvs(npj+npsp(ip)*11)=az_p(list(i))
#ifdef COOL
          tdvs(npj+npsp(ip)*12)=ram_p(list(i))
#endif
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(np*nval),nprp(ip)*nval &
           ,MPI_DOUBLE_PRECISION,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
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

        np=np+nprp(ip)
      enddo
      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          rho_p(ord(npj))=tdvr(i)
          omgh_p(ord(npj))=tdvr(i+nprp(ip))
          zetah_p(ord(npj))=tdvr(i+nprp(ip)*2)
          arot_v_p(ord(npj))=tdvr(i+nprp(ip)*3)
          d2u_p(ord(npj))=tdvr(i+nprp(ip)*4)
          dtr_p(ord(npj))=tdvr(i+nprp(ip)*5)
          myu_p(ord(npj))=tdvr(i+nprp(ip)*6)
          dtg_p(ord(npj))=tdvr(i+nprp(ip)*7)
          ltg_p(ord(npj))=tdvr(i+nprp(ip)*8)
          ax_p(ord(npj))=tdvr(i+nprp(ip)*9)
          ay_p(ord(npj))=tdvr(i+nprp(ip)*10)
          az_p(ord(npj))=tdvr(i+nprp(ip)*11)
#ifdef COOL
          ram_p(ord(npj))=tdvr(i+nprp(ip)*12)
#endif
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tdvr)

! for stars

#ifdef STAR
#ifdef SF_EFD
      nval=2
#else
      nval=1
#endif

      allocate(tdvr(0:np*nval))

      np=0
      do ip = 0,nprocs-1

        allocate(tdvs(0:npsp(ip)*nval))

        npj=0
        do i=jstap(ip),jendp(ip)
          tdvs(npj)=ts_p(list(i))
#ifdef SF_EFD
          tdvs(npj+npsp(ip))=Gsn_p(list(i))
#endif
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(np*nval),nprp(ip)*nval,MPI_DOUBLE_PRECISION &
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

        np=np+nprp(ip)
      enddo
      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          ts_p(ord(npj))=tdvr(i)
#ifdef SF_EFD
          Gsn_p(ord(npj))=tdvr(i+nprp(ip))
#endif
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tdvr)

#endif


#ifdef METAL
#ifdef METALDIFF
      nval=10
#else
      nval=9
#endif

      allocate(tdvr(0:np*nval))

      np=0
      do ip = 0,nprocs-1

        allocate(tdvs(0:npsp(ip)*nval))

        npj=0
        do i=jstap(ip),jendp(ip)
          tdvs(npj)=mzHe_p(list(i))
          tdvs(npj+npsp(ip))=mzC_p(list(i))
          tdvs(npj+npsp(ip)*2)=mzN_p(list(i))
          tdvs(npj+npsp(ip)*3)=mzO_p(list(i))
          tdvs(npj+npsp(ip)*4)=mzNe_p(list(i))
          tdvs(npj+npsp(ip)*5)=mzMg_p(list(i))
          tdvs(npj+npsp(ip)*6)=mzSi_p(list(i))
          tdvs(npj+npsp(ip)*7)=mzFe_p(list(i))
          tdvs(npj+npsp(ip)*8)=mzZ_p(list(i))
#ifdef METALDIFF
          tdvs(npj+npsp(ip)*9)=vsig_p(list(i))
#endif
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(np*nval),nprp(ip)*nval,MPI_DOUBLE_PRECISION &
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

        np=np+nprp(ip)
      enddo
      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          mzHe_p(ord(npj))=tdvr(i)
          mzC_p(ord(npj))=tdvr(i+nprp(ip))
          mzN_p(ord(npj))=tdvr(i+nprp(ip)*2)
          mzO_p(ord(npj))=tdvr(i+nprp(ip)*3)
          mzNe_p(ord(npj))=tdvr(i+nprp(ip)*4)
          mzMg_p(ord(npj))=tdvr(i+nprp(ip)*5)
          mzSi_p(ord(npj))=tdvr(i+nprp(ip)*6)
          mzFe_p(ord(npj))=tdvr(i+nprp(ip)*7)
          mzZ_p(ord(npj))=tdvr(i+nprp(ip)*8)
#ifdef METALDIFF
          vsig_p(ord(npj))=tdvr(i+nprp(ip)*9)
#endif
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tdvr)

#endif

#ifdef METAL
! *** for original metal ***
      nval=9

      allocate(tdvr(0:np*nval))

      np=0
      do ip = 0,nprocs-1

        allocate(tdvs(0:npsp(ip)*nval))

        npj=0
        do i=jstap(ip),jendp(ip)
          tdvs(npj)=mzHe0_p(list(i))
          tdvs(npj+npsp(ip))=mzC0_p(list(i))
          tdvs(npj+npsp(ip)*2)=mzN0_p(list(i))
          tdvs(npj+npsp(ip)*3)=mzO0_p(list(i))
          tdvs(npj+npsp(ip)*4)=mzNe0_p(list(i))
          tdvs(npj+npsp(ip)*5)=mzMg0_p(list(i))
          tdvs(npj+npsp(ip)*6)=mzSi0_p(list(i))
          tdvs(npj+npsp(ip)*7)=mzFe0_p(list(i))
          tdvs(npj+npsp(ip)*8)=mzZ0_p(list(i))
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(np*nval),nprp(ip)*nval,MPI_DOUBLE_PRECISION &
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

        np=np+nprp(ip)
      enddo
      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          mzHe0_p(ord(npj))=tdvr(i)
          mzC0_p(ord(npj))=tdvr(i+nprp(ip))
          mzN0_p(ord(npj))=tdvr(i+nprp(ip)*2)
          mzO0_p(ord(npj))=tdvr(i+nprp(ip)*3)
          mzNe0_p(ord(npj))=tdvr(i+nprp(ip)*4)
          mzMg0_p(ord(npj))=tdvr(i+nprp(ip)*5)
          mzSi0_p(ord(npj))=tdvr(i+nprp(ip)*6)
          mzFe0_p(ord(npj))=tdvr(i+nprp(ip)*7)
          mzZ0_p(ord(npj))=tdvr(i+nprp(ip)*8)
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tdvr)

#endif

      deallocate(srank)
      deallocate(rrank)
      deallocate(jstap)
      deallocate(jendp)
      deallocate(idisp)
      deallocate(jjlen)
      deallocate(npsp)
      deallocate(nprp)

      deallocate(istatus)
      deallocate(istatusr)

! ** after care **

      deallocate(ord)
      deallocate(list)
      allocate(list(0:np))

! ** make active particle list **
      ns=0
      ng=0
! ** devide gas and star **
      do i=0,np-1
        if(flagc_p(i).le.0) then
          list_ap(ng)=i
          ng=ng+1
        else
          list(ns)=i
          ns=ns+1
        endif
      enddo
      do i=0,ns-1
        list_ap(ng+i)=list(i)
      enddo

      deallocate(list)

      if(step.ge.0) then
! *** find active particles ***
        lwt =TM_lt+MGTL*TM_dt
        upt=TM_lt+MGTU*TM_dt

        allocate(talist(0:np))
        allocate(list(0:np))
        allocate(nalist(0:np))

        do i=0,np-1
          talist(i)=list_ap(i)
        enddo
! ** for Gas Particles **
        na=0
        nanp=0
        nagravg=0
        do i=0,ng-1
          pn=talist(i)		 
          if(lwt.lt.ltg_p(pn)+dtg_p(pn).and.ltg_p(pn)+dtg_p(pn).lt.upt) then
            list_ap(nagravg)=pn
            nagravg=nagravg+1
          else if(lwt.lt.lt_p(pn)+dt_p(pn).and.lt_p(pn)+dt_p(pn).lt.upt) then
            list(na)=pn
            na=na+1
          else
            nalist(nanp)=pn
            nanp=nanp+1
          endif
        enddo
! ** list_ap for non-active particles **
        do i=0,na-1
          list_ap(nagravg+i)=list(i)
        enddo
        nag=na+nagravg
        do i=0,nanp-1
          list_ap(nag+i)=nalist(i)
        enddo

! ** for Star Particles **
        if(ns.gt.0) then
          na=0
          nanp=0
          do i=ng,np-1
            pn=talist(i)		   
            if(lwt.lt.lt_p(pn)+dt_p(pn).and. lt_p(pn)+dt_p(pn).lt.upt) then
              list_ap(ng+na)=pn
              na=na+1
            else
              nalist(nanp)=pn
              nanp=nanp+1
            endif
          enddo
! ** list_ap for non-active particles **
          j=ng+na
          nas = na
          do i=0,nanp-1
            list_ap(j)=nalist(i)
            j=j+1
          enddo
        else
          nas=0
        endif

        deallocate(talist)
        deallocate(list)
        deallocate(nalist)

! ** setup virtual timestep **
! ** for Gas and Star Particles **
        do i=0,np-1
#ifdef FIXEDP
        if(id_p(i).lt.SI_fid(0).or.id_p(i).gt.SI_fid(1)) then
#endif
          vdt_p(i)=TM_lt+TM_dt-lt_p(i)
#ifdef FIXEDP
        else
          vdt_p(i)=0.0d0
        endif
#endif
        enddo
! *** for gas properties ***
! *** since pv30.5 need this ***
        do i=0,ng-1
          pn=list_ap(i)
          p_p(pn)=(gam-1.0d0)*rho_p(pn)*u_p(pn)
          cs_p(pn)=dsqrt(gam*p_p(pn)/rho_p(pn))
          as_p(pn)=p_p(pn)/(rho_p(pn)**GAM)
        enddo
      else
        nas=ns
        nag=ng
        nagravg=ng
        do i=0,ng-1
          pn=list_ap(i)
          x_p(pn)=xc_p(pn)
          y_p(pn)=yc_p(pn)
          z_p(pn)=zc_p(pn)
          vvnx_p(pn)=vnx_p(pn)
          vx_p(pn)=vnx_p(pn)
          vvny_p(pn)=vny_p(pn)
          vy_p(pn)=vny_p(pn)
          vvnz_p(pn)=vnz_p(pn)
          vz_p(pn) = vnz_p(pn)
! need to set p for PINI
          p_p(pn)=(gam-1.0d0)*rho_p(pn)*u_p(pn)
          cs_p(pn)=dsqrt(gam*p_p(pn)/rho_p(pn))
          as_p(pn)=p_p(pn)/(rho_p(pn)**GAM)
        enddo
        do i=ng,ng+ns-1
          pn=list_ap(i)
          x_p(pn)=xc_p(pn)
          y_p(pn)=yc_p(pn)
          z_p(pn)=zc_p(pn)
          vvnx_p(pn)=vnx_p(pn)
          vx_p(pn)=vnx_p(pn)
          vvny_p(pn)=vny_p(pn)
          vy_p(pn)=vny_p(pn)
          vvnz_p(pn)=vnz_p(pn)
          vz_p(pn) = vnz_p(pn)
        enddo

#ifdef FIXEDP
! *** change nag and nagravg

        allocate(nalist(0:ng))
        allocate(list(0:ng))

        nag=0
        nc=0
        do i=0,ng-1
          pn=list_ap(i)
          if(id_p(pn).ge.SI_fid(0).and.id_p(pn).le.SI_fid(1)) then
            nalist(nc)=pn
            nc=nc+1
          else
            list(nag)=pn
            nag=nag+1
          endif
        enddo
        nagravg=nag
! *** change list_ap putting fixedp at the end
        do i=0,nag-1
          list_ap(i)=list(i)
        enddo
        do i=0,nc-1
          list_ap(nag+i)=nalist(i)
        enddo

        deallocate(nalist)
        deallocate(list)

#endif
      endif

#ifdef FIXEDP
      do i=0,ng-1
        pn=list_ap(i)
        if(id_p(pn).ge.SI_fid(0).and.id_p(pn).le.SI_fid(1)) then
          x_p(pn)=xc_p(pn)
          y_p(pn)=yc_p(pn)
          z_p(pn)=zc_p(pn)
          vvnx_p(pn)=vnx_p(pn)
          vx_p(pn)=vnx_p(pn)
          vvny_p(pn)=vny_p(pn)
          vy_p(pn)=vny_p(pn)
          vvnz_p(pn)=vnz_p(pn)
          vz_p(pn) = vnz_p(pn)
! *** set dvx and etc.
          dvx_p(pn)=0.0d0
          dvy_p(pn)=0.0d0
          dvz_p(pn)=0.0d0
          ax_p(pn)=0.0d0
          ay_p(pn)=0.0d0
          az_p(pn)=0.0d0
          ndu_p(pn)=0.0d0
          pdu_p(pn)=0.0d0
        endif
      enddo
#endif

#if defined(SF_EFD) || defined(SF_ZFD)
! *** set total yield ***

      do i=0,np-1
! ** total **
        tnsn_p(i)=0.0d0
        tmej_p(i)=0.0d0
        tmzHe_p(i)=0.0d0
        tmzC_p(i)=0.0d0
        tmzN_p(i)=0.0d0
        tmzO_p(i)=0.0d0
        tmzNe_p(i)=0.0d0
        tmzMg_p(i)=0.0d0
        tmzSi_p(i)=0.0d0
        tmzFe_p(i)=0.0d0
        tmzZ_p(i)=0.0d0
      enddo

      allocate(tmej_sn(0:np))
      allocate(tnsn_sn(0:np))
      allocate(tmzHe_sn(0:np))
      allocate(tmzC_sn(0:np))
      allocate(tmzN_sn(0:np))
      allocate(tmzO_sn(0:np))
      allocate(tmzNe_sn(0:np))
      allocate(tmzMg_sn(0:np))
      allocate(tmzSi_sn(0:np))
      allocate(tmzFe_sn(0:np))
      allocate(tmzZ_sn(0:np))
      allocate(nsw_sn(0:np))

! *** use TM_dt, allow TMsf_dt is different from TM_dt ***
      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)

      do i=0,ng+ns-1
        pn = list_ap(i)
        if(flagrfd_p(pn).eq.1.and.i.lt.ng) then
          flagrfd_p(pn)=0
        endif
!  update tm*_p() 
!  Unit MUSM 
        tnsn_p(pn)=tnsn_sn(pn)
        tmej_p(pn)=tmej_sn(pn)
!  Unit Solar Mass 
        tmzHe_p(pn)=tmzHe_sn(pn)
        tmzC_p(pn)=tmzC_sn(pn)
        tmzN_p(pn)=tmzN_sn(pn)
        tmzO_p(pn)=tmzO_sn(pn)
        tmzNe_p(pn)=tmzNe_sn(pn)
        tmzMg_p(pn)=tmzMg_sn(pn)
        tmzSi_p(pn)=tmzSi_sn(pn)
        tmzFe_p(pn)=tmzFe_sn(pn)
        tmzZ_p(pn)=tmzZ_sn(pn)
      enddo

      deallocate(tmej_sn)
      deallocate(tnsn_sn)
      deallocate(tmzHe_sn)
      deallocate(tmzC_sn)
      deallocate(tmzN_sn)
      deallocate(tmzO_sn)
      deallocate(tmzNe_sn)
      deallocate(tmzMg_sn)
      deallocate(tmzSi_sn)
      deallocate(tmzFe_sn)
      deallocate(tmzZ_sn)
      deallocate(nsw_sn)

      if(flagcont.eq.1.or.step.le.0) then
        TM_dt=INF
      endif
#endif

!      if(step.gt.0) then
!        write(fileo,'(a4,i3.3,i3.3)') 'decb',myrank,step
!      else
!        write(fileo,'(a4,i3.3,i3.3)') 'decb',myrank,0
!      endif
!       write(fileo,'(a4,i3.3)') 'decb',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np-1
!        write(60,'(13(1pE13.5),2I10)') xc_p(i),yc_p(i),zc_p(i) &
!        ,vnx_p(i),vny_p(i),vnz_p(i),dvx_p(i),dvy_p(i),dvz_p(i) &
!        ,rho_p(i),vdt_p(i),dt_p(i),lt_p(i),id_p(i),i
!      enddo
!      close(60)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(6,*) 'after decb ng,nag,myrank=',ng,nag,ns,nas,nagravg,myrank
!      stop


end subroutine ddecb
#endif

