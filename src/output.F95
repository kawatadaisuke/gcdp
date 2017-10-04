#include "gcdp.def"
! ***********************************************
!    output.f for gcd+ ver. f03.0
!  27  Jan. 2013    written by D.KAWATA
! ***********************************************

!*****   Definition of analysis()   *****
subroutine output(npt,ndmt,np,ng,ndm,ns,step,totstep)
      use gcdp_const
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#ifdef METAL
      use gcdp_metal
#endif
#endif
#ifdef DM
      use gcdp_dm
#endif
#ifdef MHD  
      use gcdp_mhd
#endif
      use gcdp_system
    
      implicit none
      include 'mpif.h'

      integer,intent(in) :: npt,ndmt,np,ng,ndm,ns,step,totstep
      integer nag,nadm,nas,nagravg
      integer flagascii
      integer i,pn
      integer ngt,nst
      double precision mg,ms,mgt,mst,mmet,tnsn,mfdg,mfdgt,maxrhog,mdm
      integer nfdg,nfdgt
#ifdef TESTP
      integer idtps,pntps
#endif
! * for file *
      character sname*6,filename*60
      integer logs,slogs
      double precision fs
      integer uni,count
! *** for communication ***
      integer ipsta,ipend,ipd
      integer srank,npj
      integer ip,iv,nivalb,nivaldm,ndvaldm,nisupb,nisupdm
      integer ndbval,ndbhyd,ndbsf,ndbmet
! * for system value *
      double precision gkmm
      double precision iax,z
      double precision rhobm,dtu,dtdv
#ifdef MHD
      integer ndbbat,ndbmhd
      double precision valfip,pbp
#endif
#ifdef OUT_TEST
      integer nval
      double precision ene(0:1),etot(0:1),momt(0:2)
#endif
#ifdef TESTP
      double precision rp,rmp,eccp
#endif
! *** for work ***
      integer ierr,npjrmax
      integer,allocatable :: istatus(:)
      integer,allocatable :: npjr(:),ireqs(:),ireqr(:)
      integer,allocatable :: tivs(:),tivr(:)
      double precision,allocatable :: tdvs(:),tdvr(:)

      gkmm=(GAM-1.0d0)/((KCGS/(MP*MYU))/K_MU)
      iax = 1.0d0/SI_a
      z = iax-1.0d0

#ifdef OUT_TEST
! *** binary output
      flagascii = 1
! *** ascii output
!      flagascii = 0
#else
! *** binary output: no other choice ***
      flagascii = 1
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
      nivalb=4
#else
      nivalb=3
#endif
! *** number of output values ***
#ifdef GAS
      ndbval=9
      ndbhyd=8
#else
      ndbval=0
      ndbhyd=0
#endif
#ifdef STAR
      ndbsf=1
#else
      ndbsf=0
#endif
#ifdef METAL
#ifdef METALDIFF
      ndbmet=19
#else
      ndbmet=18
#endif
#else
      ndbmet=0
#endif
#ifdef GAS
      nisupb=1
#else
      nisupb=0
#endif
#ifdef DM
      nisupdm=1
#endif

#if defined(COSM) && defined(BOXSIM)
! convert vnx_ to physical values
#if defined(GAS) || defined(STAR)
      do i=0,np-1
        vnx_p(i)=SI_hub*xc_p(i)+vnx_p(i)/SI_a
        vny_p(i)=SI_hub*yc_p(i)+vny_p(i)/SI_a
        vnz_p(i)=SI_hub*zc_p(i)+vnz_p(i)/SI_a
      enddo
#endif
#ifdef DM
      do i=0,ndm-1
        vnx_dm(i)=SI_hub*xc_dm(i)+vnx_dm(i)/SI_a
        vny_dm(i)=SI_hub*yc_dm(i)+vny_dm(i)/SI_a
        vnz_dm(i)=SI_hub*zc_dm(i)+vnz_dm(i)/SI_a
      enddo
#endif
#endif

! *** for DM ***
      nivaldm=2
#ifdef TESTP
      ndvaldm=13
#else
      ndvaldm=10
#endif
! #ifdef OUT_TEST
!      ndvaldm=12
! #endif

      nag=ng
      nadm=ndm
      nas=ns
      nagravg=ng
      nfdg=0
      if(myrank.lt.SI_nof) then
        ipsta=myrank+SI_nof
        ipend=nprocs-1
        ipd=SI_nof
! dummy
        srank=myrank
      else
        ipsta=myrank
        ipend=myrank
        ipd=1
        srank=mod(myrank,SI_nof)
      endif
! * compute mass of gas and stars *
! *** for check of mass conservation ***
      mmet=0.0d0
      mg = 0.0d0 
      mfdg=0.0d0
      maxrhog=0.0d0
#ifdef GAS     
      do i=0,ng-1
        pn=list_ap(i)
        if(rho_p(pn).gt.maxrhog) then  
          maxrhog=rho_p(pn)
        endif
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).lt.0) then
          nfdg=nfdg+1
          mfdg=mfdg+m_p(pn)
        else
#endif
          mg=mg+m_p(pn)
#if defined(SF_EFD) || defined(SF_ZFD)  
        endif
#endif

#ifdef METAL
        mmet=mmet+mzZ_p(pn)
#endif
      enddo
#endif
      ms = 0.0d0      
#ifdef STAR
      do i=ng,ng+ns-1
        pn=list_ap(i)
        ms = ms+m_p(pn)
#ifdef METAL
        mmet=mmet+mzZ_p(pn)
#endif
      enddo
#endif
      mdm=0.0d0
#ifdef DM
      do i=0,ndm-1
        mdm=mdm+m_dm(i)
      enddo
#endif
! allocate work space
      allocate(tdvs(0:4))
      allocate(tdvr(0:4))
! *** get maxrhog

      tdvs(0)=maxrhog
      tdvr(0)=0.0d0
      call MPI_ALLREDUCE(tdvs,tdvr,1,MPI_DOUBLE_PRECISION &
       ,MPI_MAX,MPI_COMM_WORLD,ierr)
      maxrhog=tdvr(0)
! *** get total mass
      tdvs(0)=mg
      tdvs(1)=ms
      tdvs(2)=mmet
      tdvs(3)=mfdg
      tdvs(4)=mdm
      tdvr(0)=0.0d0
      tdvr(1)=0.0d0
      tdvr(2)=0.0d0
      tdvr(3)=0.0d0
      tdvr(4)=0.0d0
      call MPI_ALLREDUCE(tdvs,tdvr,5,MPI_DOUBLE_PRECISION &
       ,MPI_SUM,MPI_COMM_WORLD,ierr)
      mgt=tdvr(0)
      mst=tdvr(1)
      mmet=tdvr(2)
      mfdgt=tdvr(3)
      mdm=tdvr(4)
! *** for total number of star particles ***
      allocate(tivs(0:2))
      allocate(tivr(0:2))
      tivs(0)=ng
      tivs(1)=ns
      tivs(2)=nfdg
      tivr(0)=0
      tivr(1)=0
      tivr(2)=0
      tivr(3)=0
      call MPI_ALLREDUCE(tivs,tivr,3,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      ngt=tivr(0)
      nst=tivr(1)-tivr(2)
      nfdgt=tivr(2)

      deallocate(tivs)
      deallocate(tivr)

#ifdef OUT_TEST
      ene(0)=0.0d0
      ene(1)=0.0d0
      momt(0)=0.0d0
      momt(1)=0.0d0
      momt(2)=0.0d0
#ifdef GAS
      do i=0,ng-1
        pn=list_ap(i)
        ene(0)=ene(0)+0.5d0*m_p(pn)*(vnx_p(pn)**2+vny_p(pn)**2+vnz_p(pn)**2)
        ene(1)=ene(1)+m_p(pn)*u_p(pn)
        momt(0)=momt(0)+m_p(pn)*vnx_p(pn)
        momt(1)=momt(1)+m_p(pn)*vny_p(pn)
        momt(2)=momt(2)+m_p(pn)*vnz_p(pn)
      enddo
#endif
#ifdef DM
      do i=0,ndm-1
        ene(0)=ene(0)+0.5d0*m_dm(i)*(vnx_dm(i)**2+vny_dm(i)**2+vnz_dm(i)**2)
        momt(0)=momt(0)+m_dm(i)*vnx_dm(i)
        momt(1)=momt(1)+m_dm(i)*vny_dm(i)
        momt(2)=momt(2)+m_dm(i)*vnz_dm(i)
      enddo
#endif
      if(nprocs.gt.1) then
        tdvs(0)=ene(0)
        tdvs(1)=ene(1)
        tdvs(2)=momt(0)
        tdvs(3)=momt(1)
        tdvs(4)=momt(2)
        nval=5
        do i=0,nval-1  
          tdvr(0)=0.0d0       
        enddo
        call MPI_ALLREDUCE(tdvs,tdvr,nval,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
        etot(0)=tdvr(0)
        etot(1)=tdvr(1)
        momt(0)=tdvr(2)
        momt(1)=tdvr(3)
        momt(2)=tdvr(4)
      else
        etot(0)=ene(0)
        etot(1)=ene(1)
      endif
#endif
      deallocate(tdvs)
      deallocate(tdvr)

! open system.dat
      if(myrank.eq.0) then
        if(step.eq.0) then
          open(60,file='./output/ana/system.dat',status='unknown')
          write(60,*) &
           't tstep step ns ng nfdg ns flagrad a mg mfdg ms mb mZ mZdiff'
        else
          open(60,file='./output/ana/system.dat',status='old')
          do i=1,step+1
            read(60,*)
          enddo
        endif
! ***  Calculate System Value ***
! *** Total ***
! *** estimate number of supernovae ***
#ifdef STAR
        if(SNEU.gt.0.0d0) then
          tnsn=SI_esns/SNEU
        else
          tnsn=0.0d0
        endif
#endif

        write(60,169) TM_tot,totstep,step,ngt,nfdgt,nst &
         ,SI_flagrad,SI_a,mgt,mfdgt,mst,mgt+mfdgt+mst,mmet,SI_mms &
! 13: mmet, 14: SI_mms
#ifdef OUT_TEST
#ifdef STAR
! 15, 16, 17-25
         ,SI_esns,tnsn,SI_zsns(0),SI_zsns(1),SI_zsns(2),SI_zsns(3) &
         ,SI_zsns(4),SI_zsns(5),SI_zsns(6),SI_zsns(7),SI_zsns(8) &
! 26-35
         ,SI_tnsf,0,0,TM_tot*TMUGYR &
         ,etot(0)+etot(1),etot(0),etot(1),momt(0),momt(1),momt(2) &
         ,maxrhog,TM_dt
! *** assuming NYEL=9
  169   format((1pE13.5),6I10,7(1pE13.5),11(1pE13.5),3I10,9(1pE13.5))
#else
         ,etot(0)+etot(1),etot(0),etot(1),momt(0),momt(1),momt(2),TM_dt
  169   format((1pE13.5),6I10,14(1pE13.5))
#endif
#else
#ifdef STAR
         ,SI_esns,tnsn,SI_zsns(0),SI_zsns(1),SI_zsns(2),SI_zsns(3) &
         ,SI_zsns(4),SI_zsns(5),SI_zsns(6),SI_zsns(7),SI_zsns(8) &
         ,SI_tnsf,TM_tot*TMUGYR,maxrhog,TM_dt
  169   format((1pE13.5),6I10,7(1pE13.5),11(1pE13.5),I10,3(1pE13.5))
#else
         ,TM_dt
  169   format((1pE13.5),6I10,8(1pE13.5))
#endif
#endif
        close(60)
#ifdef DEL_DMPART
! output mass and particle number history
        open(60,file='./output/ana/del_dmpart.dat',status='old')
        do i=0,step+4-1
          read(60,*)
        enddo
        write(60,170) TM_tot,mdm,SI_ndmtdc,ndmt
  170   format(2(1pE13.5),2I10)
        close(60)
#endif
#ifdef ADD_DMPART
        open(60,file='./output/ana/add_dmpart.dat',status='old')
        do i=0,step+6-1
          read(60,172)

        enddo
        write(60,172) TM_tot,mdm,SI_ndmtadded,ndmt
  172   format(2(1pE13.5),2I10)
#endif

      endif

#ifdef TESTP

! update rperi, rapo and zmax
      do i=ndm-SI_ntpdm,ndm-1
        rp=dsqrt(x_dm(i)**2+y_dm(i)**2)
        if(rp.lt.rperi_dm(i)) then
          rperi_dm(i)=rp
        endif
        if(rp.gt.rapo_dm(i)) then
          rapo_dm(i)=rp 
        endif
        if(dabs(z_dm(i)).gt.zmax_dm(i)) then
          zmax_dm(i)=dabs(z_dm(i))
        endif
      enddo

! starting (unique) ID of test particle
      idtps=ndmt-SI_nttpdm
      do i=ndm-SI_ntpdm,ndm-1
        if(id_dm(i).ge.idtps.and.id_dm(i).lt.idtps+OF_nporb) then
          write(filename,'(a16,i6.6)') './output/data/id',id_dm(i)-idtps
          if(step.eq.0) then
            open(70,file=filename,status='unknown')
          else
            open(70,file=filename,status='old',access='append')
          endif
          rp=dsqrt(x_dm(i)**2+y_dm(i)**2)
          rmp=0.5d0*(rperi_dm(i)+rapo_dm(i))
          eccp=(rapo_dm(i)-rperi_dm(i))/(rperi_dm(i)+rapo_dm(i))
          write(70,'(10(1pE13.5))') x_dm(i),y_dm(i),z_dm(i) &
           ,rperi_dm(i),rapo_dm(i),zmax_dm(i),rmp,eccp,rp,TM_tot
          close(70)
        endif
      enddo
!      write(6,*) ' id output range=',idtps,idtps+OF_nporb
!      stop

#endif

      if(step.eq.0.or.SI_flagout.ne.0) then
        if(myrank.eq.0) then
          if(step.eq.0) then
            open(70,file='./output/ana/ostep.dat',status='unknown')
          else
            open(70,file='./output/ana/ostep.dat',status='old',access='append')
          endif
          write(70,'(I8)') step
          close(70)
        endif

#ifdef TESTP

        allocate(istatus(MPI_STATUS_SIZE))

! *** for binary output ***
! allocate array for MPI
        allocate(npjr(0:nprocs))
        allocate(ireqs(0:nprocs))
        allocate(ireqr(0:nprocs))

        do i=0,nprocs-1
          npjr(i)=0.0d0
        enddo
        ip=myrank

! starting (unique) ID of test particle
        idtps=ndmt-SI_nttpdm
! starting position in array (pn) in each proc
        pntps=ndm-SI_ntpdm
! *** output DM data ***
        do i=0,nprocs-1
          npjr(i)=0
        enddo
        ip=myrank
        if(myrank.lt.SI_nof) then
          write(filename,'(a20,i6.6,a1,i4.4)') &
           './output/data/bdvals',step,'n',myrank
          open(62,file=filename,status='unknown' &
            ,form='unformatted')
          write(62) npt,SI_nttpdm,SI_nttpdm,SI_a,TM_tot
          write(62) nprocs,SI_nof,nivaldm,ndvaldm
        endif
          
        allocate(tivs(0:6))
        allocate(tivr(0:6))

! output only test particles
        tivs(0)=ng
        tivs(1)=SI_ntpdm
        tivs(2)=ns
        tivs(3)=SI_ntpdm
        tivs(4)=nag
        tivs(5)=SI_ntpdm
        tivs(6)=nas
        npj=tivs(1)
        if(myrank.lt.SI_nof) then
          write(62) (tivs(i),i=0,6)
          npjr(myrank)=tivs(1)
          do ip=ipsta,ipend,ipd
            call MPI_IRECV(tivr,7,MPI_INTEGER,ip,1,MPI_COMM_WORLD &
              ,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
            write(62) (tivr(i),i=0,6)
            npjr(ip)=tivr(1)
          enddo
        else
          call MPI_ISEND(tivs,7,MPI_INTEGER,srank,1 &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif

! *** output integer values ***
        npjrmax=0
        do ip=0,nprocs-1
          if(npjrmax.lt.npjr(ip)) then
            npjrmax=npjr(ip)
          endif
        enddo

        deallocate(tivs)
        deallocate(tivr)
        allocate(tivs(0:npj-1))
        allocate(tivr(0:npjrmax-1))

        do iv=0,nivaldm-1
          if(iv.eq.0) then
            do i=0,npj-1
! ID=0 starts from test particles
              tivs(i)=id_dm(i+pntps)-idtps
            enddo
          else if(iv.eq.1) then
            do i=0,npj-1
              tivs(i)=list_adm(i)-pntps
            enddo
          endif
          if(myrank.lt.SI_nof) then
            write(62) (tivs(i),i=0,npj-1)
            do ip=ipsta,ipend,ipd
              call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              write(62) (tivr(i),i=0,npjr(ip)-1)
            enddo
          else
            call MPI_ISEND(tivs,npj,MPI_INTEGER,srank,iv &
              ,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        enddo

! *** output double values ***
        allocate(tdvs(0:npj-1))
        allocate(tdvr(0:npjrmax-1))

        do iv=0,ndvaldm-1
          if(iv.eq.0) then
            do i=0,npj-1
              tdvs(i)=xc_dm(i+pntps)
            enddo
          else if(iv.eq.1) then
            do i=0,npj-1 
              tdvs(i)=yc_dm(i+pntps)
            enddo
          else if(iv.eq.2) then
            do i=0,npj-1 
              tdvs(i)=zc_dm(i+pntps)
            enddo
          else if(iv.eq.3) then
            do i=0,npj-1 
              tdvs(i)=vnx_dm(i+pntps)
            enddo
          else if(iv.eq.4) then
            do i=0,npj-1 
              tdvs(i)=vny_dm(i+pntps)
            enddo
          else if(iv.eq.5) then
            do i=0,npj-1 
              tdvs(i)=vnz_dm(i+pntps)
            enddo
          else if(iv.eq.6) then
            do i=0,npj-1 
              tdvs(i)=m_dm(i+pntps)
            enddo
          else if(iv.eq.7) then
            do i=0,npj-1 
              tdvs(i)=rho_dm(i+pntps)
            enddo
          else if(iv.eq.8) then
            do i=0,npj-1 
              tdvs(i)=h_dm(i+pntps)
            enddo
          else if(iv.eq.9) then
            do i=0,npj-1 
              tdvs(i)=tadd_dm(i+pntps)
            enddo
          else if(iv.eq.10) then
            do i=0,npj-1 
              tdvs(i)=rperi_dm(i+pntps)
            enddo
          else if(iv.eq.11) then
            do i=0,npj-1 
              tdvs(i)=rapo_dm(i+pntps)
            enddo
          else if(iv.eq.12) then
            do i=0,npj-1 
              tdvs(i)=zmax_dm(i+pntps)
            enddo
          endif
          if(myrank.lt.SI_nof) then
            write(62) (tdvs(i),i=0,npj-1)
            do ip=ipsta,ipend,ipd
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
               ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              write(62) (tdvr(i),i=0,npjr(ip)-1)
            enddo
          else
            call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION,srank &
              ,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        enddo
        if(myrank.lt.SI_nof) then     
          close(62)
        endif

       deallocate(npjr)
       deallocate(ireqs)
       deallocate(ireqr)
       deallocate(istatus)

#else

! *** Open Output File ***
        if(flagascii.eq.0) then
#ifdef DM
          if(ndm.gt.0) then
            write(filename,'(a16,i6.6,a1,i4.4)') &
             './output/data/ds',step,'n',myrank
            open(61,file=filename,status='unknown')
            do i=0,ndm-1
              write(61,161) xc_dm(i),yc_dm(i),zc_dm(i) &
               ,vnx_dm(i),vny_dm(i),vnz_dm(i),m_dm(i),nnb_dm(i) &
               ,id_dm(i),rho_dm(i),h_dm(i),omgh_dm(i),-zetah_dm(i),div_v_dm(i) &
! 15
               ,dsqrt(xc_dm(i)**2+yc_dm(i)**2+zc_dm(i)**2) &
! 16
               ,dabs(dvx_dm(i)*xc_dm(i)+dvy_dm(i)*yc_dm(i) &
               +dvz_dm(i)*zc_dm(i))/dsqrt(xc_dm(i)**2+yc_dm(i)**2+zc_dm(i)**2) &
! 17-20
               ,dvx_dm(i),dvy_dm(i),dvz_dm(i),dt_dm(i),tadd_dm(i)
  161         format(7(1pE13.5),2I10,12(1pE13.5))
            enddo
            close(61)
          endif
#endif
#ifdef GAS
          if(ng.gt.0.or.ns.gt.0) then

#ifdef MHD
            write(filename,'(a17,i6.6,a1,i4.4)') &
           './output/data/mhd',step,'n',myrank
            open(68,file=filename,status='unknown')
            do i=0,ng-1
              pn=list_ap(i)
              valfip=(bx_p(pn)+by_p(pn)+bz_p(pn))/dsqrt(rho_p(pn))
              pbp=0.5d0*(bx_p(pn)**2+by_p(pn)**2+bz_p(pn)**2)
              write(68,168) xc_p(pn),yc_p(pn),zc_p(pn),vnx_p(pn) &
                   ,vny_p(pn),vnz_p(pn),m_p(pn),rho_p(pn),u_p(pn) &
                   ,p_p(pn),bx_p(pn),by_p(pn),bz_p(pn),divb_p(pn) &
                   ,arotb_p(pn),agradb_p(pn),alpb_p(pn),valfip,pbp &
#ifdef DCLEAN
                   ,psi_p(pn),maxvmhd_p(pn) &
#endif
                   ,divberr_p(pn) !&
!                   ,dvx_p(pn),dvy_p(pn),dvz_p(pn),dbx_p(pn),dby_p(pn) &
!                   ,dbz_p(pn),pdbx_p(pn),pdby_p(pn),pdbz_p(pn)
#ifdef DCLEAN
  168          format(22(1pE13.5))
#else
  168          format(20(1pE13.5))
#endif
!  168          format(29(1pE13.5))
            enddo
            close(68)

#ifdef BBAT
            write(filename,'(a18,i6.6,a1,i4.4)') &
             './output/data/bbat',step,'n',myrank
            open(71,file=filename,status='unknown')
            do i=0,ng-1
              pn=list_ap(i)
              write(71,171) chi_p(pn),lognh_p(pn),logne_p(pn)
  171         format(3(1pE13.5))
            enddo
            close(71)
#endif
#endif

            write(filename,'(a16,i6.6,a1,i4.4)') &
             './output/data/gs',step,'n',myrank
            open(61,file=filename,status='unknown')
#ifdef COSM
            rhobm=SI_omgb*3.0d0*HUB0*SI_h0*HUB0*SI_h0/(8.0d0*M_PI*G)/(SI_a**3)
#endif
            do i=0,ng-1
              pn=list_ap(i)
! *** output dt requirements ***
              dtdv=dsqrt(dvx_p(pn)**2+dvy_p(pn)**2+dvz_p(pn)**2)
              if(dtdv.gt.0.0d0) then
                dtdv=CDYN*dsqrt(h_p(pn)*0.5d0/dtdv)
              else
                dtdv=0.0d0
              endif
#ifdef SF_EFD
            if(flagfd_p(pn).lt.0.and.Gsn_p(pn).gt.0.0d0) then
              dtu=LDUFD*u_p(pn)/(Gsn_p(pn)/TMsf_dt)
            else        
#endif
              dtu=0.0d0
#ifdef SF_EFD
            endif
#endif
              write(61,162) xc_p(pn),yc_p(pn),zc_p(pn) &
               ,vnx_p(pn),vny_p(pn),vnz_p(pn),m_p(pn),rho_p(pn),u_p(pn) &
! 10-15
#if defined(SF_EFD) || defined(SF_ZFD)
               ,p_p(pn),h_p(pn),nnb_p(pn),flagfd_p(pn),id_p(pn),div_v_p(pn) &
#else
               ,p_p(pn),h_p(pn),nnb_p(pn),pn,id_p(pn),div_v_p(pn) &
#endif
! 16-20
#ifdef METALDIFF
! *** diffusion coefficient
              ,2.0d0*vsig_p(pn)*VUKMS*1.0e5*h_p(pn)*LUCM,vsig_p(pn) &
#else
              ,arot_v_p(pn),-zetah_p(pn) &
#endif
              ,CCFL*hvsigdt_p(pn),dtdv,dtu &
! 21-24
              ,cs_p(pn),d2u_p(pn),alpv_p(pn),alpu_p(pn) &
! 25-26
              ,alpu_p(pn)*V_ALPHAMIN/alpv_p(pn),omgh_p(pn) &
! 27
              ,dsqrt(xc_p(pn)**2+yc_p(pn)**2+zc_p(pn)**2) &
! 28 T
              ,(GAM-1.0d0)*u_p(pn)*myu_p(pn)/(TPRHO*MYU)*TUK &
! 29
              ,dabs(ax_p(pn)*xc_p(pn)+ay_p(pn)*yc_p(pn) &
               +az_p(pn)*zc_p(pn))/dsqrt(xc_p(pn)**2+yc_p(pn)**2+zc_p(pn)**2) &
! 30-32
              ,dvx_p(pn),dvy_p(pn),dvz_p(pn) &
!     &   ,dsqrt(ax_p(pn)**2+ay_p(pn)**2+az_p(pn)**2)
! 33 (19 since i)
              ,ndu_p(pn) &
! 34
#ifdef METAL
! 34:nh, 35, 36
              ,((m_p(pn)-((mzZ_p(pn)+mzHe_p(pn))/MUSM))/m_p(pn)) &
               *rho_p(pn)*(DU/MP),mzZ_p(pn)/(m_p(pn)*MUSM) &
#if defined(METALDIFF) && defined(CSZDIFF)
! *** for test 
              ,XZSOL*dexp(-(xc_p(pn)**2+yc_p(pn)**2+zc_p(pn)**2) &
               /(4.0d0*cs_p(pn)*h_p(pn)*TM_tot)) &
               /((4.0d0*M_PI*cs_p(pn)*h_p(pn)*TM_tot)**(3.0d0/2.0d0)) &
#else
              ,mzZ0_p(pn)/(m_p(pn)*MUSM) &
#endif
#endif
#ifdef COOL
! 37,38
              ,ram_p(pn),dabs(u_p(pn)/ram_p(pn)) &
#endif
#ifdef COSM
! 39
              ,rho_p(pn)/rhobm-1.0d0 &
#endif
! end of the line
              ,0
! format
#ifdef METAL
#ifdef COOL
#ifdef COSM
  162         format(11(1pE13.5),3I10,25(1pE13.5),I2)
#else
! no COSM
  162         format(11(1pE13.5),3I10,24(1pE13.5),I2)
#endif
#else
! no COOL with COSM
#ifdef COSM
  162         format(11(1pE13.5),3I10,24(1pE13.5),I2)
#else
! no COSM no COOL
  162         format(11(1pE13.5),3I10,23(1pE13.5),I2)
#endif
#endif

#else
! no metal
#ifdef COOL
#ifdef COSM
! no metal + cool + cosm
  162         format(11(1pE13.5),3I10,22(1pE13.5),I2)
#else
! no metal no cosm + cool
  162         format(11(1pE13.5),3I10,21(1pE13.5),I2)
#endif
#else
! no metal no COOL + cosm
#ifdef COSM
  162         format(11(1pE13.5),3I10,20(1pE13.5),I2)
#else
! no metal, no cosm, no cool
  162         format(11(1pE13.5),3I10,19(1pE13.5),I2)
#endif
#endif
#endif
            enddo
            close(61)
#ifdef STAR
            if(ns.gt.0) then
              write(filename,'(a16,i6.6,a1,i4.4)') &
               './output/data/ss',step,'n',myrank
              open(61,file=filename,status='unknown')
              do i=ng,ng+ns-1
                pn=list_ap(i) 
              write(61,163) xc_p(pn),yc_p(pn),zc_p(pn) &
              ,vnx_p(pn),vny_p(pn),vnz_p(pn),m_p(pn),rho_p(pn) &
! 9-15
              ,u_p(pn),ts_p(pn),h_p(pn),div_v_p(pn),nnb_p(pn),pn,id_p(pn) &
! 16-19
              ,dt_p(pn),dvx_p(pn),dvy_p(pn),dvz_p(pn) &
#if defined(SF_EFD) || defined(SF_ZFD)
              ,omgh_p(pn),zetah_p(pn) &
#ifdef METAL
              ,mzZ_p(pn)/(m_p(pn)*MUSM) &
#else
              ,0.0d0 &
#endif
             ,dsqrt(xc_p(pn)**2+yc_p(pn)**2+zc_p(pn)**2),flagfd_p(pn)
  163        format(12(1pE13.5),3I10,8(1pE13.5),I8)
#else
             ,dsqrt(xc_p(pn)**2+yc_p(pn)**2+zc_p(pn)**2)
  163        format(12(1pE13.5),3I10,5(1pE13.5))
#endif
              enddo 
              close(61)
              close(62)
            endif
#endif
          endif
#endif
        else

          allocate(istatus(MPI_STATUS_SIZE))

! *** for binary output ***
! allocate array for MPI
          allocate(npjr(0:nprocs))
          allocate(ireqs(0:nprocs))
          allocate(ireqr(0:nprocs))
          do i=0,nprocs-1
            npjr(i)=0.0d0
          enddo
          ip=myrank
#ifdef GAS
          if(npt.gt.0) then
            if(myrank.lt.SI_nof) then
! *** file for baryon basic values ***
              write(filename,'(a20,i6.6,a1,i4.4)') &
               './output/data/bbvals',step,'n',myrank
              open(61,file=filename,status='unknown',form='unformatted')
              write(61) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(61) nprocs,SI_nof,nivalb,ndbval
            endif

            allocate(tivs(0:7))
            allocate(tivr(0:7))

            tivs(0)=ng
            tivs(1)=ndm
            tivs(2)=ns
            tivs(3)=SI_ndm1
            tivs(4)=nag
            tivs(5)=nadm
            tivs(6)=nas
            tivs(7)=nagravg
! *** number of sending particles ****
            npj=tivs(0)+tivs(2)
            if(myrank.lt.SI_nof) then
              write(61) (tivs(i),i=0,7)
              npjr(myrank)=tivs(0)+tivs(2)
              do ip=ipsta,ipend,ipd
                call MPI_IRECV(tivr,8,MPI_INTEGER,ip,1 &
                 ,MPI_COMM_WORLD,ireqr(ip),ierr)
                call MPI_WAIT(ireqr(ip),istatus,ierr)
                write(61) (tivr(i),i=0,7)
! *** number of recieving particles ****
                npjr(ip)=tivr(0)+tivr(2)
              enddo
            else
              call MPI_ISEND(tivs,8,MPI_INTEGER,srank,1 &
               ,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif

            npjrmax=0
            do ip=0,nprocs-1
              if(npjrmax.lt.npjr(ip)) then
                npjrmax=npjr(ip)
              endif
            enddo      
            deallocate(tivs)
            deallocate(tivr)
            allocate(tivs(0:npj))  
            allocate(tivr(0:npjrmax))

! *** output integer values ***
            do iv=0,nivalb-1
              if(iv.eq.0) then
                do i=0,npj-1
                  tivs(i)=id_p(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npj-1
                  tivs(i)=flagc_p(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npj-1
                 tivs(i)=list_ap(i)
                enddo
#if defined(SF_EFD) || defined(SF_ZFD)
              else if(iv.eq.3) then                  
                do i=0,npj-1
                  tivs(i)=flagfd_p(i)
                enddo
#endif
              endif
              if(myrank.lt.SI_nof) then
                write(61) (tivs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                   ,MPI_COMM_WORLD,ireqr(ip),ierr)
                   call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(61) (tivr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tivs,npj,MPI_INTEGER,srank,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo

! *** output double values ***
            allocate(tdvs(0:npj))  
            allocate(tdvr(0:npjrmax))
! *** output baryon values ***
            do iv=0,ndbval-1
! *** set data to send ***
              if(iv.eq.0) then
                do i=0,npj-1 
                  tdvs(i)=xc_p(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npj-1 
                  tdvs(i)=yc_p(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npj-1 
                  tdvs(i)=zc_p(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npj-1 
                  tdvs(i)=vnx_p(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npj-1 
                  tdvs(i)=vny_p(i)
                enddo
              else if(iv.eq.5) then
                do i=0,npj-1 
                  tdvs(i)=vnz_p(i)
                enddo
              else if(iv.eq.6) then
                do i=0,npj-1 
                  tdvs(i)=m_p(i)
                enddo
              else if(iv.eq.7) then
                do i=0,npj-1 
                  tdvs(i)=rho_p(i)
                enddo
              else if(iv.eq.8) then
                do i=0,npj-1 
                  tdvs(i)=u_p(i)
                enddo
              endif
              if(myrank.lt.SI_nof) then
                write(61) (tdvs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                   ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(61) (tdvr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
                 ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
              close(61)
            endif

! *** file for baryon hydro values ***
            if(myrank.lt.SI_nof) then
              write(filename,'(a20,i6.6,a1,i4.4)') &
               './output/data/bbhyds',step,'n',myrank
              open(63,file=filename,status='unknown',form='unformatted')
              write(63) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(63) nprocs,SI_nof,ndbhyd
            endif
! *** output hydro vaules ***
            do iv=0,ndbhyd-1
              if(iv.eq.0) then
                do i=0,npj-1 
                  tdvs(i)=h_p(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npj-1 
                  tdvs(i)=div_v_p(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npj-1 
                  tdvs(i)=alpv_p(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npj-1 
                  tdvs(i)=alpu_p(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npj-1 
                  tdvs(i)=myu_p(i)
                enddo
              else if(iv.eq.5) then
                do i=0,npj-1 
                  tdvs(i)=dt_p(i)
                enddo
              else if(iv.eq.6) then
                do i=0,npj-1 
                  tdvs(i)=hvsigdt_p(i)
                enddo
              else if(iv.eq.7) then
                do i=0,npj-1 
                  tdvs(i)=arot_v_p(i)
                enddo
              endif
              if(myrank.lt.SI_nof) then
                write(63) (tdvs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                   ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(63) (tdvr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
                 ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
              close(63)
            endif
! *** file for baryon H,He cooling values ***
#ifdef STAR
            if(myrank.lt.SI_nof) then
              write(filename,'(a20,i6.6,a1,i4.4)') &
               './output/data/bbsfis',step,'n',myrank
              open(65,file=filename,status='unknown',form='unformatted')
              write(65) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(65) nprocs,SI_nof,ndbsf
            endif
! *** output star values ***
            do iv=0,ndbsf-1
              if(iv.eq.0) then
                do i=0,npj-1 
                  tdvs(i)=ts_p(i)
                enddo
              endif
              if(myrank.lt.SI_nof) then
                write(65) (tdvs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                   ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(65) (tdvr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
                 ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
              close(65)
            endif
#endif
! *** file for baryon metal values ***
#ifdef METAL
            if(myrank.lt.SI_nof) then
              write(filename,'(a20,i6.6,a1,i4.4)') &
               './output/data/bbmets',step,'n',myrank
              open(67,file=filename,status='unknown',form='unformatted')
              write(67) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(67) nprocs,SI_nof,ndbmet
            endif
            do iv=0,ndbmet-1
              if(iv.eq.0) then
                do i=0,npj-1 
                  tdvs(i)=mzHe_p(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npj-1 
                  tdvs(i)=mzZ_p(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npj-1 
                  tdvs(i)=mzC_p(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npj-1 
                  tdvs(i)=mzN_p(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npj-1 
                  tdvs(i)=mzO_p(i)
                enddo
              else if(iv.eq.5) then
                do i=0,npj-1 
                  tdvs(i)=mzNe_p(i)
                enddo
              else if(iv.eq.6) then
                do i=0,npj-1 
                  tdvs(i)=mzMg_p(i)
                enddo
              else if(iv.eq.7) then
                do i=0,npj-1 
                  tdvs(i)=mzSi_p(i)
                enddo
              else if(iv.eq.8) then
                do i=0,npj-1 
                  tdvs(i)=mzFe_p(i)
                enddo
              else if(iv.eq.9) then
                do i=0,npj-1 
                  tdvs(i)=mzHe0_p(i)
                enddo
              else if(iv.eq.10) then
                do i=0,npj-1 
                  tdvs(i)=mzZ0_p(i)
                enddo
              else if(iv.eq.11) then
                do i=0,npj-1 
                  tdvs(i)=mzC0_p(i)
                enddo
              else if(iv.eq.12) then
                do i=0,npj-1 
                  tdvs(i)=mzN0_p(i)
                enddo
              else if(iv.eq.13) then
                do i=0,npj-1 
                  tdvs(i)=mzO0_p(i)
                enddo
              else if(iv.eq.14) then
                do i=0,npj-1 
                  tdvs(i)=mzNe0_p(i)
                enddo
              else if(iv.eq.15) then
                do i=0,npj-1 
                  tdvs(i)=mzMg0_p(i)
                enddo
              else if(iv.eq.16) then
                do i=0,npj-1 
                  tdvs(i)=mzSi0_p(i)
                enddo
              else if(iv.eq.17) then
                do i=0,npj-1 
                  tdvs(i)=mzFe0_p(i)
                enddo
#ifdef METALDIFF
              else if(iv.eq.18) then
                do i=0,npj-1 
                  tdvs(i)=vsig_p(i)
                enddo
#endif
              endif
              if(myrank.lt.SI_nof) then
                write(67) (tdvs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                   ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(67) (tdvr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
                 ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
              close(67)
            endif
#endif
! *** file for supplementary information ***
            if(myrank.lt.SI_nof) then
              write(filename,'(a20,i6.6,a1,i4.4)') &
               './output/data/bbsups',step,'n',myrank
              open(68,file=filename,status='unknown',form='unformatted')
              write(68) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(68) nprocs,SI_nof,nisupb
            endif
! *** for supplementary information ***
            do iv=0,nisupb-1
              if(iv.eq.0) then
                do i=0,npj-1 
                  tivs(i)=nnb_p(i)
                enddo
              endif
              if(myrank.lt.SI_nof) then
                write(68) (tivs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                   ,MPI_COMM_WORLD,ireqr(ip),ierr)
                   call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(68) (tivr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tivs,npj,MPI_INTEGER,srank,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
              close(68)
            endif

            deallocate(tivs)
            deallocate(tivr)

#ifdef MHD
            ndbmhd=10
#ifdef DCLEAN
            ndbmhd=ndbmhd+2
#endif            
            if (myrank.lt.SI_nof) then
              write(filename, '(a20,i6.6,a1,i4.4)') &
                './output/data/bbmhds',step,'n',myrank
              open(67,file=filename,status='unknown',form='unformatted')
              write(67) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(67) nprocs,SI_nof,ndbmhd
            endif
            do iv=0,ndbmhd-1
              if (iv.eq.0)  then
                  do i=0,npj-1
                     tdvs(i)=bx_p(i)
                  enddo
              else if (iv.eq.1) then
                  do i=0,npj-1
                     tdvs(i)=by_p(i)
                  enddo
              else if (iv.eq.2) then
                  do i=0,npj-1
                     tdvs(i)=bz_p(i)
                  enddo
              else if (iv.eq.3) then
                  do i=0,npj-1
                     tdvs(i)=divb_p(i)
                  enddo
              else if (iv.eq.4) then
                  do i=0,npj-1
                     tdvs(i)=arotb_p(i)
                  enddo
              else if (iv.eq.5) then
                  do i=0,npj-1
                     tdvs(i)=agradb_p(i)
                  enddo
              else if (iv.eq.6) then
                  do i=0,npj-1
                     tdvs(i)=alpb_p(i)
                  enddo
              else if (iv.eq.7) then
                  do i=0,npj-1
                    valfip=(bx_p(i)+by_p(i)+bz_p(i))/dsqrt(rho_p(i))
                    tdvs(i)=valfip
                  enddo
              else if (iv.eq.8) then
                  do i=0,npj-1
                    pbp=0.5d0*(bx_p(i)**2+by_p(i)**2+bz_p(i)**2)
                    tdvs(i)=pbp
                  enddo
              else if (iv.eq.9) then
                  do i=0,npj-1
                     tdvs(i)=divberr_p(i)
                  enddo
#ifdef DCLEAN
              else if (iv.eq.10) then
                  do i=0,npj-1
                     tdvs(i)=psi_p(i)
                  enddo
              else if (iv.eq.11) then
                  do i=0,npj-1
                     tdvs(i)=maxvmhd_p(i)
                  enddo
#endif
              endif
              if (myrank.lt.SI_nof) then
                write(67) (tdvs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                    ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(67) (tdvr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
                  ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
              close(67)
            endif
#ifdef BBAT
            ndbbat=3
            if (myrank.lt.SI_nof) then
              write(filename, '(a20,i6.6,a1,i4.4)') &
               './output/data/bbbats',step,'n',myrank
              open(72,file=filename,status='unknown',form='unformatted')
              write(72) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(72) nprocs,SI_nof,ndbmhd
            endif
            do iv=0,ndbbat-1
              if (iv.eq.0)  then
                  do i=0,npj-1
                     tdvs(i)=bx_p(i)
                  enddo
              else if (iv.eq.1) then
                  do i=0,npj-1
                     tdvs(i)=by_p(i)
                  enddo
              else if (iv.eq.2) then
                  do i=0,npj-1
                     tdvs(i)=bz_p(i)
                  enddo
              endif
              if (myrank.lt.SI_nof) then
                write(72) (tdvs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                   call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                    ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                   call MPI_WAIT(ireqr(ip),istatus,ierr)
                   write(72) (tdvr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
                 ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                 call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
               close(72)
            endif
#endif
#endif

            deallocate(tdvs)
            deallocate(tdvr)

! endif for if(npt.gt.0)
          endif



! end for GAS
#endif
! end for GAS

#ifdef DM
! *** output DM data ***
          if(ndmt.gt.0) then
            do i=0,nprocs-1
              npjr(i)=0
            enddo
            ip=myrank
            if(myrank.lt.SI_nof) then
              write(filename,'(a20,i6.6,a1,i4.4)') &
               './output/data/bdvals',step,'n',myrank
              open(62,file=filename,status='unknown' &
               ,form='unformatted')
              write(62) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(62) nprocs,SI_nof,nivaldm,ndvaldm
            endif
          
            allocate(tivs(0:6))
            allocate(tivr(0:6))

            tivs(0)=ng
            tivs(1)=ndm
            tivs(2)=ns
            tivs(3)=SI_ndm1
            tivs(4)=nag
            tivs(5)=nadm
            tivs(6)=nas
            npj=tivs(1)
            if(myrank.lt.SI_nof) then
              write(62) (tivs(i),i=0,6)
              npjr(myrank)=tivs(1)
              do ip=ipsta,ipend,ipd
                call MPI_IRECV(tivr,7,MPI_INTEGER,ip,1,MPI_COMM_WORLD &
                 ,ireqr(ip),ierr)
                call MPI_WAIT(ireqr(ip),istatus,ierr)
                write(62) (tivr(i),i=0,6)
                npjr(ip)=tivr(1)
              enddo
            else
              call MPI_ISEND(tivs,7,MPI_INTEGER,srank,1 &
               ,MPI_COMM_WORLD,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif

! *** output integer values ***

            npjrmax=0
            do ip=0,nprocs-1
              if(npjrmax.lt.npjr(ip)) then
                npjrmax=npjr(ip)
              endif
            enddo
            deallocate(tivs)
            deallocate(tivr)
            allocate(tivs(0:npj))
            allocate(tivr(0:npjrmax))

            do iv=0,nivaldm-1
              if(iv.eq.0) then
                do i=0,npj-1
                  tivs(i)=id_dm(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npj-1
                  tivs(i)=list_adm(i)
                enddo
              endif
              if(myrank.lt.SI_nof) then
                write(62) (tivs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                   ,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(62) (tivr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tivs,npj,MPI_INTEGER,srank,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo

! *** output double values ***
            allocate(tdvs(0:npj))
            allocate(tdvr(0:npjrmax))

            do iv=0,ndvaldm-1
              if(iv.eq.0) then
                do i=0,npj-1 
                  tdvs(i)=xc_dm(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npj-1 
                  tdvs(i)=yc_dm(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npj-1 
                  tdvs(i)=zc_dm(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npj-1 
                  tdvs(i)=vnx_dm(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npj-1 
                  tdvs(i)=vny_dm(i)
                enddo
              else if(iv.eq.5) then
                do i=0,npj-1 
                  tdvs(i)=vnz_dm(i)
                enddo
              else if(iv.eq.6) then
                do i=0,npj-1 
                  tdvs(i)=m_dm(i)
                enddo
              else if(iv.eq.7) then
                do i=0,npj-1 
                  tdvs(i)=rho_dm(i)
                enddo
              else if(iv.eq.8) then
                do i=0,npj-1 
                  tdvs(i)=h_dm(i)
                enddo
              else if(iv.eq.9) then
                do i=0,npj-1 
                  tdvs(i)=tadd_dm(i)
                enddo
              endif
              if(myrank.lt.SI_nof) then
                write(62) (tdvs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                   ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(62) (tdvr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION,srank &
                 ,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then     
              close(62)
            endif

! *** file for supplementary information ***
            if(myrank.lt.SI_nof) then     
              write(filename,'(a20,i6.6,a1,i4.4)') &
               './output/data/bdsups',step,'n',myrank
              open(69,file=filename,status='unknown',form='unformatted')
              write(69) npt,ndmt,SI_ndm1t,SI_a,TM_tot
              write(69) nprocs,SI_nof,nisupdm
            endif
! *** output supplementary integer values ***
            do iv=0,nisupdm-1
              if(iv.eq.0) then
                do i=0,npj-1
                  tivs(i)=nnb_dm(i)
                enddo
              endif
              if(myrank.lt.SI_nof) then
                write(69) (tivs(i),i=0,npj-1)
                do ip=ipsta,ipend,ipd
                  call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER &
                   ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                  call MPI_WAIT(ireqr(ip),istatus,ierr)
                  write(69) (tivr(i),i=0,npjr(ip)-1)
                enddo
              else
                call MPI_ISEND(tivs,npj,MPI_INTEGER,srank,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            enddo
            if(myrank.lt.SI_nof) then
              close(69)
            endif

            deallocate(tivs)
            deallocate(tivr)
            deallocate(tdvs)
            deallocate(tdvr)

          endif


! *** endif for DM ***
#endif

          deallocate(npjr)
          deallocate(ireqs)
          deallocate(ireqr)
          deallocate(istatus)

! *** end for if(flagascii.eq.0) then
        endif

! endif for TESTP
#endif

        if(step.ne.0) then
#ifdef OFTIMESPEC
          OF_iof=OF_iof+1
          if(OF_iof.ge.OF_nof.and.TM_tot.lt.TM_t0-0.0000001d0) then
            if(myrank.eq.0) then
              write(6,*) ' Error in output(): OF_iof>=OF_nof'
              write(6,*) ' Ttot,OF_nof,iof=',TM_tot,OF_nof,OF_iof
            endif
          else if(OF_iof.lt.OF_nof) then
            OF_todt=OF_time(OF_iof)
            if(OF_todt.lt.TM_tot) then
              if(myrank.eq.0) then
                write(6,*) ' Error in output(): OF_todt<TM_tot'
                write(6,*) ' Ttot,OF_nof,iof,todt=',TM_tot,OF_nof,OF_iof,OF_todt
              endif
              call MPI_FINALIZE()
              stop
            endif
          endif
#else       
!          if(z.lt.1.0d0) then
!            OF_todt=OF_todt+4.0d0*OF_odt
!          else
            OF_todt=OF_todt+OF_odt
!          endif
#endif
          SI_flagout = 0
        endif


      endif

#if defined(COSM) && defined(BOXSIM)
! convert vnx_ back to a^2 vpec unit
#if defined(GAS) || defined(STAR)
      do i=0,np-1
        vnx_p(i)=(vnx_p(i)-SI_hub*xc_p(i))*SI_a
        vny_p(i)=(vny_p(i)-SI_hub*yc_p(i))*SI_a
        vnz_p(i)=(vnz_p(i)-SI_hub*zc_p(i))*SI_a
      enddo
#endif
#ifdef DM
      do i=0,ndm-1
        vnx_dm(i)=(vnx_dm(i)-SI_hub*xc_dm(i))*SI_a
        vny_dm(i)=(vny_dm(i)-SI_hub*yc_dm(i))*SI_a
        vnz_dm(i)=(vnz_dm(i)-SI_hub*zc_dm(i))*SI_a
      enddo
#endif
#endif


end subroutine

