#include "gcdp.def"
! ***********************************************
!   forcont.F95 for GCD+ ver. f03.0
!   output data for continue calculation
!  3  Feb. 2013   written by D.KAWATA
! ***********************************************

subroutine forcont(npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg,step,totstep,lstep)
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
   
      integer,intent(in) :: npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg,step &
       ,totstep,lstep 
      integer i,ip,iv,nivalb,nivaldm,ndvaldm,pn
      integer ipsta,ipend,ipd
      integer npj,iwdummy
      integer ndbval,ndbhyd,ndbsf,ndbmet,ndbcon,nddmcon,nibcon,nidmcon
#ifdef MHD
      integer ndbmhd
#endif
      integer srank
      character fileo*60
! *** for work ***
      integer ierr,npjrmax
      integer,allocatable :: npjr(:),ireqs(:),ireqr(:),istatus(:)
      integer,allocatable :: tivs(:),tivr(:)
      double precision,allocatable :: tdvs(:),tdvr(:)

! allocate array for MPI
      allocate(npjr(0:nprocs))
      allocate(ireqs(0:nprocs))
      allocate(ireqr(0:nprocs))

      if(myrank.eq.0) then
        open(60,file='./output/data/cont.dat',status='unknown' &
         ,form='unformatted')
         iwdummy=0
! /*** Output continue data ***/
#ifdef OFTIMESPEC
         write(60) step,totstep,lstep,SI_flagout,OF_todt &
#else
         write(60) step,totstep,lstep,SI_flagout,OF_odt,OF_todt &
#endif
#ifdef BOXSIM
           ,SI_lbox(0),SI_lbox(1),SI_lbox(2) &
#ifdef COSM
           ,SI_lbox0(0),SI_lbox0(1),SI_lbox0(2) &
#endif
#endif
         ,TM_tot,TM_sdt,TM_lt,TM_ldtmin,TM_k,TM_kn,SI_a,SI_eps,SI_mms &
#ifdef COSM
         ,SI_omg0,SI_lam0,SI_h0,SI_eps0,SI_epsa &
#endif
#ifdef COOL
         ,SI_flagrad &
#endif
#ifdef STAR
        ,SI_esns,(SI_zsns(i),i=0,NYEL-1),TMsf_t,TMsf_dt,SI_tnsf &
#endif
        ,iwdummy
#ifdef OFTIMESPEC
        write(60) OF_nof,OF_iof
        write(60) OF_time
#endif
        close(60)
      endif

#if defined(SF_EFD) || defined(SF_ZFD)
      nivalb=4
#else
      nivalb=3
#endif
      nivaldm=2
#ifdef TESTP_HYBRID
      ndvaldm=13
#else
      ndvaldm=10
#endif
#if defined(COSM) && defined(BOXSIM)
      nddmcon=14
#else
      nddmcon=12
#endif
      ndbcon=25
      nibcon=1
      nidmcon=1
#if defined(COSM) && defined(BOXSIM)
      ndbcon=ndbcon+4
#endif

! *** number of output values ***
#ifdef GAS
      ndbval=9
      ndbhyd=5
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

      if(myrank.lt.SI_nof) then
        ipsta=myrank+SI_nof
        ipend=nprocs-1
        ipd=SI_nof
        srank=myrank
      else
        ipsta=myrank
        ipend=myrank
        ipd=1
        srank=mod(myrank,SI_nof)
      endif
      do i=0,nprocs-1
        npjr(i)=0
      enddo
      ip=myrank

      allocate(istatus(MPI_STATUS_SIZE))

#ifdef GAS
      if(npt.gt.0) then

        allocate(tivs(0:7))
        allocate(tivr(0:7))

        if(myrank.lt.SI_nof) then
          write(fileo,'(a17,i4.4,a4)') 'output/data/pbval',myrank,'.dat'
          open(61,file=fileo,status='unknown',form='unformatted')
          write(61) npt,ndmt,SI_ndm1t,SI_a,TM_tot
          write(61) nprocs,SI_nof,nivalb,ndbval
        endif
        tivs(0)=ng
        tivs(1)=ndm
        tivs(2)=ns
        tivs(3)=SI_ndm1
        tivs(4)=nag
        tivs(5)=nadm
        tivs(6)=nas
        tivs(7)=nagravg
! *** number of sending particles ***
        npj=tivs(0)+tivs(2)        
        if(myrank.lt.SI_nof) then
          write(61) (tivs(i),i=0,7)
          do ip=ipsta,ipend,ipd
            call MPI_IRECV(tivr,8,MPI_INTEGER,ip,1,MPI_COMM_WORLD &
             ,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
            write(61) (tivr(i),i=0,7)
            npjr(ip)=tivr(0)+tivr(2)
          enddo
        else
          call MPI_ISEND(tivs,8,MPI_INTEGER,srank,1,MPI_COMM_WORLD &
           ,ireqs(ip),ierr)
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

! *** output integer variables ***
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
! *** output double variables ***

        allocate(tdvs(0:npj))  
        allocate(tdvr(0:npjrmax))

! *** output baryon values ***
        do iv=0,ndbval-1
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
! *** hydro values
        if(myrank.lt.SI_nof) then
          write(fileo,'(a17,i4.4,a4)') 'output/data/pbhyd',myrank,'.dat'
          open(63,file=fileo,status='unknown',form='unformatted')
          write(63) npt,ndmt,SI_ndm1t,SI_a,TM_tot
          write(63) nprocs,SI_nof,ndbhyd
        endif
! *** for hydro values ***
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

#ifdef MHD
! *** file for MHD values ***
      ndbmhd=10
#ifdef DCLEAN
      ndbmhd=ndbmhd+1
#endif
      if(myrank.lt.SI_nof) then
        write(fileo,'(a17,i4.4,a4)') 'output/data/pbmhd',myrank,'.dat'
        open(68,file=fileo,status='unknown',form='unformatted')
        write(68) npt,ndmt,SI_ndm1t,SI_a,TM_tot
        write(68) nprocs,SI_nof,ndbmhd
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
            tdvs(i)=valfi_p(i)
          enddo
        else if (iv.eq.8) then
          do i=0,npj-1
            tdvs(i)=pb_p(i)
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
#endif
        endif
        if(myrank.lt.SI_nof) then
          write(68) (tdvs(i),i=0,npj-1)
          do ip=ipsta,ipend,ipd
            call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
              ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
            write(68) (tdvr(i),i=0,npjr(ip)-1)
          enddo
        else
          call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
           ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif
      enddo
      if(myrank.lt.SI_nof) then
        close(68)
      endif
#endif

! *** file for baryon H,He cooling values ***
#ifdef STAR
        if(myrank.lt.SI_nof) then
          write(fileo,'(a17,i4.4,a4)') 'output/data/pbsfi',myrank,'.dat'
          open(65,file=fileo,status='unknown',form='unformatted')
          write(65) npt,ndmt,SI_ndm1t,SI_a,TM_tot
          write(65) nprocs,SI_nof,ndbsf
        endif
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
          write(fileo,'(a17,i4.4,a4)') 'output/data/pbmet',myrank,'.dat'
          open(67,file=fileo,status='unknown',form='unformatted')
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

        if(myrank.lt.SI_nof) then        
          write(fileo,'(a17,i4.4,a4)') 'output/data/pbcon',myrank,'.dat'
          open(71,file=fileo,status='unknown',form='unformatted')
          write(71) npt,ndmt,SI_ndm1t,SI_a,TM_tot
          write(71) nprocs,SI_nof,nibcon,ndbcon
        endif
! *** output cont info ****
! *** output integer variables ***
        do iv=0,nibcon-1
          if(iv.eq.0) then
            do i=0,npj-1
              tivs(i)=flagt_p(i)
            enddo
          endif
          if(myrank.lt.SI_nof) then
            write(71) (tivs(i),i=0,npj-1)
            do ip=ipsta,ipend,ipd
              call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,ip,iv &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              write(71) (tivr(i),i=0,npjr(ip)-1)
            enddo
          else
            call MPI_ISEND(tivs,npj,MPI_INTEGER,srank,iv &
             ,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        enddo
! *** double precision ***
        do iv=0,ndbcon-1
          if(iv.eq.0) then
            do i=0,npj-1 
              tdvs(i)=lt_p(i)
            enddo
          else if(iv.eq.1) then
            do i=0,npj-1 
              tdvs(i)=dt_p(i)
            enddo
          else if(iv.eq.2) then
            do i=0,npj-1 
              tdvs(i)=x_p(i)
            enddo
          else if(iv.eq.3) then
            do i=0,npj-1 
              tdvs(i)=y_p(i)
            enddo
          else if(iv.eq.4) then
            do i=0,npj-1 
              tdvs(i)=z_p(i)
            enddo
          else if(iv.eq.5) then
            do i=0,npj-1 
              tdvs(i)=vvnx_p(i)
            enddo
          else if(iv.eq.6) then
            do i=0,npj-1 
              tdvs(i)=vvny_p(i)
            enddo
          else if(iv.eq.7) then
            do i=0,npj-1 
              tdvs(i)=vvnz_p(i)
            enddo
          else if(iv.eq.8) then
            do i=0,npj-1 
              tdvs(i)=dvx_p(i)
            enddo
          else if(iv.eq.9) then
            do i=0,npj-1 
              tdvs(i)=dvy_p(i)
            enddo
          else if(iv.eq.10) then
            do i=0,npj-1 
              tdvs(i)=dvz_p(i)
            enddo
          else if(iv.eq.11) then
            do i=0,npj-1 
              tdvs(i)=zetah_p(i)
            enddo
          else if(iv.eq.12) then
            do i=0,npj-1 
              tdvs(i)=omgh_p(i)
            enddo
          else if(iv.eq.13) then
            do i=0,npj-1 
              tdvs(i)=arot_v_p(i)
            enddo
          else if(iv.eq.14) then
            do i=0,npj-1 
              tdvs(i)=d2u_p(i)
            enddo
          else if(iv.eq.15) then
            do i=0,npj-1 
              tdvs(i)=dtr_p(i)
            enddo
          else if(iv.eq.16) then
            do i=0,npj-1 
              tdvs(i)=ndu_p(i)
            enddo
          else if(iv.eq.17) then
            do i=0,npj-1 
              tdvs(i)=pdu_p(i)
            enddo
          else if(iv.eq.18) then
            do i=0,npj-1 
              tdvs(i)=dtg_p(i)
            enddo
          else if(iv.eq.19) then
            do i=0,npj-1 
              tdvs(i)=ltg_p(i)
            enddo
          else if(iv.eq.20) then
            do i=0,npj-1 
              tdvs(i)=ax_p(i)
            enddo
          else if(iv.eq.21) then
            do i=0,npj-1 
              tdvs(i)=ay_p(i)
            enddo
          else if(iv.eq.22) then
            do i=0,npj-1 
              tdvs(i)=az_p(i)
            enddo
          else if(iv.eq.23) then
            do i=0,npj-1 
#ifdef COOL
              tdvs(i)=ram_p(i)
#else
! dummy
              tdvs(i)=0.0d0
#endif
            enddo
          else if(iv.eq.24) then

            do i=0,npj-1 
#ifdef SF_EFD
              tdvs(i)=Gsn_p(i)
#else
! dummy
              tdvs(i)=0.0d0
#endif
            enddo
#if defined(COSM) && defined(BOXSIM)
          else if(iv.eq.25) then
            do i=0,npj-1 
              tdvs(i)=asc_p(i)
            enddo
          else if(iv.eq.26) then
            do i=0,npj-1 
              tdvs(i)=tc_p(i)
            enddo
          else if(iv.eq.27) then
            do i=0,npj-1 
              tdvs(i)=ascg_p(i)
            enddo
          else if(iv.eq.28) then
            do i=0,npj-1 
              tdvs(i)=tcg_p(i)
            enddo
#endif
          endif
          if(myrank.lt.SI_nof) then
            write(71) (tdvs(i),i=0,npj-1)
            do ip=ipsta,ipend,ipd
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
               ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              write(71) (tdvr(i),i=0,npjr(ip)-1)
            enddo
          else
            call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION &
             ,srank,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        enddo
        if(myrank.lt.SI_nof) then
          close(71)
        endif

        deallocate(tivs)
        deallocate(tivr)
        deallocate(tdvs)
        deallocate(tdvr)

! *** endif for if(npt.gt.0) then
      endif

! *** endif for if(npt.gt.0) then
! *** end GAS ***
#endif
! *** end GAS ***
#ifdef DM
! *** output DM data ***
      if(ndmt.gt.0) then
        do i=0,nprocs-1
          npjr(i)=0
        enddo
        ip=myrank
        if(myrank.lt.SI_nof) then
          write(fileo,'(a17,i4.4,a4)') 'output/data/pdval',myrank,'.dat'
          open(62,file=fileo,status='unknown',form='unformatted')
          write(62) npt,ndmt,SI_ndm1t,SI_a,TM_tot          
          write(62) nprocs,SI_nof,nivaldm,ndvaldm
        endif

        allocate(tivs(0:6))
        allocate(tivr(0:6))

        tivs(0)=ng
        tivs(1)=ndm
        tivs(2)=ns
#ifdef TESTP_HYBRID        
        tivs(3)=SI_nttpdm
#else
        tivs(3)=SI_ndm1
#endif        
        tivs(4)=nag
        tivs(5)=nadm
        tivs(6)=nas

! *** number of particles sending ****
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
! *** output integer variables ***

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

! *** output double variables ***
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
#ifdef TESTP_HYBRID               
          else if(iv.eq.10) then
            do i=0,npj-1 
              tdvs(i)=rperi_dm(i)
            enddo
          else if(iv.eq.11) then
            do i=0,npj-1 
              tdvs(i)=rapo_dm(i)
            enddo
          else if(iv.eq.12) then
            do i=0,npj-1 
              tdvs(i)=zmax_dm(i)
            enddo
#endif
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

! *** for continue
        if(myrank.lt.SI_nof) then     
          write(fileo,'(a17,i4.4,a4)') 'output/data/pdcon',myrank,'.dat'
          open(64,file=fileo,status='unknown',form='unformatted')
          write(64) npt,ndmt,SI_ndm1t,SI_a,TM_tot
          write(64) nprocs,SI_nof,nidmcon,nddmcon
        endif
! *** for continue 
! *** output integer variables ***
        do iv=0,nidmcon-1
          if(iv.eq.0) then
            do i=0,npj-1
              tivs(i)=flagt_dm(i)
            enddo
          endif
          if(myrank.lt.SI_nof) then
            write(64) (tivs(i),i=0,npj-1)
            do ip=ipsta,ipend,ipd
              call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER &
               ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              write(64) (tivr(i),i=0,npjr(ip)-1)
            enddo
          else
            call MPI_ISEND(tivs,npj,MPI_INTEGER,srank,iv &
             ,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        enddo
! *** double precision ***
        do iv=0,nddmcon-1
          if(iv.eq.0) then
            do i=0,npj-1 
              tdvs(i)=lt_dm(i)
            enddo
          else if(iv.eq.1) then
            do i=0,npj-1 
              tdvs(i)=dt_dm(i)
            enddo
          else if(iv.eq.2) then
            do i=0,npj-1 
              tdvs(i)=x_dm(i)
            enddo
          else if(iv.eq.3) then
            do i=0,npj-1 
              tdvs(i)=y_dm(i)
            enddo
          else if(iv.eq.4) then
            do i=0,npj-1 
              tdvs(i)=z_dm(i)
            enddo
          else if(iv.eq.5) then
            do i=0,npj-1 
              tdvs(i)=dvx_dm(i)
            enddo
          else if(iv.eq.6) then
            do i=0,npj-1 
              tdvs(i)=dvy_dm(i)
            enddo
          else if(iv.eq.7) then
            do i=0,npj-1 
              tdvs(i)=dvz_dm(i)
            enddo
          else if(iv.eq.8) then
            do i=0,npj-1 
              tdvs(i)=zetah_dm(i)
            enddo
          else if(iv.eq.9) then
            do i=0,npj-1 
              tdvs(i)=omgh_dm(i)
            enddo
          else if(iv.eq.10) then
            do i=0,npj-1 
              tdvs(i)=div_v_dm(i)
            enddo
          else if(iv.eq.11) then
            do i=0,npj-1 
              tdvs(i)=dtr_dm(i)
            enddo
#if defined(COSM) && defined(BOXSIM)
          else if(iv.eq.12) then
            do i=0,npj-1 
              tdvs(i)=asc_dm(i)
            enddo
          else if(iv.eq.13) then
            do i=0,npj-1 
              tdvs(i)=tc_dm(i)
            enddo
#endif
          endif
          if(myrank.lt.SI_nof) then
            write(64) (tdvs(i),i=0,npj-1)
            do ip=ipsta,ipend,ipd
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
               ,ip,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              write(64) (tdvr(i),i=0,npjr(ip)-1)
            enddo
          else
            call MPI_ISEND(tdvs,npj,MPI_DOUBLE_PRECISION,srank &
             ,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
            call MPI_WAIT(ireqs(ip),istatus,ierr)
          endif
        enddo
        if(myrank.lt.SI_nof) then     
          close(64)
        endif

        deallocate(tivs)
        deallocate(tivr)
        deallocate(tdvs)
        deallocate(tdvr)

      endif

#ifdef DEL_DMPART
      if(myrank.eq.0) then
        open(60,file='output/data/del_dmpartcon.dat',status='unknown' &
          ,form='unformatted')   
        write(60) SI_dmdelci(0),SI_dmdelci(1)
        write(60) SI_deldmdt
        write(60) SI_mdmp0
        write(60) SI_mtdmdc0
        write(60) SI_ndmtdc
        write(60) SI_ndm1t
        close(60)
      endif
#endif

! end for DM
#endif

       deallocate(npjr)
       deallocate(ireqs)
       deallocate(ireqr)
       deallocate(istatus)

end subroutine
