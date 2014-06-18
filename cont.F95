#include "gcdp.def"
! ***************************************
!      cont.F95 for GCD+ ver. f03.0
!  18 Apr. 2014  written by D.KAWATA
! *************************************** 

subroutine set_ini(npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
      ,step,totstep,lstep,flagcont)
      use gcdp_const
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#ifdef METAL
      use gcdp_metal
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
      use gcdp_yields
#endif
#endif
#ifdef DM
      use gcdp_dm
#endif
#ifdef MHD
      use gcdp_mhd
#endif
      use gcdp_system
#ifdef TREEPM
      use gcdp_pm
#endif
 
      implicit none
      include 'mpif.h'

      integer,intent(out) :: npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg, &
        step,totstep,lstep,flagcont

      integer i,pn,ifl,ip,iv,nst
      integer inprocs,inof,inpt,indmt,indm1t
      integer npc,rrank,ipsta,ipend,ipd
      integer nivalb,nivaldm,ndvaldm,ndvalc,nivalc
      integer ndbval,ndbhyd,ndbsf,ndbmet,ndbcon,nddmcon,nibcon,nidmcon
      integer invali,invald
      integer ngt
      double precision dtmp(0:3)
      double precision mgt,mst,mbt,mms,mmet,tnsn
      character tchar*60
! * om0 (km s^-1 kpc^-1) *      
      double precision Ini_OMGg
! * (10^12 M@), (100 kpc),(1.0e4 K),(4.71e8 yr) * 
      double precision Ini_T0,Ini_ze,Ini_zi
      integer Ini_mul,Ini_nd
      integer flagtime,Ini_flagcop
      character filei*60
      character fileo*60
#if defined(COSM) && defined(BOXSIM)
      double precision iai
#endif
#ifdef MHD
      integer ndbmhd
#endif
! for work
      integer ierr,irdummy
      integer,allocatable :: istatus(:)
      integer,allocatable :: npjr(:),ireqs(:),ireqr(:)
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:)
#if defined(SF_EFD) || defined(SF_ZFD)
! work
      double precision,allocatable :: mfd(:)
! * for feedback  common *
! * Total *
      double precision,allocatable :: tmej_sn(:),tnsn_sn(:),tmzHe_sn(:) &
       ,tmzC_sn(:),tmzN_sn(:),tmzO_sn(:),tmzNe_sn(:),tmzMg_sn(:) &
       ,tmzSi_sn(:),tmzFe_sn(:),tmzZ_sn(:),nsw_sn(:)
#endif

#ifdef BCX_FIX  
      if(myrank.eq.0) then
        write(6,*) ' Error: BCX_FIX is not yet avairable for gcdpc'
      endif
      stop
#endif

! allocate array for MPI
      allocate(istatus(MPI_STATUS_SIZE))
      allocate(npjr(0:nprocs))
      allocate(ireqs(0:nprocs))
      allocate(ireqr(0:nprocs))

#ifdef COOL
      SI_flagrad=0
#else
      SI_flagrad=1
#endif

      if(myrank.eq.0) then
        open(50,file='ini/input.dat',status='old')
! *** Read Initial Parameter ***
! * Number of Particle etc. *  
        read(50,*) Ini_mul
! * Interval of Output Data, Distribution *
        read(50,*) Ini_nd
        read(50,*) SI_nof
#ifdef COSM
        read(50,*) Ini_OMGg
        read(50,*) Ini_ze
        TM_t0=0.0d0
#else
! * Initial Values unit is our Unit *
        read(50,*) Ini_zi
        read(50,*) flagtime
        read(50,*) TM_tot,TM_t0
        read(50,*) 
        Ini_ze=Ini_zi
! *** change time unit ***
        if(flagtime.ne.0) then
          TM_t0 = TM_t0/TMUGYR
        endif
#endif
#ifdef TREEPM
        read(50,*) nx_m,ny_m,nz_m
        if((mod(nx_m,2).ne.0)) then
          if(myrank.eq.0) then
            write(6,*) ' Error in set_ini(): nx_m,ny_m,nz_m=',nx_m,ny_m,nz_m
            write(6,*) ' nx_m should be able to be devided by 2.'
          endif
          call MPI_FINALIZE(ierr)
          stop
        endif
#endif
        read(50,*) Ini_flagcop
        if(Ini_flagcop.ne.0) then
          read(50,*) step
        endif
        close(50)      

! *** reading the data from inipara ***
        open(50,file='./output/ana/inipara.dat',status='old')
        read(50,'(a34,8I10)') tchar,npt,ndmt,np,ng,ns,ndm,SI_ndm1,SI_ndm1t
#ifdef COSM
        read(50,'(a19,4(1pE13.5))') tchar,SI_h0,SI_omg0,SI_lam0,Ini_OMGg
#endif
! *** need a better accuracy, changed since pv30.4 ***
        read(50,'(a5,1pE25.17)') tchar,OF_odt
#ifdef COSM
        read(50,'(a14,2(1pE25.17))') tchar,SI_eps0,SI_epsa
        SI_eps=SI_eps0
#else
        read(50,'(a8,1pE25.17)') tchar,SI_eps
#endif
#ifdef BOXSIM
        do i=1,100
          read(50,'(a50)') tchar
          if(tchar(2:4).eq.'box') then
            read(tchar(11:23),'(1pE13.5)') SI_lbox(0)
            read(tchar(24:36),'(1pE13.5)') SI_lbox(1)
            read(tchar(37:49),'(1pE13.5)') SI_lbox(2)
#ifdef COSM
            read(50,'(a22,3(1pE13.5))') tchar,SI_lbox0(0),SI_lbox0(1) &
              ,SI_lbox0(2)
#endif
            goto 70
          endif
        enddo
        write(6,*) ' Error: cannot find box size in output/ana/inipara.dat'
        stop
   70   continue
#endif        
#ifdef BCX_FIX
       write(6,*) ' Error: continuing BCX_FIX simulation is not supported yet.'
       stop 
#endif
        close(50)
        if(Ini_flagcop.eq.0) then

          open(51,file='output/data/cont.dat',status='old',form='unformatted')
! *** Read continue data ***
          read(51) step,totstep,lstep,SI_flagout,OF_odt,OF_todt &
#ifdef BOXSIM
           ,SI_lbox(0),SI_lbox(1),SI_lbox(2) &
#ifdef COSM
           ,SI_lbox0(0),SI_lbox0(1),SI_lbox0(2) &
#endif
#endif
           ,TM_tot,TM_sdt,TM_lt,TM_ldtmin,TM_k,TM_kn,SI_a,SI_eps,SI_mms &
#ifdef COSM
#if defined(COOL) || defined(STAR)
           ,SI_omg0,SI_lam0,SI_h0,SI_eps0,SI_epsa &
#else
           ,SI_omg0,SI_lam0,SI_h0,SI_eps0,SI_epsa &
#endif
#endif
#ifdef COOL
           ,SI_flagrad &
#endif
#ifdef STAR
           ,SI_esns,(SI_zsns(i),i=0,NYEL-1),TMsf_t,TMsf_dt,SI_tnsf &
#endif
           ,irdummy
          close(51)

        else

          allocate(tivr(0:2))
          allocate(tdvr(0:1))

          open(50,file='./output/ana/system.dat',status='old')
          write(6,*) ' myrank,step=',myrank,step

          do i=0,step
            read(50,*) 
          enddo
          read(50,150) TM_tot,totstep,step,ngt,tivr(0),nst &
#ifdef STAR
           ,SI_flagrad,SI_a,mgt,tdvr(0),mst,mbt,mmet,SI_mms &
           ,SI_esns,tnsn,SI_zsns(0),SI_zsns(1),SI_zsns(2),SI_zsns(3) &
           ,SI_zsns(4),SI_zsns(5),SI_zsns(6),SI_zsns(7),SI_zsns(8) &
           ,SI_tnsf,tivr(1),tivr(2) 
  150     format((1pE13.5),6I10,7(1pE13.5),11(1pE13.5),3I10)
#else
           ,SI_flagrad,SI_a,mgt,tdvr(0),mst,mbt,mmet,SI_mms 
  150     format((1pE13.5),6I10,7(1pE13.5))
#endif

          close(50)
#ifdef STAR
          TMsf_t=TM_tot
#endif    

          deallocate(tivr)
          deallocate(tdvr)

        endif

        if(Ini_flagcop.eq.0) then
          if(TM_kn.eq.TM_k.and.SI_flagout.ne.0) then
             OF_todt=OF_todt+OF_odt
             SI_flagout = 0
          endif
       else
          SI_flagout=0
          TM_dt=0.0d0
          TM_sdt=0.0d0
          TM_ldtmin=0.0d0
          lstep=0
          TM_k=0
          TM_kn=0
! *** This may need to be revised.
          OF_todt = OF_odt+TM_tot
        endif
      endif

! *** for integer values ***
#ifdef STAR
      nivalc=13
#else 
      nivalc=12
#endif

      allocate(tivr(0:nivalc))

      if(myrank.eq.0) then
        tivr(0)=step
        tivr(1)=totstep
        tivr(2)=lstep
        tivr(3)=SI_flagout
        tivr(4)=TM_k
        tivr(5)=TM_kn
        tivr(6)=SI_flagrad
        tivr(7)=SI_nof
        tivr(8)=npt
        tivr(9)=ndmt
        tivr(10)=SI_ndm1t
        tivr(11)=Ini_flagcop
#ifdef STAR
        tivr(12)=SI_tnsf
#endif
      endif
      call MPI_BCAST(tivr,nivalc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      step=tivr(0)
      totstep=tivr(1)
      lstep=tivr(2)
      SI_flagout=tivr(3)
      TM_k=tivr(4)
      TM_kn=tivr(5)
      SI_flagrad=tivr(6)
      SI_nof=tivr(7)
      npt=tivr(8)
      ndmt=tivr(9)
      SI_ndm1t=tivr(10)
      Ini_flagcop=tivr(11)
#ifdef STAR
      SI_tnsf=tivr(12)
#endif

      deallocate(tivr)


#ifdef TREEPM
! send grid information for TREEPM calculation
      nivalc=3

      allocate(tivr(0:nivalc))

      if(myrank.eq.0) then
        tivr(0)=nx_m
        tivr(1)=ny_m
        tivr(2)=nz_m
      endif
      call MPI_BCAST(tivr,nivalc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nx_m=tivr(0)
      ny_m=tivr(1)
      nz_m=tivr(2)

      deallocate(tivr)
#endif

#ifdef COSM
#ifdef STAR
      ndvalc=17+NYEL
#else
      ndvalc=15
#endif
#else
#ifdef STAR
      ndvalc=13+NYEL
#else
      ndvalc=11
#endif
#endif

      allocate(tdvr(0:ndvalc))

! *** for double values ***
      if(myrank.eq.0) then
        tdvr(0)=OF_odt
        tdvr(1)=OF_todt
        tdvr(2)=TM_tot
        tdvr(3)=TM_sdt
        tdvr(4)=TM_lt
        tdvr(5)=TM_ldtmin
        tdvr(6)=SI_a
        tdvr(7)=Ini_ze
        tdvr(8)=SI_eps
        tdvr(9)=SI_mms
! *** for cosmological simulation
#ifdef COSM
        tdvr(10)=SI_omg0
        tdvr(11)=SI_lam0
        tdvr(12)=SI_h0
        tdvr(13)=SI_eps0
        tdvr(14)=SI_epsa
#ifdef STAR
        tdvr(15)=SI_esns
        tdvr(16)=TMsf_t
        tdvr(17)=TMsf_dt
        do i=0,NYEL-1
          tdvr(18+i)=SI_zsns(i)
        enddo
#endif
#else
! *** for non-cosmological simulation with STAR
        tdvr(10)=TM_t0
#ifdef STAR
        tdvr(11)=SI_esns
        tdvr(12)=TMsf_t
        tdvr(13)=TMsf_dt
        do i=0,NYEL-1
          tdvr(14+i)=SI_zsns(i)
        enddo
#endif
#endif
      endif

      call MPI_BCAST(tdvr,ndvalc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! *** set the values ***
      OF_odt=tdvr(0)
      OF_todt=tdvr(1)
      TM_tot=tdvr(2)
      TM_sdt=tdvr(3)
      TM_lt=tdvr(4)
      TM_ldtmin=tdvr(5)
      SI_a=tdvr(6)
      Ini_ze=tdvr(7)
      SI_eps=tdvr(8)
      SI_mms=tdvr(9)
! *** for cosmological simulation
#ifdef COSM
      SI_omg0=tdvr(10)
      SI_lam0=tdvr(11)
      SI_h0=tdvr(12)
      SI_eps0=tdvr(13)
      SI_epsa=tdvr(14)
#ifdef STAR
      SI_esns=tdvr(15)
      TMsf_t=tdvr(16)
      TMsf_dt=tdvr(17)
      do i=0,NYEL-1
        SI_zsns(i)=tdvr(18+i)
      enddo
#endif
#else
      TM_t0=tdvr(10)
! *** for non-cosmological simulation with STAR
#ifdef STAR
      SI_esns=tdvr(11)
      TMsf_t=tdvr(12)
      TMsf_dt=tdvr(13)
      do i=0,NYEL-1
        SI_zsns(i)=tdvr(14+i)
      enddo
#endif
#endif
#if defined(COSM) && defined(BOXSIM)
! *** Set Hubble parameter ***
      iai = 1.0d0/SI_a
      SI_hub = SI_h0*HUB0*dsqrt(SI_omg0*iai*iai*iai &
       -(SI_omg0+SI_lam0-1.0)*iai*iai+SI_lam0)
#endif

      deallocate(tdvr)

! *** set TM_mindtsfd
#if defined(SF) || defined(SF_EFD) || defined(SF_ZFD)
! *** minimum time step for SF and feedback
      if(TM_ldtmin.gt.0.0d0) then
        i=int(dlog(MINTMSFD/TM_ldtmin)/dlog(2.0d0))
        TM_mindtsfd=TM_ldtmin*(2.0d0**dble(i))
   90   if(TM_mindtsfd.gt.MINTMSFD) then
          TM_mindtsfd=TM_mindtsfd*0.5d0
        endif
        if(TM_mindtsfd.gt.MINTMSFD) then 
          i=i-1
          if(myrank.eq.0) then
            write(6,*) ' TM_mindtsfd devided by 2',i
          endif
          goto 90
        endif
        if(myrank.eq.0) then
          write(6,*) ' Min dt for SF and FD=',i,TM_mindtsfd,TM_mindtsfd*TMUGYR
        endif
      else
        TM_mindtsfd=0.0d0
      endif
#endif

#ifdef BOXSIM

#ifdef COSM
      ndvalc=6
#else
      ndvalc=3
#endif

      allocate(tdvr(0:ndvalc-1))

      if(myrank.eq.0) then
        tdvr(0)=SI_lbox(0)
        tdvr(1)=SI_lbox(1)
        tdvr(2)=SI_lbox(2)
#ifdef COSM
        tdvr(3)=SI_lbox0(0)
        tdvr(4)=SI_lbox0(1)
        tdvr(5)=SI_lbox0(2)
#endif       
      endif
      call MPI_BCAST(tdvr,ndvalc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      SI_lbox(0)=tdvr(0)
      SI_lbox(1)=tdvr(1)
      SI_lbox(2)=tdvr(2)
#ifdef COSM
      SI_lbox0(0)=tdvr(3)
      SI_lbox0(1)=tdvr(4)
      SI_lbox0(2)=tdvr(5)
#endif       

      deallocate(tdvr)

#endif

! *** set number of output values ***
#if defined(SF_EFD) || defined(SF_ZFD)
      nivalb=4
#else
      nivalb=3
#endif
      nivaldm=2
      ndvaldm=9
#ifdef GAS
      ndbval=9
! since pv33.17.6
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
#ifdef MHD
      ndbmhd=10
#ifdef DCLEAN
      ndbmhd=ndbmhd+1
#endif
#endif

      if(Ini_flagcop.eq.0) then
        if(myrank.eq.0) then
          write(6,*) ' data read from output*cont.dat'
        endif
        flagcont=1
        ndbcon=25
#if defined(COSM) && defined(BOXSIM)
        ndbcon=ndbcon+4
#endif
        nibcon=1
#if defined(COSM) && defined(BOXSIM)
        nddmcon=14
#else
        nddmcon=12
#endif
        nidmcon=1
      else
        if(myrank.eq.0) then
          write(6,*) ' data read from output of step=',step
        endif
        flagcont=2
        ndbcon=0
        nibcon=0
        nddmcon=0
        nidmcon=0
      endif

! *** Change Lambda for time ***
#ifdef COSM
      SI_lam0 = SI_lam0/(HUB0*HUB0*SI_h0*SI_h0)             
!*********************************
!  Make Time tabele and Set TM_t0
! ******************************** 
      call makettable(1.0d0/(1.0d0+Ini_ze))
#endif
      if(TM_tot.ge.TM_t0-0.000001d0) then
        write(6,*) ' *** This task is already finished !! ***',myrank
        write(6,*) ' Ttot,Tt0=',TM_tot,TM_t0
        call MPI_FINALIZE()
        stop
      endif

      if(myrank.eq.0) then
        write(6,*) '##### Initial  Condition Information #####'      
#ifdef COSM
        write(6,*) ' h, Omg0, lamda0 = ',SI_h0,SI_omg0,SI_lam0
#endif
        write(6,*) ' gamma       = ',GAM
        write(6,*) ' simulate z  = ',1.0d0/SI_a-1.0d0,'->',Ini_ze
        write(6,*) '          t  = ',TM_tot,'->',TM_t0
        write(6,*) '##### Continue  Condition Information #####'
        write(6,*) ' step, t, l  = ',step,totstep,lstep
        write(6,*) ' SI_flagrad=',SI_flagrad
        write(6,*) ' # of output file per step =',SI_nof
        if(Ini_flagcop.eq.0) then
          write(6,*) ' k, kn=',TM_k,TM_kn &
          ,' sdt=',TM_sdt,' lt=',TM_lt,' ldt=',TM_ldtmin
        endif
#ifdef BOXSIM
        write(6,*) ' box size=',SI_lbox(0),SI_lbox(1),SI_lbox(2)
#ifdef COSM
        write(6,*) ' box size (comoving)=',SI_lbox0(0),SI_lbox0(1) &
          ,SI_lbox0(2)
#endif
#endif
      endif

! *** Lambda for acceleration ***
#ifdef COSM
      SI_lam0 = SI_lam0*HUB0*HUB0*SI_h0*SI_h0
#endif
      if(myrank.lt.SI_nof) then
        ipsta=myrank
        ipend=nprocs-1
        ipd=SI_nof
! dummy
        rrank=myrank
      else
        ipsta=myrank
        ipend=myrank
        ipd=1
        rrank=mod(myrank,SI_nof)
      endif
      do i=0,nprocs-1
        npjr(i)=0
      enddo
#ifdef GAS
      if(npt.gt.0) then
! *** reading baryon particle data **      
        if(myrank.lt.SI_nof) then
          if(Ini_flagcop.eq.0) then
            write(fileo,'(a17,i4.4,a4)') 'output/data/pbval',myrank,'.dat'
          else
            write(fileo,'(a20,i6.6,a1,i4.4)') &
             './output/data/bbvals',step,'n',myrank
          endif
          open(51,file=fileo,status='old',form='unformatted')
          read(51) inpt,indmt,indm1t,SI_a,TM_tot
          read(51) inprocs,inof,invali,invald
          if(nivalb.ne.invali) then
            write(6,*) myrank,' Error: nivalb is inconsistent'
            write(6,*) ' in file, cont.F=',invali,nivalb
            stop
          endif
          if(invald.ne.ndbval) then
            write(6,*) myrank,' Error: ndbval is inconsistent'
            write(6,*) ' in file, cont.F=',invald,ndbval
            stop
          endif
! check consistency
          if(inpt.ne.npt.or.(indmt.ne.ndmt.or.indm1t.ne.SI_ndm1t)) then
            write(6,*) ' Error in cont.f: # of particles is inconsistent'
            write(6,*) ' npt,ndmt(inipara,output/data)',npt,ndmt,inpt,indmt
            stop
          endif
          if(inprocs.ne.nprocs.or.inof.ne.SI_nof) then
            write(6,*) ' Error in cont.f :'
            write(6,*) '  nprocs or SI_nof should not be changed.'
            write(6,*) '  previous nprocs and SI_nof=',inprocs,inof
            write(6,*) '  current nprocs and SI_nof=',nprocs,SI_nof
            stop
          endif
        endif
! *** reading the data ***

        allocate(tivr(0:7))

        do ip=ipsta,ipend,ipd
          if(myrank.lt.SI_nof) then
            read(51) (tivr(i),i=0,7)
            if(ip.ne.myrank) then
              call MPI_ISEND(tivr,8,MPI_INTEGER,ip,1,MPI_COMM_WORLD &
               ,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif
          else
            call MPI_IRECV(tivr,8,MPI_INTEGER,rrank,1,MPI_COMM_WORLD &
             ,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(ip.eq.myrank) then
            ng=tivr(0)
            ndm=tivr(1)
            ns=tivr(2)
            SI_ndm1=tivr(3)
            nag=tivr(4)
            nadm=tivr(5)
            nas=tivr(6)
            nagravg=tivr(7)
          endif
          npjr(ip)=tivr(0)+tivr(2)
        enddo

        deallocate(tivr)

! allocate baryon data
#if defined(GAS) || defined(STAR)
        call allocate_baryon_all(npjr(myrank))
#endif

! *** read integer variables ***
        do iv=0,nivalb-1
          do ip=ipsta,ipend,ipd

            allocate(tivr(0:npjr(ip)))
    
            if(myrank.lt.SI_nof) then
              read(51) (tivr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                  ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,rrank,iv &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
            if(ip.eq.myrank) then
              if(iv.eq.0) then
                do i=0,npjr(ip)-1
                  id_p(i)=tivr(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npjr(ip)-1
                  flagc_p(i)=tivr(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npjr(ip)-1
                  list_ap(i)=tivr(i)
                enddo
#if defined(SF_EFD) || defined(SF_ZFD)
              else if(iv.eq.3) then                  
                do i=0,npjr(ip)-1
                  flagfd_p(i)=tivr(i)
                enddo
#endif
              endif
            endif

            deallocate(tivr)

          enddo
        enddo

! *** read baryon basic values ***
        do iv=0,ndbval-1
          do ip=ipsta,ipend,ipd

            allocate(tdvr(0:npjr(ip)))

            if(myrank.lt.SI_nof) then
              read(51) (tdvr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                 ,ip,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
               ,rrank,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
            if(ip.eq.myrank) then
              if(iv.eq.0) then
                do i=0,npjr(ip)-1 
                  xc_p(i)=tdvr(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npjr(ip)-1 
                  yc_p(i)=tdvr(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npjr(ip)-1 
                  zc_p(i)=tdvr(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npjr(ip)-1 
                  vnx_p(i)=tdvr(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npjr(ip)-1 
                  vny_p(i)=tdvr(i)
                enddo
              else if(iv.eq.5) then
                do i=0,npjr(ip)-1 
                  vnz_p(i)=tdvr(i)
                enddo
              else if(iv.eq.6) then
                do i=0,npjr(ip)-1 
                  m_p(i)=tdvr(i)
                enddo
              else if(iv.eq.7) then
                do i=0,npjr(ip)-1 
                  rho_p(i)=tdvr(i)
                enddo
              else if(iv.eq.8) then
                do i=0,npjr(ip)-1 
                  u_p(i)=tdvr(i)
                enddo
              endif
            endif

            deallocate(tdvr)

          enddo
        enddo
        if(myrank.lt.SI_nof) then
          close(51)
        endif
! *** file for baryon hydro values ***
        if(myrank.lt.SI_nof) then
          if(Ini_flagcop.eq.0) then
            write(fileo,'(a17,i4.4,a4)') 'output/data/pbhyd',myrank,'.dat'
          else
            write(fileo,'(a20,i6.6,a1,i4.4)') &
              './output/data/bbhyds',step,'n',myrank
          endif
          open(53,file=fileo,status='unknown',form='unformatted')
          read(53) inpt,indmt,indm1t,dtmp(0),dtmp(1)
          read(53) inprocs,inof,invald
          if(ndbhyd.ne.invald.and.Ini_flagcop.eq.0) then
            write(6,*) myrank,' Error: ndbhyd is inconsistent'
            write(6,*) ' in file, cont.F=',invald,ndbhyd
            stop
          endif
        endif
! *** read hydro values ***
        do iv=0,ndbhyd-1
          do ip=ipsta,ipend,ipd

            allocate(tdvr(0:npjr(ip)))

            if(myrank.lt.SI_nof) then
              read(53) (tdvr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,ip,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,rrank,iv &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
            if(ip.eq.myrank) then
              if(iv.eq.0) then
                do i=0,npjr(ip)-1 
                  h_p(i)=tdvr(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npjr(ip)-1 
                  div_v_p(i)=tdvr(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npjr(ip)-1 
                  alpv_p(i)=tdvr(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npjr(ip)-1 
                  alpu_p(i)=tdvr(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npjr(ip)-1 
                  myu_p(i)=tdvr(i)
                enddo
              endif
            endif

            deallocate(tdvr)

          enddo
        enddo
        if(myrank.lt.SI_nof) then
          close(53)
        endif
! *** Read MHD values ***
#ifdef MHD 
        call allocate_mhd(npjr(myrank))
      
        if(myrank.lt.SI_nof) then
          if(Ini_flagcop.eq.0) then
            write(fileo,'(a17,i4.4,a4)') 'output/data/pbmhd',myrank,'.dat'
          else
            write(fileo,'(a20,i6.6,a1,i4.4)') &
             './output/data/bbmhds',step,'n',myrank
          endif
          open(58,file=fileo,status='old',form='unformatted')
          read(58) inpt,indmt,indm1t,dtmp(0),dtmp(1)
          read(58) inprocs,inof,invald
        endif 
        do iv=0,ndbmhd-1
          do ip=ipsta,ipend,ipd

            allocate(tdvr(0:npjr(ip)))

            if(myrank.lt.SI_nof) then
              read(58) (tdvr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                 ,ip,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
               ,rrank,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif

            if(ip.eq.myrank) then
              if(iv.eq.0)  then
                do i=0,npjr(ip)-1
                  bx_p(i)=tdvr(i)
                enddo
              else if (iv.eq.1) then
                do i=0,npjr(ip)-1
                  by_p(i)=tdvr(i)
                enddo
              else if (iv.eq.2) then
                do i=0,npjr(ip)-1
                  bz_p(i)=tdvr(i)
                enddo
              else if (iv.eq.3) then
                do i=0,npjr(ip)-1
                  divb_p(i)=tdvr(i)
                enddo
              else if (iv.eq.4) then
                do i=0,npjr(ip)-1
                  arotb_p(i)=tdvr(i)
                enddo
              else if (iv.eq.5) then
                do i=0,npjr(ip)-1
                  agradb_p(i)=tdvr(i)
                enddo
              else if (iv.eq.6) then
                do i=0,npjr(ip)-1
                  alpb_p(i)=tdvr(i)
                enddo
              else if (iv.eq.7) then
                do i=0,npjr(ip)-1
                  valfi_p(i)=tdvr(i)
                enddo
              else if (iv.eq.8) then
                do i=0,npjr(ip)-1
                  pb_p(i)=tdvr(i)
                enddo
              else if (iv.eq.9) then
                do i=0,npjr(ip)-1
                  divberr_p(i)=tdvr(i)
                enddo
#ifdef DCLEAN
              else if (iv.eq.10) then
                do i=0,npjr(ip)-1
                  psi_p(i)=tdvr(i)
                enddo
#endif
              endif
            endif

            deallocate(tdvr)

          enddo
        enddo
        if(myrank.lt.SI_nof) then
          close(58)
        endif
#endif

#ifdef STAR
        if(myrank.lt.SI_nof) then
          if(Ini_flagcop.eq.0) then
            write(fileo,'(a17,i4.4,a4)') 'output/data/pbsfi',myrank,'.dat'
          else
            write(fileo,'(a20,i6.6,a1,i4.4)') &
             './output/data/bbsfis',step,'n',myrank
          endif
          open(55,file=fileo,status='unknown',form='unformatted')
          read(55) inpt,indmt,indm1t,dtmp(0),dtmp(1)
          read(55) inprocs,inof,invald
          if(ndbsf.ne.invald) then
            write(6,*) myrank,' Error: ndbsf is inconsistent'
            write(6,*) ' in file, cont.F=',invald,ndbsf
            stop
          endif
        endif
! *** read star formatioin info ****
        do iv=0,ndbsf-1
          do ip=ipsta,ipend,ipd

            allocate(tdvr(0:npjr(ip)))

            if(myrank.lt.SI_nof) then
              read(55) (tdvr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,ip,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,rrank,iv &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
            if(ip.eq.myrank) then
              if(iv.eq.0) then
                do i=0,npjr(ip)-1 
                  ts_p(i)=tdvr(i)
                enddo
              endif
            endif

            deallocate(tdvr)

          enddo
        enddo
        if(myrank.lt.SI_nof) then
          close(55)
        endif
#endif
! *** file for baryon metal values ***
#ifdef METAL
        if(myrank.lt.SI_nof) then
          if(Ini_flagcop.eq.0) then
            write(fileo,'(a17,i4.4,a4)') 'output/data/pbmet',myrank,'.dat'
          else
            write(fileo,'(a20,i6.6,a1,i4.4)') &
             './output/data/bbmets',step,'n',myrank
          endif
          open(57,file=fileo,status='unknown',form='unformatted')
          read(57) inpt,indmt,indm1t,dtmp(0),dtmp(1)
          read(57) inprocs,inof,invald
          if(ndbmet.ne.invald) then
            write(6,*) myrank,' Error: ndbmet is inconsistent'
            write(6,*) ' in file, cont.F=',invald,ndbmet
            stop
          endif
        endif
! *** read metal info ****
        do iv=0,ndbmet-1
          do ip=ipsta,ipend,ipd

            allocate(tdvr(0:npjr(ip)))

            if(myrank.lt.SI_nof) then
              read(57) (tdvr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,ip,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,rrank,iv &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
            if(ip.eq.myrank) then
              if(iv.eq.0) then
                do i=0,npjr(ip)-1 
                  mzHe_p(i)=tdvr(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npjr(ip)-1 
                  mzZ_p(i)=tdvr(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npjr(ip)-1 
                  mzC_p(i)=tdvr(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npjr(ip)-1 
                  mzN_p(i)=tdvr(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npjr(ip)-1 
                  mzO_p(i)=tdvr(i)
                enddo
              else if(iv.eq.5) then
                do i=0,npjr(ip)-1 
                  mzNe_p(i)=tdvr(i)
                enddo
              else if(iv.eq.6) then
                do i=0,npjr(ip)-1 
                  mzMg_p(i)=tdvr(i)
                enddo
              else if(iv.eq.7) then
                do i=0,npjr(ip)-1 
                  mzSi_p(i)=tdvr(i)
                enddo
              else if(iv.eq.8) then
                do i=0,npjr(ip)-1 
                  mzFe_p(i)=tdvr(i)
                enddo
              else if(iv.eq.9) then
                do i=0,npjr(ip)-1 
                  mzHe0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.10) then
                do i=0,npjr(ip)-1 
                  mzZ0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.11) then
                do i=0,npjr(ip)-1 
                  mzC0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.12) then
                do i=0,npjr(ip)-1 
                  mzN0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.13) then
                do i=0,npjr(ip)-1 
                  mzO0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.14) then
                do i=0,npjr(ip)-1 
                  mzNe0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.15) then
                do i=0,npjr(ip)-1 
                  mzMg0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.16) then
                do i=0,npjr(ip)-1 
                  mzSi0_p(i)=tdvr(i)
                enddo
              else if(iv.eq.17) then
                do i=0,npjr(ip)-1 
                  mzFe0_p(i)=tdvr(i)
                enddo
#ifdef METALDIFF
              else if(iv.eq.18) then
                do i=0,npjr(ip)-1 
                  vsig_p(i)=tdvr(i)
                enddo
#endif
              endif
            endif

            deallocate(tdvr)

          enddo
        enddo
        if(myrank.lt.SI_nof) then
          close(57)
        endif
#endif
        if(Ini_flagcop.eq.0) then
          if(myrank.lt.SI_nof) then
            write(fileo,'(a17,i4.4,a4)') 'output/data/pbcon',myrank,'.dat'
            open(41,file=fileo,status='unknown',form='unformatted')
            read(41) inpt,indmt,indm1t,dtmp(0),dtmp(1)
            read(41) inprocs,inof,invali,invald
            if(ndbcon.ne.invald) then
              write(6,*) myrank,' Error: ndbcon is inconsistent'
              write(6,*) ' in file, cont.F=',invald,ndbcon
              stop
            endif
            if(nibcon.ne.invali) then
              write(6,*) myrank,' Error: nibcon is inconsistent'
              write(6,*) ' in file, cont.F=',invali,nibcon
              stop
            endif
          endif
! *** read continue files ***
! *** read integer variables ***
          do iv=0,nibcon-1
            do ip=ipsta,ipend,ipd

              allocate(tivr(0:npjr(ip)))

              if(myrank.lt.SI_nof) then
                read(41) (tivr(i),i=0,npjr(ip)-1)
                if(ip.ne.myrank) then
                  call MPI_ISEND(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                   ,MPI_COMM_WORLD,ireqs(ip),ierr)
                  call MPI_WAIT(ireqs(ip),istatus,ierr)
                endif
              else
                call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,rrank,iv &
                 ,MPI_COMM_WORLD,ireqr(ip),ierr)
                call MPI_WAIT(ireqr(ip),istatus,ierr)
              endif
              if(ip.eq.myrank) then
                if(iv.eq.0) then
                  do i=0,npjr(ip)-1
                    flagt_p(i)=tivr(i)
                  enddo
                endif
              endif

              deallocate(tivr)

            enddo
          enddo
! *** read double precision ***
          do iv=0,ndbcon-1
            do ip=ipsta,ipend,ipd

              allocate(tdvr(0:npjr(ip)))

              if(myrank.lt.SI_nof) then
                read(41) (tdvr(i),i=0,npjr(ip)-1)
                if(ip.ne.myrank) then
                  call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                   ,ip,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                  call MPI_WAIT(ireqs(ip),istatus,ierr)
                endif
              else
                call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                 ,rrank,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                call MPI_WAIT(ireqr(ip),istatus,ierr)
              endif
              if(ip.eq.myrank) then
                if(iv.eq.0) then
                  do i=0,npjr(ip)-1 
                    lt_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.1) then
                  do i=0,npjr(ip)-1 
                    dt_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.2) then
                  do i=0,npjr(ip)-1 
                    x_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.3) then
                  do i=0,npjr(ip)-1 
                    y_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.4) then
                  do i=0,npjr(ip)-1 
                    z_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.5) then
                  do i=0,npjr(ip)-1 
                    vvnx_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.6) then
                  do i=0,npjr(ip)-1 
                    vvny_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.7) then
                  do i=0,npjr(ip)-1 
                    vvnz_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.8) then
                  do i=0,npjr(ip)-1 
                    dvx_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.9) then
                  do i=0,npjr(ip)-1 
                    dvy_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.10) then
                  do i=0,npjr(ip)-1 
                    dvz_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.11) then
                  do i=0,npjr(ip)-1 
                    zetah_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.12) then
                  do i=0,npjr(ip)-1 
                    omgh_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.13) then
                  do i=0,npjr(ip)-1 
                    arot_v_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.14) then
                  do i=0,npjr(ip)-1 
                    d2u_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.15) then
                  do i=0,npjr(ip)-1 
                    dtr_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.16) then
                  do i=0,npjr(ip)-1 
                    ndu_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.17) then
                  do i=0,npjr(ip)-1 
                    pdu_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.18) then
                  do i=0,npjr(ip)-1 
                    dtg_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.19) then
                  do i=0,npjr(ip)-1 
                    ltg_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.20) then
                  do i=0,npjr(ip)-1 
                    ax_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.21) then
                  do i=0,npjr(ip)-1 
                    ay_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.22) then
                  do i=0,npjr(ip)-1 
                    az_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.23) then
#ifdef COOL 
                  do i=0,npjr(ip)-1
                    ram_p(i)=tdvr(i)
                  enddo
#endif
#ifdef SF_EFD
                else if(iv.eq.24) then
                  do i=0,npjr(ip)-1 
                    Gsn_p(i)=tdvr(i)
                  enddo
#endif
#if defined(COSM) && defined(BOXSIM)
                else if(iv.eq.25) then
                  do i=0,npjr(ip)-1 
                    asc_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.26) then
                  do i=0,npjr(ip)-1 
                    tc_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.27) then
                  do i=0,npjr(ip)-1 
                    ascg_p(i)=tdvr(i)
                  enddo
                else if(iv.eq.28) then
                  do i=0,npjr(ip)-1 
                    tcg_p(i)=tdvr(i)
                  enddo
#endif
                endif
              endif

              deallocate(tdvr)

            enddo
          enddo
          if(myrank.lt.SI_nof) then
            close(41)
          endif
        endif
      endif
! *** end GAS ***
#endif
! *** end GAS ***

#ifdef DM
! *** for DM data ***
      if(ndmt.gt.0) then
! *** reading DM particle data **      
        if(myrank.lt.SI_nof) then
          if(Ini_flagcop.eq.0) then
            write(fileo,'(a17,i4.4,a4)') 'output/data/pdval',myrank,'.dat'
          else
            write(fileo,'(a20,i6.6,a1,i4.4)') &
             './output/data/bdvals',step,'n',myrank
          endif
          open(52,file=fileo,status='old',form='unformatted')
          read(52) inpt,indmt,indm1t,SI_a,TM_tot
          read(52) inprocs,inof,invali,invald
          if(nivaldm.ne.invali) then
            write(6,*) myrank,' Error: nivaldm is inconsistent'
            write(6,*) ' in file, cont.F=',invali,nivaldm
            stop
          endif
          if(invald.ne.ndvaldm) then
            write(6,*) myrank,' Error: ndvaldm is inconsistent'
            write(6,*) ' in file, cont.F=',invald,ndvaldm
            stop
          endif
          if(inpt.ne.npt.or.(indmt.ne.ndmt.or.indm1t.ne.SI_ndm1t)) then
            write(6,*) ' Error in cont.f: # of particles is inconsistent'
            write(6,*) ' npt,ndmt(inipara,output/data/dm)',npt,ndmt,inpt,indmt
            stop
          endif
          if(inprocs.ne.nprocs.or.inof.ne.SI_nof) then
            write(6,*) ' Error in cont.f :'
            write(6,*) '  nprocs or SI_nof should not be changed.'
            write(6,*) '  previous nprocs and SI_nof=',inprocs,inof
            write(6,*) '  current nprocs and SI_nof from dm file=',nprocs,SI_nof
            stop
          endif
        endif

        allocate(tivr(0:6))

        do ip=ipsta,ipend,ipd
          if(myrank.lt.SI_nof) then
            read(52) (tivr(i),i=0,6)
            if(ip.ne.myrank) then
              call MPI_ISEND(tivr,7,MPI_INTEGER,ip,1,MPI_COMM_WORLD &
               ,ireqs(ip),ierr)
              call MPI_WAIT(ireqs(ip),istatus,ierr)
            endif
          else
            call MPI_IRECV(tivr,7,MPI_INTEGER,rrank,1,MPI_COMM_WORLD &
             ,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
          if(ip.eq.myrank) then
            ng=tivr(0)
            ndm=tivr(1)
            ns=tivr(2)
            SI_ndm1=tivr(3)
            nag=tivr(4)
            nadm=tivr(5)
            nas=tivr(6)
          endif
          npjr(ip)=tivr(1)
        enddo

        deallocate(tivr)

! allocate DM values
        call allocate_dm_all(ndm)

! *** read integer variables ***
        do iv=0,nivaldm-1
          do ip=ipsta,ipend,ipd

            allocate(tivr(0:npjr(ip)))

            if(myrank.lt.SI_nof) then
              read(52) (tivr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,rrank,iv &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
            if(ip.eq.myrank) then
              if(iv.eq.0) then
                do i=0,npjr(ip)-1
                  id_dm(i)=tivr(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npjr(ip)-1
                  list_adm(i)=tivr(i)
                enddo
              endif
            endif

            deallocate(tivr)

          enddo
        enddo
! *** read DM double variables ***
!          write(6,*) ' ndvaldm=',ndvaldm,
        do iv=0,ndvaldm-1
          do ip=ipsta,ipend,ipd

            allocate(tdvr(0:npjr(ip)))

            if(myrank.lt.SI_nof) then
              read(52) (tdvr(i),i=0,npjr(ip)-1)
              if(ip.ne.myrank) then
                call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,ip,iv &
                 ,MPI_COMM_WORLD,ireqs(ip),ierr)
                call MPI_WAIT(ireqs(ip),istatus,ierr)
              endif
            else
              call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION,rrank,iv &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
            if(ip.eq.myrank) then
              if(iv.eq.0) then
                do i=0,npjr(ip)-1 
                  xc_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.1) then
                do i=0,npjr(ip)-1 
                  yc_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.2) then
                do i=0,npjr(ip)-1 
                  zc_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.3) then
                do i=0,npjr(ip)-1 
                  vnx_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.4) then
                do i=0,npjr(ip)-1 
                  vny_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.5) then
                do i=0,npjr(ip)-1 
                  vnz_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.6) then
                do i=0,npjr(ip)-1 
                  m_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.7) then
                do i=0,npjr(ip)-1 
                  rho_dm(i)=tdvr(i)
                enddo
              else if(iv.eq.8) then
                do i=0,npjr(ip)-1 
                  h_dm(i)=tdvr(i)
                enddo
              endif
            endif

            deallocate(tdvr)

          enddo
        enddo
        if(myrank.lt.SI_nof) then
          close(52) 
        endif
        if(Ini_flagcop.eq.0) then
          if(myrank.lt.SI_nof) then
! *** for continue
            write(fileo,'(a17,i4.4,a4)') 'output/data/pdcon',myrank,'.dat'
            open(54,file=fileo,status='old',form='unformatted') 
            read(54) inpt,indmt,indm1t,dtmp(0),dtmp(1)
            read(54) inprocs,inof,invali,invald
            if(invald.ne.nddmcon) then
              write(6,*) myrank,' Error: nddmcon is inconsistent'
              write(6,*) ' in file, cont.F=',invald,nddmcon
              stop
            endif
            if(invali.ne.nidmcon) then
              write(6,*) myrank,' Error: nidmcon is inconsistent'
              write(6,*) ' in file, cont.F=',invali,nidmcon
              stop
            endif
          endif
! *** read integer variables ***
          do iv=0,nidmcon-1
            do ip=ipsta,ipend,ipd

              allocate(tivr(0:npjr(ip)))

              if(myrank.lt.SI_nof) then
                read(54) (tivr(i),i=0,npjr(ip)-1)
                if(ip.ne.myrank) then
                  call MPI_ISEND(tivr,npjr(ip),MPI_INTEGER,ip,iv &
                   ,MPI_COMM_WORLD,ireqs(ip),ierr)
                  call MPI_WAIT(ireqs(ip),istatus,ierr)
                endif
              else
                call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,rrank,iv &
                 ,MPI_COMM_WORLD,ireqr(ip),ierr)
                call MPI_WAIT(ireqr(ip),istatus,ierr)
              endif
              if(ip.eq.myrank) then
                if(iv.eq.0) then
                  do i=0,npjr(ip)-1
                    flagt_dm(i)=tivr(i)
                  enddo
                endif
              endif

              deallocate(tivr)

            enddo
          enddo
! *** for continue variables ***
          do iv=0,nddmcon-1
            do ip=ipsta,ipend,ipd

              allocate(tdvr(0:npjr(ip)))

              if(myrank.lt.SI_nof) then
                read(54) (tdvr(i),i=0,npjr(ip)-1)
                if(ip.ne.myrank) then
                  call MPI_ISEND(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                   ,ip,iv,MPI_COMM_WORLD,ireqs(ip),ierr)
                  call MPI_WAIT(ireqs(ip),istatus,ierr)
                endif
              else
                call MPI_IRECV(tdvr,npjr(ip),MPI_DOUBLE_PRECISION &
                 ,rrank,iv,MPI_COMM_WORLD,ireqr(ip),ierr)
                call MPI_WAIT(ireqr(ip),istatus,ierr)
              endif
              if(ip.eq.myrank) then
                if(iv.eq.0) then
                  do i=0,npjr(ip)-1 
                    lt_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.1) then
                  do i=0,npjr(ip)-1 
                    dt_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.2) then
                  do i=0,npjr(ip)-1 
                    x_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.3) then
                  do i=0,npjr(ip)-1 
                    y_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.4) then
                  do i=0,npjr(ip)-1 
                    z_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.5) then
                  do i=0,npjr(ip)-1 
                    dvx_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.6) then
                  do i=0,npjr(ip)-1 
                    dvy_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.7) then
                  do i=0,npjr(ip)-1 
                    dvz_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.8) then
                  do i=0,npjr(ip)-1 
                    zetah_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.9) then
                  do i=0,npjr(ip)-1 
                    omgh_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.10) then
                  do i=0,npjr(ip)-1 
                    div_v_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.11) then
                  do i=0,npjr(ip)-1 
                    dtr_dm(i)=tdvr(i)
                  enddo
#if defined(COSM) && defined(BOXSIM)
                else if(iv.eq.12) then
                  do i=0,npjr(ip)-1 
                    asc_dm(i)=tdvr(i)
                  enddo
                else if(iv.eq.13) then
                  do i=0,npjr(ip)-1 
                    tc_dm(i)=tdvr(i)
                  enddo
#endif
                endif
              endif

              deallocate(tdvr)

            enddo
          enddo
          if(myrank.lt.SI_nof) then
            close(54)
          endif
        endif
      endif
! *** end for DM ***
#endif
! *** end for DM ***
!      write(6,*) ' myrank=',myrank,' DM data read h_dm=',h_dm(1)
! *** set the various values ***
! *** send TM_tot,SI_a in case reading from system.dat***

      allocate(tdvr(0:1))

      if(myrank.eq.0) then
        tdvr(0)=TM_tot
        tdvr(1)=SI_a
      endif
      call MPI_BCAST(tdvr,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! *** set the values ***
      TM_tot=tdvr(0)
      SI_a=tdvr(1)

      deallocate(tdvr)

! reset
      if(Ini_flagcop.ne.0) then
        OF_todt = OF_odt+TM_tot
      endif
      np=ng+ns
#ifdef GAS
      do i = 0,ng-1
! *** Set Velocity and Position ***
        pn = list_ap(i)
! *** Set Hydro Value ***
        p_p(pn)=(GAM-1.0)*rho_p(pn)*u_p(pn)
        cs_p(pn)=dsqrt(gam*p_p(pn)/rho_p(pn))
        as_p(pn)=p_p(pn)/(rho_p(pn)**GAM)
        if(flagc_p(pn).gt.0) then
          write(6,*) ' Error in cont.f: stars in gas list'
          write(6,*) ' i,pn,id,flagc=',i,pn,id_p(pn),flagc_p(pn)
        endif
        hvsigdt_p(i)=INF
      enddo
#endif
! *** update position and velocity ***	
      if(flagcont.eq.2) then
#if defined(GAS) || defined(STAR)
        do i=0,np-1
#if defined(COSM) && defined(BOXSIM)
          asc_p(i)=SI_a
          tc_p(i)=TM_tot
          ascg_p(i)=SI_a
          tcg_p(i)=TM_tot
          vnx_p(i)=SI_a*(vnx_p(i)-SI_hub*xc_p(i))
          vny_p(i)=SI_a*(vny_p(i)-SI_hub*yc_p(i))
          vnz_p(i)=SI_a*(vnz_p(i)-SI_hub*zc_p(i))
          vvnx_p(i)=vnx_p(i)/SI_a
          vvny_p(i)=vny_p(i)/SI_a
          vvnz_p(i)=vnz_p(i)/SI_a
#else
          vvnx_p(i)=vnx_p(i)
          vvny_p(i)=vny_p(i)
          vvnz_p(i)=vnz_p(i)
#endif
          x_p(i)=xc_p(i)
          y_p(i)=yc_p(i)
          z_p(i)=zc_p(i)
          flagt_p(i)=0
          dtr_p(i)=INF
          dt_p(i)=INF
          dtg_p(i)=INF
          ltg_p(i)=INF
          dtmnb_p(i)=INF
          omgh_p(i)=1.0d0
          zetah_p(i)=0.0d0
          ax_p(i)=0.0d0
          ay_p(i)=0.0d0
          az_p(i)=0.0d0
        enddo
#endif
#ifdef DM
        do i=0,ndm-1
#if defined(COSM) && defined(BOXSIM)
          asc_dm(i)=SI_a
          tc_dm(i)=TM_tot
          vnx_dm(i)=SI_a*(vnx_dm(i)-SI_hub*xc_dm(i))
          vny_dm(i)=SI_a*(vny_dm(i)-SI_hub*yc_dm(i))
          vnz_dm(i)=SI_a*(vnz_dm(i)-SI_hub*zc_dm(i))
#endif
          x_dm(i)=xc_dm(i)
          y_dm(i)=yc_dm(i)
          z_dm(i)=zc_dm(i)
          flagt_dm(i)=0
          dtr_dm(i)=INF
          dtmnb_dm(i)=INF
          dt_dm(i)=INF
          omgh_dm(i)=1.0d0
          zetah_dm(i)=0.0d0
        enddo
#endif
#ifdef STAR
! *** for stars ***
        do i=ng,ng+ns-1
          pn=list_ap(i)
          p_p(pn)=(GAM-1.0d0)*rho_p(pn)*u_p(pn)
#ifdef COOL
          ram_p(pn)=0.0d0
#endif
        enddo
#endif
      endif

#ifdef MHD
! for MHD
        do i=0,np-1
          pbx_p(i)=bx_p(i)
          pby_p(i)=by_p(i)
          pbz_p(i)=bz_p(i)
          dbx_p(i)=0.0d0
          dby_p(i)=0.0d0
          dbz_p(i)=0.0d0
          pdbx_p(i)=0.0d0
          pdby_p(i)=0.0d0
          pdbz_p(i)=0.0d0
#ifdef DCLEAN
          ppsi_p(i)=psi_p(i)
          dpsi_p(i)=0.0d0
          pdpsi_p(i)=0.0d0
#endif
        enddo
#endif

! *** setting the total yields ***
#if defined(SF_EFD) || defined(SF_ZFD)
      do i=0,np-1
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
#endif
! *** get total ***
#if defined(SF_EFD) || defined(SF_ZFD)
      call setyield()

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

      dtmp(0)=TMsf_dt
      dtmp(1)=TM_dt
      dtmp(2)=TM_tot
      TMsf_dt = 0.0d0
      TM_dt=0.0d0
      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)
! for both gas and stars 
      do i = 0,ng+ns-1
        pn=i
! *** re-set flagrfd_p to be 0
        if(flagfd_p(pn).eq.SI_snii.or. &
          (flagfd_p(pn).gt.SI_snii.and.flagrfd_p(pn).eq.1)) then
          flagrfd_p(pn)=0
        endif
! * update tm*_p() *
! * Unit MUSM *
        tnsn_p(pn)=tnsn_sn(pn)
        tmej_p(pn)=tmej_sn(pn)
! * Unit Solar Mass *
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
! check feedback energy
      TMsf_dt=DTSF
      TM_dt=TMsf_dt
      TM_tot=TM_tot-TMsf_dt
! *** dtyield at T-DTSF, to keep previous time step t*_p
! *** for new feedback particles
      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)

      allocate(mfd(0:np))

      do i=0,ng-1
        pn=list_ap(i)
#ifdef SF_EFD
! *** set Gsn_p for time step for feedback particles
        if(flagfd_p(pn).ne.0) then
! *** number of supernovae ***
          mfd(pn)=dabs(tnsn_sn(pn)-tnsn_p(pn))
          if(mfd(pn).lt.0.0d0 &
!         .and.dabs(mfd(pn)).lt.MININF*tnsn_p(pn)) then
           ) then
            mfd(pn)=0.0d0
          endif
          if(nsw_sn(pn).lt.0.0d0) then
            nsw_sn(pn)=0.0d0
          endif
          Gsn_p(pn)=mfd(pn)*SNEU+SWEU*nsw_sn(pn)
#ifdef CONSTFDE
! *** erg/s unit 5.752e43
          Gsn_p(pn)=(1.739e-6)*TMsf_dt
#endif
          Gsn_p(pn)=Gsn_p(pn)/m_p(pn)
        endif
#endif
! *** new feedback particles ***
        if(flagfd_p(pn).gt.0.and.i.lt.nag) then
          flagfd_p(pn)=-flagfd_p(pn)
! * keep tm*_p() at previous timestep *
! * Unit MUSM *
          tnsn_p(pn)=tnsn_sn(pn)
          tmej_p(pn)=tmej_sn(pn)
! * Unit Solar Mass *
          tmzHe_p(pn)=tmzHe_sn(pn)
          tmzC_p(pn)=tmzC_sn(pn)
          tmzN_p(pn)=tmzN_sn(pn)
          tmzO_p(pn)=tmzO_sn(pn)
          tmzNe_p(pn)=tmzNe_sn(pn)
          tmzMg_p(pn)=tmzMg_sn(pn)
          tmzSi_p(pn)=tmzSi_sn(pn)
          tmzFe_p(pn)=tmzFe_sn(pn)
          tmzZ_p(pn)=tmzZ_sn(pn)
        endif
      enddo
! TMsf_dt is used at *time.F
!        TMsf_dt=dtmp(0) 
      TM_dt=dtmp(1)
      TM_tot=dtmp(2)

      deallocate(mfd)
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

#endif

      if(myrank.eq.0) then
        write(6,*) ' output dt, next output t = ',OF_odt,OF_todt
      endif

!      if(ng.gt.0) then
!      write(fileo,'(a4,i3.3,a4)') 'cont',myrank,'.dat'
!      open(60,file=fileo,status='unknown')
!      do i=0,ng-1
!        pn=list_ap(i)
!        write(60,'(10(1pE13.5),4I10,4(1pE13.5))')
!     &    x_p(pn),y_p(pn),z_p(pn)
!     &   ,xc_p(pn),yc_p(pn),zc_p(pn),m_p(pn),rho_p(pn)
!     &   ,h_p(pn),alpv_p(pn),id_p(pn),flagc_p(pn),list_ap(i)
!     &   ,flagfd_p(pn),u_p(pn),ram_p(pn),myu_p(pn),mzZ_p(pn)
!      enddo
!      close(60)
!      endif

!      if(ns.gt.0) then
!      write(fileo,'(a5,i3.3,a4)') 'conts',myrank,'.dat'
!      open(60,file=fileo,status='unknown')
!      do i=ng,ng+ns-1
!        pn=list_ap(i)
!        write(60,'(10(1pE13.5),3I10,4(1pE13.5))') 
!     &    x_p(pn),y_p(pn),z_p(pn)
!     &   ,xc_p(pn),yc_p(pn),zc_p(pn),m_p(pn),rho_p(pn)
!     &   ,h_p(pn),alpv_p(pn),id_p(pn),flagc_p(pn),list_ap(i)
!     &   ,ts_p(pn),omgh_p(pn),zetah_p(pn),div_v_p(pn)
!      enddo
!      close(60)
!      endif

!      write(fileo,'(a6,i3.3,a4)') 'dmcont',myrank,'.dat'
!      open(60,file=fileo,status='unknown')
!      do i=0,ndm-1
!        write(60,'(13(1pE13.5),3I10,3(1pE13.5))') 
!     &    xc_dm(i),yc_dm(i),zc_dm(i)
!     &   ,vnx_dm(i),vny_dm(i),vnz_dm(i),m_p(i),lt_dm(i)
!     &   ,dt_dm(i),eps_dm(i),xc_dm(i),yc_dm(i),zc_dm(i)
!     &   ,i,id_dm(i),list_adm(i),dvx_dm(i),dvy_dm(i),dvz_dm(i)
!      enddo 
!      close(60)
!      stop

!      write(fileo,'(a7,i3.3,a4)') 'syscont',myrank,'.dat'
!      open(60,file=fileo,status='unknown')
!      write(60,'(13(1pE13.5),16I8)') OF_odt,OF_todt,TM_tot,TM_sdt
!     & ,TM_lt,TM_ldtmin,SI_a,SI_omg0,SI_lam0,SI_h0,SI_esng
!     & ,SI_zsns,SI_zsng
!     & ,step,totstep,lstep,SI_flagout,TM_k,TM_kn,SI_flagc,SI_flags
!     & ,SI_flagfd,SI_flagrad,SI_flagbh,SI_nof
!     & ,npt,ndmt,SI_ndm1t,SI_flagcosm
!      close(60) 

!      if(ns.gt.0) then
!      write(fileo,'(a5,i3.3)') 'contp',myrank
!      open(60,file=fileo,status='unknown')
!      do i=ng,ng+ns-1
!        pn=list_ap(i)
!        write(60,'(12(1pE13.5),2I7)') x_p(pn),y_p(pn),z_p(pn)
!     &   ,xc_p(pn),yc_p(pn),zc_p(pn),m_p(pn),rho_p(pn)
!     &   ,h_p(pn),dt_p(i),lt_p(i),dtr_p(i) 
!     &   ,flagt_p(i),list_ap(i)
!      enddo
!      close(60)
!      endif
!      stop

! *** random sequence cannot be preserved anyway ***
      idum = -111
#ifdef SGRAV
      if(myrank.eq.0) then
        write(6,*) ' softening limit = mp^(1/3) x',SI_eps
      endif
#ifdef COSM
      if(myrank.eq.0) then
        write(6,*) ' SI_epsa=',SI_epsa
      endif
      if(SI_a.lt.SI_epsa) then
        SI_eps=(1.0d0-SI_epsa+SI_a)*SI_eps
      endif
#endif
#else
      SI_eps=0.0d0
#endif

#ifdef FIXEDP
      open(50,file='ini/fixed.dat',status='old')      
      read(50,*) SI_fid(0),SI_fid(1)
      if(myrank.eq.0) then
        write(6,*) ' Fixed particles id between ',SI_fid(0),SI_fid(1)
      endif
      if(SI_fid(0).gt.SI_fid(1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1.'
        endif
        call MPI_FINALIZE()
        stop
      endif
! *** change nag and nagravg
      nag=0
      nc=0
      do i=0,ng-1
        pn=list_ap(i)
        if(id_p(pn).ge.SI_fid(0).and.id_p(pn).le.SI_fid(1)) then
          nalist(nc)=pn
          nc=nc+1
! *** set dvx and etc.
          dvx_p(pn)=0.0d0
          dvy_p(pn)=0.0d0
          dvz_p(pn)=0.0d0
          ndu_p(pn)=0.0d0
          pdu_p(pn)=0.0d0
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
#endif

#ifdef FIXEDDMP
      open(50,file='ini/fixeddm.dat',status='old')      
      read(50,*) SI_fiddm(0),SI_fiddm(1)
      if(myrank.eq.0) then
        write(6,*) ' Fixed DM particles id between ',SI_fiddm(0),SI_fiddm(1)
      endif
      if(SI_fiddm(0).gt.SI_fiddm(1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: DM id0 should be bigger than id1.'
        endif
        call MPI_FINALIZE()
        stop
      endif
! *** change nag and nagravg
      nadm=0
      nc=0
      do i=0,ndm-1
        if(id_dm(i).ge.SI_fiddm(0).and.id_dm(i).le.SI_fiddm(1)) then
          nalist(nc)=i
          nc=nc+1
! *** set dvx and etc.
          dvx_dm(i)=0.0d0
          dvy_dm(i)=0.0d0
          dvz_dm(i)=0.0d0
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

#ifdef RANAZ
      open(50,file='ini/ranaz.dat',status='old')      
      read(50,*) SI_ranazi(0),SI_ranazi(1),SI_ranazi(2)
      read(50,*) SI_ranazt
      if(myrank.eq.0) then
        if(SI_ranazi(0).eq.0) then
          write(6,*) ' RANAZ baryon particles id between ' &
          ,SI_ranazi(1),SI_ranazi(2)
       else
          write(6,*) ' RANAZ DM particles id between ' &
          ,SI_ranazi(1),SI_ranazi(2)
        endif
        write(6,*) ' until t (Gyr)=',SI_ranazt
      endif
      SI_ranazt=SI_ranazt/TMUGYR
      if(SI_ranazi(1).gt.SI_ranazi(2)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1 RANAZ.'
        endif
        call MPI_FINALIZE()
        stop
      endif
#endif

#ifdef ARTFDM
      open(50,file='ini/artfdm.dat',status='old')      
      read(50,*) SI_artfdm(0),SI_artfdm(1)
      if(myrank.eq.0) then
        write(6,*) ' ARTF DM particles id between ',SI_artfdm(0),SI_artfdm(1)
      endif
      if(SI_artfdm(0).gt.SI_artfdm(1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1 for DM.'
        endif
        call MPI_FINALIZE()
        stop
      endif
#endif

#endif

      deallocate(istatus)
      deallocate(npjr)
      deallocate(ireqs)
      deallocate(ireqr)

      return
      end
