#include "gcdp.def"
! ***************************************
!      set_ini.F95 for GCD+ ver. f03.0
!  04 Feb., 2013    produced by D.KAWATA
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

      integer,intent(out) :: npt,ndmt,np,ng,ndm,ns,nag,nadm,nas,nagravg &
       ,step,totstep,lstep,flagcont
      integer MNL
      parameter (MNL=3)
      character fin*20,tcha*1,fileo*60,ctmp*60,mname*20
      integer i,j,il,nc,pn,itmp
      integer flagrf
      integer nif,idf,ndmf,ngf,nsf,ngt,nst,ngi
      integer flagfdm,flagfg,flagfs
      double precision gamf,gkmm
! *** nd0, nd1 for number of particles for level 0 and 1 ***
      integer nd(0:MNL),ndl,ndt(0:MNL)
! * om0 (km s^-1 kpc^-1) *      
      double precision Ini_OMGg,dxi(0:MNL)
! * (10^12 M@), (100 kpc),(1.0e4 K),(4.71e8 yr) * /     
      integer flagtime
      double precision Ini_T0,Ini_ze,Ini_zi
! *** to read recfast.out **
      double precision zrec(2),xrec(4,2),xe(4),hnf,henf
      double precision Crec(6)
      integer Ini_nd
! *** physical value ***
      double precision rhob,rhob0,hubi,iai,rhoba
! * Mass of SPH, DM, total *      
      double precision m,mdm
      double precision kmmp
      double precision hmin,hmint,etahmax
! *** Initial parameter ***
! *** Ini_hyd: 0: N-body only, 1: include gas in high reso region ***
! *** Ini_mul: 0: only level0, 1: level 0 and 1 ***
      integer Ini_hyd,Ini_mul
! work
      integer jsta,jend,ierr
      double precision,allocatable :: tx(:),ty(:),tz(:),adv(:)
#if defined(SF_EFD) || defined(SF_ZFD)
! work
      double precision,allocatable :: mfd(:)
! * for feedback  common *
! * Total *
      double precision,allocatable :: tmej_sn(:),tnsn_sn(:),tmzHe_sn(:) &
       ,tmzC_sn(:),tmzN_sn(:),tmzO_sn(:),tmzNe_sn(:),tmzMg_sn(:) &
       ,tmzSi_sn(:),tmzFe_sn(:),tmzZ_sn(:),nsw_sn(:)
#endif


#ifdef COSM
      if(myrank.eq.0) then
        write(6,*) ' H0,HUB0 = ',H0_1,HUB0
      endif
#endif

#ifdef COOL
      SI_flagrad=0
#else
      SI_flagrad=1
#endif
!*****   Open Initial Data File ./ini/input.dat
      open(50,file='ini/input.dat',status='old')      
! *** Read Initial Parameter ***
! * Number of Particle etc. *  
      read(50,*) Ini_mul
! * Interval of Output Data, Distribution *
      read(50,*) Ini_nd
      read(50,*) SI_nof
! * Initial Values unit is our Unit *
#ifdef COSM
      read(50,*) Ini_OMGg
      read(50,*) Ini_ze
#else
      read(50,*) Ini_zi
      read(50,*) flagtime
      read(50,*) TM_tot,TM_t0
      read(50,*) flagfg,flagfdm,flagfs
      Ini_ze=Ini_zi
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

      close(50)
! *** initialization for number of particles ***   
      npt=0
      ndmt=0
      nst=0
      ngt=0
      np=0
      ng=0
      ndm=0
      ns=0
#ifndef MULTI_LEV
      if(Ini_mul.gt.0) then
        if(myrank.eq.0) then
          write(6,*) ' Error: MULTI_LEV is off in gcdp.def, but Ini_mul=' &
           ,Ini_mul
        endif
        stop
      endif
#endif

#ifdef COSM
#ifdef STAR
      SI_tnsf=0
#endif
! *** cosmological simulation ***
      SI_omgb=Ini_OMGg
! *** Check Input Parameter ***
      if(Ini_OMGg.lt.0.0d0.or.Ini_OMGg.gt.1.0d0) then
        if(myrank.eq.0) then
          write(6,*) ' Error in set_ini() : OMGg = ',Ini_OMGg
        endif
        stop
      endif
! *** reading initial condision for cosmological simulation
! *** Open Initial Data File for level 1 ***
      ndmt = 0
      ndm = 0
#ifdef DM
! get ndm for each rank
      do il=Ini_mul,0,-1
        write(fin,'(a9,i1,a5)') 'ini/level',il,'d.dat'
        open(51,file=fin,status='old',form='unformatted')
        read(51) nif,idf,ndl,nd(il),flagrf
! *** read the data ***
        call para_range(0,ndl-1,nprocs,myrank,jsta,jend)
        ndm=ndm+jend-jsta+1
        close(51)
      enddo

! allocate dm variables
      call allocate_dm_all(ndm)      

      allocate(tx(0:ndm))
      allocate(ty(0:ndm))
      allocate(tz(0:ndm))

      ndm=0
      do il=Ini_mul,0,-1
        write(fin,'(a9,i1,a5)') 'ini/level',il,'d.dat'
        open(51,file=fin,status='old',form='unformatted')
        read(51) nif,idf,ndl,nd(il),flagrf
! *** flagrf: 0: read only position and velocity, otherwise read also mass
! *** dx (Mpc, Comoving), zstart,omegam,lambda ***
        read(51) dxi(il),SI_a,SI_omg0,SI_lam0,SI_h0
        dxi(il) = dxi(il)*10.0d0
        if(myrank.eq.0) then
          write(6,*) il,'level ndm=',ndl
        endif
!        call MPI_BCAST(ndl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! *** read the data ***
        call para_range(0,ndl-1,nprocs,myrank,jsta,jend)
!       call MPI_BCAST(flagrf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        ndt(il)=0
        nd(il)=0
        do i=0,ndl-1
          if(i.ge.jsta.and.i.le.jend) then
            if(flagrf.eq.0) then
              read(51) x_dm(ndm),y_dm(ndm),z_dm(ndm) &
               ,tx(ndm),ty(ndm),tz(ndm)
            else
              read(51) x_dm(ndm),y_dm(ndm),z_dm(ndm) &
               ,tx(ndm),ty(ndm),tz(ndm),m_dm(ndm)
            endif
            id_dm(ndm)=ndmt
            ndm = ndm+1
            nd(il)=nd(il)+1
          else
            read(51)
          endif
          ndmt = ndmt+1
          ndt(il)=ndt(il)+1
        enddo
#ifdef BOXSIM
        if(il.eq.0) then
          read(51) SI_lbox0(0),SI_lbox0(1),SI_lbox0(2)
! Mpc -> 100 kpc
          SI_lbox0(0)=SI_lbox0(0)*10.0d0
          SI_lbox0(1)=SI_lbox0(1)*10.0d0
          SI_lbox0(2)=SI_lbox0(2)*10.0d0
! set box size at the initial redshift
          SI_lbox(0)=SI_a*SI_lbox0(0)
          SI_lbox(1)=SI_a*SI_lbox0(1)
          SI_lbox(2)=SI_a*SI_lbox0(2)
          if(myrank.eq.0) then
            write(6,*) ' box size (comoving)=',SI_lbox0(0),SI_lbox0(1) &
              ,SI_lbox0(2)
            write(6,*) ' box size (initial) =',SI_lbox(0),SI_lbox(1) &
              ,SI_lbox(2)
          endif
        endif
#endif
        close(51)

      enddo
#endif
  151 format(5(1pE13.5))
  152 format(6(1pE13.5))
  153 format(7(1pE13.5))
! *** for Gas particles
! *** Open Initial Data File for level 1 ***
      npt = 0
      np = 0
#ifdef GAS
      il=Ini_mul
      write(fin,'(a9,i1,a5)') 'ini/level',il,'p.dat'
      open(51,file=fin,status='old',form='unformatted')
      read(51) nif,idf,npt,itmp,flagrf
! *** dx (Mpc, Comoving), zstart,omegam,lambda ***
      read(51) dxi(il),SI_a,SI_omg0,SI_lam0,SI_h0
      dxi(il) = dxi(il)*10.0d0

#ifdef MHD
      open(57,file='ini/gas-mhd.dat',status='old',form='unformatted')
      if(myrank.eq.0) then
         write(6,*) 'reading ini/gas-mhd.dat'
      endif
#endif

! *** read the data ***
      call para_range(0,npt-1,nprocs,myrank,jsta,jend)

      allocate(adv(0:jend-jsta+1))
      np=jend-jsta+1
      call allocate_baryon_all(np)
#ifdef MHD
      call allocate_mhd(np)
#endif
      np=0

! *** Read and  Set Position ***
      do i=0,npt-1
! * x, y, z, a(t)b(t)p(q) (Mpc) *		           
        if(i.ge.jsta.and.i.le.jend) then
          if(flagrf.eq.0) then
            read(51) xc_p(np),yc_p(np),zc_p(np) &
             ,vnx_p(np),vny_p(np),vnz_p(np),adv(np)
          else
            read(51) xc_p(np),yc_p(np),zc_p(np) &
             ,vnx_p(np),vny_p(np),vnz_p(np),adv(np),m_p(np)
          endif
#if defined(SF_EFD) || defined(SF_ZFD)
          flagfd_p(np)=0
          flagrfd_p(np)=0
          ts_p(np)=0.0d0
#endif

#ifdef MHD
          read(57) bx_p(np),by_p(np),bz_p(np)
          pbx_p(np)=bx_p(np)
          pby_p(np)=by_p(np)
          pbz_p(np)=bz_p(np)
          dbx_p(np)=0.0d0
          dby_p(np)=0.0d0
          dbz_p(np)=0.0d0
          pdbx_p(np)=0.0d0
          pdby_p(np)=0.0d0
          pdbz_p(np)=0.0d0
          divb_p(np)=0.0d0
          arotb_p(np)=0.0d0
          alpb_p(np)=ALPB_MIN
#endif

          id_p(np)=i
          np = np+1
        else
          read(51)
#ifdef MHD
          read(57)
#endif
        endif
      enddo
      close(51)
#ifdef MHD   
      close(57)
#endif
!        write(6,*) myrank,' rank ',il,' level np = ',np
#endif
      SI_h0 = SI_h0/100.0d0
#ifdef GAS
      Ini_OMGg = Ini_OMGg
#else
      Ini_OMGg = 0.0d0
#endif
!*********************************
!  Make Time tabele and Set TM_t0
! ******************************** 
!      write(6,*) ' Ini_ze = ',Ini_ze
      call makettable(1.0d0/(1.0d0+Ini_ze))
! *****   Set Initial Parameter   *****
! *** Initial Time, Set TM_tot ***
      call tfroma(SI_a,TM_tot)
      Ini_zi = 1.0d0/SI_a-1.0d0
      if(myrank.eq.0) then
        write(6,*) ' cosmology: Omg0,lam0,h,omb=',SI_omg0,SI_lam0,SI_h0,SI_omgb
      endif
#ifdef GAS
      open(50,file='./ini/recfast.out',status='old')
      read(50,*)
      read(50,'(a1,6(1pE13.5))') tcha,Crec(1),Crec(2),Crec(3) &
       ,Crec(4),Crec(5),Crec(6)
      read(50,*)
      if(myrank.eq.0) then
        write(6,*) ' Ob, Odm, lam, H0, T0, Yp'
        write(6,*) Crec(1),Crec(2),Crec(3),Crec(4),Crec(5),Crec(6)      
      endif
! *** check ***
      if(dabs(Crec(1)-Ini_OMGg).gt.0.001d0) then
        if(myrank.eq.0) then
          write(6,*) ' Error in reading recfast.out: Omgb'
          write(6,*) ' input OMGg and in recfast=',Ini_OMGg,Crec(1)
        endif
        stop
      endif
      if(dabs((Crec(2)-(SI_omg0-Ini_OMGg))).gt.0.001d0) then
        if(myrank.eq.0) then
        write(6,*) ' Error in reading recfast.out: Omgdm'
        write(6,*) ' Omgdm (sim)',SI_omg0-Ini_OMGg
        endif
        stop
      endif
      if(dabs(Crec(3)-SI_lam0).gt.0.001d0) then
        if(myrank.eq.0) then
        write(6,*) ' Error in reading recfast.out: lam0'
        endif
        stop
      endif
      if(dabs(Crec(4)-100.0d0*SI_h0).gt.0.001d0) then
        if(myrank.eq.0) then
        write(6,*) ' Error in reading recfast.out: h0'
        endif
        stop
      endif
      zrec(1) = 10000.0d0
      do i=1,4
        xrec(i,1)=0.0d0
      enddo
      do i=1,2000
        do j=1,4
          xrec(j,2)=xrec(j,1)
        enddo  
        zrec(2) = zrec(1)
! *** read z xe xh xhe Tmat ***
        read(50,'(5(1pE13.5))') zrec(1),xrec(1,1),xrec(2,1),xrec(3,1),xrec(4,1)
        if(zrec(1).lt.Ini_zi) then
          goto 71
        endif
      enddo
      if(myrank.eq.0) then
        write(6,*) ' Error in reading recfast: z(sim)=',Ini_zi 
      endif
   71 if(myrank.eq.0) then
        write(6,*) ' zregrid(1,2), zi=',zrec(1),zrec(2),Ini_zi
      endif
      do i=1,4
        xe(i)=xrec(i,2)+(xrec(i,1)-xrec(i,2))*(Ini_zi-zrec(2)) &
         /(zrec(1)-zrec(2))
      enddo
! *** set initial temperature ***
      Ini_T0 = xe(4)/TUK
      if(myrank.eq.0) then
        write(6,*) ' Temperature   = ',Ini_T0,' (10^4 K)'
        write(6,*) ' xe, xH, xHe = ',xe(1),xe(2),xe(3)
      endif
#endif
! *** Set Coordinate and Data of particle **
      nag=np
      ng=np
      ns=0	
! *** Set Hubble constant ***
      iai = 1.0d0/SI_a
      hubi = SI_h0*HUB0*dsqrt(SI_omg0*iai*iai*iai &
       -(SI_omg0+SI_lam0-1.0)*iai*iai+SI_lam0)
#ifdef BOXSIM
      SI_hub=hubi
#endif
! *** Calculate Initial Value ***
      rhob0=3.0d0*HUB0*SI_h0*HUB0*SI_h0/(8.0d0*M_PI*G)
      rhoba=rhob0*iai*iai*iai
! *** Calculate Initial Baryon Density ***
      rhob=Ini_OMGg*rhob0*iai*iai*iai
! *** Lambda ***
      SI_lam0 = SI_lam0*HUB0*HUB0*SI_h0*SI_h0

! **********************************************************
!    Set Coordinate and Velocity and Density of particle
! **********************************************************
! *** Invert to this program Unit ***
#ifdef DM
      do i=0,ndm-1
        x_dm(i)=x_dm(i)*(SI_a*10.0d0)
        y_dm(i)=y_dm(i)*(SI_a*10.0d0)
        z_dm(i)=z_dm(i)*(SI_a*10.0d0)
#ifdef BOXSIM
        asc_dm(i)=SI_a
        tc_dm(i)=TM_tot
#endif
      enddo
! *** mass of DM for level 1 particles ***
      nc=0
      do il=Ini_mul,0,-1
        if(il.eq.Ini_mul) then
          mdm = (SI_omg0-Ini_OMGg)*rhob0*(dxi(il)**3)
        else
          mdm = SI_omg0*rhob0*(dxi(il)**3)
        endif
! *****   Set Initial Mass for DM *****
        do i=nc,nc+nd(il)-1
          if(flagrf.eq.0) then
            m_dm(i)=mdm
          endif
! assume density is uniform because tough to get the info
          if(il.eq.Ini_mul) then
            rho_dm(i)=(SI_omg0-Ini_OMGg)*rhoba
          else
            rho_dm(i)=SI_omg0*rhoba
          endif 
        enddo        
        nc = nc+nd(il)
        if(myrank.eq.0) then

          if(flagrf.eq.0) then
            write(6,*) ' mdm (level ',il,') = ',mdm
          else
            write(6,*) ' mdm (level ',il,') = ',m_dm(nc-1)
          endif
        endif
      enddo
! *** Set Initial Velocity ***
! *** Local velocity + Hubble Flow ***
! *** for DM Particles ***
      do i=0,ndm-1      
#ifdef BOXSIM
! v = a^2 vpec, tx = a x vpec
        vnx_dm(i)=(SI_a)*tx(i)/VUKMS
        vny_dm(i)=(SI_a)*ty(i)/VUKMS
        vnz_dm(i)=(SI_a)*tz(i)/VUKMS
#else
        vnx_dm(i)=hubi*x_dm(i)+tx(i)/VUKMS
        vny_dm(i)=hubi*y_dm(i)+ty(i)/VUKMS
        vnz_dm(i)=hubi*z_dm(i)+tz(i)/VUKMS
#endif
      enddo

      deallocate(tx)
      deallocate(ty)
      deallocate(tz)

#endif
#ifdef GAS
! *** mass of Gas ***
      m = Ini_OMGg*rhob0*(dxi(Ini_mul)**3)
! * Calculate physical Value *
      kmmp=(KCGS/(MP*MYU))/K_MU

      if(myrank.eq.0) then      
        write(6,*) ' Rho_b,b0   =',rhob,rhob0
        write(6,*) ' Hub_i = ',hubi
        if(flagrf.eq.0) then
          write(6,*) ' mgas       =',m
        else
          write(6,*) ' mgas       =',m_p(0)
        endif
        write(6,*) ' dxi(0)     =',dxi(0)*SI_a
        write(6,*) ' a (zi)     = ',SI_a
      endif

! *** for gas particle ***
      do i=0,np-1
! *** positions ***
        xc_p(i)=xc_p(i)*(SI_a*10.0d0)
        yc_p(i)=yc_p(i)*(SI_a*10.0d0)
        zc_p(i)=zc_p(i)*(SI_a*10.0d0)
! *** velocities ***
#ifdef BOXSIM
        vnx_p(i)=(SI_a)*vnx_p(i)/VUKMS
        vny_p(i)=(SI_a)*vny_p(i)/VUKMS
        vnz_p(i)=(SI_a)*vnz_p(i)/VUKMS
        asc_p(i)=SI_a
        tc_p(i)=TM_tot
        ascg_p(i)=SI_a
        tcg_p(i)=TM_tot
#else
        vnx_p(i)=hubi*xc_p(i)+vnx_p(i)/VUKMS
        vny_p(i)=hubi*yc_p(i)+vny_p(i)/VUKMS
        vnz_p(i)=hubi*zc_p(i)+vnz_p(i)/VUKMS
#endif
      enddo
! *****   Set Initial Mass and Pressure
!   and Sound Velocity  and Thermal Energy*****
      do i=0,np-1
        if(flagrf.eq.0) then
          m_p(i)=m
        endif
        rho_p(i)=rhob*(1.0d0+adv(i))
        p_p(i)=kmmp*rho_p(i)*Ini_T0
! just for test
!        p_p(i)=kmmp*rho_p(i)*(Ini_T0+10.0d0**(dble(i)*1.5e-3-1.0d0))
        cs_p(i)=dsqrt(GAM*p_p(i)/rho_p(i))
        u_p(i)=p_p(i)/((GAM-1.0d0)*rho_p(i))
        myu_p(i)=MYU
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**THIRD)
! * Initialize for star *
#ifdef SF_EFD
        Gsn_p(i) = 0.0d0
#endif
#ifdef STAR
        ts_p(i)=0.0d0
#endif
        list_ap(i)=i
      enddo

      deallocate(adv)

! *** Set Initial Metallicity ***
! *** H and He mass fraction w.r.t. w.r.t. the total number of p *
      hnf = 1.0d0-Crec(6)
      henf = Crec(6)
      do i=0,np-1
#ifdef METAL
        mzC_p(i)=0.0d0
        mzN_p(i)=0.0d0
        mzO_p(i)=0.0d0
        mzNe_p(i)=0.0d0
        mzMg_p(i)=0.0d0
        mzSi_p(i)=0.0d0
        mzFe_p(i)=0.0d0
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** total ***
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
#endif
      enddo
#endif
! ***** end of input for cosmological simulation *****
#else
      Ini_mul=-1
! ***** non-cosmological simulation *****
      SI_a = 1.0d0/(Ini_zi+1.0d0)
! ***  Time bin ***
      if(flagtime.ne.0) then
        TM_tot = TM_tot/TMUGYR
        TM_t0 = TM_t0/TMUGYR
      endif
! * Set Time Parameter *      
      TM_dt = 0.0d0
      Ini_OMGg=0.0d0
#ifdef DM
      if(flagfdm.eq.0) then
        if(myrank.eq.0) then
          write(6,*) ' reading dm.dat'
        endif
! *** Open DM data file ***

        open(51,file='ini/dm.dat',status='old',form='unformatted')
        read(51) nif,idf,ndmt,ndm
        if(myrank.eq.0) then
          write(6,*) ' ndmt,ndm=',ndmt,ndm
        endif
        call para_range(0,ndmt-1,nprocs,myrank,jsta,jend)

! allocate DM variables
        ndm=jend-jsta+1
        call allocate_dm_all(ndm)

        nc=0
        ndm=0
! *** Read and  Set Position ***
        do i=0,ndmt-1
          if(i.ge.jsta.and.i.le.jend) then
            read(51) x_dm(ndm),y_dm(ndm),z_dm(ndm) &
            ,vnx_dm(ndm),vny_dm(ndm),vnz_dm(ndm),m_dm(ndm) &
            ,rho_dm(ndm)
            id_dm(ndm)=nc
            ndm=ndm+1
          else
            read(51)
          endif
          nc=nc+1
        enddo
        read(51) mname
        if(myrank.eq.0) then
          write(6,*) mname
        endif
#ifdef BOXSIM
        read(51) SI_lbox(0),SI_lbox(1),SI_lbox(2)
        if(myrank.eq.0) then
          write(6,*) ' box size=',SI_lbox(0),SI_lbox(1),SI_lbox(2)
        endif
#endif
        close(51)
        if(myrank.eq.0) then
          write(6,*) ' mdm(0)=',m_dm(0)
        endif
#ifdef DM
      else if(myrank.eq.0) then
        write(6,*) '##### Warning in set_ini()#####'
        write(6,*) ' DM data not read, but DM is on in gcdp.def.'
        write(6,*) ' This likely causes a problem.'
#endif
      endif
#endif
! *** for Gas particles
! *** Open Initial Data File for level 1 ***
      ng=0
      ns=0
      np=0
#ifdef GAS
! get ng for each rank
      if(flagfg.le.0) then
        if(myrank.eq.0) then
          write(6,*) ' reading gas.dat. flagfg=',flagfg
        endif
        open(51,file='ini/gas.dat',status='old',form='unformatted')
        read(51) nif,idf,ngt,ng,flagrf
        close(51) 
        if(myrank.eq.0) then
          write(6,*) ' nif,idf,ngt,ng,flagrf=',nif,idf,ngt,ng,flagrf
        endif    
        call para_range(0,ngt-1,nprocs,myrank,jsta,jend)
        ng=jend-jsta+1
      endif
#endif
! get ns for each rank
#ifdef STAR
      SI_tnsf=0
      if(flagfs.eq.0) then
        open(51,file='ini/star.dat',status='old',form='unformatted')
        if(myrank.eq.0) then
          write(6,*) ' reading star.dat. flagfs=',flagfs
        endif
        read(51) nif,idf,nst,ns
        close(51)

        if(myrank.eq.0) then
          write(6,*) ' nif,idf,nst,ns=',nif,idf,nst,ns
        endif      
        call para_range(0,nst-1,nprocs,myrank,jsta,jend)
      endif
#endif
#if defined(GAS) || defined(STAR)
      np=ng+ns

! initial allocate memory space
! all baryon
      call allocate_baryon_all(np)
#ifdef MHD
      call allocate_mhd(np)
#endif

      ng=0
      ns=0
      np=0
#endif
#ifdef GAS
! *** Read and Set Position and basic values ***
      if(flagfg.le.0) then
        open(51,file='ini/gas.dat',status='old',form='unformatted')
!        if(myrank.eq.0) then
!          write(6,*) ' reading gas.dat. flagfg=',flagfg
!        endif
        read(51) nif,idf,ngt,ng,flagrf
! flagrf is not 0: read flagfd_p info from gas.dat
! since pv31.19
#ifndef STAR
        if(myrank.eq.0.and.flagrf.ne.0) then
          write(6,*) ' set flagrf=0'
        endif    
        flagrf=0
#endif    
        call para_range(0,ngt-1,nprocs,myrank,jsta,jend)
        nc=0
        ng=0
        if(flagrf.eq.0) then
          do i=0,ngt-1
            if(i.ge.jsta.and.i.le.jend) then
              read(51) xc_p(ng),yc_p(ng),zc_p(ng),vnx_p(ng),vny_p(ng) &
               ,vnz_p(ng),m_p(ng),rho_p(ng),u_p(ng)
              id_p(ng)=nc
#if defined(SF_EFD) || defined(SF_ZFD)
              flagfd_p(ng)=0
              flagrfd_p(ng)=0
              ts_p(ng)=0.0d0
#endif
              ng=ng+1
            else
              read(51)
            endif
            nc=nc+1
          enddo
        else
#ifdef STAR
          do i=0,ngt-1
            if(i.ge.jsta.and.i.le.jend) then
              read(51) xc_p(ng),yc_p(ng),zc_p(ng),vnx_p(ng),vny_p(ng) &
               ,vnz_p(ng),m_p(ng),rho_p(ng),u_p(ng),ts_p(ng),itmp
#if defined(SF_EFD) || defined(SF_ZFD)
              flagfd_p(ng)=itmp
              if(flagfd_p(ng).lt.0) then
                flagrfd_p(ng)=1
              else
                flagrfd_p(ng)=0
              endif
#endif
              id_p(ng)=nc
              ng=ng+1
            else
              read(51)
            endif
            nc=nc+1
          enddo
#else
            if(myrank.eq.0) then
              write(6,*) &
               ' Error in set_ini.F: STAR is off, but reading ts_p'
            endif
            stop
#endif
        endif

        read(51) mname
        if(myrank.eq.0) then
          write(6,*) mname
        endif
        read(51) gamf
        if(dabs(GAM-gamf)/GAM.gt.0.1d0) then
          if(myrank.eq.0) then
            write(6,*) 
            if(myrank.eq.0) then
              write(6,*) ' Error(): gamma is incosnsitent'
              write(6,*) ' ini/input.dat, ini/gas.dat',GAM,gamf
            endif
            call MPI_FINALIZE()
            stop
          endif
        endif    
#ifdef BOXSIM
        read(51) SI_lbox(0),SI_lbox(1),SI_lbox(2)
        if(myrank.eq.0) then
          write(6,*) ' box size=',SI_lbox(0),SI_lbox(1),SI_lbox(2)
        endif
#ifdef BCX_FIX
        read(51) BCX_ndx
#ifdef MHD
        read(51) BCX_rhox(0),BCX_px(0),BCX_vx(0),BCX_vy(0),BCX_vz(0)
        read(51) BCX_rhox(1),BCX_px(1),BCX_vx(1),BCX_vy(1),BCX_vz(1)
        read(51) BCX_bx(0),BCX_by(0),BCX_bz(0)
        read(51) BCX_bx(1),BCX_by(1),BCX_bz(1)
        BCX_ndx=5
#else
        read(51) BCX_rhox(0),BCX_px(0),BCX_vx(0)
        read(51) BCX_rhox(1),BCX_px(1),BCX_vx(1)
#endif
        if(myrank.eq.0) then
          write(6,*) ' boundary is fixed from the edge+-dx=',BCX_ndx
        endif
#endif
#endif
        close(51)
        if(myrank.eq.0) then
          write(6,*) ' mgas(0)=',m_p(0)
        endif
      else
        if(myrank.eq.0) then
          write(6,*) &
           ' Error: not reading gas.dat, but GAS is defined.'
          write(6,*) ' input flagfg=',flagfg
        endif
        stop
      endif
#ifdef METAL
! *** reading metal data ***
      if(flagfg.eq.-1.or.flagfg.le.-3) then
        open(51,file='ini/gas-metal.dat',status='old',form='unformatted')
        if(myrank.eq.0) then
          write(6,*) ' reading ini/gas-metal.dat'
        endif 
        ngi=0
        do i=0,ngt-1
          if(i.ge.jsta.and.i.le.jend) then
            read(51) mzHe_p(ngi),mzC_p(ngi),mzN_p(ngi),mzO_p(ngi) &
             ,mzNe_p(ngi),mzMg_p(ngi),mzSi_p(ngi),mzFe_p(ngi),mzZ_p(ngi)
            ngi=ngi+1
          else
            read(51)
          endif
        enddo
        close(51)
        if(ngi.ne.ng) then
          write(6,*) ' Error in reading gas-metal.dat at rank',myrank
          write(6,*) ' ng in gas.dat gas-metail.dat=',ng,ngi
          stop
        endif
      endif
#else
      if(flagfg.eq.-1.or.flagfg.le.-3) then
        if(myrank.eq.0) then
          write(6,*) ' Error: flagfg is inconsistent=',flagfg
          write(6,*) ' METAL is NOT defined, but trying to read gas-metal.dat'
        endif
        stop
      endif
#endif
      np=ng  
#ifdef METAL
      if(flagfg.ne.-1.and.flagfg.gt.-3) then
! *** assume primordial abundance ***
        do i=0,ng-1
          mzHe_p(i)=0.24d0*MUSM*m_p(i)
          mzC_p(i)=0.0d0
          mzN_p(i)=0.0d0
          mzO_p(i)=0.0d0
          mzNe_p(i)=0.0d0
          mzMg_p(i)=0.0d0
          mzSi_p(i)=0.0d0
          mzFe_p(i)=0.0d0
          mzZ_p(i)=0.0d0
        enddo
      endif
#endif

#ifdef MHD

! *** reading B data ***
      open(51,file='ini/gas-mhd.dat',status='old',form='unformatted')
      if(myrank.eq.0) then
         write(6,*) 'reading ini/gas-mhd.dat'
      endif
      ngi=0
      do i=0,ngt-1
         if(i.ge.jsta.and.i.le.jend) then
            read(51) bx_p(ngi),by_p(ngi),bz_p(ngi)
            pbx_p(ngi)=bx_p(ngi)
            pby_p(ngi)=by_p(ngi)
            pbz_p(ngi)=bz_p(ngi)
            dbx_p(ngi)=0.0d0
            dby_p(ngi)=0.0d0
            dbz_p(ngi)=0.0d0
            pdbx_p(ngi)=0.0d0
            pdby_p(ngi)=0.0d0
            pdbz_p(ngi)=0.0d0
            divb_p(ngi)=0.0d0
            arotb_p(ngi)=0.0d0
            alpb_p(ngi)=ALPB_MIN
            ngi=ngi+1
         else
            read(51)
         endif
      enddo
      close(51)
      if(ngi.ne.ng) then
         write(6,*) 'Error in reading gad-mhd.dat at rank',myrank
         write(6,*) 'ng in gat.dat gas-mhd.dat=',ng,ngi
         stop
      endif
#endif

! ***** end of reading files *****
! *** the other properties **
      do i=0,ng-1
#ifdef SF_EFD
        Gsn_p(i) = 0.0d0
#endif
        flagc_p(i) = 0
        p_p(i) = (GAM-1.0d0)*rho_p(i)*u_p(i)
        cs_p(i)=dsqrt(GAM*p_p(i)/rho_p(i))
        myu_p(i)=MYU
#ifdef STAR
        ts_p(i)=0.0d0
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! *** total ***
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
#endif
! *** particle list ***
        list_ap(i) = i
      enddo
#ifdef STAR
      SI_tnsf=0
      if(flagfs.eq.0) then
        open(51,file='ini/star.dat',status='old',form='unformatted')
#ifdef METAL
        open(52,file='ini/star-metal.dat',status='old',form='unformatted')
#endif
        if(myrank.eq.0) then
          write(6,*) ' reading star.dat. flagfs=',flagfs
        endif
!        read(51,'(a24,4I10)') ctmp,nif,idf,ngt,ng
        read(51) nif,idf,nst,ns
        if(myrank.eq.0) then
          write(6,*) ' nif,idf,nst,ns=',nif,idf,nst,ns
        endif      
        call para_range(0,nst-1,nprocs,myrank,jsta,jend)
! keep the total number of formed star 
        SI_tnsf=nst

        nc=ngt
        ns=0
        np=ng
        do i=0,nst-1
          if(i.ge.jsta.and.i.le.jend) then
            read(51) xc_p(np),yc_p(np),zc_p(np),vnx_p(np),vny_p(np) &
             ,vnz_p(np),m_p(np),m_p(np),rho_p(np),ts_p(np),itmp
#if defined(SF_EFD) || defined(SF_ZFD)
           flagfd_p(np)=itmp
#endif
! *** age -> formation time ***
            ts_p(np)=TM_tot-ts_p(np)            
#ifdef METAL
            read(52) mzHe_p(np),mzC_p(np),mzN_p(np),mzO_p(np),mzNe_p(np) &
             ,mzMg_p(np),mzSi_p(np),mzFe_p(np),mzZ_p(np)
#endif
            id_p(np)=nc
            ns=ns+1
            np=np+1
          else
            read(51)
#ifdef METAL
            read(52)
#endif
          endif
          nc=nc+1
        enddo
        read(51) mname
        close(51)
        if(myrank.eq.0) then
          write(6,*) ' mstar(0)=',m_p(ng)
        endif
#ifdef METAL
        close(52)
#endif
        do i=ng,np-1
! *** the other properties **
#ifdef SF_EFD
          Gsn_p(i) = 0.0d0
#endif
          flagc_p(i) = 1
          flagt_p(i) = 0
! *** set zero for unused data ***
          dtr_p(i)=INF
          p_p(i) = 0.0d0
          cs_p(i) = 0.0d0
          u_p(i) = 0.0d0
          myu_p(i)=0.0d0
! *** particle list ***
          list_ap(i) = i
#if defined(SF_EFD) || defined(SF_ZFD)
          flagrfd_p(i)=0
#endif
        enddo
      endif

#ifdef METAL
! initialise all the metal 
      do i=0,np-1
        mzHe0_p(i)=mzHe_p(i)
        mzC0_p(i)=mzC_p(i)
        mzN0_p(i)=mzN_p(i)
        mzO0_p(i)=mzO_p(i)
        mzNe0_p(i)=mzNe_p(i)
        mzMg0_p(i)=mzMg_p(i)
        mzSi0_p(i)=mzSi_p(i)
        mzFe0_p(i)=mzFe_p(i)
        mzZ0_p(i)=mzZ_p(i)
      enddo
#endif

#if defined(SF_EFD) || defined(SF_ZFD)
! *** get yield even if there is no star, but for feedback gas ***
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

      TM_dt=0.0d0
      TMsf_dt=0.0d0

      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)
! *** set tmz*_p for both star and gas ***
      do i = 0,ng+ns-1
        if((flagfd_p(i).ge.1.and.flagfd_p(pn).le.SI_snii) &
         .or.(flagfd_p(i).gt.SI_snii.and.flagrfd_p(i).eq.1)) then
          flagrfd_p(i)=0
! since pv32.9
!          flagfd_p(i)=-flagfd_p(i)
        endif
! * update tm*_p() *
! * Unit MUSM *
        tmej_p(i)=tmej_p(i)
        tnsn_p(i)=tnsn_sn(i)
! * Unit Solar Mass *
        tmzHe_p(i)=tmzHe_sn(i)
        tmzC_p(i)=tmzC_sn(i)
        tmzN_p(i)=tmzN_sn(i)
        tmzO_p(i)=tmzO_sn(i)
        tmzNe_p(i)=tmzNe_sn(i)
        tmzMg_p(i)=tmzMg_sn(i)
        tmzSi_p(i)=tmzSi_sn(i)
        tmzFe_p(i)=tmzFe_sn(i)
        tmzZ_p(i)=tmzZ_sn(i)
! * check mass *
      enddo
! check feedback energy
#ifdef SF_EFD
      TMsf_dt=DTSF
      TM_dt=TMsf_dt
      call dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)

      allocate(mfd(0:np))

      do i=0,ng+ns-1
        pn=i
        if(flagfd_p(pn).lt.0) then
! *** number of supernovae ***
          mfd(pn)=tnsn_sn(pn)-tnsn_p(pn)
          if(mfd(pn).lt.0.0d0.and.dabs(mfd(pn)).lt.MININF*tnsn_p(pn)) then
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
      enddo

      deallocate(mfd)

#endif

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
      TM_dt=0.0d0
! *** end for STAR
#endif
! *** end for GAS
#endif
      npt=ngt+nst
! *** end of non-cosmological simulation input ***
#endif  
! *** end of data setting ***

! *** Set the data in common for cosmological and non-cosmological mode ***
   92 nag=ng
      nadm=ndm
      nas=ns
      nagravg=ng
      np=ng+ns
      step=0
      lstep=0
      totstep = 0
! ***  Time bin ***
! * Set Time Parameter. This is used at ddecb() for setting tm*_p *
      TM_dt=0.0d0 
! ***   about Output   ***
      OF_odt = (TM_t0-TM_tot)/dble(Ini_nd)
      OF_todt = OF_odt+TM_tot

      if(myrank.eq.0) then
        write(6,*) '#####   Initial Condition Information   #####'
#ifdef COSM
        write(6,*) ' h, Omg0, lamda0 = ',SI_h0,SI_omg0,SI_lam0
        write(6,*) ' OMGg        = ',Ini_OMGg
#endif
        write(6,*) ' total Nb Ng, Ns Ndm = ',npt,ngt,nst,ndmt
#ifdef MULTI_LEV
        do i=0,Ini_mul
          write(6,*) i,' level n,dx  = ',ndt(i),dxi(i)
        enddo
#endif
        write(6,*) ' gamma       = ',GAM
        write(6,*) ' simulate z  = ',1.0d0/SI_a-1.0d0,'->',Ini_ze
        write(6,*) '          t  = ',TM_tot,'->',TM_t0
        write(6,*) '  output dt  = ',OF_odt
        write(6,*) ' # of output = ',Ini_nd
        write(6,*) ' # of output file for each output=',SI_nof
        if(SI_nof.gt.nprocs) then
          write(6,*) ' Error: It has to be less than nprocs=',nprocs
          stop
        endif
        if(Ini_mul.eq.0) then
          write(6,*) '   level 0 Only'
        endif
      endif        

! *** setting the other parameters ***
#ifdef STAR
! *** for star ***	
      do i=ng,ng+ns-1
        x_p(i)=xc_p(i)
        y_p(i)=yc_p(i)
        z_p(i)=zc_p(i)
        vx_p(i)=vnx_p(i)
        vy_p(i)=vny_p(i)
        vz_p(i)=vnz_p(i)
#if defined(COSM) && defined(BOXSIM)
        vvnx_p(i)=hubi*xc_p(i)+vnx_p(i)/SI_a
        vvny_p(i)=hubi*yc_p(i)+vny_p(i)/SI_a
        vvnz_p(i)=hubi*zc_p(i)+vnz_p(i)/SI_a
#else
        vvnx_p(i)=vnx_p(i)
        vvny_p(i)=vny_p(i)
        vvnz_p(i)=vnz_p(i)
#endif
#ifdef SIM1D
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i)))
#elif defined(SIM2D)
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**0.5d0)
#else
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**THIRD)
#endif
        omgh_p(i)=1.0d0
        zetah_p(i)=0.0d0
        flagt_p(i)=0
        dtr_p(i)=INF
        dt_p(i)=INF
        dtmnb_p(i)=INF
        dtg_p(i)=INF
        ax_p(i)=0.0d0
        ay_p(i)=0.0d0
        az_p(i)=0.0d0
      enddo
#endif
#ifdef DM
! *** for DM ***	
      do i=0,ndm-1
        xc_dm(i)=x_dm(i)
        yc_dm(i)=y_dm(i)
        zc_dm(i)=z_dm(i)
#ifdef SIM1D
        h_dm(i)=(1.0d0+ERRH)*ETAHDM*((m_dm(i)/rho_dm(i)))
#elif defined(SIM2D)
        h_dm(i)=(1.0d0+ERRH)*ETAHDM*((m_dm(i)/rho_dm(i))**0.5d0)
#else
        h_dm(i)=(1.0d0+ERRH)*ETAHDM*((m_dm(i)/rho_dm(i))**THIRD)
#endif
        omgh_dm(i)=1.0d0
        zetah_dm(i)=0.0d0
        list_adm(i)=i
        flagt_dm(i)=0
        dtr_dm(i)=INF
        dt_dm(i)=0.0d0
        dtmnb_dm(i)=INF
      enddo
#endif
#ifdef GAS
      do i=0,ng-1
        x_p(i)=xc_p(i)
        y_p(i)=yc_p(i)
        z_p(i)=zc_p(i)
        vx_p(i)=vnx_p(i)
        vy_p(i)=vny_p(i)
        vz_p(i) = vnz_p(i)
#if defined(COSM) && defined(BOXSIM)
        vvnx_p(i)=hubi*xc_p(i)+vnx_p(i)/SI_a
        vvny_p(i)=hubi*yc_p(i)+vny_p(i)/SI_a
        vvnz_p(i)=hubi*zc_p(i)+vnz_p(i)/SI_a
#else
        vvnx_p(i)=vnx_p(i)
        vvny_p(i)=vny_p(i)
        vvnz_p(i)=vnz_p(i)
#endif

#ifdef SIM1D
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i)))
#elif defined(SIM2D)
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**0.5d0)
#else
        h_p(i)=(1.0d0+ERRH)*ETAH*((m_p(i)/rho_p(i))**THIRD)
#endif
        omgh_p(i)=1.0d0
        zetah_p(i)=0.0d0
        as_p(i)=p_p(i)/(rho_p(i)**GAM)
      enddo
! *** For Gas, Set flagc ***
      do i=0,ng-1
        flagc_p(i)=0
        flagt_p(i)=0
        dtr_p(i)=INF
        dt_p(i)=INF
        dtmnb_p(i)=INF
!        alpv_p(i)=V_ALPHAMAX
        alpv_p(i)=V_ALPHAINI
!        alpu_p(i)=U_ALPHAMAX
        alpu_p(i)=0.0d0
        hvsigdt_p(i)=INF
        dtg_p(i)=INF
        ax_p(i)=0.0d0
        ay_p(i)=0.0d0
        az_p(i)=0.0d0
#ifdef METALDIFF
        vsig_p(i)=0.0d0
#endif
      enddo
#endif      

#ifdef METAL
! initialise all the metal 
      do i=0,np-1
        mzHe0_p(i)=mzHe_p(i)
        mzC0_p(i)=mzC_p(i)
        mzN0_p(i)=mzN_p(i)
        mzO0_p(i)=mzO_p(i)
        mzNe0_p(i)=mzNe_p(i)
        mzMg0_p(i)=mzMg_p(i)
        mzSi0_p(i)=mzSi_p(i)
        mzFe0_p(i)=mzFe_p(i)
        mzZ0_p(i)=mzZ_p(i)
      enddo
#endif

      if(Ini_mul.gt.0) then
        SI_ndm1=nd(Ini_mul)
        SI_ndm1t=ndt(Ini_mul)
      else
        SI_ndm1t=ndmt
        SI_ndm1=ndm
      endif
      flagcont=0
      SI_flagout = 0
      SI_mms = 0.0d0
#ifdef STAR
      SI_esns = 0.0d0
      SI_zsns = 0.0d0
! *** star formation time ***
      TMsf_t=TM_tot
#endif
      idum = -111      

#ifdef SGRAV
! *** set softening length limit ***
#ifdef SF
      SI_eps=(XHSOL*(DU/MP)/NSTH)**(1.0d0/3.0d0)*ETAH
#else
      SI_eps=(XHSOL*(DU/MP)/NSTH)**(1.0d0/3.0d0)*ETAH
#endif
      if(myrank.eq.0) then
        write(6,*) ' softening limit = mp^1/3 x',SI_eps
#if defined(GAS) || defined(STAR)
        write(6,*) '  for gas and/or star ',SI_eps*(m_p(0)**THIRD)
        write(6,*) '  corresponding nh(Zsolar)=' &
         ,0.706d0*(DU/MP)*((ETAH/SI_eps)**3)
#endif
#if defined(GEPSLIM_G) && defined(GAS)
        write(6,*) '  minimu h=',FH*SI_eps*(m_p(0)**THIRD)
#endif
#ifdef DM
        write(6,*) '  for DM ',SI_eps*(m_dm(0)**THIRD)
#endif
      endif
! *** check the minimum h ***
      hmin=INF
#if defined(GAS) || defined(STAR)
      do i=0,np-1
        if(h_p(i)/(m_p(i)**THIRD).lt.hmin) then 
          hmin=h_p(i)/(m_p(i)**THIRD)
        endif
      enddo
#endif
#ifdef DM
      do i=0,ndm-1
        if(h_dm(i)/(m_dm(i)**THIRD).lt.hmin) then 
          hmin=h_dm(i)/(m_dm(i)**THIRD)
       endif
      enddo
#endif
      hmint=INF
      call MPI_ALLREDUCE(hmin,hmint,1,MPI_DOUBLE_PRECISION &
       ,MPI_MIN,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0) then
        write(6,*) ' minimum h/(m^1/3)=',hmin
      endif
#ifdef COSM
      SI_eps0=SI_eps
! *** set redshift when SI_eps becomes effective ***
! *** assume mass of DM is always bigger ***
#ifdef DM
      etahmax=ETAHDM
! *** try slightly higher density ***
      SI_epsa=SI_eps*2.0d0*(m_dm(0)**THIRD)/(dxi(Ini_mul)*etahmax)
#else
      etahmax=ETAH
      SI_epsa=SI_eps*2.0d0*(m_p(0)**THIRD)/(dxi(Ini_mul)*etahmax)
#endif
      if(myrank.eq.0) then
        write(6,*) ' SI_epsa=',SI_epsa
      endif
      if(SI_a.lt.SI_epsa) then
        SI_eps=(1.0d0-SI_epsa+SI_a)*SI_eps
      endif
      if(hmint.lt.SI_eps) then
        if(myrank.eq.0) then
          write(6,*) ' Warning in set_ini(): SI_eps may be too big'
          write(6,*) ' hmint,SI_eps=',hmint,SI_eps
#ifdef GAS
          write(6,*) ' m/rho**THIRD=',(m_p(0)/rho_p(0))**THIRD
#else
          write(6,*) ' m/rho**THIRD=',(m_dm(1)/rho_dm(1))**THIRD
#endif
        endif
!        call MPI_ABORT(MPI_COMM_WORLD,ierr)
!        stop
      endif
#endif
! *** fixed softening case.
#ifdef DM
#ifdef FIXEDDMEPS
! *** set to minimum eps 
      do i=0,ndm-1
        h_dm(i)=SI_eps*(m_dm(i)**THIRD)
      enddo
#endif
#endif
#ifdef STAR
! *** set softening for stars
      do i=ng,ng+ns-1
        pn=list_ap(i)
        h_p(pn)=SI_eps*(m_p(pn)**THIRD)
        omgh_p(pn)=1.0d0
        zetah_p(pn)=0.0d0
        myu_p(pn)=MYU
      enddo
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
          write(6,*) ' Error: id0 should be bigger than id1 for DM.'
        endif
        call MPI_FINALIZE()
        stop
      endif
! *** change nadm
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
#endif

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
        write(6,*) ' ARTF DM particles id between ' &
         ,SI_artfdm(0),SI_artfdm(1)
      endif
      if(SI_artfdm(0).gt.SI_artfdm(1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: id0 should be bigger than id1 for DM.'
        endif
        call MPI_FINALIZE()
        stop
      endif
#endif

! *** output initial parameters ***
      if(myrank.eq.0) then
! ****** Initialzation for output file *****
        call system("rm ./output/ana/*")
        call system("rm ./output/data/*")           
        call system("rm ./output/*anim*")
        open(60,file='./output/ana/inipara.dat',status='unknown')
        write(60,'(a34,8I10)') 'npt,ndmt,np,ng,ns,ndm,ndm1,ndm1t=' &
         ,npt,ndmt,np,ng,ns,ndm,SI_ndm1,SI_ndm1t
#ifdef COSM
        write(60,'(a19,4(1pE13.5))') 'h,Omg0,lamda0,omgb=' &
         ,SI_h0,SI_omg0,SI_lam0,Ini_OMGg
#endif
!        write(60,'(a5,1pE13.5)') 'OFdt=',OF_odt
! since pv30.4
        write(60,'(a5,1pE25.17)') 'OFdt=',OF_odt
#ifdef COS
!                                      1234567890123456789012345678901234
        write(60,'(a14,2(1pE25.17))') ' SI_eps0,epsa=' &
         ,SI_eps0,SI_epsa
#else
        write(60,'(a8,2(1pE25.17))') ' SI_eps=',SI_eps
#endif
        write(60,'(a6,2(1pE13.5))') 'zi,ze=',1.0d0/SI_a-1.0d0,Ini_ze
        write(60,'(a6,2(1pE12.5))') 'ti,te=',TM_tot,TM_t0
        write(60,'(a4,1pE13.5)') 'GAM=',GAM
#ifdef GAS
        write(60,'(a3,(1pE13.5))') 'mb=',m_p(0)
#endif
        do i=0,Ini_mul
          write(60,'(i2,a12,I10,1pE13.5)') i,' level n,dx=',ndt(i),dxi(i)
           
        enddo
        nc = 0
#ifdef DM
        do il=Ini_mul,0,-1
          nc = nc+nd(il)
          write(60,'(a4,i2,a2,(1pE13.5))') 'mdm(l',il,')=',m_dm(nc-1)
        enddo
#endif
        write(60,'(a21,2(1pE13.5))') ' MGROOT,THETA,ALPHAT=',MGROOT,THETA,ALPHAT
        write(60,'(a10,I10)') ' NATTABLE=',NATTABLE
        write(60,'(a24,4(1pE13.5))')  &
!          12345678901234567890123456789012345
         ' CCFL,CDYN,CGRAV,FITLIM=',CCFL,CDYN,CGRAV,FITLIM
        write(60,'(a35,7(1pE13.5))')  &
!          12345678901234567890123456789012345
          ' AV ALPHAINI,MAX,MIN,BETA,SHARE,U,C' &
          ,V_ALPHAINI,V_ALPHAMAX,V_ALPHAMIN,V_BETA,SHAREV,U_ALPHAMAX,C_AV
        write(60,'(a18,3(1pE13.5))')  &
         ' ETAH,ETAHDM,ERRH=',ETAH,ETAHDM,ERRH
        write(60,'(a8,3(1pE13.5))')  &
         ' DHFLIM=',DHFLIM
        write(60,'(a38,5(1pE13.5))')  &
!          123456789012345678901234567890123456789012345678901
         ' MGTL,MGTU,CCOOLDT,CCOOLEDT,CCOOLRIDT=' &
         ,MGTL,MGTU,CCOOLDT,CCOOLEDT,CCOOLRIDT
        write(60,'(a36,5(1pE13.5))')  &
!         123456789012345678901234567890123456789012345678901
         ' TLLIMIT,THLIMIT,CLIMIT,EPSC,DU_UPU=' &
         ,TLLIMIT,THLIMIT,CLIMIT,EPSC,DU_UPU
        write(60,'(a22,2(1pE13.5))') &
!         12345678901234567890123456789012345
         ' DTSF,MINTMSFD,DLDUFD=',DTSF,MINTMSFD,LDUFD 
        write(60,'(a30,9(1pE13.5))') &
         ' NSTH,NSTHLIM,CSEFF,SNEU,SWEU=',NSTH,NSTHLIM,CSEFF,SNEU,SWEU
!         12345678901234567890123456789012345
#ifdef BOXSIM
        write(60,'(a10,3(1pE13.5))') ' box size=',SI_lbox(0) &
         ,SI_lbox(1),SI_lbox(2)
#ifdef COSM
!                                     12345678901234567890123456789012345
        write(60,'(a22,3(1pE13.5))') ' box size (comoving) =',SI_lbox0(0) &
         ,SI_lbox0(1),SI_lbox0(2)
#endif
#endif
#ifdef TREEPM
        write(60,'(a15,3I10)') ' Mesh nx,ny,nz=',nx_m,ny_m,nz_m
#endif
#ifdef BCX_FIX
        write(60,'(a9,I10)') ' BCX_ndx=',BCX_ndx
#ifdef MHD
        write(60,'(a25,5(1pE13.5))') 'Left BCX_rho,p,vx,vy,vz=' &
          ,BCX_rhox(0),BCX_px(0),BCX_vx(0),BCX_vy(0),BCX_vz(0)
        write(60,'(a13,3(1pE13.5))')'BCX_bx,by,bz=' &
          ,BCX_bx(0),BCX_by(0),BCX_bz(0)
        write(60,'(a26,5(1pE13.5))') 'Right BCX_rho,p,vx,vy,vz=' &
            ,BCX_rhox(1),BCX_px(1),BCX_vx(1),BCX_vy(1),BCX_vz(1)
        write(60,'(a13,3(1pE13.5))')'BCX_bx,by,bz=' &
          ,BCX_bx(1),BCX_by(1),BCX_bz(1)
#else
        write(60,'(a18,3(1pE13.5))') 'Left BCX_rho,p,vx=',BCX_rhox(0) &
         ,BCX_px(0),BCX_vx(0)
        write(60,'(a19,3(1pE13.5))') 'Right BCX_rho,p,vx=' &
         ,BCX_rhox(1),BCX_px(1),BCX_vx(1)
#endif
#endif
#ifdef MAXR
        write(60,'(a9,1pE13.5)') ' MAXRLIM=',MAXRLIM
#endif
#ifdef METALDIFF
        write(60,'(a9,1pE13.5)') ' ZDFCFAC=',ZDFCFAC
#endif
#ifdef UTHLIM_JM
        write(60,'(a9,(1pE13.5))') ' ULJMFAC=',ULJMFAC
#endif
#ifdef COSM
!                           12345678901234567890123456789012345
        write(60,'(a33)') '###  cosmological simulation  ###'
#endif
#ifdef MULTI_LEV
        write(60,'(a33)') '###  multi level  ###            '
#endif
#ifdef GAS
        write(60,'(a33)') '###  including GAS  ###          '
#endif
#ifdef DM
        write(60,'(a33)') '###  including DM  ###           '
#endif
#ifdef SGRAV
        write(60,'(a33)') '###  include self-gravity  ###   '
#endif
#ifdef COOL
        write(60,'(a33)') '###  include cooling  ###        '
#ifdef NOUVB
        write(60,'(a33)') '###  no UVB heating  ###         '
#endif
#endif
#ifdef STAR
        write(60,'(a22)') '###  influde star ### '
#endif
#ifdef SF
        write(60,'(a33)') '###  influde star formation  ### '
#endif
#ifdef SF_EFD
        write(60,'(a33)') '   ###  include energy feedback  ###'
#ifdef DTFD
        write(60,'(a18)') '   ### DTFD on ###'
#endif
#endif
#ifdef SF_ZFD
        write(60,'(a33)') '   ###  include metal feedback  ###'
#endif
#ifdef COREDW
        write(60,'(a33)') '   ###  core dw/ds for P grad   ###'
#endif
#ifdef GEPSLIM_G
        write(60,'(a31,1pE13.5)') &
         '###  gas softening limit on FH=',FH
#endif
#ifdef GEPSLIM_VH
        write(60,'(a36)') '   ###  gas no lim for smoothing ###'
#endif
#ifdef UTHLIM_JM
        write(60,'(a33)') '   ###  UTHLIM (JM) on  ###'
#endif
#ifdef METALDIFF
        write(60,'(a33)') '   ###  METALDIFF on ###'
#endif
#ifdef MAXR
        write(60,'(a33)') '   ###  MAXR on ###'
#endif
#ifdef MEXT
        write(60,'(a33)') '   ###  MEXT on ###'
#endif
#ifdef FIXEDDMEPS
!                          12345678901234567890123456789012345
        write(60,'(a30)') '   ###  Fixed DM softening ###'
#endif
#ifdef RANAZ
!                          1234567890123456789012345678901234567890
        write(60,'(a38)') '   ###  move particles azimuthally ###'
        write(60,'(a10,3I10,1pE13.5)') '   i012,t=',SI_ranazi(0) &
         ,SI_ranazi(1),SI_ranazi(2),SI_ranazt*TMUGYR
#endif
#ifdef QPTREEF
        write(60,'(a24)') '   ###  QP TREE F on ###'
#endif
#ifdef FASTG
        write(60,'(a20)') '   ###  FASTG on ###'
#endif
#ifdef DTDADT
!                          123456789012345678901234567890123456
        write(60,'(a21)') '   ###  DTDADT on ###'
#endif
#ifdef MINDT_SFD
!                          123456789012345678901234567890123456
        write(60,'(a24)') '   ###  MINDT_SFD on ###'
#endif
#ifdef NOAC
        write(60,'(a18)') '   ###  AC off ###'
#endif
#ifdef NOAVSWITCH
        write(60,'(a25)') '   ###  AV switch off ###'
#endif
#ifdef RHOINI
        write(60,'(a27)') '   ###  use initial RHO ###'
#endif
#ifdef PINI
        write(60,'(a25)') '   ###  use initial P ###'
#endif
#ifdef GLASS
!                          123456789012345678901234567890123456
        write(60,'(a24,F7.4,a4)') '   ###  GLASS on FGLASS=' &
        ,FGLASS,' ###'
#endif
#ifdef FIXEDP
        write(60,'(a28)') '   ###  fixed particlese ###'
        write(60,'(a12,2I10)') '   ID range=',SI_fid(0),SI_fid(1)
#endif
#ifdef FIXEDDMP
!                          123456789012345678901234567890123456
        write(60,'(a31)') '   ###  fixed DM particlese ###'
        write(60,'(a12,2I10)') '   ID range=',SI_fiddm(0),SI_fiddm(1)
#endif
#ifdef ARTFDM
!                          123456789012345678901234567890123456
        write(60,'(a30)') '   ###  artF DM particlese ###'
        write(60,'(a12,2I10)') '   ID range=',SI_artfdm(0),SI_artfdm(1)
#endif

        close(60)
        write(6,*) ' end writing inipara.dat'
      endif

#ifdef GAS
!      write(fileo,'(a4,i3.3)') 'inip',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,ng-1
!        pn=list_ap(i) 
!        write(60,'(9(1pE13.5),2I10)') xc_p(pn),yc_p(pn),zc_p(pn) &
!         ,vnx_p(pn),vny_p(pn),vnz_p(pn),m_p(pn),rho_p(pn) &
!         ,h_p(pn),id_p(pn),flagc_p(pn)
!      enddo
!      close(60)
#endif

end subroutine
