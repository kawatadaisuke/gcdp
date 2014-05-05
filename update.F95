#include "gcdp.def"
! *****************************************************
!    update.F for gcd+ ver f03.0
!  27 Jan. 2013  written by D.KAWATA
! ***************************************************** 

! ***********  Definition of updatev()   **********
! * Setup Vn+1/2 for active particle *
! this is used in updatev_pos and updatevn.
subroutine updatev(np,ndm)
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif      
#ifdef DM
      use gcdp_dm
#endif
      implicit none

      integer,intent(in) :: np,ndm
      integer i,pn
#if defined(COSM) && defined(BOXSIM)
      double precision,allocatable :: int1oa(:),t0(:),t1(:)
#endif

#if defined(GAS) || defined(STAR)
#if defined(COSM) && defined(BOXSIM)
! time at n+1/3
      allocate(t0(0:np))
      allocate(t1(0:np))
      allocate(int1oa(0:np))
! assign time range 
      do i=0,np-1
        t0(i)=tc_p(i)
        t1(i)=tc_p(i)+0.5d0*vdt_p(i)
      enddo
! integrate 1/a from n to n+1/2 for each particle)
      call integrate1oa(np,t0,t1,int1oa)
#endif
! *** baryon Particle ***
      do i=0,np-1
#if defined(COSM) && defined(BOXSIM)
! store a^2 vpec at n+1/2
        vx_p(i)=vnx_p(i)+(asc_p(i)**2)*dvx_p(i)*int1oa(i)
        vy_p(i)=vny_p(i)+(asc_p(i)**2)*dvy_p(i)*int1oa(i)
        vz_p(i)=vnz_p(i)+(asc_p(i)**2)*dvz_p(i)*int1oa(i)
#else
        vx_p(i)=vnx_p(i)+(vdt_p(i)*dvx_p(i))*0.5d0
        vy_p(i)=vny_p(i)+(vdt_p(i)*dvy_p(i))*0.5d0
        vz_p(i)=vnz_p(i)+(vdt_p(i)*dvz_p(i))*0.5d0
#endif
      enddo

#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa)
#endif

#endif
#ifdef DM
#if defined(COSM) && defined(BOXSIM)
! time at n+1/3
      allocate(t0(0:ndm))
      allocate(t1(0:ndm))
      allocate(int1oa(0:ndm))
! assign time range 
      do i=0,ndm-1
        t0(i)=tc_dm(i)
        t1(i)=tc_dm(i)+0.5d0*vdt_dm(i)
      enddo
! integrate 1/a from n to n+1/2 for each particle)
      call integrate1oa(ndm,t0,t1,int1oa)
#endif
! *** DM Particle ***
      do i=0,ndm-1
#if defined(COSM) && defined(BOXSIM)
! store a^2 vpec at n+1/2
        vx_dm(i)=vnx_dm(i)+(asc_dm(i)**2)*dvx_dm(i)*int1oa(i)
        vy_dm(i)=vny_dm(i)+(asc_dm(i)**2)*dvy_dm(i)*int1oa(i)
        vz_dm(i)=vnz_dm(i)+(asc_dm(i)**2)*dvz_dm(i)*int1oa(i)
#else
        vx_dm(i)=vnx_dm(i)+(vdt_dm(i)*dvx_dm(i))*0.5d0
        vy_dm(i)=vny_dm(i)+(vdt_dm(i)*dvy_dm(i))*0.5d0
        vz_dm(i)=vnz_dm(i)+(vdt_dm(i)*dvz_dm(i))*0.5d0
#endif

      enddo
#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa)
#endif
#endif
      return
end subroutine

#ifdef GAS
! **********   Definition of updatev_vn()   **********
! * Define Vn+1 *
subroutine updatev_vn(ng)
      use gcdp_baryon
#if defined(COSM) && defined(BOXSIM)   
      use gcdp_system
#endif
      implicit none

      integer,intent(in) :: ng
      integer i,pn
#if defined(COSM) && defined(BOXSIM)
      double precision,allocatable :: int1oa(:),t0(:),t1(:)
#endif

#if defined(COSM) && defined(BOXSIM)
! time at n+1/3
      allocate(t0(0:ng))
      allocate(t1(0:ng))
      allocate(int1oa(0:ng))
! assign time range 
      do i=0,ng-1
        pn=list_ap(i)
        t0(i)=tc_p(pn)
        t1(i)=tc_p(pn)+vdt_p(pn)
      enddo
! integrate 1/a from n to n+1/2 for each particle)
      call integrate1oa(ng,t0,t1,int1oa)
#endif

      do i=0,ng-1
        pn=list_ap(i)          
#if defined(COSM) && defined(BOXSIM)
! vvnx_p a*vpec,  note that vnx_p and vx_p are a^2*vpec
        vvnx_p(pn)=(vnx_p(pn)+(asc_p(pn)**2)*dvx_p(i)*int1oa(i))/SI_a
        vvny_p(pn)=(vny_p(pn)+(asc_p(pn)**2)*dvy_p(i)*int1oa(i))/SI_a
        vvnz_p(pn)=(vnz_p(pn)+(asc_p(pn)**2)*dvz_p(i)*int1oa(i))/SI_a
#else
        vvnx_p(pn)=vnx_p(pn)+(vdt_p(pn)*dvx_p(pn))
        vvny_p(pn)=vny_p(pn)+(vdt_p(pn)*dvy_p(pn))
        vvnz_p(pn)=vnz_p(pn)+(vdt_p(pn)*dvz_p(pn))
#endif
      enddo

#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa)
#endif

end subroutine
#endif      

! *****   Definition of updatevn()   *****
subroutine updatevn(ng,nag,nadm,nas,nagravg)
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif      
#ifdef DM
      use gcdp_dm
#endif
#if defined(COSM) && defined(BOXSIM)   
      use gcdp_system
#endif
      implicit none

      integer,intent(in) :: ng,nag,nadm,nas,nagravg
      integer i,pn
#if defined(COSM) && defined(BOXSIM)   
      double precision,allocatable :: t0(:),t1(:),int1oa(:)
      double precision,allocatable :: t0g(:),t1g(:),int1oag(:)
#endif

#ifdef GAS
! *** Gas Particles *** 
#if defined(COSM) && defined(BOXSIM)
      allocate(t0g(0:nagravg))
      allocate(t1g(0:nagravg))
      allocate(int1oag(0:nagravg))
! set time range for each particle
      do i=0,nagravg-1
        pn=list_ap(i)
        t0g(i)=tcg_p(pn)+0.5d0*dtg_p(pn)
        t1g(i)=TM_tot
      enddo
! integrate(1/a) from n+1/2 to n+1 for each particle
      call integrate1oa(nagravg,t0g,t1g,int1oag)
! for integration for hydro value
      allocate(t0(0:nag))
      allocate(t1(0:nag))
      allocate(int1oa(0:nag))
! set time range for each particle
      do i=0,nag-1
        pn=list_ap(i)
        t0(i)=tc_p(pn)+0.5d0*dt_p(pn)
        t1(i)=TM_tot
      enddo
! integrate(1/a) from n+1/2 to n+1 for each particle
      call integrate1oa(nag,t0,t1,int1oa)
#endif
! *** for gravity acceleration ***
      do i=0,nagravg-1
        pn=list_ap(i)
#if defined(COSM) && defined(BOXSIM)
!  physical velocity
        vnx_p(pn)=vx_p(pn)+(SI_a**2)*dvx_p(pn)*int1oa(i) &
          +(SI_a**2)*ax_p(pn)*int1oag(i)
        vny_p(pn)=vy_p(pn)+(SI_a**2)*dvy_p(pn)*int1oa(i) &
          +(SI_a**2)*ay_p(pn)*int1oag(i)
        vnz_p(pn)=vz_p(pn)+(SI_a**2)*dvz_p(pn)*int1oa(i) &
          +(SI_a**2)*az_p(pn)*int1oag(i)
! update H, a and t
        ascg_p(pn)=SI_a
        tcg_p(pn)=TM_tot
        asc_p(pn)=SI_a
        tc_p(pn)=TM_tot
#else
        vnx_p(pn)=vx_p(pn)+(dtg_p(pn)*ax_p(pn))*0.5d0+(dt_p(pn)*dvx_p(pn))*0.5d0
        vny_p(pn)=vy_p(pn)+(dtg_p(pn)*ay_p(pn))*0.5d0+(dt_p(pn)*dvy_p(pn))*0.5d0
        vnz_p(pn)=vz_p(pn)+(dtg_p(pn)*az_p(pn))*0.5d0+(dt_p(pn)*dvz_p(pn))*0.5d0
#endif
      enddo

! *** for hydro acceleration ***
      do i=nagravg,nag-1
        pn=list_ap(i)
#if defined(COSM) && defined(BOXSIM)
!  physical velocity
        vnx_p(pn)=vx_p(pn)+(SI_a**2)*dvx_p(pn)*int1oa(i)
        vny_p(pn)=vy_p(pn)+(SI_a**2)*dvy_p(pn)*int1oa(i)
        vnz_p(pn)=vz_p(pn)+(SI_a**2)*dvz_p(pn)*int1oa(i)
! update H, a and t
        asc_p(pn)=SI_a
        tc_p(pn)=TM_tot
#else
        vnx_p(pn)=vx_p(pn)+(dt_p(pn)*dvx_p(pn))*0.5d0
        vny_p(pn)=vy_p(pn)+(dt_p(pn)*dvy_p(pn))*0.5d0
        vnz_p(pn)=vz_p(pn)+(dt_p(pn)*dvz_p(pn))*0.5d0
#endif
      enddo

#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa)
      deallocate(t0g)
      deallocate(t1g)
      deallocate(int1oag)
#endif

#endif

#ifdef DM
! *** DM Particles *** 
#if defined(COSM) && defined(BOXSIM)
      allocate(t0(0:nadm))
      allocate(t1(0:nadm))
      allocate(int1oa(0:nadm))
! set time range for each particle
      do i=0,nadm-1
        pn=list_adm(i)
        t0(i)=tc_dm(pn)+0.5d0*dt_dm(pn)
        t1(i)=TM_tot
      enddo
! integrate(1/a) from n+1/2 to n+1 for each particle
      call integrate1oa(nadm,t0,t1,int1oa)
#endif
      do i=0,nadm-1
        pn=list_adm(i)         
#if defined(COSM) && defined(BOXSIM)
! x  physical velocity
        vnx_dm(pn)=vx_dm(pn)+(SI_a**2)*dvx_dm(pn)*int1oa(i)
        vny_dm(pn)=vy_dm(pn)+(SI_a**2)*dvy_dm(pn)*int1oa(i)
        vnz_dm(pn)=vz_dm(pn)+(SI_a**2)*dvz_dm(pn)*int1oa(i)

!        if(id_dm(pn).eq.0) then
!          write(6,*) xc_dm(pn)/SI_a,yc_dm(pn)/SI_a,zc_dm(pn)/SI_a,vnx_dm(pn)-SI_hub*xc_dm(pn)
!        endif
! update H, a and t
        asc_dm(pn)=SI_a
        tc_dm(pn)=TM_tot
#else
        vnx_dm(pn)=vx_dm(pn)+(dt_dm(pn)*dvx_dm(pn))*0.5d0
        vny_dm(pn)=vy_dm(pn)+(dt_dm(pn)*dvy_dm(pn))*0.5d0
        vnz_dm(pn)=vz_dm(pn)+(dt_dm(pn)*dvz_dm(pn))*0.5d0
#endif

!        if(id_dm(pn).eq.0) then
!          write(6,*) xc_dm(pn)/SI_a,yc_dm(pn)/SI_a,zc_dm(pn)/,vx_dm(pn),dvx_dm(pn)
!        endif

      enddo

#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa)
#endif

#endif
#ifdef STAR
! *** Star Particles *** 
#if defined(COSM) && defined(BOXSIM)
      allocate(t0(0:nas))
      allocate(t1(0:nas))
      allocate(int1oa(0:nas))
! set time range for each particle
      do i=0,nas-1
        pn=list_ap(i+ng)
        t0(i)=tcg_p(pn)+0.5d0*dtg_p(pn)
        t1(i)=TM_tot
      enddo
! integrate(1/a) from n+1/2 to n+1 for each particle
      call integrate1oa(nas,t0,t1,int1oa)
#endif

      do i=0,nas-1
        pn=list_ap(i+ng)         
#if defined(COSM) && defined(BOXSIM)
!  physical velocity
        vnx_p(pn)=vx_p(pn)+(SI_a**2)*dvx_p(pn)*int1oa(i)
        vny_p(pn)=vy_p(pn)+(SI_a**2)*dvy_p(pn)*int1oa(i)
        vnz_p(pn)=vz_p(pn)+(SI_a**2)*dvz_p(pn)*int1oa(i)
! update H, a and t
        ascg_p(pn)=SI_a
        tcg_p(pn)=TM_tot
        asc_p(pn)=SI_a
        tc_p(pn)=TM_tot
#else
! *** use dtg_p, because new stars may have different dtg from dt
! *** and dvx_p is already only from gravity 
        vnx_p(pn)=vx_p(pn)+(dtg_p(pn)*dvx_p(pn))*0.5d0
        vny_p(pn)=vy_p(pn)+(dtg_p(pn)*dvy_p(pn))*0.5d0
        vnz_p(pn)=vz_p(pn)+(dtg_p(pn)*dvz_p(pn))*0.5d0
#endif
      enddo

#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa)
#endif

#endif
end subroutine

! **********   Definition of updatepos()   **********
! * Define Virtual Xn+1 *
subroutine updatev_pos(np,ndm)
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif      
#ifdef DM
      use gcdp_dm
#endif
#if defined(BOXSIM) || defined(FIXEDP) || defined(FIXEDDMP)
      use gcdp_system
#endif
#ifdef MAXR
      use gcdp_const
#endif
      implicit none

      integer,intent(in) :: np,ndm      
      integer i
#ifdef MAXR
      double precision rp2,rm2,dvr,rp,vr,vt,th,ph,rm,rpxy
#endif
#if defined(COSM) && defined(BOXSIM)
      double precision,allocatable :: t0(:),t1(:),int1oa2(:)
#endif
#ifdef BOXSIM
      double precision L2(0:2)
      L2(0)=0.5d0*SI_lbox(0)
      L2(1)=0.5d0*SI_lbox(1)
      L2(2)=0.5d0*SI_lbox(2)
#endif

#ifdef MAXR
! *** maximum radius (100 kpc unit)
      rm2=MAXRLIM
      rm2=rm2**2
#endif

#ifdef GAS	
! *** Gas and Star Particle ***
#if defined(COSM) && defined(BOXSIM)
! get integrate(1/a^2) at n+1
      allocate(t0(0:np))
      allocate(t1(0:np))
      allocate(int1oa2(0:np))
! set time range
      do i=0,np-1
        t0(i)=tc_p(i)
        t1(i)=TM_tot+TM_dt
      enddo
! integrate(1/a^2) from n to n+1 for each particle
      call integrate1oa2(np,t0,t1,int1oa2)
#endif

      do i=0,np-1
#ifdef FIXEDP
        if(id_p(i).lt.SI_fid(0).or.id_p(i).gt.SI_fid(1)) then
#endif
#if defined(COSM) && defined(BOXSIM)
! vx_p a^2 vpec at n+1/2
! SI_a at n+1 step
        x_p(i)=SI_a*((xc_p(i)/asc_p(i))+vx_p(i)*int1oa2(i))
#else
        x_p(i)=xc_p(i)+vdt_p(i)*(vnx_p(i)+0.5d0*vdt_p(i)*dvx_p(i))
#endif
#ifdef BCX_PERIOD
        if(x_p(i).gt.L2(0)) then
          x_p(i)=x_p(i)-SI_lbox(0)
        else if(x_p(i).lt.-L2(0)) then
          x_p(i)=x_p(i)+SI_lbox(0)
        endif
#endif
#ifdef SIM1D
        y_p(i)=0.0d0
        z_p(i)=0.0d0
#else
#if defined(COSM) && defined(BOXSIM)
! vx_p a^2 vpec at n+1/2
! SI_a at n+1 step
        y_p(i)=SI_a*((yc_p(i)/asc_p(i))+vy_p(i)*int1oa2(i))
#else
        y_p(i)=yc_p(i)+vdt_p(i)*(vny_p(i)+0.5d0*vdt_p(i)*dvy_p(i))
#endif
#ifdef BCY_PERIOD
        if(y_p(i).gt.L2(1)) then
          y_p(i)=y_p(i)-SI_lbox(1)
        else if(y_p(i).lt.-L2(1)) then
          y_p(i)=y_p(i)+SI_lbox(1)
        endif
#endif
#ifdef SIM2D
        z_p(i)=0.0d0
#else
#if defined(COSM) && defined(BOXSIM)
! vz_p a^2 vpec at n+1/2
! SI_a at n+1 step
        z_p(i)=SI_a*((zc_p(i)/asc_p(i))+vz_p(i)*int1oa2(i))
#else
        z_p(i)=zc_p(i)+vdt_p(i)*(vnz_p(i)+0.5d0*vdt_p(i)*dvz_p(i))
#endif
#ifdef BCZ_PERIOD
        if(z_p(i).gt.L2(2)) then
          z_p(i)=z_p(i)-SI_lbox(2)
        else if(z_p(i).lt.-L2(2)) then
          z_p(i)=z_p(i)+SI_lbox(2)
        endif
#endif
#endif
#endif
#ifdef FIXEDP  
        endif
#endif
      enddo

#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa2)
#endif

#ifndef COSM
#ifdef MAXR 
      do i=0,np-1
        rp2=x_p(i)**2+y_p(i)**2+z_p(i)**2
        if(rp2.gt.rm2) then
!          x_p(i)=xc_p(i)
!          y_p(i)=yc_p(i)
!          z_p(i)=zc_p(i)
          rp2=xc_p(i)**2+yc_p(i)**2+zc_p(i)**2
          rp=dsqrt(rp2)
          vr=(vnx_p(i)*xc_p(i)+vny_p(i)*yc_p(i)+vnz_p(i)*zc_p(i))/rp
          vt=dsqrt(vnx_p(i)**2+vny_p(i)**2+vnz_p(i)**2-vr**2)
          if(vr.gt.0.0d0) then
            vnx_p(i)=vnx_p(i)-vr*xc_p(i)/rp
            vny_p(i)=vny_p(i)-vr*yc_p(i)/rp
            vnz_p(i)=vnz_p(i)-vr*zc_p(i)/rp
          endif
          dvr=((dvx_p(i)+ax_p(i))*xc_p(i)+(dvy_p(i)+ay_p(i))*yc_p(i) &
           +(dvz_p(i)+az_p(i))*zc_p(i))/rp+(vt**2)/rp
          if(dvr.gt.0.0d0) then
            dvx_p(i)=dvx_p(i)-dvr*xc_p(i)/rp
            dvy_p(i)=dvy_p(i)-dvr*yc_p(i)/rp
            dvz_p(i)=dvz_p(i)-dvr*zc_p(i)/rp
          endif
          pdu_p(i)=0.0d0
          ndu_p(i)=0.0d0
!          vnx_p(i)=0.0d0
!          vny_p(i)=0.0d0
!          vnz_p(i)=0.0d0
!          dvx_p(i)=0.0d0
!          dvy_p(i)=0.0d0
!          dvz_p(i)=0.0d0
          x_p(i)=xc_p(i)+vdt_p(i)*(vnx_p(i)+0.5d0*vdt_p(i)*dvx_p(i))
          y_p(i)=yc_p(i)+vdt_p(i)*(vny_p(i)+0.5d0*vdt_p(i)*dvy_p(i))
          z_p(i)=zc_p(i)+vdt_p(i)*(vnz_p(i)+0.5d0*vdt_p(i)*dvz_p(i))
          rp2=x_p(i)**2+y_p(i)**2+z_p(i)**2
          if(rp2.gt.rm2) then
            rp=dsqrt(rp2)
            rpxy=dsqrt(x_p(i)**2+y_p(i)**2)
            th=dasin(z_p(i)/rp)
            ph=dacos(x_p(i)/rpxy)
            if(y_p(i).lt.0.0d0) then  
              ph=2.0d0*M_PI-ph
            endif
            z_p(i)=MAXRLIM*dsin(th)
            rpxy=MAXRLIM*dcos(th)
            x_p(i)=rpxy*dcos(ph)
            y_p(i)=rpxy*dsin(ph)
          endif
        endif
      enddo
#endif
#endif
#endif

#ifdef DM      

! *** DM Particle ***
#if defined(COSM) && defined(BOXSIM)
! get integrate(1/a^2) at n+1
      allocate(t0(0:ndm))
      allocate(t1(0:ndm))
      allocate(int1oa2(0:ndm))
! set time range
      do i=0,ndm-1
        t0(i)=tc_dm(i)
        t1(i)=TM_tot+TM_dt
      enddo
! integrate(1/a^2) from n to n+1 for each particle
      call integrate1oa2(ndm,t0,t1,int1oa2)
#endif

      do i=0,ndm-1
#ifdef FIXEDDMP
        if(id_dm(i).lt.SI_fiddm(0).or.id_dm(i).gt.SI_fiddm(1)) then
#endif
#if defined(COSM) && defined(BOXSIM)
! vx_dm a^2 vpec at n+1/2
! SI_a at n+1 step
        x_dm(i)=SI_a*((xc_dm(i)/asc_dm(i))+vx_dm(i)*int1oa2(i))
#else
        x_dm(i)=xc_dm(i)+vdt_dm(i)*vx_dm(i)
#endif
#ifdef BCX_PERIOD
        if(x_dm(i).gt.L2(0)) then
          x_dm(i)=x_dm(i)-SI_lbox(0)
        else if(x_dm(i).lt.-L2(0)) then
          x_dm(i)=x_dm(i)+SI_lbox(0)
        endif
#endif
#ifdef SIM1D
        y_dm(i)=0.0d0
        z_dm(i)=0.0d0
#else
#if defined(COSM) && defined(BOXSIM)
! vy_dm a^2 vpec at n+1/2
! SI_a at n+1 step
        y_dm(i)=SI_a*((yc_dm(i)/asc_dm(i))+vy_dm(i)*int1oa2(i))
#else
        y_dm(i)=yc_dm(i)+vdt_dm(i)*vy_dm(i)
#endif
#ifdef BCY_PERIOD
        if(y_dm(i).gt.L2(1)) then
          y_dm(i)=y_dm(i)-SI_lbox(1)
        else if(y_dm(i).lt.-L2(1)) then
          y_dm(i)=y_dm(i)+SI_lbox(1)
        endif
#endif
#ifdef SIM2D
        z_dm(i)=0.0d0
#else
#if defined(COSM) && defined(BOXSIM)
! vy_dm a^2 vpec at n+1/2
! SI_a at n+1 step
        z_dm(i)=SI_a*((zc_dm(i)/asc_dm(i))+vz_dm(i)*int1oa2(i))
#else
        z_dm(i)=zc_dm(i)+vdt_dm(i)*vz_dm(i)
#endif
#ifdef BCZ_PERIOD
        if(z_dm(i).gt.L2(2)) then
          z_dm(i)=z_dm(i)-SI_lbox(2)
        else if(z_dm(i).lt.-L2(2)) then
          z_dm(i)=z_dm(i)+SI_lbox(2)
        endif
#endif
#endif
#endif
#ifdef FIXEDDMP
        endif
#endif
      enddo

#if defined(COSM) && defined(BOXSIM)
      deallocate(t0)
      deallocate(t1)
      deallocate(int1oa2)
#endif

#ifndef COSM
#ifdef MAXR 
      do i=0,ndm-1
        rp2=xc_dm(i)**2+yc_dm(i)**2+zc_dm(i)**2
        if(rp2.gt.rm2) then
          x_dm(i)=xc_dm(i)
          y_dm(i)=yc_dm(i)
          z_dm(i)=zc_dm(i)
          vnx_dm(i)=0.0d0
          vny_dm(i)=0.0d0
          vnz_dm(i)=0.0d0
        endif
      enddo
#endif
#endif

#endif

end subroutine

! **********   Definition of updatepos   **********
! * Define Xn+1 for active particle *
subroutine updatepos(ng,nag,nadm,nas)
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif      
#ifdef DM
      use gcdp_dm
#endif
      implicit none

      integer,intent(in) :: ng,nag,nadm,nas
      integer i,pn
#ifdef GAS
! *** Gas Particles ***
      do i=0,nag-1
        pn=list_ap(i)         
        xc_p(pn)=x_p(pn)
        yc_p(pn)=y_p(pn)
        zc_p(pn)=z_p(pn)
      enddo
#endif
#ifdef DM
! *** DM Particles ***
      do i=0,nadm-1
        pn=list_adm(i)
        xc_dm(pn)=x_dm(pn)
        yc_dm(pn)=y_dm(pn)
        zc_dm(pn)=z_dm(pn)
      enddo
#endif
#ifdef STAR
! *** Star Particles ***
      do i=ng,ng+nas-1
        pn = list_ap(i)         
        xc_p(pn)=x_p(pn)
        yc_p(pn)=y_p(pn)
        zc_p(pn)=z_p(pn)
      enddo
#endif

end subroutine

! *********   Definition of updatesl()  since pv.29.1 **********
subroutine updatesl(ng,nag,ns,nas,ndm)
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif      
#ifdef DM
#ifndef FIXEDDMEPS
      use gcdp_dm
#endif
#endif
#if defined(GEPSLIM_G) || defined(DMEPSLIM_G)
      use gcdp_system
#endif
      use gcdp_const

      implicit none

      integer,intent(in) :: ng,ns,ndm,nag,nas
      integer i,j,pn,np
      double precision sh,ah,epsh,hp,decfac

! *** factor of decrease if hp<0 ***
      decfac=0.8d0   

#ifdef GAS
! *** for SPH ***
      do i=0,nag-1
        pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).le.0) then
#endif
#ifdef SIM1D
        hp=h_p(pn)+h_p(pn)*div_v_p(pn)*dt_p(pn)/omgh_p(pn)
#elif defined(SIM2D)
        hp=h_p(pn)+0.5d0*h_p(pn)*div_v_p(pn)*dt_p(pn)/omgh_p(pn)
#else
        hp=h_p(pn)+THIRD*h_p(pn)*div_v_p(pn)*dt_p(pn)/omgh_p(pn)
#endif
        if(hp.gt.DHFLIM*h_p(pn)) then
!          write(6,'(a11,5(1pE13.5))') 'large dh_g=',hp,h_p(pn)
!     &     ,div_v_p(pn),rho_p(pn),omgh_p(pn)
          h_p(pn)=DHFLIM*h_p(pn)
        else if(hp.gt.0.0d0) then
          h_p(pn)=hp
        else
          h_p(pn)=h_p(pn)*decfac
        endif
#ifdef SGRAV
#ifdef GEPSLIM_G
        if(h_p(pn).lt.SI_eps*(m_p(pn)**THIRD)) then
          h_p(pn)=SI_eps*(m_p(pn)**THIRD)
          omgh_p(pn)=1.0d0
          zetah_p(pn)=0.0d0
        endif
#endif
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
        endif
#endif
      enddo
#endif

#ifdef DM
#ifndef FIXEDDMEPS
! *** for DM: only active paritcles ***
      do i=0,ndm-1
        pn=list_adm(i)
#ifdef SIM1D
        hp=h_dm(pn)+h_dm(pn)*div_v_dm(pn)*dt_dm(pn)/omgh_dm(pn)
#elif defined(SIM2D)
        hp=h_dm(pn)+0.5d0*h_dm(pn)*div_v_dm(pn)*dt_dm(pn)/omgh_dm(pn)
#else
        hp=h_dm(pn)+THIRD*h_dm(pn)*div_v_dm(pn)*dt_dm(pn)/omgh_dm(pn)
#endif
        if(hp.gt.DHFLIM*h_dm(pn)) then
          h_dm(pn)=DHFLIM*h_dm(pn)
        else if(hp.gt.0.0d0) then
          h_dm(pn)=hp
        else 
          h_dm(pn)=h_dm(pn)*decfac
        endif
#ifdef DMEPSLIM_G
        if(h_dm(pn).lt.SI_eps*(m_dm(pn)**THIRD)) then
          h_dm(pn)=SI_eps*(m_dm(pn)**THIRD)
          omgh_dm(pn)=1.0d0
          zetah_dm(pn)=0.0d0
        endif
#endif
      enddo
#endif
#endif

end subroutine

#ifdef GAS
! **********   Definition of updatep_cs()   **********
subroutine updatep_cs(ng)
      use gcdp_baryon
      use gcdp_const
      implicit none

      integer,intent(in) :: ng
      integer i,pn

! *** Pressure ***
      if(GAM.eq.1.0d0) then
        do i=0,ng-1
          pn = list_ap(i)
          p_p(pn)=u_p(pn)*rho_p(pn)
        enddo
      else 
! #pragma loop novrec p_p,u_p,rho_p
        do i=0,ng-1
          pn = list_ap(i)           
#if defined(SF_EFD) || defined(SF_ZFD)
          if(flagfd_p(pn).le.0) then
#endif
          p_p(pn)=(GAM-1.0d0)*rho_p(pn)*u_p(pn)
#if defined(SF_EFD) || defined(SF_ZFD)
          endif
#endif 
        enddo
      endif        

! *** Sound Velocity ***
      do i=0,ng-1
        pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
       if(flagfd_p(pn).le.0) then
#endif
        cs_p(pn)=dsqrt(GAM*p_p(pn)/rho_p(pn))
        as_p(pn)=p_p(pn)/(rho_p(pn)**GAM)
#if defined(SF_EFD) || defined(SF_ZFD)
        endif
#endif
      enddo

end subroutine


! **********   Definition of updatalv()   **********
subroutine updatealpvu(ng)
      use gcdp_const
      use gcdp_baryon
#if defined(MHD) && defined(TIMEDIS)
      use gcdp_mhd
#endif
      implicit none

      integer,intent(in) :: ng
      double precision tau,sav,dalpv_tau,dalpu_tau,sau
      integer i,pn
#if defined(MHD) && defined(TIMEDIS)
      double precision dalpb_tau,sab
#endif

      do i=0,ng-1
        pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).le.0) then
#endif
! *** excluding new feedback gas particles ***
!        if(vsigmax_p(pn).gt.0.0d0) then
! Price (2008): looks better and faster off
!          tau=0.5d0*h_p(pn)/(C_AV*vsigmax_p(pn))

! Rosswog & Price (2007)
        tau=0.5d0*h_p(pn)/(C_AV*cs_p(pn))
        dalpv_tau=-(alpv_p(pn)-V_ALPHAMIN)/tau
        dalpu_tau=-alpu_p(pn)/tau
        sav=0.0d0
        if(-div_v_p(pn).gt.sav) then
          sav=-div_v_p(pn)
        endif    
! Rosswog & Price (2007): Price (2008) seems quick enough.
        sav=sav*(V_ALPHAMAX-alpv_p(pn)) &
! no fij factor seems better, when strong shock is involved. (6/10 2008)
         *dabs(div_v_p(pn))/(dabs(div_v_p(pn))+SHAREV*arot_v_p(pn) &
         +0.0002d0*cs_p(pn)/h_p(pn))
! *** update alp_v ***
        alpv_p(pn)=alpv_p(pn)+dt_p(pn)*(dalpv_tau+sav)
        if(alpv_p(pn).gt.V_ALPHAMAX) then
          alpv_p(pn)=V_ALPHAMAX
        endif
        if(alpv_p(pn).lt.V_ALPHAMIN) then
          alpv_p(pn)=V_ALPHAMIN
        endif 
! *** update alp_u (0.1 x 0.5 for larger h in gcd+) ***
        sau=0.05d0*h_p(pn)*dabs(d2u_p(pn))/dsqrt(u_p(pn))
        alpu_p(pn)=alpu_p(pn)+dt_p(pn)*(dalpu_tau+sau)
        if(alpu_p(pn).gt.U_ALPHAMAX) then
          alpu_p(pn)=U_ALPHAMAX
        endif
        if(alpu_p(pn).lt.0.0d0) then
          alpu_p(pn)=0.0d0
        endif

#if defined(MHD) && defined(TIMEDIS)
        dalpb_tau=-(alpb_p(pn)-ALPB_MIN)/tau
        sab=arotb_p(pn)/dsqrt(rho_p(pn))
        if(dabs(divb_p(pn))/dsqrt(rho_p(pn)).gt.sav) then
           sab=dabs(divb_p(pn))/dsqrt(rho_p(pn))
        endif
        alpb_p(pn)=alpb_p(pn)+dt_p(pn)*(dalpb_tau+sab)
! Check allowed min/max diss
        if(alpb_p(pn).lt.ALPB_MIN) then
           alpb_p(pn)=ALPB_MIN
        else if(alpb_p(pn).gt.ALPB_MAX) then
           alpb_p(pn)=ALPB_MAX

        endif
#endif

#if defined(SF_EFD) || defined(SF_ZFD)
        endif
#endif
      enddo

      return
end subroutine
! endif for GAS
#endif

#ifdef MHD
! *** B predictor ***
subroutine updatevb(ng)
      use gcdp_baryon
      use gcdp_mhd

      implicit none

      integer ng,i,pn

      do i=0,ng-1
         pn=list_ap(i)
         bx_p(pn)=pbx_p(pn)+dbx_p(pn)*vdt_p(pn)
         by_p(pn)=pby_p(pn)+dby_p(pn)*vdt_p(pn)
         bz_p(pn)=pbz_p(pn)+dbz_p(pn)*vdt_p(pn)
      enddo

      return
end subroutine

! *** B corrector only for active particles ***
subroutine updatebn(nag)
      use gcdp_baryon
      use gcdp_mhd

      implicit none

      integer nag, i, pn

      do i=0,nag-1
         pn=list_ap(i)
         bx_p(pn)=pbx_p(pn)+0.5d0*(dbx_p(pn)+pdbx_p(pn))*dt_p(pn)
         by_p(pn)=pby_p(pn)+0.5d0*(dby_p(pn)+pdby_p(pn))*dt_p(pn)
         bz_p(pn)=pbz_p(pn)+0.5d0*(dbz_p(pn)+pdbz_p(pn))*dt_p(pn)
         pbx_p(pn)=bx_p(pn)
         pby_p(pn)=by_p(pn)
         pbz_p(pn)=bz_p(pn)
      enddo

      return
end subroutine
#endif

#ifdef GAS
#ifdef METALDIFF
! *** metal diffusion for active gas particles ***
subroutine metal_diff(nag)
      use gcdp_baryon
      use gcdp_metal
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: nag   
      integer i,pn,ierr
      double precision mz0,tmms,tmmsr

      tmms=0.0d0
      do i=0,nag-1
        pn=list_ap(i)
#if defined(SF_EFD) || defined(SF_ZFD)
        if(flagfd_p(pn).le.0.and.zdA_p(pn).ne.0.0d0) then
#else
        if(zdA_p(pn).ne.0.0d0) then
#endif
! *** update metal abundance ***
! *** eq.5 of Grief et al. (2009) 
        mzHe_p(pn)=mzHe_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdHeB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzHe_p(pn).lt.0.0d0) then
          mzHe_p(pn)=0.0d0
        endif
        mzC_p(pn)=mzC_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdCB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzC_p(pn).lt.0.0d0) then
          mzC_p(pn)=0.0d0
        endif
        mzN_p(pn)=mzN_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdNB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzN_p(pn).lt.0.0d0) then
          mzN_p(pn)=0.0d0
        endif
        mzO_p(pn)=mzO_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdOB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzO_p(pn).lt.0.0d0) then
          mzO_p(pn)=0.0d0
        endif
        mzNe_p(pn)=mzNe_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdNeB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzNe_p(pn).lt.0.0d0) then
          mzNe_p(pn)=0.0d0
        endif
        mzMg_p(pn)=mzMg_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdMgB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzMg_p(pn).lt.0.0d0) then
          mzMg_p(pn)=0.0d0
        endif
        mzSi_p(pn)=mzSi_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdSiB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzSi_p(pn).lt.0.0d0) then
          mzSi_p(pn)=0.0d0
        endif
        mzFe_p(pn)=mzFe_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdFeB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
        if(mzFe_p(pn).lt.0.0d0) then
          mzFe_p(pn)=0.0d0
        endif
        mz0=mzZ_p(pn)
!        mzZ_p(pn)=(mzZ_p(pn)/m_p(pn))*dexp(zdA_p(pn)*dt_p(pn)) &
!         +(zdZB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))
!        mzZ_p(pn)=mzZ_p(pn)*m_p(pn)
!
        mzZ_p(pn)=mzZ_p(pn)*dexp(zdA_p(pn)*dt_p(pn)) &
         +(zdZB_p(pn)/zdA_p(pn))*(1.0d0-dexp(zdA_p(pn)*dt_p(pn)))

        if(mzZ_p(pn).lt.0.0d0) then
          write(6,*) ' Waring in metal_diff: mZ<0 at id',id_p(pn)
          write(6,*) ' original and new mZ (set to 0)=',mz0,mzZ_p(pn)
          mzZ_p(pn)=0.0d0
        endif
! for record
        tmms=tmms+(mzZ_p(pn)-mz0)
! check
        endif
      enddo

! *** for check ***
      tmmsr=0.0d0
      call MPI_ALLREDUCE(tmms,tmmsr,1,MPI_DOUBLE_PRECISION &
       ,MPI_SUM,MPI_COMM_WORLD,ierr)
      SI_mms=SI_mms+tmmsr

end subroutine
#endif
#endif
  

#ifdef RANAZ
! **********   Definition of ranazimuth   **********
subroutine ranazimuth(ng,ndm,ns,nag,nadm,nas,nagravg)
      use gcdp_const
#if defined(GAS) || defined(STARS)
      use gcdp_baryon
#endif
#ifdef DM
      use gcdp_dm
#endif
      use gcdp_system
      implicit none

      double precision FRANAZ
      parameter (FRANAZ=0.1d0)
      integer ng,ndm,ns,nag,nadm,nas,nagravg
      integer i,pn,nap,nnap,nagp
      double precision th,pi2,dtgc
      integer naz,listaz(0:MN-1)
      real ran1
      external ran1

      pi2=2.0d0*M_PI
 
      if(TM_tot.ge.SI_ranazt) then
        return
      endif

#ifdef GAS
! number of active gravitationally particles
      nagp=0
! number of active particles
      nap=0
! number of active particles
      nnap=0
      naz=0
! *** Gas Particles ***
      do i=0,ng-1
        pn=list_ap(i)
        if(SI_ranazi(0).eq.0.and.(id_p(pn).ge.SI_ranazi(1) &
         .and.id_p(pn).le.SI_ranazi(2))) then
          if(dble(ran1(idum)).lt.FRANAZ) then
! *** store randomly moved particle list ***
            listaz(naz)=pn
            naz=naz+1
! *** store active particle list ***
            list(nagp)=pn
            nagp=nagp+1
          else if(i.lt.nagravg) then
            list(nagp)=pn
            nagp=nagp+1
          else if(i.lt.nag) then
            talist(nap)=pn
            nap=nap+1
          else
! *** store non-active particles ***
            nalist(nnap)=pn
            nnap=nnap+1
          endif
        else if(i.lt.nagravg) then
          list(nagp)=pn
          nagp=nagp+1
        else if(i.lt.nag) then
          talist(nap)=pn
          nap=nap+1
        else
! *** store non-active particles ***
          nalist(nnap)=pn
          nnap=nnap+1
        endif
      enddo

! move particle by angle th
      do i=0,naz-1
        pn=listaz(i)
        th=pi2*ran1(idum)
! *** correction for initial kick ***
! *** overestimated dt ***
        dtgc=0.5d0*(TM_lt+TM_dt-(dtg_p(pn)+ltg_p(pn)))
! *** corrected vn?_p ***
        vnx_p(pn)=vnx_p(pn)+dtgc*ax_p(pn)
        vny_p(pn)=vny_p(pn)+dtgc*ay_p(pn)
        vnz_p(pn)=vnz_p(pn)+dtgc*az_p(pn)
        dtg_p(pn)=TM_lt+TM_dt-ltg_p(pn)
! *** rotate ***
! *** position ***
        tx(pn)=xc_p(pn)
        ty(pn)=yc_p(pn)
        xc_p(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        yc_p(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** velocity ***
        tx(pn)=vnx_p(pn)
        ty(pn)=vny_p(pn)
        vnx_p(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        vny_p(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** acceleration ***
        tx(pn)=dvx_p(pn)
        ty(pn)=dvy_p(pn)
        dvx_p(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        dvy_p(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** make it active particles ***
        dt_p(pn)=TM_lt+TM_dt-lt_p(pn)
        flagt_p(pn)=1
        dtr_p(pn)=TM_dt
      enddo
! *** update nag and list_ap(i)
      nagravg=nagp
      nag=nap+nagp
      do i=0,nagp-1
        list_ap(i)=list(i)
      enddo
      do i=0,nap-1
        list_ap(nagravg+i)=talist(i)
      enddo
      do i=0,nnap-1
        list_ap(nag+i)=nalist(i)
      enddo
#endif
#ifdef DM
      nap=0
      nnap=0
      naz=0
! *** DM Particles ***
      do i=0,ndm-1
        pn=list_adm(i)
        if(SI_ranazi(0).ne.0.and. &
         (id_dm(pn).ge.SI_ranazi(1).and.id_dm(pn).le.SI_ranazi(2))) then
          if(dble(ran1(idum)).lt.FRANAZ) then
! *** store active az list ***
            listaz(naz)=pn
            naz=naz+1
! *** store active particle list ***
            talist(nap)=pn
            nap=nap+1
          else if(i.lt.nadm) then
            talist(nap)=pn
            nap=nap+1
          else
! *** store non-active particles ***
            nalist(nnap)=pn
            nnap=nnap+1
          endif
        else if(i.lt.nadm) then
          talist(nap)=pn
          nap=nap+1
        else
! *** store non-active particles ***
          nalist(nnap)=pn
          nnap=nnap+1
        endif
      enddo

! move particle by angle th
      do i=0,naz-1
        pn=listaz(i)
        th=pi2*ran1(idum)
! *** rotate ***
! *** position ***
        tx(pn)=xc_dm(pn)
        ty(pn)=yc_dm(pn)
        xc_dm(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        yc_dm(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** velocity ***
        tx(pn)=vnx_dm(pn)
        ty(pn)=vny_dm(pn)
        vnx_dm(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        vny_dm(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** acceleration ***
        tx(pn)=dvx_dm(pn)
        ty(pn)=dvy_dm(pn)
        dvx_dm(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        dvy_dm(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** make it active particles ***
        dt_dm(pn)=TM_lt+TM_dt-lt_dm(pn)
        flagt_dm(pn)=1
        dtr_dm(pn)=TM_dt
      enddo
! *** update nadm and list_adm(i)
      nadm=nap
      do i=0,nadm-1
        list_adm(i)=talist(i)
      enddo
      do i=0,nnap-1
        list_adm(nadm+i)=nalist(i)
      enddo
#endif
#ifdef STAR
      nap=0
      nnap=0
      naz=0
! *** Star Particles ***
      do i=ng,ng+nas-1
        pn = list_ap(i)         
        pn=list_ap(i)
        if(SI_ranazi(0).eq.0.and. &
         (id_p(pn).ge.SI_ranazi(1).and.id_p(pn).le.SI_ranazi(2))) then
          if(dble(ran1(idum)).lt.FRANAZ) then
! *** store randomly moved particle list ***
            listaz(naz)=pn
            naz=naz+1
! *** store active particle list ***
            talist(nap)=pn
            nap=nap+1
          else if(i.lt.ng+nas) then
            talist(nap)=pn
            nap=nap+1
          else
! *** store non-active particles ***
            nalist(nnap)=pn
            nnap=nnap+1
          endif
        else if(i.lt.ng+nas) then
          talist(nap)=pn
          nap=nap+1
        else
! *** store non-active particles ***
          nalist(nnap)=pn
          nnap=nnap+1
        endif
      enddo
! *** random azimuthal move
      do i=0,naz-1
        pn=listaz(i)
! move particle by angle th
        th=pi2*ran1(idum)
! *** rotate ***
! *** position ***
        tx(pn)=xc_p(pn)
        ty(pn)=yc_p(pn)
        xc_p(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        yc_p(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** velocity ***
        tx(pn)=vnx_p(pn)
        ty(pn)=vny_p(pn)
        vnx_p(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        vny_p(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** acceleration ***
        tx(pn)=dvx_p(pn)
        ty(pn)=dvy_p(pn)
        dvx_p(pn)=dcos(th)*tx(pn)-dsin(th)*ty(pn)
        dvy_p(pn)=dsin(th)*tx(pn)+dcos(th)*ty(pn)
! *** make it active particles ***
        dt_p(pn)=TM_lt+TM_dt-lt_p(pn)
        flagt_p(pn)=1
        dtr_p(pn)=TM_dt
      enddo
! *** update nag and list_ap(i)
      nas=nap
      do i=0,nas-1
        list_ap(ng+i)=talist(i)
      enddo
      do i=0,nnap-1
        list_ap(ng+nas+i)=nalist(i)
      enddo
#endif

end subroutine
#endif

