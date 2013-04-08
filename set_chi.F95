#include "gcdp.def"
#if defined(MHD) && defined(BBAT)

! *********************************************
!  setcool.F95 for GCD+ ver. f03.0
!  17 April 2012    written by D. Barnes 
! *********************************************
! ***********************************************************
!  set ionization fraction for battery source term
! ************************************************************

SUBROUTINE setchi(np,plist,current_z)
      use gcdp_const
      use gcdp_baryon
#ifdef METAL
      use gcdp_metal
#endif
      use gcdp_cool
      use gcdp_mhd
      use gcdp_system

      implicit none

      integer,intent(in) :: np,plist(0:np)
      double precision,intent(in) :: current_z

      integer i
      double precision cradp,hradp,radp
      double precision lowmet
      double precision nhp,logTp,metp,logmetp,logTlimK
      double precision dnh,dmet,dt,logchi
      double precision wt1,wt2,wd1,wd2,wz1,wz2,dz
      double precision wm1,wm2
      integer pn,it,id,iz,im

! lowest metallicity
      lowmet=10.0d0**met_crtb(0)
! temperature limit at K
      logTlimK=dlog10(CLIMIT*TUK)
! set delta nH, met, T, ne
      dnh=nh_crtb(2)-nh_crtb(1)
      dmet=met_crtb(2)-met_crtb(1)
      dt=t_crtb(2)-t_crtb(1)

! find redshift 
      if(current_z.lt.z_crtb(nz_crtb-1)) then
        do i=1,nz_crtb-1
          if(z_crtb(i).gt.current_z) then
            iz = i
            goto 70
          endif
        enddo
        if(myrank.eq.0) then
          write(6,*) ' cool(): cannot find z =',current_z
          write(6,*) ' it must be smaller than ',z_crtb(nz_crtb-1)
        endif
        iz = nz_crtb-1
   70   dz = z_crtb(iz)-z_crtb(iz-1)
! *** weight for redshift ***
        wz1 = (z_crtb(iz)-current_z)/dz
        wz2 = (current_z-z_crtb(iz-1))/dz
      else
        iz = nz_crtb-1
        wz1 = 0.0d0
        wz2 = 1.0d0
      endif

      do i =0,np-1
        pn = plist(i)

! *** nH ***
#ifdef METAL
        nhp = ((m_p(pn)-((mzZ_p(pn)+mzHe_p(pn))/MUSM))/m_p(pn))*rho_p(pn)*(DU/MP)
#else
! Hydrogen density
!        nhp = 0.76d0*rho_p(pn)*(DU/MP)
! Proton Density
        nhp = rho_p(pn)*(DU/MP)
#endif
        lognh_p(pn)=dlog10(nhp)

! *** Find Metallicity ***
#ifdef METAL
        metp = (mzZ_p(pn))/(m_p(pn)*MUSM*XZSOL)
#else
        metp = 0.0d0
#endif
        if(metp.gt.lowmet) then
          logmetp=dlog10(metp)
        else
          logmetp=met_crtb(0)
        endif
        im = int((logmetp-met_crtb(0))/dmet)
        if(im.lt.0) then
! *** use the vale at met_crtb(0) ***
          im = 0
          wm1=1.0d0
          wm2=0.0d0
        else if(im.ge.nmet_crtb-1) then
! *** use the vale at met_crtb(nmet_crtb-1) ***
          im = nmet_crtb-2
          wm1=0.0d0
          wm2=1.0d0     
        else
          wm1 = (met_crtb(im+1)-logmetp)/dmet
          wm2 = (logmetp-met_crtb(im))/dmet
        endif
        id=int((lognh_p(pn)-nh_crtb(0))/dnh)
        if(id.lt.0) then
! *** use the vale at nh_crtb(0) ***
          id = 0
        else if(id.ge.nnh_crtb-1) then
! *** use the vale at met_crtb(nmet_crtb-1) ***
          id = nnh_crtb-2
        endif
        wd1 = (nh_crtb(id+1)-lognh_p(pn))/dnh
        wd2 = (lognh_p(pn)-nh_crtb(id))/dnh
! *** find density, allow extrapolation ***
        if(SI_flagrad.gt.-2) then
! *** before EoR use cooling rate lognh=0.0***
          id = int((0.0d0-nh_crtb(0))/dnh)
          wd1=(nh_crtb(id+1)-0.0)/dnh
          wd2=(0.0d0-nh_crtb(id))/dnh
        endif

! *** Temperature Weight ***
        logTp = dlog10((p_p(pn)*myu_p(pn)/(rho_p(pn)*TPRHO*MYU))*TUK)
        if(logTp.lt.logTlimK) then
          logTp=logTlimK
        endif
        it = int((logTp-t_crtb(0))/dt)
        if(it.lt.0) then
          it = 0
        else if(it.ge.nt_crtb-1) then
          it = nt_crtb-2
        endif
        wt1 = (t_crtb(it+1)-logTp)/dt
        wt2 = (logTp-t_crtb(it))/dt
! *** if low and high temperature, no extrapolation ***
        if(logTp.lt.t_crtb(0)) then
          it = 0
          wt1=1.0d0
          wt2=0.0d0
        else if(logTp.gt.t_crtb(nt_crtb-1)) then
          it = nt_crtb-2
          wt1=0.0d0
          wt2=1.0d0
        endif
! *** if low and high density, no extrapolation ***
        if(lognh_p(pn).lt.nh_crtb(0)) then
          id=0
          wd1=1.0d0
          wd2=0.0d0
        else if(lognh_p(pn).gt.nh_crtb(nnh_crtb-1)) then
          id = nnh_crtb-2
          wd1=0.0d0
          wd2=1.0d0
        endif

! *** Interpolate log(ne) from file ***
        logne_p(pn) &
! iz-1
!   id
!     im
!       it
         =wt1*wm1*wd1*wz1*(ne_crtb(it,im,id,iz-1)) &
!       it+1 
         +wt2*wm1*wd1*wz1*(ne_crtb(it+1,im,id,iz-1)) &
!     im+1
!       it
         +wt1*wm2*wd1*wz1*(ne_crtb(it,im+1,id,iz-1)) &
!       it+1
         +wt2*wm2*wd1*wz1*(ne_crtb(it+1,im+1,id,iz-1)) &
!   id+1
!     im
!       it
         +wt1*wm1*wd2*wz1*(ne_crtb(it,im,id+1,iz-1)) &
!       it+1
         +wt2*wm1*wd2*wz1*(ne_crtb(it+1,im,id+1,iz-1)) &
!     im+1
!       it
         +wt1*wm2*wd2*wz1*(ne_crtb(it,im+1,id+1,iz-1)) &
!       it+1
         +wt2*wm2*wd2*wz1*(ne_crtb(it+1,im+1,id+1,iz-1)) &
! iz
!   id
!     im
!       it
         +wt1*wm1*wd1*wz2*(ne_crtb(it,im,id,iz)) &
!       it+1
         +wt2*wm1*wd1*wz2*(ne_crtb(it+1,im,id,iz)) &
!     im+1
!       it
         +wt1*wm2*wd1*wz2*(ne_crtb(it,im+1,id,iz)) &
!       it+1
         +wt2*wm2*wd1*wz2*(ne_crtb(it+1,im+1,id,iz)) &
!   id+1
!     im
!       it
         +wt1*wm1*wd2*wz2*(ne_crtb(it,im,id+1,iz)) &
!       it+1
         +wt2*wm1*wd2*wz2*(ne_crtb(it+1,im,id+1,iz)) &
!     im+1
!       it
         +wt1*wm2*wd2*wz2*(ne_crtb(it,im+1,id+1,iz)) &
!       it+1
         +wt2*wm2*wd2*wz2*(ne_crtb(it+1,im+1,id+1,iz))

! *** Calc Chi for particle ***
        logchi=logne_p(pn)/lognh_p(pn)
        chi_p(pn)=10.0d0**logchi
      enddo
      return

end subroutine
#endif
