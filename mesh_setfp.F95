#include "gcdp.def"
! *********************************************
!  mesh_setrho.F95 for GCD+
!  13 Aug. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
subroutine  mesh_setfp(nag,nas,nadm)
      use gcdp_const
      use gcdp_system
      use gcdp_pm
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif
#ifdef DM
      use gcdp_dm
#endif

      implicit none
      include 'mpif.h'

      integer,intent(in) :: nag,nas,nadm

      integer pn,i,j,k,ip0,jp0,kp0,ip1,jp1,kp1,w1,w2
      double precision df,vm
      double precision wx0,wx1,wy0,wy1,wz0,wz1
! *** for test file *** 
      character fileo*60

! for gas
#ifdef GAS
      do i=0,nag-1
        pn=list_ap(i)
! x
        if(x_p(pn).le.x0_m) then
          ip0=nx_m-1
          ip1=0
          df=(x_m1d(0)-x_p(pn))/dx_m
          wx0=df
          wx1=1.0d0-df
        else if(x_p(pn).ge.x_m1d(nx_m-1)) then
          ip0=nx_m-1
          ip1=0
          df=(x_p(pn)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        else
          ip0=int((x_p(pn)-x0_m)/dx_m)
          ip1=ip0+1
          df=(x_p(pn)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        endif
! y
        if(y_p(pn).le.y0_m) then
          jp0=ny_m-1
          jp1=0
          df=(y_m1d(0)-y_p(pn))/dy_m
          wy0=df
          wy1=1.0d0-df
        else if(y_p(pn).ge.y_m1d(ny_m-1)) then
          jp0=ny_m-1
          jp1=0
          df=(y_p(pn)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        else
          jp0=int((y_p(pn)-y0_m)/dy_m)
          jp1=jp0+1
          df=(y_p(pn)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        endif
! z
        if(z_p(pn).le.z0_m) then
          kp0=nz_m-1
          kp1=0
          df=(z_m1d(0)-z_p(pn))/dz_m
          wz0=df
          wz1=1.0d0-df
        else if(z_p(pn).ge.z_m1d(nz_m-1)) then
          kp0=nz_m-1
          kp1=0
          df=(z_p(pn)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        else
          kp0=int((z_p(pn)-z0_m)/dz_m)
          kp1=kp0+1
          df=(z_p(pn)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        endif

! assign force to the particles
! dvx_p kp0
        dvx_p(pn)=dvx_p(pn)+wx0*wy0*wz0*fx_m(ip0,jp0,kp0) &
         +wx1*wy0*wz0*fx_m(ip1,jp0,kp0) &
         +wx1*wy1*wz0*fx_m(ip1,jp1,kp0) &
         +wx0*wy1*wz0*fx_m(ip0,jp1,kp0) &
! kp1
         +wx0*wy0*wz1*fx_m(ip0,jp0,kp1) &
         +wx1*wy0*wz1*fx_m(ip1,jp0,kp1) &
         +wx1*wy1*wz1*fx_m(ip1,jp1,kp1) &
         +wx0*wy1*wz1*fx_m(ip0,jp1,kp1)
! dvy_p kp0
        dvy_p(pn)=dvy_p(pn)+wx0*wy0*wz0*fy_m(ip0,jp0,kp0) &
         +wx1*wy0*wz0*fy_m(ip1,jp0,kp0) &
         +wx1*wy1*wz0*fy_m(ip1,jp1,kp0) &
         +wx0*wy1*wz0*fy_m(ip0,jp1,kp0) &
! kp1
         +wx0*wy0*wz1*fy_m(ip0,jp0,kp1) &
         +wx1*wy0*wz1*fy_m(ip1,jp0,kp1) &
         +wx1*wy1*wz1*fy_m(ip1,jp1,kp1) &
         +wx0*wy1*wz1*fy_m(ip0,jp1,kp1)
! dvz_p kp0
        dvz_p(pn)=dvz_p(pn)+wx0*wy0*wz0*fz_m(ip0,jp0,kp0) &
         +wx1*wy0*wz0*fz_m(ip1,jp0,kp0) &
         +wx1*wy1*wz0*fz_m(ip1,jp1,kp0) &
         +wx0*wy1*wz0*fz_m(ip0,jp1,kp0) &
! kp1
         +wx0*wy0*wz1*fz_m(ip0,jp0,kp1) &
         +wx1*wy0*wz1*fz_m(ip1,jp0,kp1) &
         +wx1*wy1*wz1*fz_m(ip1,jp1,kp1) &
         +wx0*wy1*wz1*fz_m(ip0,jp1,kp1)
         

      enddo

! test output
      write(fileo,'(a2,i3.3)') 'fp',myrank
      open(60,file=fileo,status='unknown')
      do i=0,nag-1
        if(x_p(i)**2+y_p(i)**2+z_p(i)**2.gt.0.0d0) then
        write(60,'(8(1pE13.5))') x_p(i),y_p(i),z_p(i) &
         ,dvx_p(i),dvy_p(i),dvz_p(i) &
         ,dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) &
         ,-(dvx_p(i)*x_p(i)+dvy_p(i)*y_p(i)+dvz_p(i)*z_p(i)) &
         /dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) 
        else
        write(60,'(8(1pE13.5))') x_p(i),y_p(i),z_p(i) &
         ,dvx_p(i),dvy_p(i),dvz_p(i) &
         ,dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) &
         ,-(dvx_p(i)*x_p(i)+dvy_p(i)*y_p(i)+dvz_p(i)*z_p(i))
        endif
      enddo
      close(60)
      stop

#endif


! DM
#ifdef DM
      do i=0,nadm-1
        pn=list_adm(i)
! x
        if(x_dm(pn).le.x0_m) then
          ip0=nx_m-1
          ip1=0
          df=(x_m1d(0)-x_dm(pn))/dx_m
          wx0=df
          wx1=1.0d0-df
        else if(x_dm(pn).ge.x_m1d(nx_m-1)) then
          ip0=nx_m-1
          ip1=0
          df=(x_dm(pn)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        else
          ip0=int((x_dm(pn)-x0_m)/dx_m)
          ip1=ip0+1
          df=(x_dm(pn)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        endif
! y
        if(y_dm(pn).le.y0_m) then
          jp0=ny_m-1
          jp1=0
          df=(y_m1d(0)-y_dm(pn))/dy_m
          wy0=df
          wy1=1.0d0-df
        else if(y_dm(pn).ge.y_m1d(ny_m-1)) then
          jp0=ny_m-1
          jp1=0
          df=(y_dm(pn)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        else
          jp0=int((y_dm(pn)-y0_m)/dy_m)
          jp1=jp0+1
          df=(y_dm(pn)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        endif
! z
        if(z_dm(pn).le.z0_m) then
          kp0=nz_m-1
          kp1=0
          df=(z_m1d(0)-z_dm(pn))/dz_m
          wz0=df
          wz1=1.0d0-df
        else if(z_dm(pn).ge.z_m1d(nz_m-1)) then
          kp0=nz_m-1
          kp1=0
          df=(z_dm(pn)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        else
          kp0=int((z_dm(pn)-z0_m)/dz_m)
          kp1=kp0+1
          df=(z_dm(pn)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        endif

! assign force to the particles
! dvx_p kp0
        dvx_dm(pn)=dvx_dm(pn)+wx0*wy0*wz0*fx_m(ip0,jp0,kp0) &
         +wx1*wy0*wz0*fx_m(ip1,jp0,kp0) &
         +wx1*wy1*wz0*fx_m(ip1,jp1,kp0) &
         +wx0*wy1*wz0*fx_m(ip0,jp1,kp0) &
! kp1
         +wx0*wy0*wz1*fx_m(ip0,jp0,kp1) &
         +wx1*wy0*wz1*fx_m(ip1,jp0,kp1) &
         +wx1*wy1*wz1*fx_m(ip1,jp1,kp1) &
         +wx0*wy1*wz1*fx_m(ip0,jp1,kp1)
! dvy_p kp0
        dvy_dm(pn)=dvy_dm(pn)+wx0*wy0*wz0*fy_m(ip0,jp0,kp0) &
         +wx1*wy0*wz0*fy_m(ip1,jp0,kp0) &
         +wx1*wy1*wz0*fy_m(ip1,jp1,kp0) &
         +wx0*wy1*wz0*fy_m(ip0,jp1,kp0) &
! kp1
         +wx0*wy0*wz1*fy_m(ip0,jp0,kp1) &
         +wx1*wy0*wz1*fy_m(ip1,jp0,kp1) &
         +wx1*wy1*wz1*fy_m(ip1,jp1,kp1) &
         +wx0*wy1*wz1*fy_m(ip0,jp1,kp1)
! dvz_p kp0
        dvz_dm(pn)=dvz_dm(pn)+wx0*wy0*wz0*fz_m(ip0,jp0,kp0) &
         +wx1*wy0*wz0*fz_m(ip1,jp0,kp0) &
         +wx1*wy1*wz0*fz_m(ip1,jp1,kp0) &
         +wx0*wy1*wz0*fz_m(ip0,jp1,kp0) &
! kp1
         +wx0*wy0*wz1*fz_m(ip0,jp0,kp1) &
         +wx1*wy0*wz1*fz_m(ip1,jp0,kp1) &
         +wx1*wy1*wz1*fz_m(ip1,jp1,kp1) &
         +wx0*wy1*wz1*fz_m(ip0,jp1,kp1)

      enddo


#endif


end subroutine
#endif





