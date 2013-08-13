#include "gcdp.def"
! *********************************************
!  mesh_setrho.F95 for GCD+
!  13 Aug. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
subroutine  mesh_setrho(np,ndm)
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

      integer,intent(in) :: np,ndm

      integer i,j,k,ip0,jp0,kp0,ip1,jp1,kp1,w1,w2
      double precision df,vm
      double precision wx0,wx1,wy0,wy1,wz0,wz1

! initialisation
      do k=0,nz_m-1
        do j=0,ny_m-1
          do i=0,nx_m-1
            rho_m(i,j,k)=0.0d0
          enddo
        enddo
      enddo

! assign baryon mass
#if defined(GAS) || defined(STAR)
      do i=0,np-1
! x
        ip0=int((x_p(i)-x0_m)/dx_m)
        ip1=ip0+1
        if(ip0.lt.0) then
          ip0=nx_m
          ip1=0
          df=(x_m1d(0)-x_p(i))/dx_m
          wx0=1.0d0-df
          wx1=df
        else if(ip1.ge.nx_m) then
          ip1=0
          df=(x_p(i)-x_m1d(ip0))/dx_m
          wx0=df
          wx1=1.0-df
        else
          df=(x_p(i)-x_m1d(ip0))/dx_m
          wx0=df
          wx1=1.0-df
        endif
! for check
        if((wx0.lt.0.0d0.or.wx0.gt.1.0d0) &
          .or.(wx1.lt.0.0d0.or.wx1.gt.1.0d0)) then
          write(6,*) ' Error in mesh_setrho() for baryon: wx0,wx1=',wx0,wx1
          write(6,*) ' x,ip0,ip1=',x_p(i),ip0,ip1
          stop
        endif
! y
        jp0=int((y_p(i)-y0_m)/dy_m)
        jp1=jp0+1
        if(jp0.lt.0) then
          jp0=ny_m
          jp1=0
          df=(y_m1d(0)-y_p(i))/dy_m
          wy0=1.0d0-df
          wy1=df
        else if(jp1.ge.ny_m) then
          jp1=0
          df=(y_p(i)-y_m1d(jp0))/dy_m
          wy0=df
          wy1=1.0-df
        else
          df=(y_p(i)-y_m1d(ip0))/dy_m
          wy0=df
          wy1=1.0-df
        endif
! for check
        if((wy0.lt.0.0d0.or.wy0.gt.1.0d0) &
          .or.(wy1.lt.0.0d0.or.wy1.gt.1.0d0)) then
          write(6,*) ' Error in mesh_setrho() for baryon: wy0,wy1=',wy0,wy1
          write(6,*) ' y,jp0,jp1=',y_p(i),jp0,jp1
          stop
        endif
! z
        kp0=int((z_p(i)-z0_m)/dz_m)
        kp1=kp0+1
        if(kp0.lt.0) then
          kp0=nz_m
          kp1=0
          df=(z_m1d(0)-z_p(i))/dz_m
          wz0=1.0d0-df
          wz1=df
        else if(kp1.ge.nz_m) then
          kp1=0
          df=(z_p(i)-z_m1d(kp0))/dz_m
          wz0=df
          wz1=1.0-df
        else
          df=(z_p(i)-z_m1d(kp0))/dz_m
          wz0=df
          wz1=1.0-df
        endif
! for check
        if((wz0.lt.0.0d0.or.wz0.gt.1.0d0) &
          .or.(wz1.lt.0.0d0.or.wz1.gt.1.0d0)) then
          write(6,*) ' Error in mesh_setrho() for baryon: wz0,wz1=',wz0,wz1
          write(6,*) ' z,kp0,kp1=',z_p(i),kp0,kp1
          stop
        endif

! assign the mass
! kp0
        rho_m(ip0,jp0,kp0)=rho_m(ip0,jp0,kp0) &
         +wx0*wy0*wz0*m_p(i)
        rho_m(ip1,jp0,kp0)=rho_m(ip1,jp0,kp0) &
         +wx1*wy0*wz0*m_p(i)
        rho_m(ip1,jp1,kp0)=rho_m(ip1,jp1,kp0) &
         +wx1*wy1*wz0*m_p(i)
        rho_m(ip0,jp1,kp0)=rho_m(ip0,jp1,kp0) &
         +wx0*wy1*wz0*m_p(i)
! kp1
        rho_m(ip0,jp0,kp1)=rho_m(ip0,jp0,kp1) &
         +wx0*wy0*wz1*m_p(i)
        rho_m(ip1,jp0,kp1)=rho_m(ip1,jp0,kp1) &
         +wx1*wy0*wz1*m_p(i)
        rho_m(ip1,jp1,kp1)=rho_m(ip1,jp1,kp1) &
         +wx1*wy1*wz1*m_p(i)
        rho_m(ip0,jp1,kp1)=rho_m(ip0,jp1,kp1) &
         +wx0*wy1*wz1*m_p(i)

      enddo
#endif

! assign DM mass
#ifdef DM
      do i=0,ndm-1
! x
        ip0=int((x_dm(i)-x0_m)/dx_m)
        ip1=ip0+1
        if(ip0.lt.0) then
          ip0=nx_m
          ip1=0
          df=(x_m1d(0)-x_dm(i))/dx_m
          wx0=1.0d0-df
          wx1=df
        else if(ip1.ge.nx_m) then
          ip1=0
          df=(x_dm(i)-x_m1d(ip0))/dx_m
          wx0=df
          wx1=1.0-df
        else
          df=(x_dm(i)-x_m1d(ip0))/dx_m
          wx0=df
          wx1=1.0-df
        endif
! for check
        if((wx0.lt.0.0d0.or.wx0.gt.1.0d0) &
          .or.(wx1.lt.0.0d0.or.wx1.gt.1.0d0)) then
          write(6,*) ' Error in mesh_setrho() for DM: wx0,wx1=',wx0,wx1
          write(6,*) ' x_dm,ip0,ip1=',x_dm(i),ip0,ip1
          stop
        endif
! y
        jp0=int((y_dm(i)-y0_m)/dy_m)
        jp1=jp0+1
        if(jp0.lt.0) then
          jp0=ny_m
          jp1=0
          df=(y_m1d(0)-y_dm(i))/dy_m
          wy0=1.0d0-df
          wy1=df
        else if(jp1.ge.ny_m) then
          jp1=0
          df=(y_dm(i)-y_m1d(jp0))/dy_m
          wy0=df
          wy1=1.0-df
        else
          df=(y_dm(i)-y_m1d(ip0))/dy_m
          wy0=df
          wy1=1.0-df
        endif
! for check
        if((wy0.lt.0.0d0.or.wy0.gt.1.0d0) &
          .or.(wy1.lt.0.0d0.or.wy1.gt.1.0d0)) then
          write(6,*) ' Error in mesh_setrho() for DM: wy0,wy1=',wy0,wy1
          write(6,*) ' y_dm,jp0,jp1=',y_dm(i),jp0,jp1
          stop
        endif
! z
        kp0=int((z_dm(i)-z0_m)/dz_m)
        kp1=kp0+1
        if(kp0.lt.0) then
          kp0=nz_m
          kp1=0
          df=(z_m1d(0)-z_dm(i))/dz_m
          wz0=1.0d0-df
          wz1=df
        else if(kp1.ge.nz_m) then
          kp1=0
          df=(z_dm(i)-z_m1d(kp0))/dz_m
          wz0=df
          wz1=1.0-df
        else
          df=(z_dm(i)-z_m1d(kp0))/dz_m
          wz0=df
          wz1=1.0-df
        endif
! for check
        if((wz0.lt.0.0d0.or.wz0.gt.1.0d0) &
          .or.(wz1.lt.0.0d0.or.wz1.gt.1.0d0)) then
          write(6,*) ' Error in mesh_setrho() for DM: wz0,wz1=',wz0,wz1
          write(6,*) ' z_dm,kp0,kp1=',z_dm(i),kp0,kp1
          stop
        endif

! assign the density
! kp0
        rho_m(ip0,jp0,kp0)=rho_m(ip0,jp0,kp0) &
         +wx0*wy0*wz0*m_dm(i)
        rho_m(ip1,jp0,kp0)=rho_m(ip1,jp0,kp0) &
         +wx1*wy0*wz0*m_dm(i)
        rho_m(ip1,jp1,kp0)=rho_m(ip1,jp1,kp0) &
         +wx1*wy1*wz0*m_dm(i)
        rho_m(ip0,jp1,kp0)=rho_m(ip0,jp1,kp0) &
         +wx0*wy1*wz0*m_dm(i)
! kp1
        rho_m(ip0,jp0,kp1)=rho_m(ip0,jp0,kp1) &
         +wx0*wy0*wz1*m_dm(i)
        rho_m(ip1,jp0,kp1)=rho_m(ip1,jp0,kp1) &
         +wx1*wy0*wz1*m_dm(i)
        rho_m(ip1,jp1,kp1)=rho_m(ip1,jp1,kp1) &
         +wx1*wy1*wz1*m_dm(i)
        rho_m(ip0,jp1,kp1)=rho_m(ip0,jp1,kp1) &
         +wx0*wy1*wz1*m_dm(i)

      enddo
#endif

! volume of each mesh
      vm=dx_m*dy_m*dz_m
! density
      do k=0,nz_m-1
        do j=0,ny_m-1
          do i=0,nx_m-1
            rho_m(i,j,k)=rho_m(i,j,k)/vm
          enddo
        enddo
      enddo

end subroutine
#endif





