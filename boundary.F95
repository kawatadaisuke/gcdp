#include "gcdp.def"
! ******************************************
!    boundary.F95 
!  04 Feb. 2013    produced by D.KAWATA
! ******************************************
#ifdef BCX_FIX 
SUBROUTINE boundary_xfix(np,npt)
#ifdef GAS
      use gcdp_const
      use gcdp_baryon
      use gcdp_system

      integer np,i,nfix,npt

      nfix=0
!      open(60,file='fix.dat',status='unknown')
      do i=0,np-1  
        if(id_p(i).lt.BCX_ndx) then
          nfix=nfix+1
!          x_p(i)=xc_p(i)
          y_p(i)=yc_p(i)
          z_p(i)=zc_p(i)
          vnx_p(i)=BCX_vx(0)
          vny_p(i)=0.0d0
          vnz_p(i)=0.0d0
          vx_p(i)=vnx_p(i)
          vy_p(i)=vny_p(i)
          vz_p(i)=vnz_p(i)
          dvx_p(i)=0.0d0
          dvy_p(i)=0.0d0
          dvz_p(i)=0.0d0
          ndu_p(i)=0.0d0
          pdu_p(i)=0.0d0
          rho_p(i)=BCX_rhox(0)
          p_p(i)=BCX_px(0)
          u_p(i)=p_p(i)/((gam-1.0d0)*rho_p(i))
          cs_p(i)=dsqrt(gam*p_p(i)/rho_p(i))
          as_p(i)=p_p(i)/(rho_p(i)**GAM)
!          write(60,*) xc_p(i),yc_p(i),zc_p(i)
!        else if(xc_p(i).gt.0.5d0*SI_lbox(0)-BCX_dx) then
         else if(npt-id_p(i).le.BCX_ndx) then 
          nfix=nfix+1
!          x_p(i)=xc_p(i)
          y_p(i)=yc_p(i)
          z_p(i)=zc_p(i)
          vnx_p(i)=BCX_vx(1)
          vny_p(i)=0.0d0
          vnz_p(i)=0.0d0
          vx_p(i)=vnx_p(i)
          vy_p(i)=vny_p(i)
          vz_p(i)=vnz_p(i)
          dvx_p(i)=0.0d0
          dvy_p(i)=0.0d0
          dvz_p(i)=0.0d0
          rho_p(i)=BCX_rhox(1)
          p_p(i)=BCX_px(1)
          u_p(i)=p_p(i)/((gam-1.0d0)*rho_p(i))
          cs_p(i)=dsqrt(gam*p_p(i)/rho_p(i))
          as_p(i)=p_p(i)/(rho_p(i)**GAM)
!          write(60,*) xc_p(i),yc_p(i),zc_p(i)
        endif
      enddo
!      write(6,*) ' nfix=',nfix
#endif
end subroutine
#endif
