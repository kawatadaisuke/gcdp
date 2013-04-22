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
#ifdef MHD
      use gcdp_mhd
#endif

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
#ifdef MHD
          vny_p(i)=BCX_vy(0)
          vnz_p(i)=BCX_vz(0)
#else
          vny_p(i)=0.0d0
          vnz_p(i)=0.0d0
#endif
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
#ifdef MHD
          bx_p(i)=BCX_bx(0)
          by_p(i)=BCX_by(0)
          bz_p(i)=BCX_bz(0)
          pbx_p(i)=BCX_bx(0)
          pby_p(i)=BCX_by(0)
          pbz_p(i)=BCX_bz(0)
          dbx_p(i)=0.0d0
          dby_p(i)=0.0d0
          dbz_p(i)=0.0d0
          pdbx_p(i)=0.0d0
          pdby_p(i)=0.0d0
          pdbz_p(i)=0.0d0
          divb_p(i)=0.0d0
          rotbx_p(i)=0.0d0
          rotby_p(i)=0.0d0
          rotbz_p(i)=0.0d0       
#endif
         else if(npt-id_p(i).le.BCX_ndx) then 
          nfix=nfix+1
!          x_p(i)=xc_p(i)
          y_p(i)=yc_p(i)
          z_p(i)=zc_p(i)
          vnx_p(i)=BCX_vx(1)
#ifdef MHD
          vny_p(i)=BCX_vy(1)
          vnz_p(i)=BCX_vz(1)
#else
          vny_p(i)=0.0d0
          vnz_p(i)=0.0d0
#endif
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
#ifdef MHD
          bx_p(i)=BCX_bx(1)
          by_p(i)=BCX_by(1)
          bz_p(i)=BCX_bz(1)
          pbx_p(i)=BCX_bx(1)
          pby_p(i)=BCX_by(1)
          pbz_p(i)=BCX_bz(1)
          dbx_p(i)=0.0d0
          dby_p(i)=0.0d0
          dbz_p(i)=0.0d0
          pdbx_p(i)=0.0d0
          pdby_p(i)=0.0d0
          pdbz_p(i)=0.0d0
          divb_p(i)=0.0d0
          rotbx_p(i)=0.0d0
          rotby_p(i)=0.0d0
          rotbz_p(i)=0.0d0          
#endif
!          write(60,*) xc_p(i),yc_p(i),zc_p(i)
        endif
      enddo
!      write(6,*) ' nfix=',nfix
#endif
end subroutine
#endif
