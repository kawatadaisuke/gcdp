#include "gcdp.def"
! *****************************************************
!    allocate.F95 for GCD+ ver. f03.0
! 25  Feb. 2013   written by D. Kawata
! ***************************************************** 

#if defined(GAS) || defined(STAR)
subroutine allocate_baryon_all(np)
      use gcdp_baryon
#ifdef METAL
      use gcdp_metal
#endif

      implicit none

      integer,intent(in) :: np

      allocate(id_p(0:np))
      allocate(nnb_p(0:np))

      allocate(m_p(0:np))
      allocate(x_p(0:np))
      allocate(y_p(0:np))
      allocate(z_p(0:np))
      allocate(xc_p(0:np))
      allocate(yc_p(0:np))
      allocate(zc_p(0:np))
      allocate(vx_p(0:np))
      allocate(vy_p(0:np))
      allocate(vz_p(0:np))
      allocate(vnx_p(0:np))
      allocate(vny_p(0:np))
      allocate(vnz_p(0:np))
      allocate(vvnx_p(0:np))
      allocate(vvny_p(0:np))
      allocate(vvnz_p(0:np))
      allocate(u_p(0:np))
      allocate(myu_p(0:np))
      allocate(h_p(0:np))
      allocate(rho_p(0:np))
      allocate(p_p(0:np))
      allocate(cs_p(0:np))
      allocate(as_p(0:np))
      allocate(div_v_p(0:np))
      allocate(arot_v_p(0:np))
      allocate(hvsigdt_p(0:np))
      allocate(alpv_p(0:np))
      allocate(alpu_p(0:np))
      allocate(d2u_p(0:np))
      allocate(omgh_p(0:np))
      allocate(zetah_p(0:np))
      allocate(dvx_p(0:np))
      allocate(dvy_p(0:np))
      allocate(dvz_p(0:np))
      allocate(pdu_p(0:np))
      allocate(ndu_p(0:np))
      allocate(ax_p(0:np))
      allocate(ay_p(0:np))
      allocate(az_p(0:np))
      allocate(pax_p(0:np))
      allocate(pay_p(0:np))
      allocate(paz_p(0:np))
      allocate(dt_p(0:np))
      allocate(lt_p(0:np))
      allocate(vdt_p(0:np))
      allocate(dtg_p(0:np))
      allocate(ltg_p(0:np))
      allocate(dtmnb_p(0:np))
      allocate(dtr_p(0:np))
      allocate(flagc_p(0:np))
      allocate(list_ap(0:np))
      allocate(flagt_p(0:np))
#ifdef COOL
      allocate(ram_p(0:np))
#endif
#ifdef STAR
      allocate(ts_p(0:np))
#ifdef SF_EFD
      allocate(Gsn_p(0:np))
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
      allocate(flagfd_p(0:np))
      allocate(flagrfd_p(0:np))
      allocate(tnsn_p(0:np))
      allocate(tmej_p(0:np))
      allocate(tmzHe_p(0:np))
      allocate(tmzC_p(0:np))
      allocate(tmzN_p(0:np))
      allocate(tmzO_p(0:np))
      allocate(tmzNe_p(0:np))
      allocate(tmzMg_p(0:np))
      allocate(tmzSi_p(0:np))
      allocate(tmzFe_p(0:np))
      allocate(tmzZ_p(0:np))
#endif
#endif

#ifdef METAL
      allocate(mzHe_p(0:np))
      allocate(mzC_p(0:np))
      allocate(mzN_p(0:np))
      allocate(mzO_p(0:np))
      allocate(mzNe_p(0:np))
      allocate(mzMg_p(0:np))
      allocate(mzSi_p(0:np))
      allocate(mzFe_p(0:np))
      allocate(mzZ_p(0:np))

      allocate(mzHe0_p(0:np))
      allocate(mzC0_p(0:np))
      allocate(mzN0_p(0:np))
      allocate(mzO0_p(0:np))
      allocate(mzNe0_p(0:np))
      allocate(mzMg0_p(0:np))
      allocate(mzSi0_p(0:np))
      allocate(mzFe0_p(0:np))
      allocate(mzZ0_p(0:np))
#ifdef METALDIFF
      allocate(vsig_p(0:np))
      allocate(zdA_p(0:np))
      allocate(zdHeB_p(0:np))
      allocate(zdCB_p(0:np))
      allocate(zdNB_p(0:np))
      allocate(zdOB_p(0:np))
      allocate(zdNeB_p(0:np))
      allocate(zdMgB_p(0:np))
      allocate(zdSiB_p(0:np))
      allocate(zdFeB_p(0:np))
      allocate(zdZB_p(0:np))
#endif
#endif

end subroutine

! reallocate integer baryon variables
subroutine reallocate_baryon_int(np)
      use gcdp_baryon

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(id_p)
      deallocate(nnb_p)
      deallocate(flagc_p)
      deallocate(flagt_p)
      deallocate(list_ap)
#if defined(SF_EFD) || defined(SF_ZFD)
      deallocate(flagrfd_p)
      deallocate(flagfd_p)
#endif

! allocate
      allocate(id_p(0:np))
      allocate(nnb_p(0:np))
      allocate(flagc_p(0:np))
      allocate(flagt_p(0:np))
      allocate(list_ap(0:np))
#if defined(SF_EFD) || defined(SF_ZFD)
      allocate(flagfd_p(0:np))
      allocate(flagrfd_p(0:np))
#endif

end subroutine


! reallocate first set of double precision baryon variables
subroutine reallocate_baryon_d1(np)
      use gcdp_baryon

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(xc_p)
      deallocate(yc_p)
      deallocate(zc_p)
      deallocate(vnx_p)
      deallocate(vny_p)
      deallocate(vnz_p)
      deallocate(m_p)
      deallocate(u_p)
      deallocate(h_p)
      deallocate(dvx_p)
      deallocate(dvy_p)
      deallocate(dvz_p)
      deallocate(ndu_p)
      deallocate(pdu_p)
      deallocate(div_v_p)    
      deallocate(lt_p)
      deallocate(dt_p)
      deallocate(alpv_p)
      deallocate(alpu_p)
! related variables
      deallocate(x_p)
      deallocate(y_p)
      deallocate(z_p)
      deallocate(vx_p)
      deallocate(vy_p)
      deallocate(vz_p)
      deallocate(vvnx_p)
      deallocate(vvny_p)
      deallocate(vvnz_p)
      deallocate(p_p)
      deallocate(cs_p)
      deallocate(as_p)
      deallocate(hvsigdt_p)
      deallocate(pax_p)
      deallocate(pay_p)
      deallocate(paz_p)
      deallocate(vdt_p)
      deallocate(dtmnb_p)

! allocate
      allocate(xc_p(0:np))
      allocate(yc_p(0:np))
      allocate(zc_p(0:np))
      allocate(vnx_p(0:np))
      allocate(vny_p(0:np))
      allocate(vnz_p(0:np))
      allocate(m_p(0:np))
      allocate(u_p(0:np))
      allocate(h_p(0:np))
      allocate(dvx_p(0:np))
      allocate(dvy_p(0:np))
      allocate(dvz_p(0:np))
      allocate(ndu_p(0:np))
      allocate(pdu_p(0:np))
      allocate(div_v_p(0:np))    
      allocate(lt_p(0:np))
      allocate(dt_p(0:np))
      allocate(alpv_p(0:np))
      allocate(alpu_p(0:np))
! related variables
      allocate(x_p(0:np))
      allocate(y_p(0:np))
      allocate(z_p(0:np))
      allocate(vx_p(0:np))
      allocate(vy_p(0:np))
      allocate(vz_p(0:np))
      allocate(vvnx_p(0:np))
      allocate(vvny_p(0:np))
      allocate(vvnz_p(0:np))
      allocate(p_p(0:np))
      allocate(cs_p(0:np))
      allocate(as_p(0:np))
      allocate(hvsigdt_p(0:np))
      allocate(pax_p(0:np))
      allocate(pay_p(0:np))
      allocate(paz_p(0:np))
      allocate(vdt_p(0:np))
      allocate(dtmnb_p(0:np))

end subroutine


subroutine reallocate_baryon_d2(np)
      use gcdp_baryon

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(rho_p)
      deallocate(omgh_p)
      deallocate(zetah_p)
      deallocate(arot_v_p)
      deallocate(d2u_p)
      deallocate(dtr_p)
      deallocate(myu_p)
      deallocate(dtg_p)
      deallocate(ltg_p)
      deallocate(ax_p)
      deallocate(ay_p)
      deallocate(az_p)
#ifdef COOL
      deallocate(ram_p)
#endif

! allocate
      allocate(rho_p(0:np))
      allocate(omgh_p(0:np))
      allocate(zetah_p(0:np))
      allocate(arot_v_p(0:np))
      allocate(d2u_p(0:np))
      allocate(dtr_p(0:np))
      allocate(myu_p(0:np))
      allocate(dtg_p(0:np))
      allocate(ltg_p(0:np))
      allocate(ax_p(0:np))
      allocate(ay_p(0:np))
      allocate(az_p(0:np))
#ifdef COOL
      allocate(ram_p(0:np))
#endif

end subroutine

#ifdef STAR
subroutine reallocate_starfd(np)
      use gcdp_baryon

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(ts_p)
#ifdef SF_EFD
      deallocate(Gsn_p)
#endif

! allocate
      allocate(ts_p(0:np))
#ifdef SF_EFD
      allocate(Gsn_p(0:np))
#endif

end subroutine
#endif

#if defined(SF_EFD) || defined(SF_ZFD)
subroutine reallocate_tfeed(np)
      use gcdp_baryon

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(tnsn_p)
      deallocate(tmej_p)
      deallocate(tmzHe_p)
      deallocate(tmzC_p)
      deallocate(tmzN_p)
      deallocate(tmzO_p)
      deallocate(tmzNe_p)
      deallocate(tmzMg_p)
      deallocate(tmzSi_p)
      deallocate(tmzFe_p)
      deallocate(tmzZ_p)

! allocate
      allocate(tnsn_p(0:np))
      allocate(tmej_p(0:np))
      allocate(tmzHe_p(0:np))
      allocate(tmzC_p(0:np))
      allocate(tmzN_p(0:np))
      allocate(tmzO_p(0:np))
      allocate(tmzNe_p(0:np))
      allocate(tmzMg_p(0:np))
      allocate(tmzSi_p(0:np))
      allocate(tmzFe_p(0:np))
      allocate(tmzZ_p(0:np))

end subroutine
#endif

#ifdef METAL
subroutine reallocate_metal(np)
      use gcdp_metal

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(mzHe_p)
      deallocate(mzC_p)
      deallocate(mzN_p)
      deallocate(mzO_p)
      deallocate(mzNe_p)
      deallocate(mzMg_p)
      deallocate(mzSi_p)
      deallocate(mzFe_p)
      deallocate(mzZ_p)
#ifdef METALDIFF
      deallocate(vsig_p)
      deallocate(zdA_p)
      deallocate(zdHeB_p)
      deallocate(zdCB_p)
      deallocate(zdNB_p)
      deallocate(zdOB_p)
      deallocate(zdNeB_p)
      deallocate(zdMgB_p)
      deallocate(zdSiB_p)
      deallocate(zdFeB_p)
      deallocate(zdZB_p)
#endif

! allocate
      allocate(mzHe_p(0:np))
      allocate(mzC_p(0:np))
      allocate(mzN_p(0:np))
      allocate(mzO_p(0:np))
      allocate(mzNe_p(0:np))
      allocate(mzMg_p(0:np))
      allocate(mzSi_p(0:np))
      allocate(mzFe_p(0:np))
      allocate(mzZ_p(0:np))
#ifdef METALDIFF
      allocate(vsig_p(0:np))
      allocate(zdA_p(0:np))
      allocate(zdHeB_p(0:np))
      allocate(zdCB_p(0:np))
      allocate(zdNB_p(0:np))
      allocate(zdOB_p(0:np))
      allocate(zdNeB_p(0:np))
      allocate(zdMgB_p(0:np))
      allocate(zdSiB_p(0:np))
      allocate(zdFeB_p(0:np))
      allocate(zdZB_p(0:np))
#endif

end subroutine

subroutine reallocate_metal0(np)
      use gcdp_metal

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(mzHe0_p)
      deallocate(mzC0_p)
      deallocate(mzN0_p)
      deallocate(mzO0_p)
      deallocate(mzNe0_p)
      deallocate(mzMg0_p)
      deallocate(mzSi0_p)
      deallocate(mzFe0_p)
      deallocate(mzZ0_p)

! allocate
      allocate(mzHe0_p(0:np))
      allocate(mzC0_p(0:np))
      allocate(mzN0_p(0:np))
      allocate(mzO0_p(0:np))
      allocate(mzNe0_p(0:np))
      allocate(mzMg0_p(0:np))
      allocate(mzSi0_p(0:np))
      allocate(mzFe0_p(0:np))
      allocate(mzZ0_p(0:np))

end subroutine
#endif

! endif for GAS || STAR
#endif


#ifdef DM
subroutine allocate_dm_all(ndm)
      use gcdp_dm

      implicit none

      integer,intent(in) :: ndm

      allocate(id_dm(0:ndm))
      allocate(nnb_dm(0:ndm))
      allocate(m_dm(0:ndm))
      allocate(x_dm(0:ndm))
      allocate(y_dm(0:ndm))
      allocate(z_dm(0:ndm))
      allocate(xc_dm(0:ndm))
      allocate(yc_dm(0:ndm))
      allocate(zc_dm(0:ndm))
      allocate(vx_dm(0:ndm))
      allocate(vy_dm(0:ndm))
      allocate(vz_dm(0:ndm))
      allocate(vnx_dm(0:ndm))
      allocate(vny_dm(0:ndm))
      allocate(vnz_dm(0:ndm))
      allocate(dvx_dm(0:ndm))
      allocate(dvy_dm(0:ndm))
      allocate(dvz_dm(0:ndm))
      allocate(pdvx_dm(0:ndm))
      allocate(pdvy_dm(0:ndm))
      allocate(pdvz_dm(0:ndm))
      allocate(dt_dm(0:ndm))
      allocate(lt_dm(0:ndm))
      allocate(vdt_dm(0:ndm))
      allocate(dtmnb_dm(0:ndm))
      allocate(dtr_dm(0:ndm))
      allocate(h_dm(0:ndm)) 
      allocate(omgh_dm(0:ndm))
      allocate(zetah_dm(0:ndm))
      allocate(rho_dm(0:ndm))
      allocate(div_v_dm(0:ndm))
      allocate(list_adm(0:ndm))
      allocate(flagt_dm(0:ndm)) 

end subroutine

! reallocate DM integer variable
subroutine reallocate_dm_int(ndm)
      use gcdp_dm

      implicit none
      integer,intent(in) :: ndm

! deallocate
      deallocate(id_dm)
      deallocate(nnb_dm)
      deallocate(flagt_dm) 
      deallocate(list_adm)

! allocate
      allocate(id_dm(0:ndm))
      allocate(nnb_dm(0:ndm))
      allocate(flagt_dm(0:ndm)) 
      allocate(list_adm(0:ndm))

end subroutine

! reallocate DM double variables 1
subroutine reallocate_dm_d1(ndm)
      use gcdp_dm

      implicit none
      integer,intent(in) :: ndm

! deallocate
      deallocate(x_dm)
      deallocate(y_dm)
      deallocate(z_dm)
      deallocate(xc_dm)
      deallocate(yc_dm)
      deallocate(zc_dm)
      deallocate(vx_dm)
      deallocate(vy_dm)
      deallocate(vz_dm)
      deallocate(vnx_dm)
      deallocate(vny_dm)
      deallocate(vnz_dm)

      deallocate(vdt_dm)
      deallocate(dtmnb_dm)

! allocate
      allocate(x_dm(0:ndm))
      allocate(y_dm(0:ndm))
      allocate(z_dm(0:ndm))
      allocate(xc_dm(0:ndm))
      allocate(yc_dm(0:ndm))
      allocate(zc_dm(0:ndm))
      allocate(vx_dm(0:ndm))
      allocate(vy_dm(0:ndm))
      allocate(vz_dm(0:ndm))
      allocate(vnx_dm(0:ndm))
      allocate(vny_dm(0:ndm))
      allocate(vnz_dm(0:ndm))

      allocate(vdt_dm(0:ndm))
      allocate(dtmnb_dm(0:ndm))

end subroutine

! reallocate DM double variables 2
subroutine reallocate_dm_d2(ndm)
      use gcdp_dm

      implicit none
      integer,intent(in) :: ndm

! deallocate
      deallocate(m_dm)
      deallocate(dvx_dm)
      deallocate(dvy_dm)
      deallocate(dvz_dm)
      deallocate(pdvx_dm)
      deallocate(pdvy_dm)
      deallocate(pdvz_dm)
      deallocate(dt_dm)
      deallocate(lt_dm)
      deallocate(h_dm) 

! allocate
      allocate(m_dm(0:ndm))
      allocate(dvx_dm(0:ndm))
      allocate(dvy_dm(0:ndm))
      allocate(dvz_dm(0:ndm))
      allocate(pdvx_dm(0:ndm))
      allocate(pdvy_dm(0:ndm))
      allocate(pdvz_dm(0:ndm))
      allocate(dt_dm(0:ndm))
      allocate(lt_dm(0:ndm))
      allocate(h_dm(0:ndm)) 

end subroutine

! reallocate DM double variables 3
subroutine reallocate_dm_d3(ndm)
      use gcdp_dm

      implicit none
      integer,intent(in) :: ndm

! deallocate
      deallocate(dtr_dm)
      deallocate(omgh_dm)
      deallocate(zetah_dm)
      deallocate(rho_dm)
      deallocate(div_v_dm)

! allocate
      allocate(dtr_dm(0:ndm))
      allocate(omgh_dm(0:ndm))
      allocate(zetah_dm(0:ndm))
      allocate(rho_dm(0:ndm))
      allocate(div_v_dm(0:ndm))

end subroutine

! endif for DM
#endif

#if defined(GAS) || defined(STAR)
! reallocate btree
subroutine allocate_btree(ntr)
      use gcdp_btree

      implicit none
      integer,intent(in) :: ntr

      if(allocated(np_tr)) then
        deallocate(np_tr)
        deallocate(pn_tr)
        deallocate(l_tr)
        deallocate(cx_tr)
        deallocate(cy_tr)
        deallocate(cz_tr)
        deallocate(daughter_tr)
        deallocate(next_tr)
#ifdef SGRAV
        deallocate(cmx_tr)
        deallocate(cmy_tr)
        deallocate(cmz_tr)
        deallocate(mass_tr)
        deallocate(delta_tr)
        deallocate(hm_tr)
#ifdef QPTREEF
        deallocate(mx_tr)
        deallocate(my_tr)
        deallocate(mz_tr)
        deallocate(mxx_tr)
        deallocate(myy_tr)
        deallocate(mzz_tr)
        deallocate(mxy_tr)
        deallocate(myz_tr)
        deallocate(mzx_tr)
#endif
#endif
        deallocate(proc_tr)
      endif

      allocate(np_tr(0:ntr))
      allocate(pn_tr(0:ntr))
      allocate(l_tr(0:ntr))
      allocate(cx_tr(0:ntr))
      allocate(cy_tr(0:ntr))
      allocate(cz_tr(0:ntr))
      allocate(daughter_tr(0:ntr))
      allocate(next_tr(0:ntr))
#ifdef SGRAV
      allocate(cmx_tr(0:ntr))
      allocate(cmy_tr(0:ntr))
      allocate(cmz_tr(0:ntr))
      allocate(mass_tr(0:ntr))
      allocate(delta_tr(0:ntr))
      allocate(hm_tr(0:ntr))
#ifdef QPTREEF
      allocate(mx_tr(0:ntr))
      allocate(my_tr(0:ntr))
      allocate(mz_tr(0:ntr))
      allocate(mxx_tr(0:ntr))
      allocate(myy_tr(0:ntr))
      allocate(mzz_tr(0:ntr))
      allocate(mxy_tr(0:ntr))
      allocate(myz_tr(0:ntr))
      allocate(mzx_tr(0:ntr))
#endif
#endif
      allocate(proc_tr(0:ntr))

end subroutine

! end GAS || STAR
#endif







