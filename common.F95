#include "gcdp.def"
!************************************************
!     common.F95 for GCD+
!  22 Jan. 2013   written by D.KAWATA
! ************************************************

#if defined(GAS) || defined(STAR)
! * for baryon Particle *
module gcdp_baryon

      implicit none
! * Particle ID *
      integer,allocatable :: id_p(:)
! * number neighbour *
      integer,allocatable :: nnb_p(:)
! * Mass *      
      double precision,allocatable :: m_p(:)
! * Virtual Position Xn+1,Xn *      
      double precision,allocatable :: x_p(:),y_p(:),z_p(:)
! * Correct position *      
      double precision,allocatable :: xc_p(:),yc_p(:),zc_p(:)
! * Velosity Vn-1/2,Vn+1/2 *
      double precision,allocatable :: vx_p(:),vy_p(:),vz_p(:)
! * Vn+1 *    
      double precision,allocatable :: vnx_p(:),vny_p(:),vnz_p(:)
! * Virtual Vn+1 *      
      double precision,allocatable :: vvnx_p(:),vvny_p(:),vvnz_p(:)
! * Internal Energy Un+1,Un * 
      double precision,allocatable :: u_p(:)
! * molecular weight *
      double precision,allocatable :: myu_p(:)
! * Smoothing Length hn+1,hn *      
      double precision,allocatable :: h_p(:)
! * Density rhon+1,rhon *       
      double precision,allocatable :: rho_p(:)
! * Presure Pn+1,Pn *      
      double precision,allocatable :: p_p(:)
! * Sound Velocity *      
      double precision,allocatable :: cs_p(:)
! * Entropic function *
      double precision,allocatable :: as_p(:)
! * div V *
      double precision,allocatable :: div_v_p(:)
! * |rot V| *      
      double precision,allocatable :: arot_v_p(:)
! *** signal velocity, max, for courant condition ***
      double precision,allocatable :: hvsigdt_p(:)
! *** artificial viscosity and thermal conductivity term ***
      double precision,allocatable :: alpv_p(:),alpu_p(:),d2u_p(:)
! *** variable smoothing correction ***
      double precision,allocatable :: omgh_p(:),zetah_p(:)
! * dv/dt *      
      double precision,allocatable :: dvx_p(:),dvy_p(:),dvz_p(:)
! * du/dt_n,du/dt_n+1 *      
      double precision,allocatable :: pdu_p(:),ndu_p(:)
! * Acceleration used in star.F for feedback *       
      double precision,allocatable :: ax_p(:),ay_p(:),az_p(:)
! previous time step 
      double precision,allocatable :: pax_p(:),pay_p(:),paz_p(:)
! *** For Individual time step ***
! * Individual time bin *      
      double precision,allocatable :: dt_p(:)
! * system time bin *      
      double precision,allocatable :: lt_p(:)
! * Virtual Time bin *      
      double precision,allocatable :: vdt_p(:)
! * dt, lt for gravity 
      double precision,allocatable :: dtg_p(:),ltg_p(:)
! * minimum dt for neighbours, and previous or required time step *      
      double precision,allocatable :: dtmnb_p(:),dtr_p(:)
! * flag for <=0: gas, -1 no cooling, >=1: stars *
      integer,allocatable :: flagc_p(:)
! ***** For Individual Time Step *****
! * list of active particles for gas and star *      
      integer,allocatable :: list_ap(:)
! * flag for time step changed 0: no change, otherwise changed *
      integer,allocatable :: flagt_p(:)
#ifdef COOL
! * for cooling *
      double precision,allocatable :: ram_p(:)
#endif
#ifdef STAR
! * star formtion time *
      double precision,allocatable :: ts_p(:)
! *** for feedback ***
#ifdef SF_EFD
! * energy of supernova *      
      double precision,allocatable :: Gsn_p(:)
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! flag for feedback gas
! flagfd_p 0: normal gas, >0: star, <0: feedback gas (set in main.F)
!   star with <0: particles become gas at the next step.
! flagrfd_p 0: not receiving feedback, 1: affected by feedback
!  for stars, 0: normal, 1: changing to gas in the next step.
!  after dtyields: -1 and flagfd<0: end of the mass group
!  set in calc_dv_du()
      integer,allocatable :: flagfd_p(:),flagrfd_p(:)
! * For Total Ejected mass by SNe unit Solar Mass *
      double precision,allocatable :: tnsn_p(:),tmej_p(:),tmzHe_p(:) &
       ,tmzC_p(:),tmzN_p(:),tmzO_p(:),tmzNe_p(:) &
       ,tmzMg_p(:),tmzSi_p(:),tmzFe_p(:) &
       ,tmzZ_p(:) 
#endif

#endif
!
end module gcdp_baryon


#ifdef METAL
module gcdp_metal

      implicit none
! * For Heavy Elements unit Solar Mass *
      double precision,allocatable :: mzHe_p(:),mzC_p(:),mzN_p(:) &
       ,mzO_p(:),mzNe_p(:),mzMg_p(:) &
       ,mzSi_p(:),mzFe_p(:),mzZ_p(:)
! original metals
      double precision,allocatable :: mzHe0_p(:),mzC0_p(:),mzN0_p(:) &
       ,mzO0_p(:),mzNe0_p(:),mzMg0_p(:) &
       ,mzSi0_p(:),mzFe0_p(:),mzZ0_p(:)
#ifdef METALDIFF
      double precision,allocatable :: vsig_p(:),zdA_p(:),zdHeB_p(:) &
       ,zdCB_p(:),zdNB_p(:),zdOB_p(:) &
       ,zdNeB_p(:),zdMgB_p(:),zdSiB_p(:) &
       ,zdFeB_p(:),zdZB_p(:)
#endif
end module gcdp_metal
#endif

! end ifdef GAS
#endif

#ifdef DM
module gcdp_dm
      use gcdp_const
      implicit none
! ***** For Dark Matter *****
! * Particle ID *
      integer,allocatable :: id_dm(:)
! *** number of neighbour particles ***
      integer,allocatable :: nnb_dm(:)
! * mass *      
      double precision,allocatable :: m_dm(:)
! * Virtual Position *      
      double precision,allocatable :: x_dm(:),y_dm(:),z_dm(:)
! * Correct Position *      
      double precision,allocatable :: xc_dm(:),yc_dm(:),zc_dm(:)
! * Velocity V(n+1/2) *
      double precision,allocatable :: vx_dm(:),vy_dm(:),vz_dm(:)
! * Virtual Velocity V(n+1/2) *      
      double precision,allocatable :: vnx_dm(:),vny_dm(:),vnz_dm(:)
! * acceleration *
      double precision,allocatable :: dvx_dm(:),dvy_dm(:),dvz_dm(:)
! * acceleration  at previous time step *
      double precision,allocatable :: pdvx_dm(:),pdvy_dm(:),pdvz_dm(:)
! * for Individual time steep *
! * Individual time bin *      
      double precision,allocatable :: dt_dm(:)
! * time in system time step *      
      double precision,allocatable :: lt_dm(:)
! * virtual time bin *      
      double precision,allocatable :: vdt_dm(:)
! * minimum time step of neighbours *      
      double precision,allocatable :: dtmnb_dm(:),dtr_dm(:)
! *** for adaptive softening ***
      double precision,allocatable :: h_dm(:),omgh_dm(:),zetah_dm(:) &
       ,rho_dm(:),div_v_dm(:)
! * list of active particles for DM *      
      integer,allocatable :: list_adm(:)
! * flag for time step changed 0: no change, otherwise changed *
      integer,allocatable :: flagt_dm(:)
#if defined(COSM) && defined(BOXSIM)
! Hubbble parameter, scale factor, and cosmic time
       double precision,allocatable :: hub_dm(:),asc_dm(:),tc_dm(:)
#endif

end module gcdp_dm
#endif

#if defined(GAS) || defined(STAR)
module gcdp_btree
      use gcdp_const
      implicit none
! ***** for Tree used also in ddecb.F ****
! * Number of nodes *
      integer num_tr
! * max node id sent to the other proc 
      integer nodese_tr
! * number of contained particle *      
       integer,allocatable :: np_tr(:)
! * name of Particle *      
      integer,allocatable :: pn_tr(:)
! * length of side *       
      double precision,allocatable :: l_tr(:)
! * Coordinate of center *      
      double precision,allocatable :: cx_tr(:),cy_tr(:),cz_tr(:)
! * first child node *      
      integer,allocatable :: daughter_tr(:)
! * next node *      
      integer,allocatable :: next_tr(:)
#ifdef SGRAV
! * center of mass *      
      double precision,allocatable :: cmx_tr(:),cmy_tr(:),cmz_tr(:)
! * total of mass *      
      double precision,allocatable :: mass_tr(:)
! * distance between cm and center, maximum softening *      
      double precision,allocatable :: delta_tr(:),hm_tr(:)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision,allocatable :: mx_tr(:),my_tr(:),mz_tr(:)
      double precision,allocatable :: mxx_tr(:),myy_tr(:),mzz_tr(:)
      double precision,allocatable :: mxy_tr(:),myz_tr(:),mzx_tr(:)
#endif
#endif
! *** proc id ***
      integer,allocatable :: proc_tr(:)
end module gcdp_btree
#endif

#ifdef DM
module gcdp_dmtree
      use gcdp_const
      implicit none
! ***** For Dark Matter Tree *****
! * Number of nodes *
      integer num_dmtr
! * max node id sent to the other proc 
      integer nodese_dmtr,nodess_dmtr
! * number of contained particle *       
      integer,allocatable :: np_dmtr(:)
! * name of Particle *      
      integer,allocatable :: pn_dmtr(:)
! * length of side *       
      double precision,allocatable :: l_dmtr(:)
! * first child node *      
      integer,allocatable :: daughter_dmtr(:)
! * next node *      
      integer,allocatable :: next_dmtr(:)
! * Coordinate of center *      
      double precision,allocatable :: cx_dmtr(:),cy_dmtr(:),cz_dmtr(:)
! * center of mass *      
      double precision,allocatable :: cmx_dmtr(:),cmy_dmtr(:),cmz_dmtr(:)
! * total of mass *
      double precision,allocatable :: mass_dmtr(:) 
! * distance between cm and center *      
      double precision,allocatable :: delta_dmtr(:),hm_dmtr(:)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision,allocatable :: mx_dmtr(:),my_dmtr(:),mz_dmtr(:)
      double precision,allocatable :: mxx_dmtr(:),myy_dmtr(:),mzz_dmtr(:)
      double precision,allocatable :: mxy_dmtr(:),myz_dmtr(:),mzx_dmtr(:)
#endif
! *** proc id ***
      integer,allocatable :: proc_dmtr(:)
end module gcdp_dmtree
#endif

#ifdef MULTI_LEV
module gcdp_ldmtree
      use gcdp_const
      implicit none
! ***** For low resolution Dark Matter Tree *****
! * Number of nodes *
      integer num_ldmtr
! * max node id sent to the other proc 
      integer nodese_ldmtr,nodess_ldmtr
! * number of contained particle *   
      integer,allocatable :: np_ldmtr(:)
! * name of Particle *      
      integer,allocatable :: pn_ldmtr(:)
! * length of side *       
      double precision,allocatable :: l_ldmtr(:),hm_ldmtr(:)
! * first child node *      
      integer,allocatable :: daughter_ldmtr(:)
! * next node *      
      integer,allocatable :: next_ldmtr(:)
! * Coordinate of center *      
      double precision,allocatable :: cx_ldmtr(:),cy_ldmtr(:),cz_ldmtr(:)
! * center of mass *      
      double precision,allocatable :: cmx_ldmtr(:),cmy_ldmtr(:),cmz_ldmtr(:)
! * total of mass *
      double precision,allocatable ::  mass_ldmtr(:)
! * distance between cm and center *      
      double precision,allocatable ::  delta_ldmtr(:)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision,allocatable ::  mx_ldmtr(:),my_ldmtr(:),mz_ldmtr(:)
      double precision,allocatable ::  mxx_ldmtr(:),myy_ldmtr(:),mzz_ldmtr(:)
      double precision,allocatable ::  mxy_ldmtr(:),myz_ldmtr(:),mzx_ldmtr(:)
#endif
! *** proc id ***
      integer,allocatable :: proc_ldmtr(:)
end module gcdp_ldmtree
#endif

#ifdef GAS
module gcdp_gtree
      use gcdp_const
      implicit none
! ***** Tree for SPH particles _gtr *****
      integer num_gtr
! * max node id sent to the other proc 
      integer nodese_gtr,nodess_gtr
! * number of contained particle *       
      integer,allocatable :: np_gtr(:)
! * name of Particle *            
      integer,allocatable :: pn_gtr(:)
! * length of side *      
      double precision,allocatable :: l_gtr(:),hm_gtr(:)
      double precision,allocatable :: cx_gtr(:),cy_gtr(:),cz_gtr(:)
! * first child node *      
      integer,allocatable :: daughter_gtr(:)
! * next node *      
      integer,allocatable :: next_gtr(:)
! *** proc id ***
      integer,allocatable :: proc_gtr(:)
end module gcdp_gtree
#endif

#ifdef MEXT
module mext
      use gcdp_const
      implicit none

      integer SI_nmext
      common SI_nmext
      double precision SI_dlr,SI_lri,SI_lro
      common SI_dlr,SI_lri,SI_lro
      double precision rmext,mextr
      common rmext(0:MNMEXT-1),mextr(0:MNMEXT-1)
end module mext
#endif
 
module gcdp_system
      use gcdp_const
      implicit none

! ***** Definition of Structure *****
! * TimeStep DTn DTn-1*
      double precision TM_dt
! * End time *      
      double precision TM_t0
! * Total Time *      
      double precision TM_tot
! * System time bin *      
      double precision TM_sdt
#if defined(SF) || defined(SF_EFD) || defined(SF_ZFD)
! * Minimum local time bin *      
      double precision TM_mindtsfd
#endif
! * Minimum local time bin *      
      double precision TM_ldtmin
! * local Time in system *      
      double precision TM_lt
      integer TM_k,TM_kn
#ifdef STAR
! * used for star formation and/or feedback *
! * pv31.17 only used for feedback
      double precision TMsf_dt,TMsf_t
#endif
      
! ***** For Output Parameter *****
! * This structure have Number of output file. *
! * Make file for display by step *
      double precision OF_odt,OF_todt
! ***** For System Information *****
! * number of high resolution DM ***
      integer SI_ndm1,SI_ndm1t
! * flag for making output file *
      integer SI_flagout
! * number of output files *
      integer SI_nof
! *** for softening ***
      double precision SI_eps
#ifdef COSM
      double precision SI_eps0,SI_epsa
! * table for scale factor and time *
      double precision at_ttb(0:NATTABLE),ta_ttb(0:NATTABLE)
! * Cosmological Value *
      double precision SI_h0,SI_omg0,SI_lam0,SI_omgb
#ifdef BOXSIM
! Hubble parameter at that time
      double precision SI_hub
! * table for 1/a and 1/a^2
      double precision ia_ttb(0:NATTABLE),ia2_ttb(0:NATTABLE)
#endif
#endif 
#ifdef COOL
! * EoR redshift
      double precision SI_zeor
#endif
      double precision SI_a
      integer SI_flagrad
#ifdef STAR
! number of mass group, total number of star formation 
      integer SI_nsp,SI_snii,SI_tnsf
! * for memory of feedback history *
      double precision SI_esns,SI_zsns(0:NYEL-1)
#endif
      double precision SI_mms
#ifdef BOXSIM
! *** box size ***
      double precision SI_lbox(0:2)
#if defined(COSM) && defined(BOXSIM)
      double precision SI_lbox0(0:2)
#endif
#endif
#ifdef FIXEDP
      integer SI_fid(0:1)
#endif
#ifdef FIXEDDMP
      integer SI_fiddm(0:1)
#endif
#ifdef RANAZ
      integer SI_ranazi(0:2)
      double precision SI_ranazt
#endif
#ifdef ARTFDM
      integer SI_artfdm(0:1)
#endif

#ifdef BCX_FIX
      integer BCX_ndx
      double precision BCX_rhox(0:1),BCX_px(0:1),BCX_vx(0:1)
#ifdef MHD
      double precision BCX_bx(0:1),BCX_by(0:1),BCX_bz(0:1)
      double precision BCX_vy(0:1),BCX_vz(0:1)
#endif
#endif

#ifdef TREEPM
      double precision SI_rstpm,SI_rsp2tpm,SI_rshtpm
#endif

! *** for counting walltime ***
      double precision WTM_ini,WTM_out,WTM_ddecb,WTM_ddecdm &
      ,WTM_cont,WTM_setval,WTM_calcvu,WTM_tree,WTM_starf &
      ,WTM_starfd,WTM_setvaldm,WTM_time &
      ,WTM_setvals,WTM_setvalfd,WTM_feedb,WTM_bf,WTM_dmf &
      ,WTM_update,WTM_updateu

! **** for Random Numbrer ***
      integer idum      

! ***** for MPI *****
      integer nprocs,myrank

end module gcdp_system

module gcdp_kernel
      use gcdp_const
      implicit none
! *** kernel look up table ****
      double precision ds_tb
      double precision s_tb(0:NKTAB),w_tb(0:NKTAB),dwds_s_tb(0:NKTAB) &
       ,dwdsc_tb(0:NKTAB)
#ifdef SGRAV
      double precision dphidr_r_tb(0:NKTAB),dphidh_tb(0:NKTAB)
      double precision d2phidr2_tb(0:NKTAB),d3phidr3_tb(0:NKTAB)
#endif
end module gcdp_kernel
 
#if defined(SF_EFD) || defined(SF_ZFD)
module gcdp_yields
      use gcdp_const
      implicit none
! ***** For Yields lookup table *****
      double precision z_ytb(0:NYTZ+2),t_ytb(0:NYTT+1) &
       ,mej_ytb(0:NYTT+1,0:NYTZ+2),nsn_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzHe_ytb(0:NYTT+1,0:NYTZ+2),mzC_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzN_ytb(0:NYTT+1,0:NYTZ+2),mzO_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzNe_ytb(0:NYTT+1,0:NYTZ+2),mzMg_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzSi_ytb(0:NYTT+1,0:NYTZ+2), mzFe_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzZ_ytb(0:NYTT+1,0:NYTZ+2),nsw_ytb(0:NYTT+1,0:NYTZ+2)
! *** for stars mass group ***
      integer nsp_ytb(0:NYTZ+2)
! *** starting time ***
      double precision tspi_ytb(0:MNSPYTT-1,0:NYTZ+2)
end module gcdp_yields
#endif

#if defined(COOL) || defined(BBAT)
module gcdp_cool
      use gcdp_const
      implicit none
! ***** for metal cooling and heating ***
      integer nz_crtb,nnh_crtb,nmet_crtb,nt_crtb
      double precision z_crtb(0:MNT_CRTB-1),nh_crtb(0:MNNH_CRTB) &
      ,met_crtb(0:MNMET_CRTB),t_crtb(0:MNT_CRTB)
      double precision cr_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,hr_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,myu_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,ne_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB)
end module gcdp_cool
#endif

! ***** MHD modules ***

#ifdef MHD
module gcdp_mhd
      implicit none
      double precision,allocatable :: bx_p(:),by_p(:),bz_p(:)
      double precision,allocatable :: pbx_p(:),pby_p(:),pbz_p(:)
      double precision,allocatable :: dbx_p(:),dby_p(:),dbz_p(:) &
       ,pdbx_p(:),pdby_p(:),pdbz_p(:)
      double precision,allocatable :: divb_p(:),arotb_p(:)
      double precision,allocatable :: rotbx_p(:),rotby_p(:),rotbz_p(:)
      double precision,allocatable :: alpb_p(:),divberr_p(:)
! not used
      double precision,allocatable :: valfi_p(:),pb_p(:)
! not used
#ifdef BBAT
      double precision,allocatable :: chi_p(:),lognh_p(:),logne_p(:)
#endif
end module gcdp_mhd
#endif

! ***** PM modules *****
#ifdef TREEPM
module gcdp_pm
      implicit none
! number of mesh
      integer nx_m,ny_m,nz_m,nt_m
! number/2 of mesh
      integer nx2_m,ny2_m,nz2_m
! size of mesh and box
      double precision dv_m,vt_m
! mesh spaces in real space
      double precision dx_m,dy_m,dz_m
! mesh spaces in k space
      double precision dkx_m,dky_m,dkz_m
! mesh 0 points
      double precision x0_m,y0_m,z0_m
! mesh 1D coordinates
      double precision,allocatable :: x_m1d(:),y_m1d(:),z_m1d(:)
! density
      double precision,allocatable :: rho_m(:,:,:)
! potential
      double precision,allocatable :: pot_m(:,:,:)
! force
      double precision,allocatable :: fx_m(:,:,:),fy_m(:,:,:),fz_m(:,:,:)

#ifdef FFTW3_MPI
      integer lnz_m,lzoff_m
      integer,allocatable :: lnz_mp(:),lzoff_mp(:)
#endif

end module gcdp_pm
#endif

#ifdef FFTW3
module fftw3
      use,intrinsic :: iso_c_binding
      implicit none
#ifdef FFTW3_MPI
      include 'fftw3-mpi.f03'
#else
      include 'fftw3.f03'
#endif

#ifdef FFTW3_MPI
      integer(C_INTPTR_T) :: nx_fftw,ny_fftw,nz_fftw
      integer(C_INTPTR_T) :: lnz_fftw,lkoff_fftw
      integer(C_INTPTR_T) :: lny_fftw,ljoff_fftw
#endif

! forward and backward FFT plan
      type(C_PTR) :: planf_fftw,planb_fftw
! data space
      real(C_DOUBLE),pointer :: rsd_fftw(:,:,:)
      complex(C_DOUBLE_COMPLEX),pointer :: ksd_fftw(:,:,:)
      type(C_PTR) :: poid_fftw

end module fftw3
#endif



