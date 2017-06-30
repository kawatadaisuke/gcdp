#include "gcdp.def"
!************************************************
!     common.F95 for GCD+
!  08 June 2017   written by D.KAWATA
! ************************************************

#if defined(GAS) || defined(STAR)
! * for baryon Particle *
module gcdp_baryon

      implicit none
! * Particle ID *
      integer,allocatable,save :: id_p(:)
! * number neighbour *
      integer,allocatable,save :: nnb_p(:)
! * Mass *      
      double precision,allocatable,save :: m_p(:)
! * Virtual Position Xn+1,Xn *      
      double precision,allocatable,save :: x_p(:),y_p(:),z_p(:)
! * Correct position *      
      double precision,allocatable,save :: xc_p(:),yc_p(:),zc_p(:)
! * Velosity Vn-1/2,Vn+1/2 *
      double precision,allocatable,save :: vx_p(:),vy_p(:),vz_p(:)
! * Vn+1 *    
      double precision,allocatable,save :: vnx_p(:),vny_p(:),vnz_p(:)
! * Virtual Vn+1 *      
      double precision,allocatable,save :: vvnx_p(:),vvny_p(:),vvnz_p(:)
! * Internal Energy Un+1,Un * 
      double precision,allocatable,save :: u_p(:)
! * molecular weight *
      double precision,allocatable,save :: myu_p(:)
! * Smoothing Length hn+1,hn *      
      double precision,allocatable,save :: h_p(:)
! * Density rhon+1,rhon *       
      double precision,allocatable,save :: rho_p(:)
! * Presure Pn+1,Pn *      
      double precision,allocatable,save :: p_p(:)
! * Sound Velocity *      
      double precision,allocatable,save :: cs_p(:)
! * Entropic function *
      double precision,allocatable,save :: as_p(:)
! * div V *
      double precision,allocatable,save :: div_v_p(:)
! * |rot V| *      
      double precision,allocatable,save :: arot_v_p(:)
! *** signal velocity, max, for courant condition ***
      double precision,allocatable,save :: hvsigdt_p(:)
! *** artificial viscosity and thermal conductivity term ***
      double precision,allocatable,save :: alpv_p(:),alpu_p(:),d2u_p(:)
! *** variable smoothing correction ***
      double precision,allocatable,save :: omgh_p(:),zetah_p(:)
! * dv/dt *      
      double precision,allocatable,save :: dvx_p(:),dvy_p(:),dvz_p(:)
! * du/dt_n,du/dt_n+1 *      
      double precision,allocatable,save :: pdu_p(:),ndu_p(:)
! * Acceleration used in star.F for feedback *       
      double precision,allocatable,save :: ax_p(:),ay_p(:),az_p(:)
! previous time step 
      double precision,allocatable,save :: pax_p(:),pay_p(:),paz_p(:)
! *** For Individual time step ***
! * Individual time bin *      
      double precision,allocatable,save :: dt_p(:)
! * system time bin *      
      double precision,allocatable,save :: lt_p(:)
! * Virtual Time bin *      
      double precision,allocatable,save :: vdt_p(:)
! * dt, lt for gravity 
      double precision,allocatable,save :: dtg_p(:),ltg_p(:)
! * minimum dt for neighbours, and previous or required time step *      
      double precision,allocatable,save :: dtmnb_p(:),dtr_p(:)
! * flag for <=0: gas, -1 no cooling, >=1: stars *
      integer,allocatable,save :: flagc_p(:)
! ***** For Individual Time Step *****
! * list of active particles for gas and star *      
      integer,allocatable,save :: list_ap(:)
! * flag for time step changed 0: no change, otherwise changed *
      integer,allocatable,save :: flagt_p(:)
#ifdef COOL
! * for cooling *
      double precision,allocatable,save :: ram_p(:)
#endif
#ifdef STAR
! * star formtion time *
      double precision,allocatable,save :: ts_p(:)
! *** for feedback ***
#ifdef SF_EFD
! * energy of supernova *      
      double precision,allocatable,save :: Gsn_p(:)
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! flag for feedback gas
! flagfd_p 0: normal gas, >0: star, <0: feedback gas (set in main.F)
!   star with <0: particles become gas at the next step.
! flagrfd_p 0: not receiving feedback, 1: affected by feedback
!  for stars, 0: normal, 1: changing to gas in the next step.
!  after dtyields: -1 and flagfd<0: end of the mass group
!  set in calc_dv_du()
      integer,allocatable,save :: flagfd_p(:),flagrfd_p(:)
! * For Total Ejected mass by SNe unit Solar Mass *
      double precision,allocatable,save :: tnsn_p(:),tmej_p(:),tmzHe_p(:) &
       ,tmzC_p(:),tmzN_p(:),tmzO_p(:),tmzNe_p(:) &
       ,tmzMg_p(:),tmzSi_p(:),tmzFe_p(:) &
       ,tmzZ_p(:) 
#endif
#endif
! endif for STAR

#if defined(COSM) && defined(BOXSIM)
! Hubbble parameter, scale factor, and cosmic time
       double precision,allocatable,save :: asc_p(:),tc_p(:)
! for gravity timestep
       double precision,allocatable,save :: ascg_p(:),tcg_p(:)
#endif

end module gcdp_baryon


#ifdef METAL
module gcdp_metal

      implicit none
! * For Heavy Elements unit Solar Mass *
      double precision,allocatable,save :: mzHe_p(:),mzC_p(:),mzN_p(:) &
       ,mzO_p(:),mzNe_p(:),mzMg_p(:) &
       ,mzSi_p(:),mzFe_p(:),mzZ_p(:)
! original metals
      double precision,allocatable,save :: mzHe0_p(:),mzC0_p(:),mzN0_p(:) &
       ,mzO0_p(:),mzNe0_p(:),mzMg0_p(:) &
       ,mzSi0_p(:),mzFe0_p(:),mzZ0_p(:)
#ifdef METALDIFF
      double precision,allocatable,save :: vsig_p(:),zdA_p(:),zdHeB_p(:) &
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
      integer,allocatable,save :: id_dm(:)
! *** number of neighbour particles ***
      integer,allocatable,save :: nnb_dm(:)
! * mass *      
      double precision,allocatable,save :: m_dm(:)
! * Virtual Position *      
      double precision,allocatable,save :: x_dm(:),y_dm(:),z_dm(:)
! * Correct Position *      
      double precision,allocatable,save :: xc_dm(:),yc_dm(:),zc_dm(:)
! * Velocity V(n+1/2) *
      double precision,allocatable,save :: vx_dm(:),vy_dm(:),vz_dm(:)
! * Virtual Velocity V(n+1/2) *      
      double precision,allocatable,save :: vnx_dm(:),vny_dm(:),vnz_dm(:)
! * acceleration *
      double precision,allocatable,save :: dvx_dm(:),dvy_dm(:),dvz_dm(:)
! * acceleration  at previous time step *
      double precision,allocatable,save :: pdvx_dm(:),pdvy_dm(:),pdvz_dm(:)
! * for Individual time steep *
! * Individual time bin *      
      double precision,allocatable,save :: dt_dm(:)
! * time in system time step *      
      double precision,allocatable,save :: lt_dm(:)
! * virtual time bin *      
      double precision,allocatable,save :: vdt_dm(:)
! * minimum time step of neighbours *      
      double precision,allocatable,save :: dtmnb_dm(:),dtr_dm(:)
! *** for adaptive softening ***
      double precision,allocatable,save :: h_dm(:),omgh_dm(:),zetah_dm(:) &
       ,rho_dm(:),div_v_dm(:)
! * list of active particles for DM *      
      integer,allocatable,save :: list_adm(:)
! * flag for time step changed 0: no change, otherwise changed *
      integer,allocatable,save :: flagt_dm(:)
#if defined(COSM) && defined(BOXSIM)
! Hubbble parameter, scale factor, and cosmic time
       double precision,allocatable,save :: asc_dm(:),tc_dm(:)
#endif

#ifdef TESTP
! define variables for Test particle simulation
       double precision,allocatable,save :: rperi_dm(:),rapo_dm(:),zmax_dm(:)
#endif

end module gcdp_dm

#endif

#if defined(GAS) || defined(STAR)
module gcdp_btree
      use gcdp_const
      implicit none
! ***** for Tree used also in ddecb.F ****
! * Number of nodes *
      integer,save :: num_tr
! * max node id sent to the other proc 
      integer,save :: nodese_tr
! * number of contained particle *      
       integer,allocatable,save :: np_tr(:)
! * name of Particle *      
      integer,allocatable,save :: pn_tr(:)
! * length of side *       
      double precision,allocatable,save :: l_tr(:)
! * Coordinate of center *      
      double precision,allocatable,save :: cx_tr(:),cy_tr(:),cz_tr(:)
! * first child node *      
      integer,allocatable,save :: daughter_tr(:)
! * next node *      
      integer,allocatable,save :: next_tr(:)
#ifdef SGRAV
! * center of mass *      
      double precision,allocatable,save :: cmx_tr(:),cmy_tr(:),cmz_tr(:)
! * total of mass *      
      double precision,allocatable,save :: mass_tr(:)
! * distance between cm and center, maximum softening *      
      double precision,allocatable,save :: delta_tr(:),hm_tr(:)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision,allocatable,save :: mx_tr(:),my_tr(:),mz_tr(:)
      double precision,allocatable,save :: mxx_tr(:),myy_tr(:),mzz_tr(:)
      double precision,allocatable,save :: mxy_tr(:),myz_tr(:),mzx_tr(:)
#endif
#endif
! *** proc id ***
      integer,allocatable,save :: proc_tr(:)
end module gcdp_btree
#endif

#ifdef DM
module gcdp_dmtree
      use gcdp_const
      implicit none
! ***** For Dark Matter Tree *****
! * Number of nodes *
      integer,save :: num_dmtr
! * max node id sent to the other proc 
      integer,save :: nodese_dmtr,nodess_dmtr
! * number of contained particle *       
      integer,allocatable,save :: np_dmtr(:)
! * name of Particle *      
      integer,allocatable,save :: pn_dmtr(:)
! * length of side *       
      double precision,allocatable,save :: l_dmtr(:)
! * first child node *      
      integer,allocatable,save :: daughter_dmtr(:)
! * next node *      
      integer,allocatable,save :: next_dmtr(:)
! * Coordinate of center *      
      double precision,allocatable,save :: cx_dmtr(:),cy_dmtr(:),cz_dmtr(:)
! * center of mass *      
      double precision,allocatable,save :: cmx_dmtr(:),cmy_dmtr(:),cmz_dmtr(:)
! * total of mass *
      double precision,allocatable,save :: mass_dmtr(:) 
! * distance between cm and center *      
      double precision,allocatable,save :: delta_dmtr(:),hm_dmtr(:)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision,allocatable,save :: mx_dmtr(:),my_dmtr(:),mz_dmtr(:)
      double precision,allocatable,save :: mxx_dmtr(:),myy_dmtr(:),mzz_dmtr(:)
      double precision,allocatable,save :: mxy_dmtr(:),myz_dmtr(:),mzx_dmtr(:)
#endif
! *** proc id ***
      integer,allocatable,save :: proc_dmtr(:)
end module gcdp_dmtree
#endif

#ifdef MULTI_LEV
module gcdp_ldmtree
      use gcdp_const
      implicit none
! ***** For low resolution Dark Matter Tree *****
! * Number of nodes *
      integer,save :: num_ldmtr
! * max node id sent to the other proc 
      integer,save :: nodese_ldmtr,nodess_ldmtr
! * number of contained particle *   
      integer,allocatable,save :: np_ldmtr(:)
! * name of Particle *      
      integer,allocatable,save :: pn_ldmtr(:)
! * length of side *       
      double precision,allocatable,save :: l_ldmtr(:),hm_ldmtr(:)
! * first child node *      
      integer,allocatable,save :: daughter_ldmtr(:)
! * next node *      
      integer,allocatable,save :: next_ldmtr(:)
! * Coordinate of center *      
      double precision,allocatable,save :: cx_ldmtr(:),cy_ldmtr(:),cz_ldmtr(:)
! * center of mass *      
      double precision,allocatable,save :: cmx_ldmtr(:),cmy_ldmtr(:),cmz_ldmtr(:)
! * total of mass *
      double precision,allocatable,save ::  mass_ldmtr(:)
! * distance between cm and center *      
      double precision,allocatable,save ::  delta_ldmtr(:)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision,allocatable,save ::  mx_ldmtr(:),my_ldmtr(:),mz_ldmtr(:)
      double precision,allocatable,save ::  mxx_ldmtr(:),myy_ldmtr(:),mzz_ldmtr(:)
      double precision,allocatable,save ::  mxy_ldmtr(:),myz_ldmtr(:),mzx_ldmtr(:)
#endif
! *** proc id ***
      integer,allocatable,save :: proc_ldmtr(:)
end module gcdp_ldmtree
#endif

#ifdef GAS
module gcdp_gtree
      use gcdp_const
      implicit none
! ***** Tree for SPH particles _gtr *****
      integer,save :: num_gtr
! * max node id sent to the other proc 
      integer,save :: nodese_gtr,nodess_gtr
! * number of contained particle *       
      integer,allocatable,save :: np_gtr(:)
! * name of Particle *            
      integer,allocatable,save :: pn_gtr(:)
! * length of side *      
      double precision,allocatable,save :: l_gtr(:),hm_gtr(:)
      double precision,allocatable,save :: cx_gtr(:),cy_gtr(:),cz_gtr(:)
! * first child node *      
      integer,allocatable,save :: daughter_gtr(:)
! * next node *      
      integer,allocatable,save :: next_gtr(:)
! *** proc id ***
      integer,allocatable,save :: proc_gtr(:)
end module gcdp_gtree
#endif

#ifdef MEXT
module mext
      use gcdp_const
      implicit none

      integer,save :: SI_nmext
      double precision,save :: SI_dlr,SI_lri,SI_lro
      double precision,save :: rmext(0:MNMEXT-1),mextr(0:MNMEXT-1)
end module mext
#endif
 
module gcdp_system
      use gcdp_const
      implicit none

! ***** Definition of Structure *****
! * TimeStep DTn DTn-1*
      double precision,save :: TM_dt
! * End time *      
      double precision,save :: TM_t0
! * Total Time *      
      double precision,save :: TM_tot
! * System time bin *      
      double precision,save :: TM_sdt
#if defined(SF) || defined(SF_EFD) || defined(SF_ZFD)
! * Minimum local time bin *      
      double precision,save :: TM_mindtsfd
#endif
! * Minimum local time bin *      
      double precision,save :: TM_ldtmin
! * local Time in system *      
      double precision,save :: TM_lt
      integer,save :: TM_k,TM_kn
#ifdef STAR
! * used for star formation and/or feedback *
! * pv31.17 only used for feedback
      double precision,save :: TMsf_dt,TMsf_t
#endif
      
! ***** For Output Parameter *****
! * This structure have Number of output file. *
! * Make file for display by step *
      double precision,save :: OF_todt
#ifdef OFTIMESPEC
      integer,save :: OF_nof,OF_iof
      double precision,allocatable,save :: OF_time(:)
#else
      double precision,save :: OF_odt
#endif

! ***** For System Information *****
! * number of high resolution DM ***
      integer,save :: SI_ndm1,SI_ndm1t
! * flag for making output file *
      integer,save :: SI_flagout
! * number of output files *
      integer,save :: SI_nof
! *** for softening ***
      double precision,save :: SI_eps
#ifdef COSM
      double precision,save :: SI_eps0,SI_epsa
! * table for scale factor and time *
      double precision,save :: at_ttb(0:NATTABLE),ta_ttb(0:NATTABLE)
! * Cosmological Value *
      double precision,save :: SI_h0,SI_omg0,SI_lam0,SI_omgb
#ifdef BOXSIM
! Hubble parameter at that time
      double precision,save :: SI_hub
! * table for 1/a and 1/a^2
      double precision,save :: ia_ttb(0:NATTABLE),ia2_ttb(0:NATTABLE)
#endif
#endif 
#ifdef COOL
! * EoR redshift
      double precision,save :: SI_zeor
#endif
      double precision,save :: SI_a
      integer,save :: SI_flagrad
#ifdef STAR
! number of mass group, total number of star formation 
      integer,save :: SI_nsp,SI_snii,SI_tnsf
! * for memory of feedback history *
      double precision,save :: SI_esns,SI_zsns(0:NYEL-1)
#endif
      double precision,save :: SI_mms
#if defined(BOXSIM) || defined(ISOPM)
! *** box size ***
      double precision,save :: SI_lbox(0:2)
#if defined(COSM) && defined(BOXSIM)
      double precision,save :: SI_lbox0(0:2)
#endif
#endif
#ifdef FIXEDP
      integer,save :: SI_fid(0:1)
#endif
#ifdef FIXEDDMP
      integer,save :: SI_fiddm(0:1)
#endif
#ifdef TESTP
! test particle ID range, total number and number within proc
      integer,save :: SI_nttpdm,SI_ntpdm
! number of orbit particle output
      integer,save :: OF_nporb
! constant timestep
      double precision,save :: TM_cdt
#endif
#ifdef RANAZ
      integer,save :: SI_ranazi(0:2)
      double precision,save :: SI_ranazt
#ifdef RANAZRLIM
      double precision,save :: SI_ranazrlim
#endif
#endif
#ifdef ARTFDM
      integer,save :: SI_artfdm(0:1)
#endif
#ifdef DEL_DMPART
! number of particles in the component where the particles are deleted
      integer,save :: SI_ndmtdc
! delete particle ID range
      integer,save :: SI_dmdelci(0:1)
! mass of deleting component, delete rate, and particle mass
      double precision,save :: SI_mtdmdc0,SI_deldmdt,SI_mdmp0
#endif
#ifdef ADD_DMPART
! number of particles will be added, and added already
      integer,save :: SI_ndmtadded
! total mass added component, mass of DM particles added
      double precision,save :: SI_mdmpadd
! mass addition rate
      double precision,save :: SI_dmdtadd
! adding component, radial scale length (constant) 
      double precision,save :: SI_hradd,SI_rmaxadd
! velocity disperion, scale length, sig_0, for vR, vth, vz
      double precision,save :: SI_sigr0add,SI_hsigradd
      double precision,save :: SI_sigph0add,SI_hsigphadd
      double precision,save :: SI_sigz0add,SI_hsigzadd
! particle mass and scale length
      double precision,save :: SI_mdmp0,SI_hdmp0
#endif

#ifdef BCX_FIX
      integer,save :: BCX_ndx
      double precision,save :: BCX_rhox(0:1),BCX_px(0:1),BCX_vx(0:1)
#ifdef MHD
      double precision,save :: BCX_bx(0:1),BCX_by(0:1),BCX_bz(0:1)
      double precision,save :: BCX_vy(0:1),BCX_vz(0:1)
#endif
#endif

#ifdef TREEPM
      double precision,save :: SI_rstpm,SI_rsp2tpm,SI_rshtpm
#endif

! *** for counting walltime ***
      double precision,save :: WTM_ini,WTM_out,WTM_ddecb,WTM_ddecdm &
      ,WTM_cont,WTM_setval,WTM_calcvu,WTM_tree,WTM_starf &
      ,WTM_starfd,WTM_setvaldm,WTM_time &
      ,WTM_setvals,WTM_setvalfd,WTM_feedb,WTM_bf,WTM_dmf &
      ,WTM_update,WTM_updateu
#ifdef TREEPM
      double precision,save :: WTM_pm
#endif

! **** for Random Numbrer ***
      integer,save :: idum      

! ***** for MPI *****
      integer,save :: nprocs,myrank

#ifdef EXTF
! for external potential
#ifdef EXTF_Miyamoto_Nagai
! Miyamoto Nagai parameter
      double precision,save :: MN_mdisk,MN_a,MN_b
#endif
#ifdef EXTF_NFW
! NFW halo parameters
      double precision,save :: NFW_rho0,NFW_rs
#endif
#ifdef EXTF_Hernquist
! Hernquist bulge parameter
      double precision,save :: HQ_mtot,HQ_ra
#endif
#endif

end module gcdp_system

module gcdp_kernel
      use gcdp_const
      implicit none
! *** kernel look up table ****
      double precision,save :: dnktab
      double precision,save :: s_tb(0:NKTAB),w_tb(0:NKTAB),dwds_s_tb(0:NKTAB) &
       ,dwdsc_tb(0:NKTAB)
#ifdef SGRAV
      double precision,save :: dphidr_r_tb(0:NKTAB),dphidh_tb(0:NKTAB)
      double precision,save :: d2phidr2_tb(0:NKTAB),d3phidr3_tb(0:NKTAB)
#ifdef TREEPM
      double precision,save :: dnftpmtab
      double precision,save :: fshtpm_tb(0:NFTPMTAB),rstpm_tb(0:NFTPMTAB)
#endif
#endif
end module gcdp_kernel
 
#if defined(SF_EFD) || defined(SF_ZFD)
module gcdp_yields
      use gcdp_const
      implicit none
! ***** For Yields lookup table *****
      double precision,save :: z_ytb(0:NYTZ+2),t_ytb(0:NYTT+1) &
       ,mej_ytb(0:NYTT+1,0:NYTZ+2),nsn_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzHe_ytb(0:NYTT+1,0:NYTZ+2),mzC_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzN_ytb(0:NYTT+1,0:NYTZ+2),mzO_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzNe_ytb(0:NYTT+1,0:NYTZ+2),mzMg_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzSi_ytb(0:NYTT+1,0:NYTZ+2), mzFe_ytb(0:NYTT+1,0:NYTZ+2) &
       ,mzZ_ytb(0:NYTT+1,0:NYTZ+2),nsw_ytb(0:NYTT+1,0:NYTZ+2)
! *** for stars mass group ***
      integer,save :: nsp_ytb(0:NYTZ+2)
! *** starting time ***
      double precision,save :: tspi_ytb(0:MNSPYTT-1,0:NYTZ+2)
end module gcdp_yields
#endif

#if defined(COOL) || defined(BBAT)
module gcdp_cool
      use gcdp_const
      implicit none
! ***** for metal cooling and heating ***
      integer,save :: nz_crtb,nnh_crtb,nmet_crtb,nt_crtb
      double precision,save :: z_crtb(0:MNT_CRTB-1),nh_crtb(0:MNNH_CRTB) &
      ,met_crtb(0:MNMET_CRTB),t_crtb(0:MNT_CRTB)
      double precision,save :: cr_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,hr_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,myu_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,ne_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB)
end module gcdp_cool
#endif

! look up table for mass profiles for adding particles.
#ifdef ADD_EXPDISK
module gcdp_expdisk
      use gcdp_const
      implicit none
! scale length, delta ln x, ln minimum x point
      double precision,save :: hd_expd,dlx_expd,lxmin_expd
      double precision,save :: lx_expd(0:NEXPDTAB-1),mx_expd(0:NEXPDTAB-1)
end module gcdp_expdisk
#endif

! ***** MHD modules ***

#ifdef MHD
module gcdp_mhd
      implicit none
      double precision,allocatable,save :: bx_p(:),by_p(:),bz_p(:)
      double precision,allocatable,save :: pbx_p(:),pby_p(:),pbz_p(:)
      double precision,allocatable,save :: dbx_p(:),dby_p(:),dbz_p(:) &
       ,pdbx_p(:),pdby_p(:),pdbz_p(:)
      double precision,allocatable :: divb_p(:),arotb_p(:),agradb_p(:)
      double precision,allocatable :: alpb_p(:),divberr_p(:)
! not used
      double precision,allocatable,save :: valfi_p(:),pb_p(:)
! not used
#ifdef BBAT
      double precision,allocatable,save :: chi_p(:),lognh_p(:),logne_p(:)
#endif
#ifdef DCLEAN
      double precision,allocatable :: psi_p(:),ppsi_p(:)
      double precision,allocatable :: dpsi_p(:),pdpsi_p(:)
      double precision,allocatable :: maxvmhd_p(:)
#endif
end module gcdp_mhd
#endif

! ***** PM modules *****
#ifdef TREEPM
module gcdp_pm
      implicit none
! number of mesh
      integer,save :: nx_m,ny_m,nz_m,nt_m
! number/2 of mesh
      integer,save :: nx2_m,ny2_m,nz2_m
! size of mesh and box
      double precision,save :: dv_m,vt_m
! mesh spaces in real space
      double precision,save :: dx_m,dy_m,dz_m
! mesh spaces in k space
      double precision,save :: dkx_m,dky_m,dkz_m
! mesh 0 points
      double precision,save :: x0_m,y0_m,z0_m
! mesh 1D coordinates
      double precision,allocatable,save :: x_m1d(:),y_m1d(:),z_m1d(:)
! density
      double precision,allocatable,save :: rho_m(:,:,:)
! potential
      double precision,allocatable,save :: pot_m(:,:,:)
! force
      double precision,allocatable,save :: fx_m(:,:,:),fy_m(:,:,:),fz_m(:,:,:)

#ifdef FFTW3_MPI
      integer,save :: lnz_m,lzoff_m
      integer,allocatable,save :: lnz_mp(:),lzoff_mp(:)
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
      integer(C_INTPTR_T),save :: nx_fftw,ny_fftw,nz_fftw
      integer(C_INTPTR_T),save :: lnz_fftw,lkoff_fftw
      integer(C_INTPTR_T),save :: lny_fftw,ljoff_fftw
#endif

! forward and backward FFT plan
      type(C_PTR),save :: planf_fftw,planb_fftw
! data space
      real(C_DOUBLE),pointer,save :: rsd_fftw(:,:,:)
      complex(C_DOUBLE_COMPLEX),pointer,save :: ksd_fftw(:,:,:)
      type(C_PTR),save :: poid_fftw
#ifdef ISOPM
! forward plan
      type(C_PTR),save :: planfgpot_fftw
! to store fourier component of 1/r potential kernel
! note that rsd and ksd are allocated to the same memory space (Sec. 7.2 p.79
! of FFTW3 manual)
      real(C_DOUBLE),pointer,save :: rsdgpot_fftw(:,:,:)
      complex(C_DOUBLE_COMPLEX),pointer,save :: ksdgpot_fftw(:,:,:)
      type(C_PTR),save :: poidgpot_fftw
#endif


end module fftw3
#endif



