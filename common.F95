#include "gcdp.def"
!************************************************
!     common.F95 for GCD+
!  22 Jan. 2013   written by D.KAWATA
! ************************************************

#if defined(GAS) || defined(STAR)
! * for baryon Particle *
module gcdp_baryon
      use gcdp_const
      implicit none
! * Particle ID *
      integer id_p(0:MNB-1)
! * number neighbour *
      integer nnb_p(0:MNB-1)
! * Mass *      
      double precision m_p(0:MNB-1)
! * Virtual Position Xn+1,Xn *      
      double precision x_p(0:MNB-1),y_p(0:MNB-1),z_p(0:MNB-1)
! * Correct position *      
      double precision xc_p(0:MNB-1),yc_p(0:MNB-1),zc_p(0:MNB-1)
! * Velosity Vn-1/2,Vn+1/2 *
      double precision vx_p(0:MNB-1),vy_p(0:MNB-1),vz_p(0:MNB-1)
! * Vn+1 *    
      double precision vnx_p(0:MNB-1),vny_p(0:MNB-1),vnz_p(0:MNB-1)
! * Virtual Vn+1 *      
      double precision vvnx_p(0:MNB-1),vvny_p(0:MNB-1),vvnz_p(0:MNB-1)
! * Internal Energy Un+1,Un * 
      double precision u_p(0:MNB-1)
! * molecular weight *
      double precision myu_p(MNB-1)
! * Smoothing Length hn+1,hn *      
      double precision h_p(0:MNB-1)
! * Density rhon+1,rhon *       
      double precision rho_p(0:MNB-1)
! * Presure Pn+1,Pn *      
      double precision p_p(0:MNB-1)
! * Sound Velocity *      
      double precision cs_p(0:MNB-1)
! * Entropic function *
      double precision as_p(0:MNB-1)
! * div V *
      double precision div_v_p(0:MNB-1)
! * |rot V| *      
      double precision arot_v_p(0:MNB-1)
! *** signal velocity, max, for courant condition ***
      double precision hvsigdt_p(0:MNB-1)
! *** artificial viscosity and thermal conductivity term ***
      double precision alpv_p(0:MNB-1),alpu_p(0:MNB-1),d2u_p(0:MNB-1)
! *** variable smoothing correction ***
      double precision omgh_p(0:MNB-1),zetah_p(0:MNB-1)
! * dv/dt *      
      double precision dvx_p(0:MNB-1),dvy_p(0:MNB-1),dvz_p(0:MNB-1)
! * du/dt_n,du/dt_n+1 *      
      double precision pdu_p(0:MNB-1),ndu_p(0:MNB-1)
! * Acceleration used in star.F for feedback *       
      double precision ax_p(0:MNB-1),ay_p(0:MNB-1),az_p(0:MNB-1)
! previous time step 
      double precision pax_p(0:MNB-1),pay_p(0:MNB-1),paz_p(0:MNB-1)
! *** For Individual time step ***
! * Individual time bin *      
      double precision dt_p(0:MNB-1)
! * system time bin *      
      double precision lt_p(0:MNB-1)
! * Virtual Time bin *      
      double precision vdt_p(0:MNB-1)
! * dt, lt for gravity 
      double precision dtg_p(0:MNB-1),ltg_p(0:MNB-1)
! * minimum dt for neighbours, and previous or required time step *      
      double precision dtmnb_p(0:MNB-1),dtr_p(0:MNB-1)
! * flag for <=0: gas, -1 no cooling, >=1: stars *
      integer flagc_p(0:MNB-1)
! ***** For Individual Time Step *****
! * list of active particles for gas and star *      
      integer list_ap(0:MNB-1)
! * flag for time step changed 0: no change, otherwise changed *
      integer flagt_p(0:MNB-1)
#ifdef COOL
! * for cooling *
      double precision ram_p(0:MNB-1)
#endif
#ifdef STAR
! * star formtion time *
      double precision ts_p(0:MNB-1)
! *** for feedback ***
#ifdef SF_EFD
! * energy of supernova *      
      double precision Gsn_p(0:MNB-1)
#endif
#if defined(SF_EFD) || defined(SF_ZFD)
! flag for feedback gas
! flagfd_p 0: normal gas, >0: star, <0: feedback gas (set in main.F)
!   star with <0: particles become gas at the next step.
! flagrfd_p 0: not receiving feedback, 1: affected by feedback
!  for stars, 0: normal, 1: changing to gas in the next step.
!  after dtyields: -1 and flagfd<0: end of the mass group
!  set in calc_dv_du()
      integer flagfd_p(0:MNB-1),flagrfd_p(0:MNB-1)
! * For Total Ejected mass by SNe unit Solar Mass *
      double precision tnsn_p(0:MNB-1),tmej_p(0:MNB-1),tmzHe_p(0:MNB-1) &
       ,tmzC_p(0:MNB-1),tmzN_p(0:MNB-1),tmzO_p(0:MNB-1),tmzNe_p(0:MNB-1) &
       ,tmzMg_p(0:MNB-1),tmzSi_p(0:MNB-1),tmzFe_p(0:MNB-1) &
       ,tmzZ_p(0:MNB-1) 
#endif
#endif
end module gcdp_baryon

#ifdef METAL
module gcdp_metal
      use gcdp_const
      implicit none
! * For Heavy Elements unit Solar Mass *
      double precision mzHe_p(0:MNB-1),mzC_p(0:MNB-1),mzN_p(0:MNB-1) &
       ,mzO_p(0:MNB-1),mzNe_p(0:MNB-1),mzMg_p(0:MNB-1) &
       ,mzSi_p(0:MNB-1),mzFe_p(0:MNB-1),mzZ_p(0:MNB-1)
! original metals
      double precision mzHe0_p(0:MNB-1),mzC0_p(0:MNB-1),mzN0_p(0:MNB-1) &
       ,mzO0_p(0:MNB-1),mzNe0_p(0:MNB-1),mzMg0_p(0:MNB-1) &
       ,mzSi0_p(0:MNB-1),mzFe0_p(0:MNB-1),mzZ0_p(0:MNB-1)
#ifdef METALDIFF
      double precision vsig_p(0:MNB-1),zdA_p(0:MNB-1),zdHeB_p(0:MNB-1) &
       ,zdCB_p(0:MNB-1),zdNB_p(0:MNB-1),zdOB_p(0:MNB-1) &
       ,zdNeB_p(0:MNB-1),zdMgB_p(0:MNB-1),zdSiB_p(0:MNB-1) &
       ,zdFeB_p(0:MNB-1),zdZB_p(0:MNB-1)
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
      integer id_dm(0:MNDM-1)
! *** number of neighbour particles ***
      integer nnb_dm(0:MNDM-1)
! * mass *      
      double precision m_dm(0:MNDM-1)
! * Virtual Position *      
      double precision x_dm(0:MNDM-1),y_dm(0:MNDM-1),z_dm(0:MNDM-1)
! * Correct Position *      
      double precision xc_dm(0:MNDM-1),yc_dm(0:MNDM-1),zc_dm(0:MNDM-1)
! * Velocity V(n+1/2) *
      double precision vx_dm(0:MNDM-1),vy_dm(0:MNDM-1),vz_dm(0:MNDM-1)
! * Virtual Velocity V(n+1/2) *      
      double precision vnx_dm(0:MNDM-1),vny_dm(0:MNDM-1),vnz_dm(0:MNDM-1)
! * acceleration *
      double precision dvx_dm(0:MNDM-1),dvy_dm(0:MNDM-1),dvz_dm(0:MNDM-1)
! * acceleration  at previous time step *
      double precision pdvx_dm(0:MNDM-1),pdvy_dm(0:MNDM-1) &
       ,pdvz_dm(0:MNDM-1)
! * for Individual time steep *
! * Individual time bin *      
      double precision dt_dm(0:MNDM-1)
! * time in system time step *      
      double precision lt_dm(0:MNDM-1)
! * virtual time bin *      
      double precision vdt_dm(0:MNDM-1)
! * minimum time step of neighbours *      
      double precision dtmnb_dm(0:MNDM-1),dtr_dm(0:MNDM-1)
! *** for adaptive softening ***
      double precision h_dm(0:MNDM-1),omgh_dm(0:MNDM-1),zetah_dm(0:MNDM-1) &
       ,rho_dm(0:MNDM-1),div_v_dm(0:MNDM-1)
! * list of active particles for DM *      
      integer list_adm(0:MNDM-1)
! * flag for time step changed 0: no change, otherwise changed *
      integer flagt_dm(0:MNDM-1)
end module gcdp_dm
#endif

#ifdef GAS
module gcdp_btree
      use gcdp_const
      implicit none
! ***** for Tree used also in ddecb.F ****
! * Number of nodes *
      integer num_tr
! * max node id sent to the other proc 
      integer nodese_tr
! * number of contained particle *      
       integer np_tr(0:MAXNODE-1)
! * name of Particle *      
      integer pn_tr(0:MAXNODE-1)
! * length of side *       
      double precision l_tr(0:MAXNODE-1)
! * Coordinate of center *      
      double precision cx_tr(0:MAXNODE-1),cy_tr(0:MAXNODE-1) &
       ,cz_tr(0:MAXNODE-1)
! * first child node *      
      integer daughter_tr(0:MAXNODE-1)
! * next node *      
      integer next_tr(0:MAXNODE-1)
#ifdef SGRAV
! * center of mass *      
      double precision cmx_tr(0:MAXNODE-1),cmy_tr(0:MAXNODE-1) &
       ,cmz_tr(0:MAXNODE-1)
! * total of mass *      
      double precision mass_tr(0:MAXNODE-1)
! * distance between cm and center, maximum softening *      
      double precision delta_tr(0:MAXNODE-1),hm_tr(0:MAXNODE-1)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision mx_tr(0:MAXNODE-1),my_tr(0:MAXNODE-1) &
       ,mz_tr(0:MAXNODE-1)
      double precision mxx_tr(0:MAXNODE-1),myy_tr(0:MAXNODE-1) &
       ,mzz_tr(0:MAXNODE-1)
      double precision mxy_tr(0:MAXNODE-1),myz_tr(0:MAXNODE-1) &
       ,mzx_tr(0:MAXNODE-1)
#endif
#endif
! *** proc id ***
      integer proc_tr(0:MAXNODE-1)
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
      integer np_dmtr(0:MAXNODEDM-1)
! * name of Particle *      
      integer pn_dmtr(0:MAXNODEDM-1)
! * length of side *       
      double precision l_dmtr(0:MAXNODEDM-1)
! * first child node *      
      integer daughter_dmtr(0:MAXNODEDM-1)
! * next node *      
      integer next_dmtr(0:MAXNODEDM-1)
! * Coordinate of center *      
      double precision cx_dmtr(0:MAXNODEDM-1),cy_dmtr(0:MAXNODEDM-1) &
       ,cz_dmtr(0:MAXNODEDM-1)
! * center of mass *      
      double precision cmx_dmtr(0:MAXNODEDM-1),cmy_dmtr(0:MAXNODEDM-1) &
       ,cmz_dmtr(0:MAXNODEDM-1)
! * total of mass *
      double precision mass_dmtr(0:MAXNODEDM-1) 
! * distance between cm and center *      
      double precision delta_dmtr(0:MAXNODEDM-1),hm_dmtr(0:MAXNODEDM-1)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision mx_dmtr(0:MAXNODEDM-1),my_dmtr(0:MAXNODEDM-1) &
       ,mz_dmtr(0:MAXNODEDM-1)
      double precision mxx_dmtr(0:MAXNODEDM-1),myy_dmtr(0:MAXNODEDM-1) &
       ,mzz_dmtr(0:MAXNODEDM-1)
      double precision mxy_dmtr(0:MAXNODEDM-1),myz_dmtr(0:MAXNODEDM-1) &
       ,mzx_dmtr(0:MAXNODEDM-1)
#endif
! *** proc id ***
      integer proc_dmtr(0:MAXNODEDM-1)
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
      integer np_ldmtr(0:MAXNODELDM-1)
! * name of Particle *      
      integer pn_ldmtr(0:MAXNODELDM-1)
! * length of side *       
      double precision l_ldmtr(0:MAXNODELDM-1),hm_ldmtr(0:MAXNODELDM-1)
! * first child node *      
      integer daughter_ldmtr(0:MAXNODELDM-1)
! * next node *      
      integer next_ldmtr(0:MAXNODELDM-1)
! * Coordinate of center *      
      double precision cx_ldmtr(0:MAXNODELDM-1) &
       ,cy_ldmtr(0:MAXNODELDM-1),cz_ldmtr(0:MAXNODELDM-1)
! * center of mass *      
      double precision cmx_ldmtr(0:MAXNODELDM-1) &
       ,cmy_ldmtr(0:MAXNODELDM-1),cmz_ldmtr(0:MAXNODELDM-1)
! * total of mass *
      double precision mass_ldmtr(0:MAXNODELDM-1)
! * distance between cm and center *      
      double precision delta_ldmtr(0:MAXNODELDM-1)
#ifdef QPTREEF
! * for Multipole Expansion *
      double precision mx_ldmtr(0:MAXNODELDM-1),my_ldmtr(0:MAXNODELDM-1) &
       ,mz_ldmtr(0:MAXNODELDM-1)
      double precision mxx_ldmtr(0:MAXNODELDM-1) &
       ,myy_ldmtr(0:MAXNODELDM-1),mzz_ldmtr(0:MAXNODELDM-1)
      double precision mxy_ldmtr(0:MAXNODELDM-1) &
       ,myz_ldmtr(0:MAXNODELDM-1),mzx_ldmtr(0:MAXNODELDM-1)
#endif
! *** proc id ***
      integer proc_ldmtr(0:MAXNODELDM-1)
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
      integer np_gtr(0:MAXNODE-1)
! * name of Particle *            
      integer pn_gtr(0:MAXNODE-1)
! * length of side *      
      double precision l_gtr(0:MAXNODE-1),hm_gtr(0:MAXNODE-1)
      double precision cx_gtr(0:MAXNODE-1),cy_gtr(0:MAXNODE-1) &
       ,cz_gtr(0:MAXNODE-1)
! * first child node *      
      integer daughter_gtr(0:MAXNODE-1)
! * next node *      
      integer next_gtr(0:MAXNODE-1)
! *** proc id ***
      integer proc_gtr(0:MAXNODE-1)
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

#ifdef COOL
module gcdp_cool
      use gcdp_const
      implicit none
! ***** for metal cooling and heating ***
      integer nz_crtb,nnh_crtb,nmet_crtb,nt_crtb
      double precision z_crtb(0:MNT_CRTB-1),nh_crtb(0:MNNH_CRTB) &
      ,met_crtb(0:MNMET_CRTB),t_crtb(0:MNT_CRTB)
      double precision cr_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,hr_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) &
       ,myu_crtb(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB)
end module gcdp_cool
#endif

