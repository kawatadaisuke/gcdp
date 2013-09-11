#include "gcdp.def"
! ****************************************
!   define.F03
! 21 Jan. 2013    written by D.KAWATA
! ****************************************

module gcdp_const
      implicit none      

      double precision,parameter :: M_PI=3.141592653589793238d0
! delta timestep for domain decomposition ***
      integer,parameter ::  NDSTEPDDEC=10
! maximum number of cpus
      integer,parameter :: NCPU=256
! maximum number of node to be sent
      integer,parameter :: MAXNODESEND=70

! *****   for baryon particles  *****
! Max Number of Particle
#ifdef GAS
      integer,parameter :: MNB=1000000
#else
      integer,parameter :: MNB=10
#endif
! specific heat
      double precision,parameter :: GAM=5.0d0/3.0d0
!      double precision,parameter :: GAM=2.0d0

! *** for DM particles ***
! Max Number of Total DM Particle
#ifdef DM
      integer,parameter :: MNDM=1000000
#else
      integer,parameter :: MNDM=10
#endif

! *** maximum number of points to set up Mext(<r) ***
#ifdef MEXT 
      integer,parameter :: MNMEXT=200000
#endif

! *** for Tree ***
! Margin for Size of root tree
      double precision,parameter :: MGROOT=1.001d0
! Torelance Parameter, following Gadget-2 parameters
#ifdef QPTREEF
! *** Power et al. (2003), MNRAS, 338, 14
      double precision,parameter :: THETA=0.6d0
      double precision,parameter :: ALPHAT=0.001d0
#else
      double precision,parameter :: THETA=0.5d0
      double precision,parameter :: ALPHAT=0.0025d0
#endif

#ifdef TREEPM
! *** for TreePM ***
! factor for scale length w.r.t. grid size
      double precision,parameter :: FRSTPM=1.25d0
! factor for neighbour search for short range force
      double precision,parameter :: FRSHTPM=4.5d0
!      double precision,parameter :: FRSHTPM=6.0d0
#endif

! *** for at Table ***
      integer,parameter :: NATTABLE=50000
! *** for kernel tables ***
      integer,parameter :: NKTAB=10000

! *****   For system value   *****
! *** DT limit ***
      double precision,parameter :: DTLIMIT=1.0d-9
! *** for cooling 
! h > FH*EPS * ifdef GEPSLIM_G
      double precision,parameter :: FH=0.5d0
! temperature limit low and high (10^4K)    
      double precision,parameter :: TLLIMIT=0.003d0
      double precision,parameter :: THLIMIT=1.0e10
! cooling limit (10^4K)
      double precision,parameter :: CLIMIT=1.0d-3
! H2 cooling limit (K)
      double precision,parameter :: TLH2C=10.0d0
      double precision,parameter :: TUH2C=1.0d5
! cooling error
      double precision,parameter :: EPSC=1.0d-3
! dlog U in updateu
! since pv31.1 (17/02/09)
      double precision,parameter :: DU_UPU=0.025d0
! *** maxmum simulation radius ***
#ifdef MAXR
      double precision,parameter :: MAXRLIM=10.0d0
#endif
       
! *** Solar Abundances from WW95 ***
! *** H1 ***
      double precision,parameter :: XHSOL=0.706d0
      double precision,parameter :: XHPR=0.76d0
! *** Fe56 ***
      double precision,parameter :: XFESOL=0.00117d0
      
! *** For Artificial Viscosity and Correction for variable h ***
      double precision,parameter :: V_ALPHAINI=1.0d0
! *** V_ALPHAMAX should be smaller than 2 or
!     edit source term of alpha in update.F
#ifdef NOAVSWITCH
      double precision,parameter :: V_ALPHAMAX=V_ALPHAINI
      double precision,parameter :: V_ALPHAMIN=V_ALPHAINI
#else
      double precision,parameter :: V_ALPHAMAX=2.0d0
      double precision,parameter :: V_ALPHAMIN=0.5d0
#endif
      double precision,parameter :: V_BETA=3.0d0
! since pv33.6
      double precision,parameter :: SHAREV=1.0d0
      double precision,parameter :: U_ALPHAMAX=2.0d0
      double precision,parameter :: C_AV=0.2d0
! *** note:  h_gcdp = 2 h_monaghan ***
      double precision,parameter :: ETAH=2.4d0
      double precision,parameter :: ERRH=0.001d0
! *** ETAH for DM, which requires bigger than SPH particles ***
#ifdef GAS
! *** for the case of ETAH=2.4
      double precision,parameter :: ETAHDM=4.0d0
#else
      double precision,parameter :: ETAHDM=2.4d0
#endif
! *** limit for the maximum change of h during the iteration ***
      double precision,parameter :: DHFLIM=4.0d0

! *****   for Time Step   *****
! * Mergin Time Interval lower  *      
      double precision,parameter :: MGTL=0.8d0
! *                      higher *      
      double precision,parameter :: MGTU=1.2d0
! * Courant Number *      
      double precision,parameter :: CCFL=0.2d0
      double precision,parameter :: CDYN=0.2d0
      double precision,parameter :: CGRAV=0.15d0
! * dt factor from cooling rate in *time.F, esolve.F durin reionization
      double precision,parameter :: CCOOLDT=1.0d8
      double precision,parameter :: CCOOLEDT=1.0d8
      double precision,parameter :: CCOOLRIDT=0.1d0
! *** Factor of Individual Timestep limiter
! *** 4: Saitoh & Makino (2008)
      double precision,parameter :: FITLIM=4.0d0
! *** FAST limit GRAVITY vs. GAS 2^LIMIT***
#ifdef FASTG
      integer,parameter :: FASTGLIMIT=4
#else
      integer,parameter :: FASTGLIMIT=0
#endif
! *** number of small timesteps after EoR. ***
      integer,parameter :: NITEOR=10

! *** For Arithmetic Value ***
      double precision,parameter :: THIRD=1.0d0/3.0d0
      double precision,parameter :: V34=3.0d0/4.0d0
      double precision,parameter :: INF=1.0d6
      double precision,parameter :: MININF=1.0d-20
      double precision,parameter :: MININFESN=1.0d-20
      integer,parameter :: MAXFILENAME=25
      integer,parameter :: MAXSTEP=1000000

! *** Normalization Unit ***
! * mass unit ( solar mass ) *      
      double precision,parameter :: MUSM=1.0d12
! * (solarmass)/(gram) *      
!      parameter (MUG=1.99e45)
! * length unit / pc ( 100kpc/pc ) *      
      double precision,parameter :: LUPC=1.0e5
! * length unit(100kpc)/kpc *      
      double precision,parameter :: LUKPC=1.0e2
! * (100 kpc)/(cm) *      
      double precision,parameter :: LUCM=3.086e23
! * Density *      
      double precision,parameter :: DU=6.77e-26
! * Pressure *      
      double precision,parameter :: PU=2.91e-11
! * tempature ( K ) *      
      double precision,parameter :: TUK=1.0e4
! * time ( yr ) *      
      double precision,parameter :: TMUYR=4.71e8
! * time ( Gyr ) *      
      double precision,parameter :: TMUGYR=0.471d0
! * time (s) *      
      double precision,parameter :: TMUS=1.488d16
! * km/s *      
      double precision,parameter :: VUKMS=207.4d0
! * cooling (energy) rate (erg s^-1 cm^-3) *      
      double precision,parameter :: ERU=1.96d-27
! * k/m unit /( cm^2 s^-2 K^-1) *       
      double precision,parameter :: K_MU=4.3d10
! * erg *      
!      parameter (EU=8.56e59)

! *** For Physical Value ***
! * Gravitatrional Constant *      
      double precision,parameter :: G=1.0d0
! * Boltzmann constarnt *    
      double precision,parameter :: KCGS=1.381d-16
! * average molecular weight (default) *      
      double precision,parameter :: MYU=0.6d0
! * proton mass *      
      double precision,parameter :: MP=1.67265d-24
! * G(cgs) *      
      double precision,parameter :: GCGS=6.672d-8
      double precision,parameter :: TPRHO=((KCGS/(MYU*MP))/K_MU)
! * Hubble constant 9.78x10^9xh^-1 yr h = 1.0*      
      double precision,parameter :: H0_1=(1.0d0/((100.0d0*1.0d5)/(3.086d24))) &
       /((((3600.0d0*24.0d0)*365.0d0)*4.0d0+(3600.0d0*24.0d0))/4.0d0)
      double precision,parameter :: HUB0=100.0d0/(VUKMS*10.0d0)
! * Solar metallicity *
      double precision,parameter :: XZSOL=0.019d0
! * photon velocity *
      double precision,parameter :: CVEL=(2.997925e5/VUKMS)

! *** For star formation ***
! *** limited factor of increase of thermal energy by feedback ***
! 27/12/2010 1e4 is ok from sfdnograv test.
      double precision,parameter :: LDUFD=10000.0d0
! *** star formation threshold ***
! *** nh ***
! *** NSTH: only for softening length ***
      double precision,parameter :: NSTH=10.0d0
! *** NSTHLIM: density limit
      double precision,parameter :: NSTHLIM=200.0d0
! *** Hopkins13 self-gravity criterion for SF
!      double precision,parameter :: ALPHASF=1.0d0
!      double precision,parameter :: BETASF=0.5d0
! *** Temperature threshold (1.0e4 K) ***
!      double precision,parameter :: TSTH=10.0d0
! * star formation efficiency coefficient *
      double precision,parameter :: CSEFF=1.0d0
! * supernova energy fraction for kinematic energy *
      doubleprecision,parameter :: FVSN=0.0d0
#ifdef SF_EFD
! * supernova energy (unit 8.56e59 erg) *
      double precision,parameter :: SNEU=0.1d0*1.168224d-9
! * stellar wind energy (erg/s unit 5.752e43)
! *** 1.739e-8 = 1e36 erg/s ***
! 1e36-37 erg? Gibson (1994) MNRAS, 271, L35
      double precision,parameter :: SWEU=1.739d-8
#else
      double precision,parameter :: SNEU=0.0d0
      double precision,parameter :: SWEU=0.0d0
#endif
! *** star foramtion and feedback detal t ***
! *** used in ?time.F
      double precision,parameter :: MINTMSFD=0.001d0/TMUGYR
! *** used in cont.F
      double precision,parameter :: DTSF=MINTMSFD

! *** for uthlimit ***
! *** Jeans mass factor (N_Jeans in Hopkins et al. 2011)
      double precision,parameter :: ULJMFAC=10.0d0

! *** for lookup table of yields ***
! *** number of elements ***
      integer,parameter :: NYEL=9
! *** Number of Metallicities and Time ***
      integer,parameter :: NYTZ=60
      integer,parameter :: NYTT=200
! *** MNSPYTT: max number of mass group for star particle 
! *** should be more than 1/3 of number of group 
      integer,parameter :: MNSPYTT=62
#ifdef COOL
      integer,parameter :: MNT_CRTB=100
      integer,parameter :: MNMET_CRTB=10
      integer,parameter :: MNNH_CRTB=20
      integer,parameter :: MNZ_CRTB=50
#endif
!*** For Parallel ***
! *** Peano-Hilbert curve ***
      integer*8,parameter :: NDPH=8
! *** for parallel metal diffusion 
#ifdef METALDIFF
! *** factor for metal diffusion coefficient Grief et al. = 2, but seems too big
      double precision,parameter :: ZDFCFAC=1.0d0
      integer,parameter :: MNVALZD=MNB*11
#endif

#ifdef GLASS
!cc random acceleration of a fraction of dvx and etc. in calc_dv_du
      double precision,parameter :: FGLASS=0.01d0
#endif

#ifdef ARTFDM
!cc parameters from Widrow (2008)
      double precision,parameter :: B1MAX_AF=0.24d0*TMUGYR
      double precision,parameter :: B2MAX_AF=0.08d0*TMUGYR
      double precision,parameter :: TG_AF=1.0d0/TMUGYR
      double precision,parameter :: TD_AF=1.0d0/TMUGYR
! *** time at the end of maximum force and adding force
      double precision,parameter :: TMC_AF=3.0d0/TMUGYR+TG_AF
      double precision,parameter :: TME_AF=TD_AF+TMC_AF
#endif

#ifdef MHD 
! *** constsnts for MHD ***
      double precision,parameter :: ALPB_MIN=0.0d0
      double precision,parameter :: ALPB_MAX=1.0d0
      double precision,parameter :: ALPHA_CONST=1.0e-3
#ifdef BBAT
! * Normalize Battery term (e/mp)*
      double precision,parameter :: MCH=4.1463e21
#endif
#endif

end module gcdp_const
