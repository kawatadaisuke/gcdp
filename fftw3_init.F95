#include "gcdp.def"
! *********************************************
!  fftw3_init.F95 for GCD+
!  15  Aug. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
#ifdef FFTW3
subroutine  fftw3_init()
      use gcdp_const
      use gcdp_system
      use gcdp_pm
      use fftw3

      implicit none
      include 'mpif.h'

      integer i,j,k
      double precision yp,zp
!      character filen*60
#ifdef FFTW3_MPI
      integer(C_INTPTR_T) :: ifftw,jfftw,kfftw,alloc_local

      nx_fftw=int(nx_m,C_INTPTR_T)
      ny_fftw=int(ny_m,C_INTPTR_T)
      nz_fftw=int(nz_m,C_INTPTR_T)
! get local data size
      alloc_local=fftw_mpi_local_size_3d_transposed( &
       nz_fftw,ny_fftw,nx_fftw,MPI_COMM_WORLD &
      ,lnz_fftw,lkoff_fftw,lny_fftw,ljoff_fftw)
! allocate data array pointer
      poid_fftw=fftw_alloc_complex(alloc_local)
! real and k-space data pointer
! note rsdata and ksdata starting from (1,1,1)
      call c_f_pointer(poid_fftw,rsd_fftw,[2*(nx_fftw/2+1),ny_fftw,lnz_fftw])
      call c_f_pointer(poid_fftw,ksd_fftw,[nx_fftw/2+1,nz_fftw,lny_fftw])
! note: nz,ny,nx order, because column ->. row measure 
      planf_fftw=fftw_mpi_plan_dft_r2c_3d(nz_fftw,ny_fftw,nx_fftw &
       ,rsd_fftw,ksd_fftw &
#ifdef FFTW3_ESTIMATE
       ,FFTW_ESTIMATE &
#else
       ,FFTW_MEASURE &
#endif
       ,FFTW_MPI_TRANSPOSED_OUT)
      planb_fftw=fftw_mpi_plan_dft_c2r_3d(nz_fftw,ny_fftw,nx_fftw &
        ,ksd_fftw,rsd_fftw &
#ifdef FFTW3_ESTIMATE
       ,FFTW_ESTIMATE &
#else
       ,FFTW_MEASURE &
#endif
       ,FFTW_MPI_TRANSPOSED_IN)

#else
! allocate data array pointer
      poid_fftw=fftw_alloc_complex(int((nx_m/2+1)*ny_m*nz_m,C_SIZE_T))
! real space data pointer
! note rsdata and ksdata starting from (1,1,1)
      call c_f_pointer(poid_fftw,rsd_fftw,[2*(nx_m/2+1),ny_m,nz_m])
      call c_f_pointer(poid_fftw,ksd_fftw,[nx_m/2+1,ny_m,nz_m])
! plan for real -> k space (forward FT)

! note: nz,ny,nx order, because column ->. row measure 
      planf_fftw=fftw_plan_dft_r2c_3d(nz_m,ny_m,nx_m,rsd_fftw,ksd_fftw &
#ifdef FFTW3_ESTIMATE
       ,FFTW_ESTIMATE)
#else
       ,FFTW_MEASURE)
#endif
      planb_fftw=fftw_plan_dft_c2r_3d(nz_m,ny_m,nx_m,ksd_fftw,rsd_fftw &
#ifdef FFTW3_ESTIMATE
       ,FFTW_ESTIMATE)
#else
       ,FFTW_MEASURE)
#endif
#endif

end subroutine
#endif FFTW3
#endif

