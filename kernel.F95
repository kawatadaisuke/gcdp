#include "gcdp.def"
! *****************************************************
!    kernel.F95 for GCD+ ver. f03.0
! 04  Jan. 2013   written by D. Kawata
! ***************************************************** 

subroutine setkernel()
      use gcdp_const
      use gcdp_kernel
      use gcdp_system

      implicit none

      integer i,j
      double precision f,cv,s

      ds_tb=1.0d0/dble(NKTAB)

! *** W: SPH kernel ***
! *** for constant ***
#ifdef SIM1D
      cv=4.0d0/3.0d0
#elif defined(SIM2D)
      cv=(40.0d0/(7.0d0*M_PI))
#else
      cv=8.0d0/M_PI
#endif
      do i=0,NKTAB
        s_tb(i)=ds_tb*dble(i)
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(1.0d0-6.0d0*(s**2)+6.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(2.0d0*((1.0d0-s)**3))
        else 
          f=0.0d0
        endif
        w_tb(i)=cv*f
      enddo
! *** dW(s)/ds/s: SPH kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=-12.0d0+18.0d0*s
        else if(s.le.1.0d0) then
          f=-6.0d0*((1.0d0-s)**2)/s
        else 
          f=0.0d0
        endif
        dwds_s_tb(i)=cv*f
! *** core dw/ds ***
#ifdef COREDW
        if(s.le.THIRD) then
          f=-2.0d0
        else if(s.le.0.5d0) then
#else
        if(s.le.0.5d0) then
#endif
          f=-12.0d0*s+18.0d0*s*s
        else if(s.le.1.0d0) then
          f=-6.0d0*((1.0d0-s)**2)
        else 
          f=0.0d0
        endif
        dwdsc_tb(i)=cv*f
      enddo
! *** h x dW/dh: SPH kernel ***
! #ifdef SIM1D
!      cv=1.0d0
! #elif defined(SIM2D)
!      cv=2.0d0
! #else
!      cv=3.0d0
! #endif
!      do i=0,NKTAB
!        s=s_tb(i)
!        hdwdh_tb(i)=-cv*w_tb(i)-dwds_s_tb(i)*(s**2)
!      enddo

#ifdef SGRAV
! *** dphi/dr/r: phi kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(32.0d0/15.0d0)*(5.0d0-18.0d0*(s**2)+15.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(4.0d0/15.0d0)*(80.0d0*s-180.0d0*(s**2)+144.0d0*(s**3) &
           -40.0d0*(s**4)-0.25d0/(s**2))/s
        endif
        dphidr_r_tb(i)=f
      enddo
! *** dphi/dh: phi kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=-40.0d0*(s**2)+120.0d0*(s**4)-96.0d0*(s**5)+7.0d0
        else if(s.le.1.0d0) then
          f=-80.0d0*(s**2)+160.0d0*(s**3)-120.0d0*(s**4)+32.0d0*(s**5)+8.0d0
        endif
        dphidh_tb(i)=0.4d0*f
      enddo
! *** for multipole ***
! *** d2phi/dr2 ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(32.0d0/15.0d0)*(5.0d0-54.0d0*(s**2)+60.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(4.0d0/15.0d0)*(80.0d0-360.0d0*s+432.0d0*(s**2) &
           -160.0d0*(s**3)+2.0d0/(4.0d0*(s**3)))
        endif
        d2phidr2_tb(i)=f
      enddo
! *** d3phi/dr3 ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(32.0d0/15.0d0)*(-108.0d0*s+180.0d0*(s**2))
        else if(s.le.1.0d0) then
          f=(4.0d0/15.0d0)*(-360.0d0+864.0d0*s-480.0d0*(s**2) &
           -6.0d0/(4.0d0*(s**4)))
        endif
        d3phidr3_tb(i)=f
      enddo
#endif


! *** for test
!      if(myrank.eq.0) then
!        open(60,file='kernel.dat',status='unknown')
!        do i=0,NKTAB 
!          s=s_tb(i)
!          write(60,160) s_tb(i),w_tb(i),dwds_s_tb(i),dphidr_r_tb(i) &
!           ,dphidh_tb(i),d2phidr2_tb(i),d3phidr3_tb(i)
!        enddo
!  160   format(7(1pE13.5))
!      endif

end subroutine
