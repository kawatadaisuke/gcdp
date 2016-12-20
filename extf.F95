#include "gcdp.def"
! **************************************************
!             extf.F 06/12/16
!  add force calculated from the external potential
! **************************************************

#ifdef EXTF
subroutine extf_paramset()
      use gcdp_const
      use gcdp_system
      implicit none
      include 'mpif.h'

      integer ierr

      if(myrank.eq.0) then
        write(6,*) ' call extf_paramset'
        open(60,file='./output/ana/inipara.dat',status='old',access='append')
      endif

#ifdef EXTF_Miyamoto_Nagai
      if(myrank.eq.0) then
        open(50,file='ini/input-mndisk.dat',status='old')
! disk mass in 1e12 Msun
        read(50,*) MN_mdisk
! kpc
        read(50,*) MN_a,MN_b
        close(50)
        write(6,*) ' Miyamoto-Nagai disk parameters'
        write(6,*) ' Mdisk (1e12 Msun)=',MN_mdisk
        write(6,*) ' a,b (kpc)=',MN_a,MN_b

        write(60,'(a31)') ' Miyamoto-Nagai disk parameters'
        write(60,'(a19,1pE13.5)') ' Mdisk (1e12 Msun)=',MN_mdisk
         write(60,'(a11,2(1pE13.5))') ' a,b (kpc)=',MN_a,MN_b
      endif
      call MPI_BCAST(MN_mdisk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(MN_a,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(MN_b,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! convert to GCD+ unit
      MN_a=MN_a/LUKPC
      MN_b=MN_b/LUKPC
#endif

#ifdef EXTF_NFW
      if(myrank.eq.0) then
        open(50,file='ini/input-nfwhalo.dat',status='old')
! Msun pc^-3
        read(50,*) NFW_rho0
! kpc
        read(50,*) NFW_rs
        close(50)

        write(6,*) ' NFW parameters'
        write(6,*) ' rho0 (Msun/pc^3)=',NFW_rho0
        write(6,*) ' rs (kpc)=',NFW_rs
        write(60,'(a15)') ' NFW parameters'
        write(60,'(a18,1pE13.5)') ' rho0 (Msun/pc^3)=',NFW_rho0
        write(60,'(a10,1pE13.5)') ' rs (kpc)=',NFW_rs
      endif
      call MPI_BCAST(NFW_rho0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NFW_rs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! convert to GCD+ unit
      NFW_rho0=NFW_rho0*(LUPC**3)/MUSM
      NFW_rs=NFW_rs/LUKPC
#endif

#ifdef EXTF_Hernquist
      if(myrank.eq.0) then
        open(50,file='ini/input-hqbulge.dat',status='old')
! 1e12 Msun
        read(50,*) HQ_mtot
! kpc
        read(50,*) HQ_ra
        close(50)

        write(6,*) ' Hernquist bulge parameters'
        write(6,*) ' Mtot (1e12 Msun) =',HQ_mtot
        write(6,*) ' a (kpc) =',HQ_ra
        write(60,'(a27)') ' Hernquist bulge parameters'
        write(60,'(a19,1pE13.5)') ' Mtot (1e12 Msun) =',HQ_mtot
!                                  1234567890123456789012345678901234
        write(60,'(a10,1pE13.5)') ' a (kpc) =',HQ_ra

      endif
      call MPI_BCAST(HQ_mtot,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(HQ_ra,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! convert to GCD+ unit
      HQ_ra=HQ_ra/LUKPC
#endif

      if(myrank.eq.0) then
        close(60)
      endif

end subroutine


#ifdef EXTF_Miyamoto_Nagai

subroutine Miyamoto_Nagai_extf(ndm,nadm)
      use gcdp_const
      use gcdp_system
      use gcdp_dm

      implicit none
      include 'mpif.h'

      integer,intent(in) :: ndm,nadm
      integer i,ir,pn
      double precision rxyp,frad,fx,fy,fz,abz,bz

      do i=0,nadm-1
         pn=list_adm(i)
         rxyp=dsqrt(x_dm(pn)**2+y_dm(pn)**2)
         bz=dsqrt(MN_b**2+z_dm(pn)**2)
         abz=MN_a+bz
         if(rxyp.gt.0.0d0) then
! radial force
           frad=G*MN_mdisk*rxyp*((rxyp**2+abz**2)**(-1.5d0))
           fx=frad*x_dm(pn)/rxyp
           fy=frad*y_dm(pn)/rxyp
         else
           fx=0.0d0
           fy=0.0d0
         endif
         if(z_dm(pn).gt.0.0d0) then
! z force
           fz=G*MN_mdisk*z_dm(pn)*abz &
             *(bz*((rxyp**2+abz**2)**(-1.5d0)))
         else
           fz=0.0d0
         endif
         dvx_dm(pn)=dvx_dm(pn)-fx
         dvy_dm(pn)=dvy_dm(pn)-fy
         dvz_dm(pn)=dvz_dm(pn)-fz
      enddo

end subroutine

#endif

#ifdef EXTF_NFW

! NFW external force
subroutine NFW_extf(ndm,nadm)
      use gcdp_const
      use gcdp_system
      use gcdp_dm

      implicit none
      include 'mpif.h'

      integer,intent(in) :: ndm,nadm
      integer i,ir,pn
      double precision rp,fr,fx,fy,fz,massr

      do i=0,nadm-1
         pn=list_adm(i)
         rp=dsqrt(x_dm(pn)**2+y_dm(pn)**2+z_dm(pn)**2)
         if(rp.gt.0.0d0) then
! mass(<r)
           massr=4.0d0*M_PI*NFW_rho0*(NFW_rs**3) &
             *(-rp/(NFW_rs+rp)-dlog(NFW_rs)+dlog(NFW_rs+rp))
! radial force / r
           fr=G*massr/rp**3
! adding x,y,z force
           dvx_dm(pn)=dvx_dm(pn)-x_dm(pn)*fr
           dvy_dm(pn)=dvy_dm(pn)-y_dm(pn)*fr
           dvz_dm(pn)=dvz_dm(pn)-z_dm(pn)*fr
        endif
      enddo

end subroutine

#endif

#ifdef EXTF_Hernquist

! Hernquist profile external force
subroutine Hernquist_extf(ndm,nadm)
      use gcdp_const
      use gcdp_system
      use gcdp_dm

      implicit none
      include 'mpif.h'

      integer,intent(in) :: ndm,nadm
      integer i,ir,pn
      double precision rp,fr,fx,fy,fz,massr

      do i=0,nadm-1
         pn=list_adm(i)
         rp=dsqrt(x_dm(pn)**2+y_dm(pn)**2+z_dm(pn)**2)
         if(rp.gt.0.0d0) then
! mass(<r)
           massr=HQ_mtot*(rp**2)/((rp+HQ_ra)**2)
! radial force / r
           fr=G*massr/rp**3
! adding x,y,z force
           dvx_dm(pn)=dvx_dm(pn)-x_dm(pn)*fr
           dvy_dm(pn)=dvy_dm(pn)-y_dm(pn)*fr
           dvz_dm(pn)=dvz_dm(pn)-z_dm(pn)*fr
        endif
      enddo

end subroutine

#endif
#endif
