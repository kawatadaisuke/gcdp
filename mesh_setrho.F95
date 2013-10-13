#include "gcdp.def"
! *********************************************
!  mesh_setrho.F95 for GCD+
!  07 Oct. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
subroutine  mesh_setrho(np,ndm)
      use gcdp_const
      use gcdp_system
      use gcdp_pm
#if defined(GAS) || defined(STAR)
      use gcdp_baryon
#endif
#ifdef DM
      use gcdp_dm
#endif

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np,ndm

      integer i,j,k,ip0,jp0,kp0,ip1,jp1,kp1,w1,w2,npmp
      integer npj
      double precision df,vm
      double precision wx0,wx1,wy0,wy1,wz0,wz1
      integer ierr
      double precision,allocatable :: tdvs(:),tdvr(:)
#ifdef FFTW3_MPI
      integer npjs,kp0proc,kp1proc,nval,ip,isend,irecv,pni
      integer,allocatable :: npjr(:),idisp(:),jjlen(:)
      integer,allocatable :: pncomp(:),proccomp(:),lcomp(:)
      double precision,allocatable :: trbuf(:)
#endif
! for DM particle data for partial grid data
      integer np_m
      double precision,allocatable :: xp_m(:),yp_m(:),zp_m(:),mp_m(:)
!      integer,allocatable :: procp_m(:),ndmsproc_m(:)
! *** for test file *** 
      character fileo*60

! initialisation
#ifdef FFTW3_MPI
      do k=0,lnz_m-1
#else
      do k=0,nz_m-1
#endif
        do j=0,ny_m-1
          do i=0,nx_m-1
            rho_m(i,j,k)=0.0d0
          enddo
        enddo
      enddo

! assign baryon mass
#if defined(GAS) || defined(STAR)
      do i=0,np-1
! x
        if(x_p(i).le.x0_m) then
          ip0=nx_m-1
          ip1=0
          df=(x_m1d(0)-x_p(i))/dx_m
          wx0=df
          wx1=1.0d0-df
        else if(x_p(i).ge.x_m1d(nx_m-1)) then
          ip0=nx_m-1
          ip1=0
          df=(x_p(i)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        else
          ip0=int((x_p(i)-x0_m)/dx_m)
          ip1=ip0+1
          df=(x_p(i)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        endif
! for check
!        if((wx0.lt.0.0d0.or.wx0.gt.1.0d0) &
!          .or.(wx1.lt.0.0d0.or.wx1.gt.1.0d0)) then
!          write(6,*) ' Error in mesh_setrho() for baryon: wx0,wx1=',wx0,wx1
!          write(6,*) ' x,ip0,ip1=',x_p(i),ip0,ip1
!          stop
!        endif
! y
        if(y_p(i).le.y0_m) then
          jp0=ny_m-1
          jp1=0
          df=(y_m1d(0)-y_p(i))/dy_m
          wy0=df
          wy1=1.0d0-df
        else if(y_p(i).ge.y_m1d(ny_m-1)) then
          jp0=ny_m-1
          jp1=0
          df=(y_p(i)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        else
          jp0=int((y_p(i)-y0_m)/dy_m)
          jp1=jp0+1
          df=(y_p(i)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        endif
! for check
!        if((wy0.lt.0.0d0.or.wy0.gt.1.0d0) &
!          .or.(wy1.lt.0.0d0.or.wy1.gt.1.0d0)) then
!          write(6,*) ' Error in mesh_setrho() for baryon: wy0,wy1=',wy0,wy1
!          write(6,*) ' y,jp0,jp1=',y_p(i),jp0,jp1
!          stop
!        endif
! z
        if(z_p(i).le.z0_m) then
          kp0=nz_m-1
          kp1=0
          df=(z_m1d(0)-z_p(i))/dz_m
          wz0=df
          wz1=1.0d0-df
        else if(z_p(i).ge.z_m1d(nz_m-1)) then
          kp0=nz_m-1
          kp1=0
          df=(z_p(i)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        else
          kp0=int((z_p(i)-z0_m)/dz_m)
          kp1=kp0+1
          df=(z_p(i)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        endif
! for check
!        if((wz0.lt.0.0d0.or.wz0.gt.1.0d0) &
!          .or.(wz1.lt.0.0d0.or.wz1.gt.1.0d0)) then
!          write(6,*) ' Error in mesh_setrho() for baryon: wz0,wz1=',wz0,wz1
!          write(6,*) ' z,kp0,kp1=',z_p(i),kp0,kp1
!          stop
!        endif

! assign the mass
! kp0
        rho_m(ip0,jp0,kp0)=rho_m(ip0,jp0,kp0) &
         +wx0*wy0*wz0*m_p(i)
        rho_m(ip1,jp0,kp0)=rho_m(ip1,jp0,kp0) &
         +wx1*wy0*wz0*m_p(i)
        rho_m(ip1,jp1,kp0)=rho_m(ip1,jp1,kp0) &
         +wx1*wy1*wz0*m_p(i)
        rho_m(ip0,jp1,kp0)=rho_m(ip0,jp1,kp0) &
         +wx0*wy1*wz0*m_p(i)
! kp1
        rho_m(ip0,jp0,kp1)=rho_m(ip0,jp0,kp1) &
         +wx0*wy0*wz1*m_p(i)
        rho_m(ip1,jp0,kp1)=rho_m(ip1,jp0,kp1) &
         +wx1*wy0*wz1*m_p(i)
        rho_m(ip1,jp1,kp1)=rho_m(ip1,jp1,kp1) &
         +wx1*wy1*wz1*m_p(i)
        rho_m(ip0,jp1,kp1)=rho_m(ip0,jp1,kp1) &
         +wx0*wy1*wz1*m_p(i)

      enddo
#endif

! assign DM mass

#ifdef FFTW3_MPI
! set which particle send to which proc
      allocate(pncomp(0:2*ndm))
      allocate(proccomp(0:2*ndm))

!      do ip=0,nprocs-1
!        write(6,*) ' lnz,lzoff=',lnz_mp(ip),lzoff_mp(ip),myrank
!      enddo

      npjs=0
      do i=0,ndm-1
        if(z_dm(i).le.z0_m) then
          kp0=nz_m-1
          kp1=0
        else if(z_dm(i).ge.z_m1d(nz_m-1)) then
          kp0=nz_m-1
          kp1=0
        else
          kp0=int((z_dm(i)-z0_m)/dz_m)
          kp1=kp0+1
        endif
! search which proc
        do ip=1,nprocs-1
          if(kp0.lt.lzoff_mp(ip)) then
            kp0proc=ip-1
            goto 70
          endif
        enddo
        kp0proc=nprocs-1
   70   do ip=1,nprocs-1
          if(kp1.lt.lzoff_mp(ip)) then
            kp1proc=ip-1
            goto 71
          endif
        enddo
        kp1proc=nprocs-1
   71   continue

        if(kp0proc.eq.kp1proc) then
          pncomp(npjs)=i
          proccomp(npjs)=kp0proc
          npjs=npjs+1
        else
          pncomp(npjs)=i
          proccomp(npjs)=kp0proc
          npjs=npjs+1
          pncomp(npjs)=i
          proccomp(npjs)=kp1proc
          npjs=npjs+1
        endif
      enddo
! set sending values, and idisp and jjlen
      nval=4

      allocate(idisp(0:nprocs))
      allocate(jjlen(0:nprocs))
      allocate(lcomp(0:npjs))

! this should be deallocated in mesh_fftf
!      allocate(ndmsproc_m(0:nprocs))

! store particle list in the order of sending proc
      isend=0
      do ip=0,nprocs-1
        idisp(ip)=isend
        jjlen(ip)=0
        do i=0,npjs-1
          if(proccomp(i).eq.ip) then
            lcomp(isend)=pncomp(i)
            jjlen(ip)=jjlen(ip)+1
            isend=isend+1            
          endif
        enddo
! store the number of particles sending to each proc
!        ndmsproc_m(ip)=jjlen(ip)
      enddo
      if(isend.ne.npjs) then
        write(6,*) ' Error in mesh_setrho(): isend.ne.ncomp'
        write(6,*) ' when counting N particles need communication'
        write(6,*) ' myrank,isend,npjs=',myrank,isend,npjs
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
        stop
      endif

      deallocate(pncomp)
      deallocate(proccomp)

! getting the total number of particles to recieved at each proc

      allocate(npjr(0:nprocs))

      np_m=0
      do ip=0,nprocs-1
        irecv=0
        call MPI_SCATTER(jjlen,1,MPI_INTEGER &
         ,irecv,1,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
        npjr(ip)=irecv
        np_m=np_m+irecv
      enddo

      allocate(tdvs(0:npjs*nval))

! store the data to be sent
      do i=0,npjs-1
        pni=lcomp(i)
        tdvs(nval*i)=x_dm(pni)
        tdvs(nval*i+1)=y_dm(pni)
        tdvs(nval*i+2)=z_dm(pni)
        tdvs(nval*i+3)=m_dm(pni)
      enddo
! reset sending parameters
      do i=0,nprocs-1
        idisp(i)=idisp(i)*nval
        jjlen(i)=jjlen(i)*nval
      enddo

      allocate(xp_m(0:np_m))
      allocate(yp_m(0:np_m))
      allocate(zp_m(0:np_m))
      allocate(mp_m(0:np_m))
!      allocate(procp_m(0:np_m))

! sending and receiving the data
      irecv=0
      do ip=0,nprocs-1

        allocate(trbuf(0:npjr(ip)*nval))

        call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
         ,trbuf,npjr(ip)*nval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! set data to tdvr
        do i=0,npjr(ip)-1
          xp_m(irecv)=trbuf(nval*i)
          yp_m(irecv)=trbuf(nval*i+1)
          zp_m(irecv)=trbuf(nval*i+2)
          mp_m(irecv)=trbuf(nval*i+3)
          irecv=irecv+1
        enddo

        deallocate(trbuf)

      enddo

      deallocate(lcomp)
      deallocate(npjr)
      deallocate(tdvs)
      deallocate(jjlen)
      deallocate(idisp)

!      write(fileo,'(a4,i3.3)') 'p_m',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np_m-1
!        write(60,'(4(1pE13.5))') xp_m(i),yp_m(i),zp_m(i),mp_m(i)
!      enddo
!      close(60)

#else
! deallocated at mesh_fftf
      allocate(xp_m(0:ndm))
      allocate(yp_m(0:ndm))
      allocate(zp_m(0:ndm))
      allocate(mp_m(0:ndm))

! setting all the data to *_m
      np_m=ndm
      do i=0,np_m-1
        xp_m(i)=x_dm(i)
        yp_m(i)=y_dm(i)
        zp_m(i)=z_dm(i)
        mp_m(i)=m_dm(i)
      enddo
#endif


#ifdef DM
      do i=0,np_m-1
! x
        if(xp_m(i).le.x0_m) then
          ip0=nx_m-1
          ip1=0
          df=(x_m1d(0)-xp_m(i))/dx_m
          wx0=df
          wx1=1.0d0-df
        else if(xp_m(i).ge.x_m1d(nx_m-1)) then
          ip0=nx_m-1
          ip1=0
          df=(xp_m(i)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        else
          ip0=int((xp_m(i)-x0_m)/dx_m)
          ip1=ip0+1
          df=(xp_m(i)-x_m1d(ip0))/dx_m
          wx0=1.0d0-df
          wx1=df
        endif
! for check
!        if((wx0.lt.0.0d0.or.wx0.gt.1.0d0) &
!          .or.(wx1.lt.0.0d0.or.wx1.gt.1.0d0)) then
!          write(6,*) ' Error in mesh_setrho() for DM: wx0,wx1=',wx0,wx1
!          write(6,*) ' x,ip0,ip1=',xp_m(i),ip0,ip1
!          stop
!        endif
! y
        if(yp_m(i).le.y0_m) then
          jp0=ny_m-1
          jp1=0
          df=(y_m1d(0)-yp_m(i))/dy_m
          wy0=df
          wy1=1.0d0-df
        else if(yp_m(i).ge.y_m1d(ny_m-1)) then
          jp0=ny_m-1
          jp1=0
          df=(yp_m(i)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        else
          jp0=int((yp_m(i)-y0_m)/dy_m)
          jp1=jp0+1
          df=(yp_m(i)-y_m1d(jp0))/dy_m
          wy0=1.0d0-df
          wy1=df
        endif
! for check
!        if((wy0.lt.0.0d0.or.wy0.gt.1.0d0) &
!          .or.(wy1.lt.0.0d0.or.wy1.gt.1.0d0)) then
!          write(6,*) ' Error in mesh_setrho() for DM: wy0,wy1=',wy0,wy1
!          write(6,*) ' y,jp0,jp1=',yp_m(i),jp0,jp1
!          stop
!        endif
! z
        if(zp_m(i).le.z0_m) then
          kp0=nz_m-1
          kp1=0
          df=(z_m1d(0)-zp_m(i))/dz_m
          wz0=df
          wz1=1.0d0-df
        else if(zp_m(i).ge.z_m1d(nz_m-1)) then
          kp0=nz_m-1
          kp1=0
          df=(zp_m(i)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        else
          kp0=int((zp_m(i)-z0_m)/dz_m)
          kp1=kp0+1
          df=(zp_m(i)-z_m1d(kp0))/dz_m
          wz0=1.0d0-df
          wz1=df
        endif
! for check
!        if((wz0.lt.0.0d0.or.wz0.gt.1.0d0) &
!          .or.(wz1.lt.0.0d0.or.wz1.gt.1.0d0)) then
!          write(6,*) ' Error in mesh_setrho() for DM: wz0,wz1=',wz0,wz1
!          write(6,*) ' z,kp0,kp1=',zp_m(i),kp0,kp1
!          stop
!        endif
#ifdef FFTW3_MPI
        kp0=kp0-lzoff_m
        kp1=kp1-lzoff_m

        if(kp0.lt.0) then
          kp0=0
          wz0=0.0d0
        else if(kp0.ge.lnz_m) then
          kp0=lnz_m-1
          wz0=0.0d0
        endif          
        if(kp1.lt.0) then
          kp1=0
          wz1=0.0d0
        else if(kp1.ge.lnz_m) then
          kp1=lnz_m-1
          wz1=0.0d0
        endif          

        if(kp0.ge.lnz_m.or.kp1.lt.0) then
          write(6,*) ' Error in mesh_setrho(): out of range kp0 or kp1'
          write(6,*) ' myrank,kp0,kp1,lnz,lzoff=',myrank,kp0,kp1,lnz_m,lzoff_m
          write(6,*) ' z=',zp_m(i)
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
       endif

#endif

! assign the density
! kp0
        rho_m(ip0,jp0,kp0)=rho_m(ip0,jp0,kp0) &
         +wx0*wy0*wz0*mp_m(i)
        rho_m(ip1,jp0,kp0)=rho_m(ip1,jp0,kp0) &
         +wx1*wy0*wz0*mp_m(i)
        rho_m(ip1,jp1,kp0)=rho_m(ip1,jp1,kp0) &
         +wx1*wy1*wz0*mp_m(i)
        rho_m(ip0,jp1,kp0)=rho_m(ip0,jp1,kp0) &
         +wx0*wy1*wz0*mp_m(i)
! kp1
        rho_m(ip0,jp0,kp1)=rho_m(ip0,jp0,kp1) &
         +wx0*wy0*wz1*mp_m(i)
        rho_m(ip1,jp0,kp1)=rho_m(ip1,jp0,kp1) &
         +wx1*wy0*wz1*mp_m(i)
        rho_m(ip1,jp1,kp1)=rho_m(ip1,jp1,kp1) &
         +wx1*wy1*wz1*mp_m(i)
        rho_m(ip0,jp1,kp1)=rho_m(ip0,jp1,kp1) &
         +wx0*wy1*wz1*mp_m(i)

      enddo
#endif

      deallocate(xp_m)
      deallocate(yp_m)
      deallocate(zp_m)
      deallocate(mp_m)

#ifndef FFTW3_MPI
      if(nprocs.gt.0) then

        allocate(tdvr(0:nt_m))
        allocate(tdvs(0:nt_m))
     
! sum up the density from all the cores
        npmp=0
        do k=0,nz_m-1
          do j=0,ny_m-1
            do i=0,nx_m-1
              tdvs(npmp)=rho_m(i,j,k)
              tdvr(npmp)=0.0d0
              npmp=npmp+1
            enddo
          enddo
        enddo
        call MPI_ALLREDUCE(tdvs,tdvr,npmp,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
        npmp=0
        do k=0,nz_m-1
          do j=0,ny_m-1
            do i=0,nx_m-1
              rho_m(i,j,k)=tdvr(npmp)
              npmp=npmp+1
            enddo
          enddo
        enddo

        deallocate(tdvr)
        deallocate(tdvs)

      endif

#endif


! volume of each mesh
      vm=dx_m*dy_m*dz_m
! density
#ifdef FFTW3_MPI
      do k=0,lnz_m-1
#else
      do k=0,nz_m-1
#endif
        do j=0,ny_m-1
          do i=0,nx_m-1
            rho_m(i,j,k)=rho_m(i,j,k)/vm
          enddo
        enddo
      enddo

! test output
!      write(fileo,'(a7,i3.3)') 'meshrho',myrank
!      open(60,file=fileo,status='unknown')

!#ifdef FFTW3_MPI
!      do k=0,lnz_m-1
!#else
!      do k=0,nz_m-1
!#endif
!        do j=0,ny_m-1
!          do i=0,nx_m-1
!#ifdef FFTW3_MPI
!            write(60,'(3I10,(1pE13.5))') i,j,k+lzoff_m,rho_m(i,j,k)
!#else
!            write(60,'(3I10,(1pE13.5))') i,j,k,rho_m(i,j,k)
!#endif
!          enddo
!        enddo
!      enddo
!      close(60)


end subroutine
#endif





