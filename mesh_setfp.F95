#include "gcdp.def"
! *********************************************
!  mesh_setrho.F95 for GCD+
!  13 Aug. 2013   written by D. Kawata
! *********************************************

! generating mesh points

#ifdef TREEPM
subroutine  mesh_setfp(ng,nag,nas,nadm)
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

      integer,intent(in) :: ng,nag,nas,nadm

      integer pn,i,j,k,ip0,jp0,kp0,ip1,jp1,kp1,w1,w2
      double precision df,vm
      double precision wx0,wx1,wy0,wy1,wz0,wz1
      integer ierr
      integer nmat,ic,pni
! *** for test file *** 
      character fileo*60
! particle data within proc
      integer npart
      integer,allocatable :: pnitp(:)
      double precision,allocatable :: xtp(:),ytp(:),ztp(:)
! for particle data for partial grid data
      integer np_m
      double precision,allocatable :: xp_m(:),yp_m(:),zp_m(:)
      double precision,allocatable :: fxp_m(:),fyp_m(:),fzp_m(:)
#ifdef FFTW3_MPI
      integer npjs,kp0proc,nval,ip,isend,irecv
      integer,allocatable :: npsproc(:)
      integer,allocatable :: npjr(:),idisp(:),jjlen(:)
      integer,allocatable :: pncomp(:),proccomp(:),lcomp(:)
      double precision,allocatable :: trbuf(:),tdvs(:)
#endif

! number of matter species
      nmat=0
#if defined(GAS) || defined(STAR)
      nmat=nmat+1
#endif
#ifdef DM  
      nmat=nmat+1
#endif

      do ic=0,nmat-1

! set particle data to xtp
        if(ic.eq.0) then
#if defined(GAS) || defined(STAR)
! allocate x,y,z,m data temporary
          npart=nag+nas

          allocate(xtp(0:npart))
          allocate(ytp(0:npart))
          allocate(ztp(0:npart))
          allocate(pnitp(0:npart))

          do i=0,nag-1
            pn=list_ap(i)
            pnitp(i)=pn
            xtp(i)=x_p(pn)
            ytp(i)=y_p(pn)
            ztp(i)=z_p(pn)
          enddo

          do i=0,nas-1
            pn=list_ap(ng+i)
            pnitp(i)=pn
            xtp(nag+i)=x_p(pn)
            ytp(nag+i)=y_p(pn)
            ztp(nag+i)=z_p(pn)
          enddo

#ifdef DM
        else if(ic.eq.1) then
#endif
#endif
#ifdef DM
          npart=nadm

          allocate(xtp(0:npart))
          allocate(ytp(0:npart))
          allocate(ztp(0:npart))
          allocate(pnitp(0:npart))

          do i=0,npart-1
            pn=list_adm(i)
            pnitp(i)=pn
            xtp(i)=x_dm(pn)
            ytp(i)=y_dm(pn)
            ztp(i)=z_dm(pn)
          enddo
#endif
        endif


#ifdef FFTW3_MPI
! set which particle send to which proc

        allocate(pncomp(0:npart))
        allocate(proccomp(0:npart))

!        do ip=0,nprocs-1
!          write(6,*) ' lnz,lzoff=',lnz_mp(ip),lzoff_mp(ip),myrank
!        enddo

        npjs=0
        do i=0,npart-1
! search which proc the force grid data are there, use kp0
          if(ztp(i).le.z0_m) then
            kp0=nz_m-1
          else if(ztp(i).ge.z_m1d(nz_m-1)) then
            kp0=nz_m-1
          else
            kp0=int((ztp(i)-z0_m)/dz_m)
          endif
! search which proc
          do ip=1,nprocs-1
            if(kp0.lt.lzoff_mp(ip)) then
              kp0proc=ip-1
              goto 70
            endif
          enddo
          kp0proc=nprocs-1
   70     continue

          pncomp(npjs)=i
          proccomp(npjs)=kp0proc
          npjs=npjs+1
        enddo
! set sending values, and idisp and jjlen
        nval=3

        allocate(idisp(0:nprocs))
        allocate(jjlen(0:nprocs))
        allocate(lcomp(0:npjs))

        allocate(npsproc(0:nprocs))

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
          npsproc(ip)=jjlen(ip)
        enddo

!      write(6,*) ' myrank,npjs=',myrank,npjs

        if(isend.ne.npjs) then
          write(6,*) ' Error in mesh_setfpo(): isend.ne.ncomp'
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
          tdvs(nval*i)=xtp(pni)
          tdvs(nval*i+1)=ytp(pni)
          tdvs(nval*i+2)=ztp(pni)
        enddo

        deallocate(xtp)
        deallocate(ytp)
        deallocate(ztp)

! reset sending parameters
        do i=0,nprocs-1
          idisp(i)=idisp(i)*nval
          jjlen(i)=jjlen(i)*nval
        enddo

        allocate(xp_m(0:np_m))
        allocate(yp_m(0:np_m))
        allocate(zp_m(0:np_m))

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
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo

        deallocate(tdvs)
        deallocate(jjlen)
        deallocate(idisp)

#else
! deallocated at mesh_fftf
        allocate(xp_m(0:npart))
        allocate(yp_m(0:npart))
        allocate(zp_m(0:npart))

! setting all the data to *_m
        np_m=npart
        do i=0,npart-1
          xp_m(i)=xtp(i)
          yp_m(i)=ytp(i)
          zp_m(i)=ztp(i)
        enddo

        deallocate(xtp)
        deallocate(ytp)
        deallocate(ztp)

#endif
        allocate(fxp_m(0:np_m))
        allocate(fyp_m(0:np_m))
        allocate(fzp_m(0:np_m))

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
#ifdef FFTW3_MPI
! in local grid   
          kp0=kp0-lzoff_m+1
          kp1=kp1-lzoff_m+1
          if(kp0.gt.lnz_m) then
            kp0=kp0-nz_m
          endif
          if(kp1.lt.0) then
            kp1=kp1+nz_m
          endif
          if((kp0.lt.0.or.kp0.gt.lnz_m+1).or.(kp1.lt.0.or.kp1.gt.lnz_m+1)) then
            write(6,*) ' Error in mesh_setfp(): kp0 or kp1 out of the range'
            write(6,*) ' myrank,kp0,kp1,z=',myrank,kp0,kp1,zp_m(i)
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
            stop
          endif
#endif

! assign force to the particles
! dvx_p kp0
          fxp_m(i)=wx0*wy0*wz0*fx_m(ip0,jp0,kp0) &
           +wx1*wy0*wz0*fx_m(ip1,jp0,kp0) &
           +wx1*wy1*wz0*fx_m(ip1,jp1,kp0) &
           +wx0*wy1*wz0*fx_m(ip0,jp1,kp0) &
! kp1
           +wx0*wy0*wz1*fx_m(ip0,jp0,kp1) &
           +wx1*wy0*wz1*fx_m(ip1,jp0,kp1) &
           +wx1*wy1*wz1*fx_m(ip1,jp1,kp1) &
           +wx0*wy1*wz1*fx_m(ip0,jp1,kp1)
! dvy_p kp0
          fyp_m(i)=wx0*wy0*wz0*fy_m(ip0,jp0,kp0) &
           +wx1*wy0*wz0*fy_m(ip1,jp0,kp0) &
           +wx1*wy1*wz0*fy_m(ip1,jp1,kp0) &
           +wx0*wy1*wz0*fy_m(ip0,jp1,kp0) &
! kp1
           +wx0*wy0*wz1*fy_m(ip0,jp0,kp1) &
           +wx1*wy0*wz1*fy_m(ip1,jp0,kp1) &
           +wx1*wy1*wz1*fy_m(ip1,jp1,kp1) &
           +wx0*wy1*wz1*fy_m(ip0,jp1,kp1)
! dvz_p kp0
          fzp_m(i)=wx0*wy0*wz0*fz_m(ip0,jp0,kp0) &
           +wx1*wy0*wz0*fz_m(ip1,jp0,kp0) &
           +wx1*wy1*wz0*fz_m(ip1,jp1,kp0) &
           +wx0*wy1*wz0*fz_m(ip0,jp1,kp0) &
! kp1
           +wx0*wy0*wz1*fz_m(ip0,jp0,kp1) &
           +wx1*wy0*wz1*fz_m(ip1,jp0,kp1) &
           +wx1*wy1*wz1*fz_m(ip1,jp1,kp1) &
           +wx0*wy1*wz1*fz_m(ip0,jp1,kp1)

        enddo

!      write(fileo,'(a3,i3.3)') 'p_m',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np_m-1
!        write(60,'(6(1pE13.5))') xp_m(i),yp_m(i),zp_m(i) &
!         ,fxp_m(i),fyp_m(i),fzp_m(i)
!      enddo
!      close(60)

#ifdef FFTW3_MPI
! send back force
        allocate(idisp(0:nprocs))
        allocate(jjlen(0:nprocs))

! set idisp and jjlen
        isend=0
        do ip=0,nprocs-1
          idisp(ip)=isend
          jjlen(ip)=npjr(ip)
          isend=isend+npjr(ip)
        enddo
        nval=3
        do ip=0,nprocs-1
          idisp(ip)=idisp(ip)*nval
          jjlen(ip)=jjlen(ip)*nval
        enddo

        allocate(tdvs(0:np_m*nval))

! store the data to be sent
        do i=0,np_m-1
          tdvs(nval*i)=fxp_m(i)
          tdvs(nval*i+1)=fyp_m(i)
          tdvs(nval*i+2)=fzp_m(i)
        enddo

! reallocate fxp,fyp,fzp
        deallocate(fxp_m)
        deallocate(fyp_m)
        deallocate(fzp_m)
        allocate(fxp_m(0:npart))
        allocate(fyp_m(0:npart))
        allocate(fzp_m(0:npart))

! sending and receiving the data
        irecv=0
        do ip=0,nprocs-1

          allocate(trbuf(0:npsproc(ip)*nval))

! use npsproc for number of receiving data
          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npsproc(ip)*nval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! set data to tdvr
          do i=0,npsproc(ip)-1
            fxp_m(irecv)=trbuf(nval*i)
            fyp_m(irecv)=trbuf(nval*i+1)
            fzp_m(irecv)=trbuf(nval*i+2)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo
        if(irecv.ne.npart) then
          write(6,*) 'Error in mesh_setfp(): irecv wrong when receiving fx'
          write(6,*) ' ic,myrank,irecv,npart=',ic,myrank,irecv,npart
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif

        deallocate(npjr)
        deallocate(npsproc)
        deallocate(idisp)
        deallocate(jjlen)
        deallocate(tdvs)

#endif

        if(ic.eq.0) then
#if defined(GAS) || defined(STAR)

          do i=0,npart-1
#ifdef FFTW3_MPI
            pn=lcomp(i)
            pni=pnitp(pn)
#else
            pni=pnitp(pn)
#endif
            dvx_p(pni)=dvx_p(pni)+fxp_m(i)
            dvy_p(pni)=dvy_p(pni)+fyp_m(i)
            dvz_p(pni)=dvz_p(pni)+fzp_m(i)
          enddo
#ifdef DM
        else if(ic.eq.1) then
#endif
#endif
#ifdef DM
          do i=0,nadm-1
#ifdef FFTW3_MPI
            pn=lcomp(i)
            pni=pnitp(pn)
#else
            pni=pnitp(pn)
#endif
            dvx_dm(pni)=dvx_dm(pni)+fxp_m(i)
            dvy_dm(pni)=dvy_dm(pni)+fyp_m(i)
            dvz_dm(pni)=dvz_dm(pni)+fzp_m(i)
          enddo
#endif
        endif

#ifdef FFTW3_MPI
        deallocate(lcomp)
#endif
        deallocate(pnitp)

        deallocate(xp_m)
        deallocate(yp_m)
        deallocate(zp_m)
        deallocate(fxp_m)
        deallocate(fyp_m)
        deallocate(fzp_m)
      enddo

end subroutine
#endif





