#include "gcdp.def"
! ******************************************
!   dtyields.F95 for GCD+ ver. f03.0
!  3  Feb. 2013  written by D. Kawata
! ******************************************
! ****************************************************
!  return total ejected mass until TM_tot+TM_dt
! ****************************************************

#if defined(SF_EFD) || defined(SF_ZFD)
subroutine dtyields(np,ng,ns,tmej_sn,tnsn_sn,tmzHe_sn &
       ,tmzC_sn,tmzN_sn,tmzO_sn,tmzNe_sn,tmzMg_sn &
       ,tmzSi_sn,tmzFe_sn,tmzZ_sn,nsw_sn)
      use gcdp_const
      use gcdp_baryon
      use gcdp_system
      use gcdp_metal
      use gcdp_yields

      implicit none       

      integer,intent(in) :: np,ng,ns
      double precision,intent(out) :: tmej_sn(0:np),tnsn_sn(0:np) &
       ,tmzHe_sn(0:np),tmzC_sn(0:np),tmzN_sn(0:np),tmzO_sn(0:np) &
       ,tmzNe_sn(0:np),tmzMg_sn(0:np),tmzSi_sn(0:np),tmzFe_sn(0:np) &
       ,tmzZ_sn(0:np),nsw_sn(0:np)
      integer i,j,pn,isp
      integer iz1,it1
      double precision ms0p
! *** for interpolation ***
      double precision wz1,wz2,wt1,wt2,dz,dt,dtsw,agesw
! *** for initial value for mass group ***
      double precision tsp1,tsp2,tspiz,tspez,dtsp
! work allocatable
      double precision,allocatable :: zg(:),age(:)

! *** Metallicity of Age for each star and gas ***

      allocate(zg(0:np))
      allocate(age(0:np))

      do i = 0,np-1
        pn = list_ap(i)
#ifdef METAL
! *** assuming m_p never change, use original Z ***
        zg(pn) = mzZ0_p(pn)/(m_p(pn)*MUSM)
#else
        zg(pn) = 0.0d0
#endif
! since pv32.0
        age(pn) = (TM_tot+TM_dt-ts_p(pn))*TMUGYR
      enddo

      dz=(dlog10(z_ytb(NYTZ+1))-dlog10(z_ytb(1)))/dble(NYTZ)
      dt=(dlog10(t_ytb(NYTT+1))-dlog10(t_ytb(1)))/dble(NYTT)
! *** Get the yields ***
      do i = 0,np-1
        pn = list_ap(i)
        tmej_sn(pn)=0.0d0
        tnsn_sn(pn)=0.0d0
! *** Unit Solar Mass ***
        tmzHe_sn(pn)=0.0d0
        tmzC_sn(pn)=0.0d0
        tmzN_sn(pn)=0.0d0
        tmzO_sn(pn)=0.0d0
        tmzNe_sn(pn)=0.0d0
        tmzMg_sn(pn)=0.0d0
        tmzSi_sn(pn)=0.0d0
        tmzFe_sn(pn)=0.0d0
        tmzZ_sn(pn)=0.0d0
! *** number of stars with stellar wind ***
        nsw_sn(pn)=0.0d0
        if(flagfd_p(pn).ne.0) then
          if(flagfd_p(pn).gt.0) then
            isp=flagfd_p(pn)-1
          else
            isp=-flagfd_p(pn)-1
          endif
! *** for metallicity weight ***
          if(zg(pn).ge.z_ytb(1)) then
            iz1=int((dlog10(zg(pn))-dlog10(z_ytb(1)))/dz)+1
          else 
            iz1=0
          endif
! *** only if the particles are qualified as feedback gas 
! *** in both metallicity tables of iz1 and iz1
          if(isp.lt.nsp_ytb(iz1)-1.and.isp.lt.nsp_ytb(iz1+1)-1) then
! *** weight for metallicity table ***
            wz2=(zg(pn)-z_ytb(iz1))/(z_ytb(iz1+1)-z_ytb(iz1))
            wz1=1.0d0-wz2
! *** calculate starting time of the group ***
            if(isp.lt.SI_snii) then
! *** if mass group includes SNe II, all treat as a SNe II group
              tspiz=0.0d0
! *** calculate starting time of the group ***
              tsp1=tspi_ytb(SI_snii,iz1)
              tsp2=tspi_ytb(SI_snii,iz1+1)
              tspez=tsp1*wz1+tsp2*wz2
            else
              tsp1=tspi_ytb(isp,iz1)
              tsp2=tspi_ytb(isp,iz1+1)
! *** starting time of the mass group ***
              tspiz=tsp1*wz1+tsp2*wz2
! *** end time of group mass ***
! *** calculate starting time of the group ***
              tsp1=tspi_ytb(isp+1,iz1)
              tsp2=tspi_ytb(isp+1,iz1+1)
              tspez=tsp1*wz1+tsp2*wz2
            endif
! *** feedback duration of the mass group ***
            dtsp=tspez-tspiz
            if(age(pn).gt.tspiz) then
              if(flagfd_p(pn).gt.0) then
                flagrfd_p(pn)=1    
              endif
              if(flagfd_p(pn).gt.0.and.tspiz.gt.age(pn)-TM_dt*TMUGYR) then
! *** set initial values, including the case of isp<SI_snii ***
                if(tspiz.ge.t_ytb(1)) then 
                  it1=int((dlog10(tspiz)-dlog10(t_ytb(1)))/dt)+1
                else
                  it1=0
                endif
                wt2=(tspiz-t_ytb(it1))/(t_ytb(it1+1)-t_ytb(it1))
                wt1=1.0d0-wt2
              else
                if(age(pn).gt.tspez) then
                  age(pn)=tspez
! *** flagrfd_p will be used in starfd ***
                  flagrfd_p(pn)=-1
                endif
! *** feedback gas particles ***
! *** search time ***
                if(age(pn).ge.t_ytb(1)) then 
                  it1=int((dlog10(age(pn))-dlog10(t_ytb(1)))/dt)+1
                else
                  it1=0
                endif
                wt2=(age(pn)-t_ytb(it1))/(t_ytb(it1+1)-t_ytb(it1))
                wt1=1.0d0-wt2
              endif
! *** each particles are 1/SI_nsp part ***
              if(isp.lt.SI_snii) then
                ms0p=m_p(pn)*dble(SI_nsp)/dble(SI_snii)
              else
                ms0p=m_p(pn)*dble(SI_nsp)
              endif  
! *** calculate total yield up the age ***
              tmej_sn(pn)=(wt1*wz1*mej_ytb(it1,iz1) &
               +wt1*wz2*mej_ytb(it1,iz1+1)+wt2*wz1*mej_ytb(it1+1,iz1) &
               +wt2*wz2*mej_ytb(it1+1,iz1+1))*ms0p 
              tnsn_sn(pn)=(wt1*wz1*nsn_ytb(it1,iz1) &
               +wt1*wz2*nsn_ytb(it1,iz1+1)+wt2*wz1*nsn_ytb(it1+1,iz1) &
               +wt2*wz2*nsn_ytb(it1+1,iz1+1))*(ms0p*MUSM)
! *** Unit Solar Mass ***
              tmzHe_sn(pn)=(wt1*wz1*mzHe_ytb(it1,iz1) &
               +wt1*wz2*mzHe_ytb(it1,iz1+1)+wt2*wz1*mzHe_ytb(it1+1,iz1) &
               +wt2*wz2*mzHe_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzC_sn(pn)=(wt1*wz1*mzC_ytb(it1,iz1) &
               +wt1*wz2*mzC_ytb(it1,iz1+1)+wt2*wz1*mzC_ytb(it1+1,iz1) &
               +wt2*wz2*mzC_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzN_sn(pn)=(wt1*wz1*mzN_ytb(it1,iz1) &
               +wt1*wz2*mzN_ytb(it1,iz1+1)+wt2*wz1*mzN_ytb(it1+1,iz1) &
               +wt2*wz2*mzN_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzO_sn(pn)=(wt1*wz1*mzO_ytb(it1,iz1) &
               +wt1*wz2*mzO_ytb(it1,iz1+1)+wt2*wz1*mzO_ytb(it1+1,iz1) &
               +wt2*wz2*mzO_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzNe_sn(pn)=(wt1*wz1*mzNe_ytb(it1,iz1) &
               +wt1*wz2*mzNe_ytb(it1,iz1+1)+wt2*wz1*mzNe_ytb(it1+1,iz1) &
               +wt2*wz2*mzNe_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzMg_sn(pn)=(wt1*wz1*mzMg_ytb(it1,iz1) &
               +wt1*wz2*mzMg_ytb(it1,iz1+1)+wt2*wz1*mzMg_ytb(it1+1,iz1) &
               +wt2*wz2*mzMg_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzSi_sn(pn)=(wt1*wz1*mzSi_ytb(it1,iz1) &
               +wt1*wz2*mzSi_ytb(it1,iz1+1)+wt2*wz1*mzSi_ytb(it1+1,iz1) &
               +wt2*wz2*mzSi_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzFe_sn(pn)=(wt1*wz1*mzFe_ytb(it1,iz1) &
               +wt1*wz2*mzFe_ytb(it1,iz1+1)+wt2*wz1*mzFe_ytb(it1+1,iz1) &
               +wt2*wz2*mzFe_ytb(it1+1,iz1+1))*(ms0p*MUSM)
              tmzZ_sn(pn)=(wt1*wz1*mzZ_ytb(it1,iz1) &
               +wt1*wz2*mzZ_ytb(it1,iz1+1)+wt2*wz1*mzZ_ytb(it1+1,iz1) &
               +wt2*wz2*mzZ_ytb(it1+1,iz1+1))*(ms0p*MUSM)
! *** number of stars with stellar wind x dt ***
! because age is modified, cannot calculate mean one.
! use the time at age
!              dtsw=age(pn)/TMUGYR+ts_p(pn)-TM_tot
!              agesw=age(pn)-0.5d0*dtsw*TMUGYR
!              if(agesw.ge.t_ytb(1)) then 
!                it1=int((dlog10(agesw)-dlog10(t_ytb(1)))/dt)+1
!              else
!                it1=0
!              endif
!              wt2=(agesw-t_ytb(it1))/(t_ytb(it1+1)-t_ytb(it1))
!              wt1=1.0d0-wt2
              dtsw=TM_dt
              nsw_sn(pn)=(wt1*wz1*nsw_ytb(it1,iz1) &
               +wt1*wz2*nsw_ytb(it1,iz1+1)+wt2*wz1*nsw_ytb(it1+1,iz1) &
               +wt2*wz2*nsw_ytb(it1+1,iz1+1))*(ms0p*MUSM)*dtsw
            endif
          endif
        endif
      enddo

      deallocate(zg)
      deallocate(age)

end subroutine
#endif
