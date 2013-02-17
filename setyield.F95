#include "gcdp.def"
!*********************************************
!   setyield.F for GCD+ ver f03.0
!   03 Feb. 2013    written by D.Kawata  
! *********************************************

#if defined(SF_EFD) || defined(SF_ZFD)
subroutine  setyield()
      use gcdp_const
      use gcdp_system
      use gcdp_yields

      implicit none

      integer i,j,k
      integer nzg
! *** metallicity and age ***
      double precision llz,luz,dlz,llt,lut,dlt
      double precision x,mu,ml,mlsnii
      character chead*22,cread*200

      open(50,file='./ini/tzy.cmp.dat',status='old')
      read(50,'(4(1pE13.5),I10)') x,mu,ml,mlsnii,SI_nsp
      read(50,'(a22)') chead
      if(myrank.eq.0) then
        write(6,*) 'IMF: x,mu,ml=',x,mu,ml
        write(6,*) ' SNe II mass limit, N mass group=',mlsnii,SI_nsp
        write(6,*) chead
      endif

      do j=0,NYTZ+1
        i=0
        t_ytb(0) =0.0d0
        mej_ytb(i,j)=0.0d0
        nsn_ytb(i,j)=0.0d0
        mzHe_ytb(i,j)=0.0d0
        mzC_ytb(i,j)=0.0d0
        mzN_ytb(i,j)=0.0d0
        mzO_ytb(i,j)=0.0d0
        mzNe_ytb(i,j)=0.0d0
        mzMg_ytb(i,j)=0.0d0
        mzSi_ytb(i,j)=0.0d0
        mzFe_ytb(i,j)=0.0d0
        mzZ_ytb(i,j)=0.0d0
        nsp_ytb(j)=0
! *** Unit Gyr ***
        do i=1,NYTT+1
          read(50,'(a200)') cread
          if(cread(1:1).eq.'#') then
            read(cread(2:14),'(1pE13.5)') tspi_ytb(nsp_ytb(j),j)
            nsp_ytb(j)=nsp_ytb(j)+1
            read(50,'(a200)') cread
          endif
          read(cread(1:195),150) t_ytb(i),z_ytb(j),mzHe_ytb(i,j) &
           ,mzC_ytb(i,j),mzN_ytb(i,j),mzO_ytb(i,j),mzNe_ytb(i,j) &
           ,mzMg_ytb(i,j),mzSi_ytb(i,j),mzFe_ytb(i,j),mzZ_ytb(i,j) &
           ,mej_ytb(i,j),nsn_ytb(i,j),nsw_ytb(i,j)
  150 format(14(1pE13.5))
! *** set initial number of SW stars ***
          if(i.eq.1) then
            nsw_ytb(0,j)=nsw_ytb(i,j)
          endif
        enddo
! *** set final mass group time ***
        tspi_ytb(nsp_ytb(j),j)=t_ytb(NYTT+1)
        nsp_ytb(j)=nsp_ytb(j)+1
      enddo
! *** read mgr group SNII end 
      read(50,'(a9,I10)') chead,SI_snii
      if(myrank.eq.0) then  
        write(6,*) ' mass group with SNeII<',SI_snii
      endif

      close(50)
! *** reset Metallicity ***
      llz = -6.0d0
      luz = 0.0d0
      dlz = (luz-llz)/dble(NYTZ)
      z_ytb(0)=0.0d0
      do  i=1,NYTZ+1
        z_ytb(i)=10.0d0**(llz+dlz*dble(i-1))
      enddo
! *** reset Age ***
      llt = -3.0d0
      lut = 1.2d0
      dlt = (lut-llt)/dble(NYTT)

      do  i=1,NYTT+1
        t_ytb(i)=10.0d0**(llt+dlt*dble(i-1))
      enddo

end subroutine
#endif

