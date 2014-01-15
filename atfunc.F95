#include "gcdp.def"
! *****************************************************
!    atfunc.F for GCD+ ver. f03.0
! 03 Feb. 2013   written by D. Kawata
! ***************************************************** 

#ifdef COSM
subroutine makettable(ae)
      use gcdp_const
      use gcdp_system

      implicit none

      double precision,intent(in) :: ae
      integer i,j,it
      double precision da,dt,f(0:NATTABLE),fa(0:NATTABLE)

! *** Initialization ***
      at_ttb(0)=0.0d0
      ta_ttb(0)=0.0d0
      da = (ae-at_ttb(0))/dble(NATTABLE)            
! *** Integration ***
      do i=0,NATTABLE
        if(i.gt.0) then
          at_ttb(i) = at_ttb(i-1)+da
        endif          
        f(i)=dsqrt(at_ttb(i)/(SI_omg0-(SI_omg0+SI_lam0-1.0d0)*at_ttb(i) &
         +SI_lam0*at_ttb(i)*at_ttb(i)*at_ttb(i)))
      enddo

      do i=1,NATTABLE
        ta_ttb(i)=ta_ttb(i-1)+0.5*da*(f(i-1)+f(i))*(H0_1/SI_h0)/TMUYR
      enddo

      TM_t0 = ta_ttb(NATTABLE)
! *** change the grid in t ***
      do i=1,NATTABLE
        f(i) = ta_ttb(i) 
        fa(i) = at_ttb(i) 
      enddo
      dt = (TM_t0-ta_ttb(0))/dble(NATTABLE)            
      do i=1,NATTABLE
        ta_ttb(i) = ta_ttb(i-1)+dt 
        do j=1,NATTABLE
          if(f(j).ge.ta_ttb(i)) then
            it=j
            goto 90
          endif
        enddo
        if(myrank.eq.0) then
          write(6,*) ' Error in updatesysha() :'
          write(6,*) '   Cannot find a(t) from table at ',ta_ttb(i)
          write(6,*) ' it = NATTABLE',it
        endif
        it=NATTABLE
   90   if(ta_ttb(i).gt.f(it).or.ta_ttb(i).lt.f(it-1)) then
          if(myrank.eq.0) then
            write(6,*) ' Error in makettabel() :'
            write(6,*) '  tn,ta(i),ta(i-1):',ta_ttb(i),f(it),f(it-1)
          endif
        endif
        at_ttb(i) = fa(it-1)+(fa(it)-fa(it-1)) &
         *(ta_ttb(i)-f(it-1))/(f(it)-f(it-1))
#ifdef BOXSIM
! 1/a
        ia_ttb(i)=1.0d0/at_ttb(i)
! 1/a^2
        ia2_ttb(i)=1.0d0/(at_ttb(i)**2)
#endif
      enddo

end subroutine


! *****   Definition of tfroma()   *****
subroutine tfroma(an,tn)
      use gcdp_const
      use gcdp_system

      implicit none

      double precision,intent(in) :: an
      double precision,intent(out) :: tn
      integer i,it

      do i=0,NATTABLE
        if(at_ttb(i).ge.an) then
          it=i
          goto 99
        endif
      enddo
      if(myrank.eq.0) then
        write(6,*) ' Error in updatesysha() :'
        write(6,*) '   Cannot find t(a) from table at ',an
        write(6,*) ' set it = NATTABLE'
      endif
      it=NATTABLE
   99 if(an.gt.at_ttb(it).or.an.lt.at_ttb(it-1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error in tfroma() :'
          write(6,*) '  tn,ta(i),ta(i-1):',an,at_ttb(it),at_ttb(it-1)
       endif
      endif
      tn = ta_ttb(it-1)+(ta_ttb(it)-ta_ttb(it-1)) &
       *(an-at_ttb(it-1))/(at_ttb(it)-at_ttb(it-1))
      
end subroutine

! *****   Definition of afromt()   *****
subroutine afromt(tn,an)
      use gcdp_const
      use gcdp_system

      implicit none

      double precision,intent(in) :: tn
      double precision,intent(out) :: an
      integer i,it

      it = int((tn-ta_ttb(0))/((ta_ttb(NATTABLE)-ta_ttb(0))/dble(NATTABLE)))+1
      if(it.gt.NATTABLE) then
        if(myrank.eq.0) then
        write(6,*) ' Error in afromt() :'
        write(6,*) '  Cannot find t(a) from table at t=',tn
        write(6,*) ' Set it = ',NATTABLE
        endif
        it = NATTABLE
      endif      
      if(tn.gt.ta_ttb(it).or.tn.lt.ta_ttb(it-1)) then
        if(myrank.eq.0) then
          write(6,*) ' Error in afromt() :'
          write(6,160) tn,ta_ttb(it),ta_ttb(it-1)
        endif
      endif
  160 format('  tn,ta(i),ta(i-1):',3(1pE13.5))
      an = at_ttb(it-1)+(at_ttb(it)-at_ttb(it-1)) &
        *(tn-ta_ttb(it-1))/(ta_ttb(it)-ta_ttb(it-1))

end subroutine

#ifdef BOXSIM
! *****   Definition of hubfroma()   *****
subroutine hubfroma(an,hub)
      use gcdp_const
      use gcdp_system

      implicit none

      double precision,intent(in) :: an
      double precision,intent(out) :: hub
      double precision iai

      iai = 1.0d0/an
      hub = SI_h0*HUB0*dsqrt(SI_omg0*iai*iai*iai &
       -(SI_omg0+SI_lam0-1.0)*iai*iai+SI_lam0)

end subroutine
#endif 

! 1/a integration
subroutine integrate1oa(np,t0,t1,in1oa)
      use gcdp_const
      use gcdp_system

      implicit none

      integer,intent(in) :: np
      double precision,intent(in) :: t0(0:np),t1(0:np)
      double precision,intent(out) :: in1oa(0:np)
      integer i,ip,it0,it1,nit
      double precision ia0,ia1,pian,ian,ptn,tn

      do ip=0,np-1
        it0 = int((t0(ip)-ta_ttb(0))/((ta_ttb(NATTABLE)-ta_ttb(0)) &
          /dble(NATTABLE)))+1
        if(it0.gt.NATTABLE) then
!          if(myrank.eq.0) then
!          write(6,*) ' Error in integrate1oa() :'
!          write(6,*) '  Cannot find t(a) from table at t0=',t0(i),i
!          write(6,*) ' Set it = ',NATTABLE
!          endif
          it0 = NATTABLE
        endif      
        it1 = int((t1(ip)-ta_ttb(0))/((ta_ttb(NATTABLE)-ta_ttb(0)) &
          /dble(NATTABLE)))+1
        if(it1.gt.NATTABLE) then
          it1 = NATTABLE
        endif      

        ia0 = ia_ttb(it0-1)+(ia_ttb(it0)-ia_ttb(it0-1)) &
          *(t0(ip)-ta_ttb(it0-1))/(ta_ttb(it0)-ta_ttb(it0-1))
        ia1 = ia_ttb(it1-1)+(ia_ttb(it1)-ia_ttb(it1-1)) &
          *(t1(ip)-ta_ttb(it1-1))/(ta_ttb(it1)-ta_ttb(it1-1))

        if(it0.eq.it1) then
          in1oa(ip)=(ia1-ia0)*(t1(ip)-t0(ip))
        else
          in1oa(ip)=0.0d0
          tn=t0(ip)
          ian=ia0
          do i=it0,it1
            ptn=tn
            pian=ian
            if(i.ne.it1) then
              tn=ta_ttb(i)
              ian=ia_ttb(i)
            else
              tn=t1(ip)
              ian=ia1
            endif
            in1oa(ip)=in1oa(ip)+(ian-pian)*(tn-ptn)
          enddo
        endif
      enddo
end subroutine

! 1/a^2 integration
subroutine integrate1oa2(np,t0,t1,in1oa2)
      use gcdp_const
      use gcdp_system

      implicit none

      integer,intent(in) :: np
      double precision,intent(in) :: t0(0:np),t1(0:np)
      double precision,intent(out) :: in1oa2(0:np)
      integer i,ip,it0,it1,nit
      double precision ia20,ia21,pian2,ian2,ptn,tn

      do ip=0,np-1
        it0 = int((t0(ip)-ta_ttb(0))/((ta_ttb(NATTABLE)-ta_ttb(0)) &
          /dble(NATTABLE)))+1
        if(it0.gt.NATTABLE) then
!          if(myrank.eq.0) then
!          write(6,*) ' Error in integrate1oa() :'
!          write(6,*) '  Cannot find t(a) from table at t0=',t0(i),i
!          write(6,*) ' Set it = ',NATTABLE
!          endif
          it0 = NATTABLE
        endif      
        it1 = int((t1(ip)-ta_ttb(0))/((ta_ttb(NATTABLE)-ta_ttb(0)) &
          /dble(NATTABLE)))+1
        if(it1.gt.NATTABLE) then
          it1 = NATTABLE
        endif      

        ia20 = ia2_ttb(it0-1)+(ia2_ttb(it0)-ia2_ttb(it0-1)) &
          *(t0(ip)-ta_ttb(it0-1))/(ta_ttb(it0)-ta_ttb(it0-1))
        ia21 = ia2_ttb(it1-1)+(ia2_ttb(it1)-ia2_ttb(it1-1)) &
          *(t1(ip)-ta_ttb(it1-1))/(ta_ttb(it1)-ta_ttb(it1-1))

        if(it0.eq.it1) then
          in1oa2(ip)=(ia21-ia20)*(t1(ip)-t0(ip))
        else
          in1oa2(ip)=0.0d0
          tn=t0(ip)
          ian2=ia20
          do i=it0,it1
            ptn=tn
            pian2=ian2
            if(i.ne.it1) then
              tn=ta_ttb(i)
              ian2=ia2_ttb(i)
            else
              tn=t1(ip)
              ian2=ia21
            endif
            in1oa2(ip)=in1oa2(ip)+(ian2-pian2)*(tn-ptn)
          enddo
        endif
      enddo
end subroutine


#endif
