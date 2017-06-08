      FUNCTION erfcc(x)

      double precision,intent(in) :: x

      double precision erfcc
      double precision t,z

      z=dabs(x)
      t=1.0d0/(1.0d0+0.5d0*z)
      erfcc=t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196+t &
      *(0.09678418d0+t*(-0.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t &
      *(1.48851587d0+t*(-0.82215223d0+t*.17087277d0))))))))) 
      if (x.lt.0.0d0) erfcc=2.0d0-erfcc
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
