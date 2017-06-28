! *****************************************************
!  Definition of Function for Paralell Computing
!   03  Feb.  2013    written by D.KAWATA
! *****************************************************

  subroutine para_range(n1,n2,np,irank,ista,iend)
      integer,intent(in) :: n1,n2,np,irank
      integer,intent(out) :: ista,iend	  
      integer iwork1,iwork2
	  
      iwork1 = (n2-n1+1)/np
      iwork2 = mod(n2-n1+1,np)
      ista=irank*iwork1+n1+min(irank,iwork2)
      iend=ista+iwork1-1
      if(iwork2.gt.irank) then
        iend=iend+1
      endif

end subroutine








