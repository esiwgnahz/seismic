module picking
  implicit none

contains

  subroutine coppens(in,ou,nl)
    real, dimension(:), intent(in) :: in
    real, dimension(:), intent(out) :: ou
    integer i,j,n,nl
    real :: tmp,maxin
    real, dimension(:), allocatable :: E1,E2

    maxin=maxval(abs(in))
    n=size(in)

    allocate(E1(n),E2(n))

    E2(1)=in(1)**2
    do i=2,n
       E2(i)=E2(i-1)+in(i)**2
    end do
    
    do i=nl+1,n
       tmp=0.
       do j=i-nl+1,i
          tmp=tmp+in(j)**2
       end do
       E1(i)=tmp
    end do

    ! takes care of starting window, not pretty
    do i=1,nl
       E1(i)=E1(nl+1)
    end do

    ou=E1/(E2+0.2*maxin)

    deallocate(E1,E2)

  end subroutine coppens

  subroutine find_maxindex(in,nt)
    real, dimension(:), intent(in):: in
    integer ::nt,i

    nt=1
    do i=2,size(in)
       if (in(i).gt.in(nt)) then
          nt=i
       end if
    end do
  end subroutine find_maxindex

end module picking
