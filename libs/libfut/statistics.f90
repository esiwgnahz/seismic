module statistics

  implicit none

contains

  function mean(x)
    real mean
    real, dimension(:), intent(in) :: x

    mean=sum(x)/size(x)
  end function mean

  function std(x)
    real std
    real, dimension(:), intent(in) :: x
    real xm

    xm=mean(x)
    std=sqrt(sum((x-xm)**2)/(size(x)-1))

  end function std

  ! n is window size
  ! in is input trace
  ! ou is output trace after Edge-preserving filtering
  subroutine running_median(in,ou,n)
    real, dimension(:), intent(in):: in
    real, dimension(:), intent(out):: ou

    real, dimension(:), allocatable :: stdin
    real, dimension(:), allocatable :: meain
    integer :: n,nt,i,j

    nt=size(in)

    allocate(stdin(nt))
    allocate(meain(nt))
    stdin=0.
    meain=0.

    do i=1,nt-n
       stdin(i)=std(in(i:i+n))
       meain(i)=mean(in(i:i+n))
    end do

    !Mirror the end of the standard deviation and mean
    do i=1,n
       stdin(nt-n+i)=stdin(nt-n-i)
       meain(nt-n+i)=meain(nt-n-i)
    end do

    do i=n,nt
       ou(i)=meain(i)
       do j=1,n-1
          if (stdin(i-j).lt.stdin(i-j+1)) then
             ou(i)=meain(i-j)
          end if
       end do
    end do

    ! Mirroring ends
    do i=1,n-1
       ou(i)=ou(2*n-i)
    end do

    deallocate(stdin,meain)

  end subroutine running_median

  subroutine derivative(in,ou)
    real, dimension(:), intent(in) :: in
    real, dimension(:), intent(out):: ou

    integer:: i
    ou=0
    do i=1,size(in)-1
       ou(i)=in(i+1)-in(i)
    end do

  end subroutine derivative
end module statistics
