module tcai1  
! Transient Convolution Adjoint Input 1-D. yy(m1+n1)
  use adj_mod
  implicit none
  real ,private,  dimension (:), pointer  :: bb 
  contains
  subroutine tcai1_init ( bb_in )
    real ,   dimension (:), pointer :: bb_in 
    bb => bb_in
  end subroutine
  function tcai1_lop ( adj, add, xx, yy) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: xx,yy 
    call adjnull (adj,add,xx,yy )
    call tcai1_lop2(adj,add,xx,yy )
    stat=0
  end function 
  subroutine tcai1_lop2(adj,add,xx,yy)
    logical,intent(in) :: adj,add 
    real, dimension (:)  :: xx 
    real, dimension (:)  :: yy 
    integer b, x, y
    if ( size(yy) < size (xx) + size(bb) - 1 ) then
      call erexit('tcai')
    end if
    do b= 1, size(bb)
      do x= 1, size(xx)
        y = x + b - 1
        if ( adj) then
          xx(x) =        xx(x) +  yy(y) * bb(b)
        else
          yy(y) =                 yy(y) +  xx(x) * bb(b)
        end if
      end do
    end do
  end subroutine 
  subroutine tcai1_close()
  end subroutine 
end module 
