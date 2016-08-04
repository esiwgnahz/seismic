module chain0_mod  
  implicit none
  logical, parameter, private :: T = .true., F = .false.
  interface chain0
  module procedure chain20
  module procedure chain30
  end interface
  contains
  subroutine array0(op1,op2, adj,add, m,d1,d2)
!        ARRAY
    interface  
      integer function op1(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function 
      integer function op2(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
    end interface 
    logical, intent(in) :: adj, add
    real, dimension(:)  :: m,d1,d2
    integer             :: st
    if (adj) then
      st = op1(T, add, m, d1)		! m = m0 + A' d1
      st = op2(T,   T, m, d2)		! m = m  + B' d2
    else
      st = op1(F, add, m, d1)		! d1 = d1 + A  m
      st = op2(F, add, m, d2)		! d2 = d2 + B  m
    end if
  end subroutine 
  subroutine chain20(op1,op2, adj,add, m,d,t1)
!        CHAIN 2
    interface  
      integer function op1(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function 
      integer function op2(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
    end interface 
    logical, intent(in) :: adj, add
    real, dimension(:)  :: m,d, t1
    integer             :: st
    if (adj) then
      st = op1( T,   F, t1,d )		!        t =    A' d
      st = op2( T, add, m, t1)		! m = B' t = B' A' d
    else
      st = op2( F,   F, m, t1)		!        t =    B  m
      st = op1( F, add, t1,d )		! d = A  t = A  B  m
    end if
  end subroutine 
  subroutine chain30(op1,op2,op3, adj,add, m,d,t1,t2)
! CHAIN 3
    interface  
      integer function op1(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function 
      integer function op2(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function 
      integer function op3(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
    end interface 
    logical, intent(in) :: adj, add
    real, dimension(:)  :: m,d, t1,t2
    integer             :: st
    if (adj) then
      st = op1( T,   F, t2, d )	!                 t1 =       A' d
      st = op2( T,   F, t1, t2)	!         t2 = B' t1 =    B' A' d
      st = op3( T, add, m , t1)	! m  = C' t2 =         C' B' A' d
    else
      st = op3( F,   F, m , t1)	!                 t1 =       C  m
      st = op2( F,   F, t1, t2)	!         t2 = B  t1 =    B  C  m
      st = op1( F, add, t2, d )	! d  = A  t2 =         A  B  C  m
    end if
  end subroutine 
end module 
