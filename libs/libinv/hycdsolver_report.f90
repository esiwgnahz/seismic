# 1 "<stdin>"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "<stdin>"
module hycdsolver_report_mod
  use ddot_mod
  use hpf_mod
  implicit none
  real, private :: scale
  contains
  subroutine hycdsolver_report(iter,x,g,rd,rm,eps)
    integer, intent(in) :: iter
    real, dimension(:), intent(in) :: x,g
    real, dimension(:), pointer :: rd,rm
    optional :: rm,eps
    real :: eps
    double precision :: dd,dm,dx,dg
    dd=hpf(rd)
    dx=dot_product( x, x)
    dg=dot_product( g, g)
    dm=0.
    if (present(rm)) then
      dm=hpf(rm)
    end if
    if (present(eps)) then
      dm=eps*dm
    end if
    if (iter.eq.1) then
      scale=(dd+dm)/1000.
      write(0,*) "DOT PRODUCT VALUES SCALED BY ",scale
    end if
    write (0,*) "iteration ", iter , "   res tot", int(&
      &(dd+dm)/scale) , "   res dat", int(dd/scale) &
      & , "   res mod", int(dm/scale) , "  &
      & mod x'x", int(dx/scale)     ,           "      grad g'g", int(dg/scale)
  end subroutine
end module
