!!$=head1 NAME
!!$
!!$dottest - Perform a dot product test on an operator
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call dot_test(oper,n_mod,n_dat,dot1,dot2)>
!!$
!!$=head1 PARAMETERS
!!$
!!$=over 4
!!$
!!$=item oper - func
!!$
!!$       standard operator
!!$
!!$=item n_mod - integer
!!$
!!$      number of model points
!!$
!!$=item n_dat - integer
!!$
!!$      number of data points
!!$
!!$=item dot1  - real(2)
!!$
!!$      result of dot product test
!!$
!!$=item dot2  - real(2)
!!$
!!$      result of dot product test using add variable
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Test whether or an operator passes the dot product test (necessary for
!!$convergence with solver)
!!$
!!$=head1 SEE ALSO
!!$
!!$L<solver>
!!$
!!$=head1 LIBRARY
!!$
!!$B<geef90>
!!$
!!$=cut
module dottest  
  implicit none
  logical, parameter, private ::  T = .true.,  F = .false.
  contains
  subroutine dot_test( oper, n_mod, n_dat, dot1, dot2)
    integer,             intent (in)  :: n_mod, n_dat 
    real, dimension (2), intent (out) :: dot1, dot2
    interface  
      function oper( adj, add, mod, dat) result(stat)
        integer              :: stat
        logical, intent (in) :: adj, add
        real, dimension (:)  :: mod, dat
      end function
    end interface 
    real, dimension( n_mod)           :: mod1, mod2
    real, dimension( n_dat)           :: dat1, dat2
    integer                           :: stat
    call random_number( mod1)
    call random_number( dat2)
    stat = oper( F, F, mod1, dat1)
    dot1( 1) = dot_product( dat1, dat2)
    stat = oper( T, F, mod2, dat2)
    dot1( 2) = dot_product( mod1, mod2)
    write (0,*) dot1 
    stat = oper( F, T, mod1, dat1)
    dot2( 1) = dot_product( dat1, dat2)
    stat = oper( T, T, mod2, dat2)
    dot2( 2) = dot_product( mod1, mod2)
    write (0,*) dot2
  end subroutine 
end module 
