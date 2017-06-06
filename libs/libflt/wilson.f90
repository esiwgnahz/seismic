! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
! Copyright 1999 Stanford University
!All rights reserved
!
!Author:  Sergey Fomel, Paul Sava, James Rickett, Jon Claerbout
!
module wilson_mod  
! Wilson's factorization
  use helicon_mod
  use polydiv_omp
  
  use omp_lib

  implicit none
  integer,                          private :: n
  real, dimension( :), allocatable, private :: auto, bb, cc, b, c 

  !$OMP THREADPRIVATE(n,auto,bb,cc,b,c)
   
  contains
  subroutine wilson_mod_init( nmax)
    integer, intent (in) :: nmax
    n = nmax
    allocate ( auto( 2*n-1), bb( 2*n-1), cc( 2*n-1), b(n), c(n))
  end subroutine 
  subroutine wilson_mod_factor( niter, s0, ss, a0, aa, verb)
    integer,       intent( in)  :: niter ! Newton iterations
    real,          intent( in)  :: s0    ! autocorrelation zero lag
    type( filter), intent( in)  :: ss    
    ! autocorrelation, other lags
    real,          intent( out) :: a0    ! factor, zero lag
    type( filter)               :: aa    ! factor, other lags
    logical,       intent( in)  :: verb
    optional                    :: verb
    real                        :: eps
    integer                     :: i, stat
    auto = 0.
    auto( n) = s0
    b( 1) =1.       ! initialize
    auto( n+ss%lag) =  ss%flt                   ! autocorrelation
    auto( n-ss%lag) =  ss%flt                   
    ! symmetrize input auto.
    call helicon_mod_init( aa)                      ! multiply polynoms
    call polydiv_omp_init( 2*n-1, aa)               ! divide   polynoms
    do i = 1, niter  
      stat= polydiv_omp_lop(.false.,.false., auto, bb)  ! bb = S/A
      stat= polydiv_omp_lop( .true.,.false., cc,   bb)  ! cc = S/(AA')
      b( 2:n) = 0.5*( cc( n+1:2*n-1  ) +cc( n-1:1    :-1)) / cc( n)   
! b = plusside(1+cc)
      eps = maxval( abs( b(2:n)))
! "L1 norm"
      if (present (verb)) then
        if (verb) then
          write (0,*) i, eps
        end if
      end if
      if ( eps < epsilon( a0)) then
        exit                  ! convergence
      end if
      stat= helicon_mod_lop( .false., .false., b, c)      ! c = A b
      aa%flt = c( 1+aa%lag)                           ! put on helix
    end do 
    a0 = sqrt( cc( n))
  end subroutine 
  subroutine wilson_mod_close ()
    deallocate( auto, bb, cc, b, c)
    call polydiv_omp_close()
  end subroutine 
end module 
