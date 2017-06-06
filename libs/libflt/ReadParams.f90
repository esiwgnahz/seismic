! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module ReadParams_mod
  
  use GenParam_types_flt
  use Filter_types
  use sep

  implicit none
contains
  subroutine readparams(param,nsmatch)
    type(GenPar_flt)::  param
    type(NSfilter)  ::        nsmatch
  
    allocate(nsmatch%npatch(param%ndim))
    allocate(nsmatch%nfilt(param%ndim))
    call from_param('npatch',nsmatch%npatch,(/1,1,1/))
    call from_param('nsp',param%nsp,nsmatch%npatch)
    call from_param('nfilt',nsmatch%nfilt,(/1,1,1/))
    call from_param('niter',param%niter,product(nsmatch%nfilt))

    call from_param('eps',param%eps,0.)
    call from_param('thresh_m',param%thresh_m,1.)
    call from_param('n_laplac',param%nlaplac,5)

    call from_param('sparse',param%sparse,.false.)
    if (.not.param%sparse) then
       call from_param('prec',param%prec,.false.)
       call from_param('hyperbolic',param%hyperbolic,.false.)
    else
       param%prec=.false.
       param%hyperbolic=.true.
!       if (param%eps.eq.0) call erexit('ERROR: epsilon cannot be zero with sparse=1, exit now')
    end if

    write(0,*) 'INFO:-------------------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO: sparse=',param%sparse
    write(0,*) 'INFO: preconditioning=',param%prec
    write(0,*) 'INFO: hyperbolic norm=',param%hyperbolic
    write(0,*) 'INFO:'
    write(0,*) 'INFO:-------------------------------------'

    call from_param('num_threads',param%nthreads,4)
    call omp_set_num_threads(param%nthreads)

  end subroutine readparams
end module ReadParams_mod
