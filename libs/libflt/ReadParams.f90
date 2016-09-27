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
    call from_param('nfilt',nsmatch%nfilt,(/1,1,1/))
    call from_param('niter',param%niter,product(nsmatch%nfilt))

    call from_param('eps',param%eps,0.)
    call from_param('thresh_m',param%thresh_m,1.)
    call from_param('n_laplac',param%nlaplac,5)
    call from_param('prec',param%prec,.false.)
    
  end subroutine readparams
end module ReadParams_mod
