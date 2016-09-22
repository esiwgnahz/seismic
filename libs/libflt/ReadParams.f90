module ReadParams_mod
  
  use GenParam_types_flt
  use sep

  implicit none
contains
  subroutine readparams(param)
    type(GenPar_flt)::  param
  end subroutine readparams
end module ReadParams_mod
