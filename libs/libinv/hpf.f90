! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module hpf_mod

  implicit none
  
contains
  
  function hpf(x) 
    real, dimension(:) ::x
    double precision   ::hpf
    hpf=sum(sqrt(1.d0+x**2)-1.d0)
  end function hpf
  
  function hpfp(x)  
    real, dimension(:) ::x
    real, dimension(size(x))::hpfp
    hpfp=x/sqrt(1.+x**2)
  end function hpfp
  
  function hpfs(x)  
    real, dimension(:) ::x
    real, dimension(size(x))::hpfs
    hpfs=1/(1.+x**2)**1.5
  end function hpfs
  
end module hpf_mod
