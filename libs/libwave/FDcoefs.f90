! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module FDcoefs_assign

  implicit none
  
contains

  subroutine FD_acoustic_rho_init_coefs(coef)

    use FD_types
    type(UnscaledFDcoefs) ::  coef

    coef%c0 =  0.
    coef%c1 =  1225./1024.
    coef%c2 = -1225./(1024.*15.)
    coef%c3 =  1225./(1024.*125.)
    coef%c4 = -1225./(1024.*1715.)
    
  end subroutine FD_acoustic_rho_init_coefs

  subroutine FD_acoustic_init_coefs(coef)

    use FD_types
    type(UnscaledFDcoefs) ::  coef

    coef%c0=-205./72.
    coef%c1=8./5
    coef%c2=-1./5.
    coef%c3=8./315.
    coef%c4=-1./560.
   
  end subroutine FD_acoustic_init_coefs

end module FDcoefs_assign
