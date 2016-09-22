program Matching

  use GenParam_types_flt
  use DataSpace_types_flt
  use ReadData_mod
  use sep

  implicit none

  type(cube)   :: obs ! Observed data to be matched
  type(cube)   :: mod ! Modeled data
  type(cube)   :: fmod! Modeled data times filter
  type(GenPar_flt) :: par

  par%obstag='obs'
  par%modtag='mod'
  par%ndim=3

  call ReadData_dim(par%obstag,obs,par)
  call ReadData_dim(par%modtag,mod,par)
  if (.not.hdrs_are_consistent(obs,mod)) call erexit('ERROR: obs and mod do not have same dimensions, exit now')
  call ReadData_cube(par%obstag,obs)
  call ReadData_cube(par%modtag,mod)

  ! Code works as follows
  !
  ! 1 - we read the modeled and observed data
  ! 2 - we read the NS matching filters parameters
  ! 3 - we pose the inverse problem to find matching filters such that
  !     Q(f)=|D_obs - D_mod*f|_h+eps|Rf|_h is minimum
  !     For the norm, we are going to use the hyperbolic norm
  !     For the the regularization term, we are going to use the Laplacian operator
  ! 4 - we output two elements
  !        a - Modeled data x filter, or residual (to be decided)
  !        b - The filter parameters which can be used for filtering an image later on 
  !            (Need to have a code that will do that)
  ! 
  ! Antoine Guitton, Sept 2016, Bellevue Geophysics LLC, All Rights Reserved
  !
  
  
  call cube_deallocate(obs)
  call cube_deallocate(mod)

end program Matching
