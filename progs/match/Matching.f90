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

program Matching

  use Filter_types
  use GenParam_types_flt
  use DataSpace_types_flt
  use ReadData_mod
  use ReadParams_mod
  use BuildAdaptiveFilter_mod
  use ComputeAdaptiveFilter_mod
  use sep

  implicit none

  type(cube)   :: wght! Data space weight
  type(cube)   :: obs ! Observed data to be matched
  type(cube)   :: mod ! Modeled data
  type(cube)   :: fmod! Modeled data times filter
  type(GenPar_flt) :: par
  type(NSfilter)   :: nmatch
  type(NSfilter)   :: rough
  call sep_init(SOURCE)

  par%wghtag='weight'
  par%obstag='obs'
  par%modtag='mod'
  par%fmodtag='fmod'
  par%filttag='filt'
  par%filtpchtag='filtpch'
  par%ndim=sep_dimension(par%obstag)

  if(par%ndim.ne.3) call erexit('ERROR: only working on 3D cubes, exit now')

  call ReadData_dim(par%obstag,obs,par)
  call ReadData_dim(par%modtag,mod,par)
  call InfoData_dim(par%obstag,obs,par)
  if (.not.hdrs_are_consistent(obs,mod)) call erexit('ERROR: obs and mod do not have same dimensions, exit now')
  call ReadData_cube(par%obstag,obs)
  call ReadData_cube(par%modtag,mod)
  call from_param('thresh_d',par%thresh_d,maxval(abs(obs%dat))/100)

  if (exist_file('weight')) then
     call ReadData_dim(par%wghtag,wght,par)
     if (.not.hdrs_are_consistent(obs,wght)) call erexit('ERROR: obs and mod do not have same dimensions, exit now')
     call ReadData_cube(par%wghtag,wght)
  end if

  allocate(fmod%d(size(mod%d)),fmod%n(size(mod%n)),fmod%o(size(mod%o)))
  fmod%d=mod%d; fmod%n=mod%n; fmod%o=mod%o

  call readparams(par,nmatch)

  call psize_init(obs,par,nmatch)
  call pch_init(obs,par,nmatch)
  call create_nsmatch_filter(obs,par,nmatch)

  allocate(rough%npatch(3))
  rough%npatch=nmatch%npatch
  rough%ncoef=nmatch%ncoef
  call create_lap_3d(rough,par%nlaplac)

  if (par%prec) then
     if (exist_file('weight')) then
        call ComputeAdaptiveFilterPrec_op(rough,nmatch,par,obs,mod,fmod,wght)
     else
        call ComputeAdaptiveFilterPrec_op(rough,nmatch,par,obs,mod,fmod)
     end if
  else
     if (exist_file('weight')) then
        call ComputeAdaptiveFilter_op(rough,nmatch,par,obs,mod,fmod,wght)
     else
        call ComputeAdaptiveFilter_op(rough,nmatch,par,obs,mod,fmod)
     end if
  end if

  call WriteData_dim(par%fmodtag,fmod,par)
  call WriteData_cube(par%fmodtag,fmod)

  call cube_deallocate(obs)
  call cube_deallocate(mod)
  call cube_deallocate(fmod)

  call NSfilter_deallocate(rough)
  call NSfilter_deallocate(nmatch)

end program Matching
