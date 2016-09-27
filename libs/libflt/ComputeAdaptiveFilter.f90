module ComputeAdaptiveFilter_mod

  use Filter_types
  use GenParam_types_flt
  use DataSpace_types_flt

  use hycdsolver_reg_mod
  use hycdsolver_prc_mod
  use hycdstep_mod
  use identity_mod

  use ncnhconest 
  use nhelicon_mod
  use ncnhelicon
  use npolydiv

  implicit none

contains

  subroutine ComputeAdaptiveFilter_op(rough,filter,param,obs,mod,fmod)
    type(cube)       :: obs    ! Observed data to be matched
    type(cube)       :: mod    ! Modeled data
    type(cube)       :: fmod   ! Modeled data
    type(GenPar_flt) :: param  ! Parameter file
    type(NSfilter)   :: filter ! Filters 
    type(NSfilter)   :: rough  ! Regularization
    
    real, dimension(:), allocatable :: filter_coefs,filter_coefsout

    integer :: i, j, stat

    write(0,*) 'INFO: NS Filter estimation with regularization'
    call nhelicon_mod_init(rough%nmatch)
    call ncnhconest_init(mod%dat,filter%nmatch,filter%ncoef)
  
    allocate(filter_coefs(filter%ncoef*size(filter%nmatch%hlx)))
    allocate(filter_coefsout(filter%ncoef*size(filter%nmatch%hlx)))

    call identity_init(param%thresh_d,param%thresh_m)

    filter_coefs=0.

    call hycdsolver_reg(m=filter_coefs,          &
    &                   d=obs%dat,               &
    &                   Fop=ncnhconest_lop,      &
    &                   Aop=nhelicon_mod_lop,    &
    &                   Wop=identityd_lop,       &
    &                   Wmop=identitym_lop,      &
    &                   stepper=hycdstep2,       &
    &                   nAop=size(filter_coefs), &
    &                   niter=param%niter,       &
    &                   eps=param%eps,           &
    &                   verb=.true.)

    !! Copy coefficients into NS filter type
    do i=1,size(filter%nmatch%hlx)
       filter%nmatch%hlx(i)%flt(1:filter%ncoef)=filter_coefs(1+(i-1)*filter%ncoef:i*filter%ncoef)
    end do
    
!    do i=1,size(filter%nmatch%hlx)
!       write(0,*) i,filter%nmatch%hlx(i)%flt(15),filter%nmatch%hlx(i)%lag(15)
!    end do
    deallocate(filter_coefs)
    deallocate(filter_coefsout)

    deallocate(obs%dat)
    allocate(fmod%dat(fmod%n(1)*fmod%n(2)*fmod%n(3)))
    
    call ncnhelicon_init(filter%nmatch)
    stat=ncnhelicon_lop(.false.,.false.,mod%dat,fmod%dat)
    
    call NSfilter_write_to_file(param%filttag,param%filtpchtag,filter,mod)

  end subroutine ComputeAdaptiveFilter_op

  subroutine ComputeAdaptiveFilterPrec_op(rough,filter,param,obs,mod,fmod)
    type(cube)       :: obs    ! Observed data to be matched
    type(cube)       :: mod    ! Modeled data
    type(cube)       :: fmod   ! Modeled data
    type(GenPar_flt) :: param  ! Parameter file
    type(NSfilter)   :: filter ! Filters 
    type(NSfilter)   :: rough  ! Regularization
    
    real, dimension(:), allocatable :: filter_coefs,filter_coefsout

    integer :: i, j, stat

    write(0,*) 'INFO: NS Filter estimation with preconditioning'
    call ncnhconest_init(mod%dat,filter%nmatch,filter%ncoef)  

    allocate(filter_coefs(filter%ncoef*size(filter%nmatch%hlx)))
    allocate(filter_coefsout(filter%ncoef*size(filter%nmatch%hlx)))
 
    call npolydiv_init(size(filter_coefs),rough%nmatch)
    call identity_init(param%thresh_d,param%thresh_m)

    filter_coefs=0.

    call hycdsolver_prc(m=filter_coefs,          &
    &                   d=obs%dat,               &
    &                   Fop=ncnhconest_lop,      &
    &                   Sop=npolydiv_lop,        &
    &                   Wop=identityd_lop,       &
    &                   Wmop=identitym_lop,      &
    &                   stepper=hycdstep2,       &
    &                   nSop=size(filter_coefs), &
    &                   niter=param%niter,       &
    &                   eps=param%eps,           &
    &                   verb=.true.)

    !! Copy coefficients into NS filter type
    do i=1,size(filter%nmatch%hlx)
       filter%nmatch%hlx(i)%flt(1:filter%ncoef)=filter_coefs(1+(i-1)*filter%ncoef:i*filter%ncoef)
    end do

    deallocate(filter_coefs)
    deallocate(filter_coefsout)

    deallocate(obs%dat)
    allocate(fmod%dat(fmod%n(1)*fmod%n(2)*fmod%n(3)))

    call ncnhelicon_init(filter%nmatch)
    call NSfilter_write_to_file(param%filttag,param%filtpchtag,filter,mod)
    
  end subroutine ComputeAdaptiveFilterPrec_op

end module ComputeAdaptiveFilter_mod
