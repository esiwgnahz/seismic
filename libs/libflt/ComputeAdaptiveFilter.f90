module ComputeAdaptiveFilter_mod

  use Filter_types
  use ReadData_mod
  use GenParam_types_flt
  use DataSpace_types_flt

  use hycdsolver_reg_mod
  use hycdsolver_prc_mod
  use hycdstep_mod
  use identity_mod

  use solver_reg_mod
  use solver_prc_mod
  use cgstep_mod

  use ncnhconest 
  use nhelicon_mod
  use ncnhelicon
  use npolydiv

  implicit none

contains

  subroutine ComputeAdaptiveFilter_op(rough,filter,param,obs,mod,fmod,wght)
    optional         :: wght
    type(cube)       :: wght
    type(cube)       :: obs    ! Observed data to be matched
    type(cube)       :: mod    ! Modeled data
    type(cube)       :: fmod   ! Modeled data
    type(GenPar_flt) :: param  ! Parameter file
    type(NSfilter)   :: filter ! Filters 
    type(NSfilter)   :: rough  ! Regularization
    type(cube)       :: illum  ! Modeled illumination
    
    real, dimension(:), allocatable :: filter_coefs,filter_coefsout

    integer :: i, j, stat
    double precision :: memory_needed

    write(0,*) 'INFO: NS Filter estimation with regularization'
    call nhelicon_mod_init(rough%nmatch)
    call ncnhconest_init(mod%dat,filter%nmatch,filter%ncoef)
  
    write(0,*) 'INFO: Done with initialization'
    allocate(filter_coefs(filter%ncoef*size(filter%nmatch%hlx)))
    allocate(filter_coefsout(filter%ncoef*size(filter%nmatch%hlx)))
    write(0,*) 'INFO: Done with filter allocation'

    if (present(wght).and.param%hyperbolic) then
       wght%dat=wght%dat/param%thresh_d
    endif

    if (param%hyperbolic) call identity_init(param%thresh_d,param%thresh_m)
    if (present(wght))    call weightd_init(wght%dat)


    filter_coefs=0.

    memory_needed=(dble(size(filter_coefs)*5)+ &  ! Model (in,out,solver)
    &              dble(size(obs%dat)*2)     + &  ! Data  (mod,obs)
    &              dble((size(obs%dat)+size(filter_coefs))*5))*1e-9*4 ! Data+Model (solver)
    write(0,*) 'INFO:'
    write(0,*) 'INFO: Memory needed for inversion =',memory_needed,'Gb'

    if (param%hyperbolic) then
       if (present(wght)) then
          write(0,*) 'INFO: Starting inversion with hyperbolic conjugate direction solver with weight'
          call hycdsolver_reg(m=filter_coefs,          &
          &                   d=obs%dat,               &
          &                   Fop=ncnhconest_lop,      &
          &                   Aop=nhelicon_mod_lop,    &
          &                   Wop=weightd_lop,       &
          &                   Wmop=identitym_lop,      &
          &                   stepper=hycdstep2,       &
          &                   nAop=size(filter_coefs), &
          &                   niter=param%niter,       &
          &                   eps=param%eps,           &
          &                   verb=.true.)
          call weightd_close()
       else
          write(0,*) 'INFO: Starting inversion with hyperbolic conjugate direction solver'
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
       end if
    else
       if (present(wght)) then
          write(0,*) 'INFO: Starting inversion with conjugate direction solver with weight'
          call solver_reg(m=filter_coefs,          &
          &               d=obs%dat,               &
          &               Fop=ncnhconest_lop,      &
          &               Aop=nhelicon_mod_lop,    &
          &               Wop=weightd_lop,         &
          &               stepper=cgstep,          &
          &               nAop=size(filter_coefs), &
          &               niter=param%niter,       &
          &               eps=param%eps,           &
          &               verb=.true.)
          call weightd_close()
       else
          write(0,*) 'INFO: Starting inversion with conjugate direction solver'
          call solver_reg(m=filter_coefs,          &
          &               d=obs%dat,               &
          &               Fop=ncnhconest_lop,      &
          &               Aop=nhelicon_mod_lop,    &
          &               stepper=cgstep,          &
          &               nAop=size(filter_coefs), &
          &               niter=param%niter,       &
          &               eps=param%eps,           &
          &               verb=.true.)
       end if

       call cgstep_close()
    end if
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
    
    allocate(illum%dat(size(mod%dat)))
    allocate(illum%d(size(mod%d)),illum%n(size(mod%n)),illum%o(size(mod%o)))
    illum%d=mod%d; illum%n=mod%n; illum%o=mod%o
    illum%dat=0.

    call PutSpikesInData(mod,param)
    stat=ncnhelicon_lop(.false.,.false.,mod%dat,illum%dat)
!    illum%dat=mod%dat
    call WriteData_dim('impulse',illum,param)
    call WriteData_cube('impulse',illum)

    ! Illumination
    do i=1,size(filter%nmatch%hlx)
       do j=1,size(filter%nmatch%hlx(i)%lag)
          if (filter%nmatch%hlx(i)%lag(j).ne.0) filter%nmatch%hlx(i)%flt(j)=0
       end do
    end do

    illum%dat=0.
    call ncnhelicon_init(filter%nmatch)
    mod%dat=1.
    stat=ncnhelicon_lop(.false.,.false.,mod%dat,illum%dat)
    
    call WriteData_dim('illum',illum,param)
    call WriteData_cube('illum',illum)

    call cube_deallocate(illum)
    call NSfilter_write_to_file(param%filttag,param%filtpchtag,filter,mod)

  end subroutine ComputeAdaptiveFilter_op

  subroutine ComputeAdaptiveFilterPrec_op(rough,filter,param,obs,mod,fmod,wght)
    optional         :: wght
    type(cube)       :: wght
    type(cube)       :: obs    ! Observed data to be matched
    type(cube)       :: mod    ! Modeled data
    type(cube)       :: fmod   ! Modeled data
    type(GenPar_flt) :: param  ! Parameter file
    type(NSfilter)   :: filter ! Filters 
    type(NSfilter)   :: rough  ! Regularization
    type(cube)       :: illum  ! Modeled illumination
    
    real, dimension(:), allocatable :: filter_coefs,filter_coefsout

    integer :: i, j, stat
    double precision :: memory_needed

    write(0,*) 'INFO: NS Filter estimation with preconditioning'
    call ncnhconest_init(mod%dat,filter%nmatch,filter%ncoef)  

    allocate(filter_coefs(filter%ncoef*size(filter%nmatch%hlx)))
    allocate(filter_coefsout(filter%ncoef*size(filter%nmatch%hlx)))
 
    call npolydiv_init(size(filter_coefs),rough%nmatch)
    
    if (present(wght).and.param%hyperbolic) then
       wght%dat=wght%dat/param%thresh_d
    endif

    if (param%hyperbolic) call identity_init(param%thresh_d,param%thresh_m)
    if (present(wght))    call weightd_init(wght%dat)


    filter_coefs=0.

    memory_needed=(dble(size(filter_coefs)*5)+ &  ! Model (in,out,solver)
    &              dble(size(obs%dat)*2)     + &  ! Data  (mod,obs)
    &              dble((size(obs%dat)+size(filter_coefs))*5))*1e-9*4 ! Data+Model (solver)

    write(0,*) 'INFO:'
    write(0,*) 'INFO: Memory needed for inversion =',memory_needed,'Gb'

    if (param%hyperbolic) then
       if (present(wght)) then
          write(0,*) 'INFO: Starting inversion with hyperbolic conjugate direction solver with weight'
          call hycdsolver_prc(m=filter_coefs,          &
          &                   d=obs%dat,               &
          &                   Fop=ncnhconest_lop,      &
          &                   Sop=npolydiv_lop,        &
          &                   Wop=weightd_lop,       &
          &                   Wmop=identitym_lop,      &
          &                   stepper=hycdstep2,       &
          &                   nSop=size(filter_coefs), &
          &                   niter=param%niter,       &
          &                   eps=param%eps,           &
          &                   verb=.true.)
          call weightd_close()
       else
          write(0,*) 'INFO: Starting inversion with hyperbolic conjugate direction solver'
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
       end if
    else
       if (present(wght)) then
          write(0,*) 'INFO: Starting inversion with conjugate direction solver with weight'
          call solver_prc(m=filter_coefs,          &
          &               d=obs%dat,               &
          &               Fop=ncnhconest_lop,      &
          &               Wop=weightd_lop,         &
          &               Sop=npolydiv_lop,        &
          &               stepper=cgstep,       &
          &               nSop=size(filter_coefs), &
          &               niter=param%niter,       &
          &               eps=param%eps,           &
          &               verb=.true.)
          call weightd_close()
       else
          write(0,*) 'INFO: Starting inversion with conjugate direction solver'
          call solver_prc(m=filter_coefs,          &
          &               d=obs%dat,               &
          &               Fop=ncnhconest_lop,      &
          &               Sop=npolydiv_lop,        &
          &               stepper=cgstep,       &
          &               nSop=size(filter_coefs), &
          &               niter=param%niter,       &
          &               eps=param%eps,           &
          &               verb=.true.)
       end if
       call cgstep_close()
    end if
       

    !! Copy coefficients into NS filter type
    do i=1,size(filter%nmatch%hlx)
       filter%nmatch%hlx(i)%flt(1:filter%ncoef)=filter_coefs(1+(i-1)*filter%ncoef:i*filter%ncoef)
    end do

    deallocate(filter_coefs)
    deallocate(filter_coefsout)

    deallocate(obs%dat)
    allocate(fmod%dat(fmod%n(1)*fmod%n(2)*fmod%n(3)))

    call ncnhelicon_init(filter%nmatch)
    stat=ncnhelicon_lop(.false.,.false.,mod%dat,fmod%dat)
    
    allocate(illum%dat(size(mod%dat)))
    allocate(illum%d(size(mod%d)),illum%n(size(mod%n)),illum%o(size(mod%o)))
    illum%d=mod%d; illum%n=mod%n; illum%o=mod%o
    illum%dat=0.

    call PutSpikesInData(mod,param)
    stat=ncnhelicon_lop(.false.,.false.,mod%dat,illum%dat)
!    illum%dat=mod%dat
    call WriteData_dim('impulse',illum,param)
    call WriteData_cube('impulse',illum)

    ! Illumination
    do i=1,size(filter%nmatch%hlx)
       do j=1,size(filter%nmatch%hlx(i)%lag)
          if (filter%nmatch%hlx(i)%lag(j).ne.0) filter%nmatch%hlx(i)%flt(j)=0
       end do
    end do

    illum%dat=0.
    call ncnhelicon_init(filter%nmatch)
    mod%dat=1.
    stat=ncnhelicon_lop(.false.,.false.,mod%dat,illum%dat)
    
    call WriteData_dim('illum',illum,param)
    call WriteData_cube('illum',illum)

    call cube_deallocate(illum)
    call NSfilter_write_to_file(param%filttag,param%filtpchtag,filter,mod)
    
  end subroutine ComputeAdaptiveFilterPrec_op

  subroutine PutSpikesInData(dat,param)
    type(cube)       ::      dat
    type(GenPar_flt) ::          param  ! Parameter file
    integer :: i,j,k


    dat%dat=0. ! First set everything to zero

    do k=dat%n(3)/(2*param%nsp(3)),(dat%n(3)/param%nsp(3))*(0.5+(param%nsp(3)-1)),dat%n(3)/param%nsp(3)
       do j=dat%n(2)/(2*param%nsp(2)),(dat%n(2)/param%nsp(2))*(0.5+(param%nsp(2)-1)),dat%n(2)/param%nsp(2)
          do i=dat%n(1)/(2*param%nsp(1)),(dat%n(1)/param%nsp(1))*(0.5+(param%nsp(1)-1)),dat%n(1)/param%nsp(1)
             dat%dat(i+(j-1)*dat%n(1)+(k-1)*dat%n(1)*dat%n(2))=1.
          end do
       end do
    end do

  end subroutine PutSpikesInData

end module ComputeAdaptiveFilter_mod
