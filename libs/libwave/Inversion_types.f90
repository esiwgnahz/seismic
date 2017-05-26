module Inversion_types

  use sep
  use ModelSpace_types

  implicit none

  type InversionParam
     integer:: iter       ! iteration number
     integer:: niter      ! Maximum number of iterations
     integer:: neval      ! Maximum number of fct/gdt eval.
     integer:: eval       ! number of fct/gdt eval.
     integer:: const_type ! constrained velocity at (1) each Fct/Gdt eval
                          !                         (2) each iteration
     logical:: freeze_soft! (0) velocity is preserved in mask area
                          ! (1) velocity is not strictly'0, enforced in masking area

     real, allocatable   :: parmin(:)
     real, allocatable   :: parmax(:)

     real, allocatable   :: sigma(:) ! Coefficients for gradient of second model vector
     
     logical:: invert_rho ! (0) Update rho
                          ! (1) Doesn't update rho

     integer:: vprho_param! (0) Vp/rho parameterization 
                          ! (1) Vp/Ip  parameterization 

     integer:: nparam     ! Number of parameters to invert

     logical:: wantreg    ! regularization parameter
     logical:: wantlog    ! regularization parameter, logistic
     real   :: illupow
     character(len=6) :: dat_nrm_type_char
     integer          :: dat_nrm_type 
     real             :: dat_thresh
     
     double precision, allocatable:: modinit(:,:)
     real,             allocatable:: modmask(:,:)

     integer :: n1
     integer :: ntotaltraces

  end type InversionParam

contains

  subroutine read_inv_params(invparam)
    type(InversionParam) ::  invparam
    logical              ::  withRho

    invparam%iter=0
    invparam%eval=0
    call from_param('niter',invparam%niter,10)
    call from_param('neval',invparam%neval,3*invparam%niter)
    
    call from_param('freeze_soft',invparam%freeze_soft,.true.)
    call from_param('bound',invparam%const_type,2)
    call from_param('wantreg',invparam%wantreg,.false.)
    call from_param('wantlog',invparam%wantlog,.false.)
    call from_param('illu_pow',invparam%illupow,1.)
    call from_param('data_nrm_type',invparam%dat_nrm_type_char,'L2norm')
    call from_param('data_threshold',invparam%dat_thresh,0.)

    if(invparam%dat_nrm_type_char(1:5).eq.'Cauch') invparam%dat_nrm_type=3
    if(invparam%dat_nrm_type_char(1:5).eq.'Huber') invparam%dat_nrm_type=12
    if(invparam%dat_nrm_type_char(1:5).eq.'L1nor') invparam%dat_nrm_type=1
    if(invparam%dat_nrm_type_char(1:5).eq.'L2nor') invparam%dat_nrm_type=2

    call from_param('withRho',withRho,.false.)
    invparam%nparam=1

    if (withRho) then
       call from_param('invert_rho',invparam%invert_rho,.false.)
       if (invparam%invert_rho) then
          invparam%nparam=2
          call from_param('vprho_param',invparam%vprho_param,0)
       end if
    end if

    if ((.not.withRho).and.(invparam%invert_rho)) call erexit('ERROR: withRho and invert_rho are incompatible, please reparameterize, exit now')

    allocate(invparam%parmin(invparam%nparam))
    allocate(invparam%parmax(invparam%nparam))
    allocate(invparam%sigma(invparam%nparam))
    invparam%sigma=1.
    call from_param('sigma1',invparam%sigma(1),1.)

    call from_param('vpmin',invparam%parmin(1),1500.)
    call from_param('vpmax',invparam%parmax(1),4500.)
    if (invparam%nparam.eq.2) then
    call from_param('sigma2',invparam%sigma(2),1.)
       call from_param('par2min',invparam%parmin(2),1000.)
       call from_param('par2max',invparam%parmax(2),3000.)
    end if

    write(0,*) 'INFO: ----------------------------'
    write(0,*) 'INFO:   Inversion Parameters      '
    write(0,*) 'INFO: ----------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO:   nrm_type   = ',invparam%dat_nrm_type_char
    if ((invparam%dat_nrm_type.eq.3).or.(invparam%dat_nrm_type.eq.12)) then
       write(0,*) 'INFO:  data thresh = ',invparam%dat_thresh
    end if
    write(0,*) 'INFO:   niter      = ',invparam%niter
    write(0,*) 'INFO:   neval      = ',invparam%neval
    write(0,*) 'INFO:   vpmin      = ',invparam%parmin(1)
    write(0,*) 'INFO:   vpmax      = ',invparam%parmax(1)
    if (withRho) then
       write(0,*) 'INFO:'
       write(0,*) 'INFO:  Inversion with density '
       write(0,*) 'INFO:'
       if (invparam%invert_rho) then
          write(0,*) 'INFO:   par2min    = ',invparam%parmin(2)
          write(0,*) 'INFO:   par2max    = ',invparam%parmax(2)
          write(0,*) 'INFO:   sigma(1,2) = ',invparam%sigma
          if (invparam%vprho_param.eq.0) then
             write(0,*) 'INFO: Inversion with Vp/Rho parameterization'
          else if (invparam%vprho_param.eq.1) then
             write(0,*) 'INFO: Inversion with Vp/Ip parameterization'
          end if
       else
          write(0,*) 'INFO:   Density is not changed'
       end if
       write(0,*) 'INFO:'
    end if
             
    write(0,*) 'INFO:   illu_pow   = ',invparam%illupow
    write(0,*) 'INFO:'
    if (invparam%const_type.eq.1) then
       write(0,*) 'INFO:   bound = 1: Model constrained at each fct/gdt eval'
    else if (invparam%const_type.eq.2) then
       write(0,*) 'INFO:   bound = 2: Model constrained at each iteration'
    end if 
    if (invparam%freeze_soft) then
       write(0,*) 'INFO:   freeze_soft = .true. : Velocity not strictly preserved in mask area'
       write(0,*) 'INFO:                          x=max(xmin,x)'
       write(0,*) 'INFO:                          x=min(xmax,x)'
       write(0,*) 'INFO:'
    else 
       write(0,*) 'INFO:   freeze_soft = .false.: Velocity strictly preservd in mask area'
       write(0,*) 'INFO:                          x=max(xmin,x), x=min(xmax,x)'
       write(0,*) 'INFO:                          x=xinit where mask=0 and mask=2'
       write(0,*) 'INFO:'
    end if
    write(0,*) 'INFO:   wantreg        =',invparam%wantreg
    write(0,*) 'INFO:   wantlog        =',invparam%wantlog
    write(0,*) 'INFO:'
    write(0,*) 'INFO: ----------------------------'

  end subroutine read_inv_params

  subroutine Init_Inversion_Array(mod,invparam)
    type(ModelSpace) ::           mod
    type(InversionParam) ::           invparam
    integer::i,j,k,n1,n2,n3

    allocate(invparam%modinit(mod%nz*mod%nx*mod%ny,invparam%nparam))
    allocate(invparam%modmask(mod%nz*mod%nx*mod%ny,invparam%nparam))

    do k=1,mod%ny
       do j=1,mod%nx
          do i=1,mod%nz
             invparam%modinit(i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx,1)=mod%vel(i,j,k)
          end do
       end do
    end do

    if (exist_file('vpmask')) then
       call from_aux('vpmask','n1',n1)
       call from_aux('vpmask','n2',n2)
       if (n1.ne.mod%nz) call erexit('Error: n1 and nz mask/vel different, exit now')
       if (n2.ne.mod%nx) call erexit('Error: n2 and nx mask/vel different, exit now')
       call sreed('vpmask',invparam%modmask(:,1),4*mod%nz*mod%nx*mod%ny)
       call auxclose('vpmask')
    else
       invparam%modmask(:,1)=1.
    end if

    if (invparam%nparam.eq.2) then  
       if (invparam%vprho_param.eq.0) then
          ! Density parameterization
          do k=1,mod%ny
             do j=1,mod%nx
                do i=1,mod%nz
                   invparam%modinit(i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx,2)=mod%rho(i,j,k)
                end do
             end do
          end do
       else
          allocate(mod%imp(mod%nz,mod%nx,mod%ny))
          ! Impedance Ip=Vp*Rho parameterization
          do k=1,mod%ny
             do j=1,mod%nx
                do i=1,mod%nz
                   invparam%modinit(i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx,2)=mod%vel(i,j,k)*mod%rho(i,j,k)
                end do
             end do
          end do
          mod%imp=mod%vel*mod%rho
       end if

       if (exist_file('par2mask')) then
          call from_aux('par2mask','n1',n1)
          call from_aux('par2mask','n2',n2)
          if (n1.ne.mod%nz) call erexit('Error: n1 and nz par2mask/vel different, exit now')
          if (n2.ne.mod%nx) call erexit('Error: n2 and nx par2mask/vel different, exit now')
          call sreed('par2mask',invparam%modmask(:,2),4*mod%nz*mod%nx*mod%ny)
          call auxclose('par2mask')
       else
          invparam%modmask(:,2)=1.
       end if      
    end if

  end subroutine Init_Inversion_Array

end module Inversion_types
