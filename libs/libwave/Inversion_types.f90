module Inversion_types

  use sep
  use ModelSpace_types

  implicit none

  type InversionParam
     integer:: iter       ! iteration number
     integer:: niter      ! Maximum number of iterations
     integer:: neval      ! Maximum number of fct/gdt eval.
     integer:: const_type ! constrained velocity at (1) each Fct/Gdt eval
                          !                         (2) each iteration
     logical:: freeze_soft! (0) velocity is preserved in mask area
                          ! (1) velocity is not strictly enforced in masking area
     real   :: vpmin
     real   :: vpmax

     real   :: eps        ! regularization parameter
     
     double precision, allocatable:: vpinit(:)
     real,             allocatable:: vpmask(:)
     real,             allocatable:: datmask(:,:)

  end type InversionParam

contains

  subroutine Init_Inversion_Array(mod,invparam)
    type(ModelSpace) ::           mod
    type(InversionParam) ::           invparam
    integer::i,j,k,n1,n2,n3

    allocate(invparam%vpinit(mod%nz*mod%nx*mod%ny))
    allocate(invparam%vpmask(mod%nz*mod%nx*mod%ny))

    do k=1,mod%ny
       do j=1,mod%nx
          do i=1,mod%nz
             invparam%vpinit(i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)=mod%vel(i,j,k)
          end do
       end do
    end do

    if (exist_file('vpmask')) then
       call from_aux('vpmask','n1',n1)
       call from_aux('vpmask','n2',n2)
       if (n1.ne.mod%nz) call erexit('Error: n1 and nz mask/vel different, exit now')
       if (n2.ne.mod%nx) call erexit('Error: n2 and nx mask/vel different, exit now')
       call sreed('vpmask',invparam%vpmask,4*mod%nz*mod%nx*mod%ny)
    else
       invparam%vpmask=1.
    end if

  end subroutine Init_Inversion_Array

end module Inversion_types
