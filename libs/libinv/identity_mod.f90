! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module identity_mod
  use adj_mod
  implicit none 

  real, private :: threshd,threshm
  real, dimension(:), private, pointer :: weightd
  real, dimension(:), private, pointer :: weightm
  
contains
  
  subroutine weightm_init(weight_init)
    real, dimension(:), target :: weight_init
    weightm=>weight_init
  end subroutine weightm_init

  subroutine weightm_close()
    if (associated(weightm)) nullify(weightm)
  end subroutine weightm_close

  integer function weightm_lop( adj, add, model, data)
    logical,intent(in) ::        adj, add
    real, dimension(:) ::                  model, data
    integer ::i
    call adjnull (adj,add,model,data)
    if (adj) then   
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)
          model(i) = model(i) + data(i)*weightm(i)  
       end do
       !$OMP END PARALLEL DO
    else       
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)
          data(i)  = data(i) + weightm(i)*model(i) 
       end do
       !$OMP END PARALLEL DO
    end if
    weightm_lop = 0
  end function weightm_lop

  subroutine weightd_init(weight_init)
    real, dimension(:), target :: weight_init
    weightd=>weight_init
  end subroutine weightd_init

  subroutine weightd_close()
    if (associated(weightd)) nullify(weightd)
  end subroutine weightd_close

  integer function weightd_lop( adj, add, model, data)
    logical,intent(in) ::        adj, add
    real, dimension(:) ::                  model, data
    integer ::i
    call adjnull (adj,add,model,data)
    if (adj) then   
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)
          model(i) = model(i) + data(i)*weightd(i)  
       end do
       !$OMP END PARALLEL DO
    else       
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)
          data(i)  = data(i) + weightd(i)*model(i) 
       end do
       !$OMP END PARALLEL DO
    end if
    weightd_lop = 0
  end function weightd_lop

  subroutine identity_init(threshd_init,threshm_init)
    real :: threshd_init
    real :: threshm_init
    threshd=threshd_init
    threshm=threshm_init
  end subroutine identity_init

  integer function identityd_lop( adj, add, model, data)
    logical,intent(in) ::        adj, add
    real, dimension(:) ::                  model, data
    call adjnull (adj,add,model,data)
    if (adj) then       
       model = model + data/threshd    
    else       
       data  = data + model/threshd      
    end if
    identityd_lop = 0
  end function identityd_lop

  integer function identitym_lop( adj, add, model, data)
    logical,intent(in) ::         adj, add
    real, dimension(:) ::                   model, data
    call adjnull (adj,add,model,data)
    if (adj) then       
       model = model + data/threshm  
    else       
       data  = data + model/threshm     
    end if
    identitym_lop = 0
  end function identitym_lop

  integer function identity_lop( adj, add, model, data)
    logical,intent(in) ::         adj, add
    real, dimension(:) ::                   model, data
    call adjnull (adj,add,model,data)
    if (adj) then       
       model = model + data
    else       
       data  = data + model     
    end if
    identity_lop = 0
  end function identity_lop

end module identity_mod

