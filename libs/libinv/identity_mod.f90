module identity_mod
  use adj_mod
  implicit none 

  real, private :: threshd,threshm
  
contains
  
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

