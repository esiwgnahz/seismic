module symmetry_mod
  
  use adj_mod
  
  implicit none
  
  integer, private :: lagn,nt2,nt
  real,allocatable,dimension(:) :: w
  
contains
  
  subroutine symmetry_init(lagn_init,nt2_init,nt_init)
    integer ::             lagn_init,nt2_init,nt_init
    integer :: it
    lagn=lagn_init
    nt2=nt2_init
    nt=nt_init

    allocate(w(nt))

     w=0.
     do it=2,lagn
        w(it)=cos(.5*3.14159*(it-1)/(lagn-1))**2
     end do

  end subroutine symmetry_init
  
  subroutine symmetry_close()
    deallocate(w)
  end subroutine symmetry_close

  function symmetry_lop(adj,add,m,d) result(stat)
    integer :: stat
    real, dimension(:)  :: m,d
    logical, intent(in) :: adj,add
    integer             :: it
    
    call adjnull(adj,add,m,d)
    
    if (adj) then
       do it=2,lagn
          m(it)        =m(it)        + w(it)*(d(it)-d(2*nt2-it+2))
          m(2*nt2-it+2)=m(2*nt2-it+2)+ w(it)*(d(2*nt2-it+2)-d(it))
       end do
       m(1)=0.     ! Zero freq
       m(nt2+1)=0. ! Nyquist
       
    else
       do it=2,lagn
          d(it)        =d(it)        + w(it)*(m(it)-m(2*nt2-it+2))
          d(2*nt2-it+2)=d(2*nt2-it+2)+ w(it)*(m(2*nt2-it+2)-m(it))
       end do
       d(1)=0.     ! Zero freq
       d(nt2+1)=0. ! Nyquist
       
    end if
    
    stat=0 
    
  end function symmetry_lop
  
end module symmetry_mod
