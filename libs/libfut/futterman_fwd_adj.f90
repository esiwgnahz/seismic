module futterman_fwd_adj

  use omp_lib
  use adj_mod
  implicit none

  real, private, dimension(:,:), pointer :: futs,futsadj
  integer, private :: n2prod,ntmax

contains

  subroutine futterman_fwd_adj_init(futs_init,futsadj_init,n2prod_init,ntmax_init)
    real, dimension(:,:), pointer :: futs_init,futsadj_init
    integer :: n2prod_init,ntmax_init
    futs=>futs_init
    futsadj=>futsadj_init
    ntmax=ntmax_init
    n2prod=n2prod_init
  end subroutine futterman_fwd_adj_init

  function futterman_fwd_adj_lop(adj,add,mm,dd) result(stat)
    integer :: stat    
    logical,intent(in) :: adj,add
    real,dimension(:) :: mm,dd
    call adjnull (adj,add,mm,dd )
    call futterman_fwd_adj_lop2(adj,add,mm,dd )
    stat=0
  end function futterman_fwd_adj_lop

  subroutine futterman_fwd_adj_lop2(adj,add,mm,dd)
    real, dimension(ntmax,n2prod) :: mm,dd
    logical,intent(in) :: adj,add
    integer :: i
   
    if (adj) then
       !$OMP PARALLEL DO
       do i=1,n2prod
!          mm(:,i)=mm(:,i)+matmul(transpose(futs),dd(:,i))
          mm(:,i)=mm(:,i)+matmul(futsadj,dd(:,i))
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO
       do i=1,n2prod
          dd(:,i)=dd(:,i)+matmul(futs,mm(:,i))
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine futterman_fwd_adj_lop2

end module futterman_fwd_adj
