module ncnmisinput  
  ! find a mask of missing filter inputs
  use ncnhelicon
  use nhelix
  use Filter_types

  implicit none

contains

  subroutine find_ncmask( known, aa)
    real, dimension(:) :: known
    type(NSfilter)            :: aa
    real, dimension(:), allocatable :: rr
    integer                         :: stat, i
 
    allocate(rr(size(known)))
    call ncnhelicon_init(aa%nmatch)

    do i = 1, size( aa%nmatch%hlx)
       aa%nmatch%hlx( i)%flt = 1.
    end do

    stat = ncnhelicon_lop( .false., .false., known, rr)

    where ( rr > 0.)
       aa%nmatch%mis = .true.        
    end where

    do i = 1, size( aa%nmatch%hlx)
       aa%nmatch%hlx( i)%flt = 0. 
    end do

    deallocate(rr)

  end subroutine find_ncmask

end module ncnmisinput
