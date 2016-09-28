! Copyright 1999 Stanford University
!All rights reserved
!
!Author: Sergey Fomel and Robert Clapp
module npef_mod
 
 ! Estimate non-stationary PEF
!  use nhconest_mod
!  use nhelicon_mod
!  use solver2_mod
!  use cgstep_omp_mod
!  use nhelix
!  use helix

  implicit none

contains

!  subroutine find_pef_mod( dd, aa, rr, niter, eps, nh, wt, verb)
!    optional                         :: wt
!    integer                          :: niter, nh      ! number of iterations
!    real                             :: eps	       ! epsilon
!    type( nfilter)                   :: aa             ! estimated PEF output.
!    type( nfilter)                   :: rr             ! roughening filter.
!    real, dimension(:)               :: dd             ! input data
!    integer                          :: ip, ih, np, nr ! filter lengths
!    real, dimension (:), allocatable :: flt	       ! np*na filter coefs
!    real, dimension (:)              :: wt
!    
!    real, dimension(:), allocatable  ::tmp
!
!
!    integer :: stat
!
!    logical :: verb
!
!    np = size( aa%hlx)			! data   length
!    nr = np*nh
!    allocate( flt( nr))
!    flt=0.
!    
!    call nhconest_mod_init(dd,aa,nh)
!    call nhelicon_mod_init(rr)
!    
!    allocate(tmp(size(dd)))
!    tmp=-dd
!
!    if (present(wt)) then
!       call solver_reg(nhconest_mod_lop,cgstep_omp,nhelicon_mod_lop,nreg=nr, &
!       & x=flt,dat=tmp,niter=niter,eps=eps,verb=verb,wt=wt)
!    else
!       call solver_reg(nhconest_mod_lop,cgstep_omp,nhelicon_mod_lop,nreg=nr, &
!       & x=flt,dat=tmp,niter=niter,eps=eps,verb=verb)
!    end if
!    
!    call cgstep_omp_close()
!    call nhelicon_mod_close() 
!    call nhconest_mod_close()
!
!    do ip = 1, np 
!       do ih = 1, nh
!          aa%hlx( ip)%flt( ih) = flt( (ip-1)*nh + ih)
!       end do
!    end do
!
!    deallocate(flt,tmp)
!    return
!
!  end subroutine find_pef_mod

end module npef_mod














