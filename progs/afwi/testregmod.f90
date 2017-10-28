! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module testreg_mod

  use Norm_mod
  use helicon_mod
  use helicon2_mod

  implicit none

  integer, private :: counter,ncount
  real, private :: epsz,epszp,epszm,epsx
  real, private :: threshp,threshm
  integer, private :: regtypep,regtypem
  real, dimension(:), pointer, private ::w

contains
  
  subroutine testreg_init(regtypep_init,regtypem_init,epszp_init,epszm_init,epsx_init,threshp_init,threshm_init,w_init)
    real :: epsx_init,epszp_init,epszm_init,threshp_init,threshm_init
    integer:: regtypep_init,regtypem_init
    real, dimension(:), target :: w_init
    
    epsz=epszp_init
    epszm=epszm_init
    epszp=epszp_init
    epsx=epsx_init
    threshp=threshp_init
    threshm=threshm_init
    regtypep=regtypep_init
    regtypem=regtypem_init
    w=>w_init

    counter=0
    ncount=0

  end subroutine testreg_init

  subroutine testreg_count(count)
    integer :: count

    count=ncount
  end subroutine testreg_count

  function simple_test(sizex,d,x,g,f) result (stat)
    integer :: stat
    integer :: sizex
    real, dimension(sizex) :: d
    real, dimension(sizex) :: x
    real, dimension(sizex) :: g
    double precision   :: f

    real, dimension(:), allocatable :: gtmpx
    real, dimension(:), allocatable :: gtmpz
    real, dimension(:), allocatable :: gx
    real, dimension(:), allocatable :: gz
    real, dimension(:), allocatable :: r

    allocate(gtmpz(sizex))
    allocate(gtmpx(sizex))
    allocate(gz(sizex))
    allocate(gx(sizex))
    allocate(r(sizex))

    counter=counter+1 

    f=0.
    g=0.
    gx=0.
    gz=0.
    gtmpx=0.
    gtmpz=0.
    r=d-x
    f=f+fct_compute(2,r,sizex)
    stat=gdt_compute(2,r,sizex)

    stat= helicon_mod_lop(.false.,.false.,x,gtmpz)
    stat=helicon2_mod_lop(.false.,.false.,x,gtmpx)
    gtmpz=w*gtmpz
    gtmpx=w*gtmpx

!    if (mod(counter,10).eq.0) then
!       ncount=ncount+1
!       call srite('grad',epsz*gtmpz+epsx*gtmpx,4*sizex)
!    end if

    f=f+epsz*fct_compute(regtypep,gtmpz,sizex,threshp)
    f=f+epsx*fct_compute(regtypep,gtmpx,sizex,threshp)
    
    stat=gdt_compute(regtypep,gtmpz,sizex,threshp)
    stat=gdt_compute(regtypep,gtmpx,sizex,threshp)

    gtmpz=w*gtmpz
    gtmpx=w*gtmpx
    
    stat= helicon_mod_lop(.true.,.false.,gz,gtmpz)
    stat=helicon2_mod_lop(.true.,.false.,gx,gtmpx)

    g=(-r+epsx*gx+epsz*gz)/sizex
    f=f/sizex
    deallocate(gtmpz,gtmpx)
    deallocate(gz,gx,r)

    stat=0

  end function simple_test

  function simple_test_hinge(sizex,d,x,g,f) result (stat)
    integer :: stat
    integer :: sizex
    real, dimension(sizex) :: d
    real, dimension(sizex) :: x
    real, dimension(sizex) :: g
    double precision   :: f

    real, dimension(:), allocatable :: gtmpx
    real, dimension(:), allocatable :: gtmpz
    real, dimension(:), allocatable :: gtmpzp
    real, dimension(:), allocatable :: gtmpzm
    real, dimension(:), allocatable :: mgtmpx
    real, dimension(:), allocatable :: mgtmpzp
    real, dimension(:), allocatable :: mgtmpzm
    real, dimension(:), allocatable :: gx
    real, dimension(:), allocatable :: gzp,gzm
    real, dimension(:), allocatable :: r

    allocate(gtmpz(sizex))
    allocate(gtmpzp(sizex))
    allocate(gtmpzm(sizex))
    allocate(gtmpx(sizex))
    allocate(mgtmpzp(sizex))
    allocate(mgtmpzm(sizex))
    allocate(mgtmpx(sizex))
    allocate(gzp(sizex))
    allocate(gzm(sizex))
    allocate(gx(sizex))
    allocate(r(sizex))

    counter=counter+1

    f=0.
    g=0.
    gx=0.
    gzp=0.
    gzm=0.
    gtmpx=0.
    gtmpz=0.
    mgtmpx=0.
    mgtmpzp=0.
    mgtmpzm=0.
    gtmpzp=0.
    gtmpzm=0.
    r=d-x
    f=f+fct_compute(2,r,sizex)
    stat=gdt_compute(2,r,sizex)

    stat= helicon_mod_lop(.false.,.false.,x,gtmpz)
    stat=helicon2_mod_lop(.false.,.false.,x,gtmpx)

    where (gtmpz.le.0) 
       mgtmpzm=1.
       mgtmpzp=0.
    elsewhere 
       mgtmpzm=0.
       mgtmpzp=1.
    end where

    gtmpzp=w*gtmpz
    gtmpzm=w*gtmpz*mgtmpzm

    gtmpx=w*gtmpx

!    if (mod(counter,10).eq.0) then
!       ncount=ncount+1
!       call srite('grad',epszp*gtmpzp+epszm*gtmpzm+epsx*gtmpx,4*sizex)
!    end if

    f=f+epszm*fct_compute(regtypem,gtmpzm,sizex,threshm)
    f=f+epszp*fct_compute(regtypep,gtmpzp,sizex,threshp)
    f=f+ epsx*fct_compute(regtypep,gtmpx, sizex,threshp)
    
    stat=gdt_compute(regtypem,gtmpzm,sizex,threshm)
    stat=gdt_compute(regtypep,gtmpzp,sizex,threshp)
    stat=gdt_compute(regtypep,gtmpx ,sizex,threshp)

    gtmpzp=w*gtmpzp
    gtmpzm=w*gtmpzm*mgtmpzm
    gtmpx=w*gtmpx

    stat= helicon_mod_lop(.true.,.false.,gzp,gtmpzp)
    stat= helicon_mod_lop(.true.,.false.,gzm,gtmpzm)
    stat=helicon2_mod_lop(.true.,.false.,gx,gtmpx)

    g=(-r+epsx*gx+epszp*gzp+epszm*gzm)/sizex
    f=f/sizex**2

    deallocate(gtmpz,gtmpx)
    deallocate(gtmpzm,gtmpzp)
    deallocate(mgtmpzp,mgtmpzm,mgtmpx)
    deallocate(gzp,gzm,gx,r)

    stat=0

  end function simple_test_hinge

end module testreg_mod
