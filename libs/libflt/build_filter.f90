! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module build_filters_mod
  
  use createnhelixmod_smp
  use lapfac3_mod
  use lapfac_mod

  use omp_lib
  
  implicit none
  
  integer,private   :: n(3)            !data size
  integer,private   :: npatch(3)       !number of patches
  integer,private   :: center(3),na(3) !box size and zero lag
  real,   private   :: o(3),d(3)       !sampling of data

contains

  subroutine init_build_filters(n_in,center_in,na_in, &
  & o_in,d_in,npatch_in)
    
    integer :: n_in(3),npatch_in(3),center_in(3),na_in(3)
    real    :: o_in(3),d_in(3)
 
    n=n_in;npatch=npatch_in;center=center_in
    na=na_in;o=o_in;d=d_in
    
  end subroutine init_build_filters
  
  subroutine build_filters(filt,rough,ncoef,mis)
    
    type(nfilter)                      :: filt,rough  ! filter and rough
    integer                            :: ncoef
    logical, dimension(:), optional    :: mis

    if(present(mis)) then
       call build_cartesian_space(filt,ncoef,mis)
    else
       call build_cartesian_space(filt,ncoef)
    end if
 
    if (na(3).ne.1) call create_lap_3d_rough(rough,ncoef)
    if (na(3).eq.1) call create_lap_2d_rough(rough,ncoef)
    
  end subroutine build_filters
 
  subroutine build_cartesian_space(filt,ncoef,mis)
    
    type(nfilter)                   :: filt
    integer                         :: ncoef
    integer                         :: i1,i2,i3,i,iref,iout(3),j(3),gap(3)
    integer, allocatable            :: map(:)
    logical, dimension(:), optional :: mis
    
    allocate(map(product(n)))
    
    j(1) = ceiling(real(n(1))/real(npatch(1))) ! number of coefficients
    j(2) = ceiling(real(n(2))/real(npatch(2))) ! per filter in the
    j(3) = ceiling(real(n(3))/real(npatch(3))) ! 3 directions
    
    do i3=1,n(3)
       iout(3)=ceiling(real(i3)/real(j(3)))
       do i2=1,n(2)
          iout(2)=ceiling(real(i2)/real(j(2)))
          do i1=1,n(1)
             iout(1)=ceiling(real(i1)/real(j(1)))
             i=i1+n(1)*(i2-1)+n(1)*n(2)*(i3-1)
             iref=iout(1)+npatch(1)*(iout(2)-1)+ &
                  &npatch(2)*npatch(1)*(iout(3)-1)
             map(i)=iref
          end do
       end do
    end do

    gap=0

    if (present(mis)) then
       filt = createnhelix_smp(n,center,gap,na,map,ncoef,mis)
    else
       filt = createnhelix_smp(n,center,gap,na,map,ncoef)
    end if    

    deallocate(map)
    
  end subroutine build_cartesian_space

  subroutine create_lap_3d_rough(rough,ncoef)
    
    type(nfilter)        :: rough
    integer              :: ncoef
    integer              :: npts(3)
    integer, allocatable :: nh(:),map(:)
    
    npts=(/npatch(1),npatch(2),npatch(3)/)
    allocate(nh(product(npts)))
    allocate(map(ncoef*product(npts)))
    
    nh=ncoef
    map=1
    
    nh=min(6,npts(1))*3
    call nallocate(rough,nh,map)
    rough%hlx(1)=lapfac3d_mod(.01,npts(1),npts(2),nh(1)/3)
        
    deallocate(nh,map)
    
  end subroutine create_lap_3d_rough
 
  subroutine create_lap_2d_rough(rough,ncoef)
    
    type(nfilter)        :: rough
    integer              :: ncoef
    integer              :: npts(2)
    integer, allocatable :: nh(:),map(:)
    
    npts=(/npatch(1),npatch(2)/)
    allocate(nh(product(npts)))
    allocate(map(ncoef*product(npts)))
    
    nh=ncoef
    map=1
    
    nh=min(6,npts(1))*2
    call nallocate(rough,nh,map)
    rough%hlx(1)=lapfac2_mod(.01,npts(1),nh(1)/2)
    
    deallocate(nh,map)
    
  end subroutine create_lap_2d_rough
  
end module build_filters_mod














