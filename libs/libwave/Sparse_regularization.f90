module Sparse_regularization_mod

  use sep
  use helix
  use print
  use helicon_mod
  use helicon2_mod
  use helicon3_mod
  use createhelixmod
  use ModelSpace_types
  
  implicit none

  type SparseRegParam
     character(len=3) :: flt_type
     character(len=6) :: nrm_type 
     type(filter)     :: xx,zz,yy
     integer          :: ntaperz
     integer          :: ntaperx
     integer          :: ntapery
     real             :: eps
     integer          :: nx,nz,ny
     integer          :: nx_tap,nz_tap,ny_tap
  end type SparseRegParam

contains

  subroutine Init_SparseRegularization(sparseparam,mod)
    type(SparseRegParam) ::            sparseparam
    type(ModelSpace)     ::                        mod

    call from_param('flt_type',sparseparam%flt_type,'ZZZ')
    call from_param('nrm_type',sparseparam%nrm_type,'Cauchy')
    call from_param('eps',sparseparam%eps,0.)
    call from_param('sparse_taperx',sparseparam%ntaperx,10)
    call from_param('sparse_tapery',sparseparam%ntapery,0)
    call from_param('sparse_taperz',sparseparam%ntaperz,10)

    sparseparam%nx_tap=mod%nx+2*sparseparam%ntaperx
    sparseparam%nz_tap=mod%nz+2*sparseparam%ntaperz

    if (mod%ny.ne.1) then
       sparseparam%ny_tap=mod%ny+2*sparseparam%ntapery
    else
       sparseparam%ntapery=0
       sparseparam%ny_tap=mod%ny
    end if
    sparseparam%nx=mod%nx
    sparseparam%ny=mod%ny
    sparseparam%nz=mod%nz

    write(0,*) 'INFO:--------------------------'
    write(0,*) 'INFO: Regularization paramters '
    write(0,*) 'INFO:--------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO:  flt_type = ',sparseparam%flt_type
    write(0,*) 'INFO:  nrm_type = ',sparseparam%nrm_type
    write(0,*) 'INFO:  eps      = ',sparseparam%eps
    write(0,*) 'INFO:  taperz    = ',sparseparam%ntaperz
    write(0,*) 'INFO:  taperx    = ',sparseparam%ntaperx
    write(0,*) 'INFO:  tapery    = ',sparseparam%ntapery
    write(0,*) 'INFO:'
    write(0,*) 'INFO:-------------------------'
    write(0,*) 'INFO:'
    call SparseRegularization_filter_init(sparseparam)
    
  end subroutine Init_SparseRegularization

  subroutine SparseRegularization_filter_init(sparseparam)
    type(SparseRegParam) ::                   sparseparam
    
    integer, dimension(:), allocatable :: n0,x,y,z,center,gap,npef
    integer :: nzt,nyt,nxt

    nxt=sparseparam%nx_tap
    nyt=sparseparam%ny_tap
    nzt=sparseparam%nz_tap

    write(0,*) 'INFO: ----------------------------'
    write(0,*) 'INFO: Built regularization filters'
    if (nyt.ne.1) then
       allocate(x(3),y(3),z(3),center(3),gap(3),n0(3),npef(3))
       center=1; gap=0; npef=(/nzt,nxt,nyt/); n0=npef
       z=(/2,1,1/)
       x=(/1,2,1/)
       y=(/1,1,2/)   
       sparseparam%yy=createhelix(npef,center,gap,y)
       sparseparam%yy%flt=-1    
       write(0,*) 'INFO: Filter y:'
       call printn(n0,center,y,sparseparam%yy)    
    else
       allocate(x(2),y(2),z(2),center(2),gap(2),n0(2),npef(2))
       center=1; gap=0; npef=(/nzt,nxt/); n0=npef
       z=(/2,1/)
       x=(/1,2/)   
    end if
    sparseparam%zz=createhelix(npef,center,gap,z)
    sparseparam%xx=createhelix(npef,center,gap,x)            
    sparseparam%zz%flt=-1
    sparseparam%xx%flt=-1

    write(0,*) 'INFO: Filter z:'
    call printn(n0,center,z,sparseparam%zz) 
    write(0,*) 'INFO: Filter x:' 
    call printn(n0,center,x,sparseparam%xx)   
    write(0,*) 'INFO: ----------------------------'

    deallocate(x,y,z,center,gap,n0,npef)

  end subroutine SparseRegularization_filter_init
  
  subroutine SparseRegularization_apply(sparseparam,mod,g,f)
    type(SparseRegParam) ::             sparseparam
    type(ModelSpace)     ::                         mod
    real, dimension(:)   ::                             g
    double precision     ::                               f
    
    real, dimension(:), allocatable :: rm   
    real, dimension(:), allocatable :: rmtap   
    real, dimension(:), allocatable :: rtmptap

    integer :: nxt,nyt,nzt
    integer :: nz,ny,nx

    nz=sparseparam%nz
    ny=sparseparam%ny
    nx=sparseparam%nx
    nxt=sparseparam%nx_tap
    nyt=sparseparam%ny_tap
    nzt=sparseparam%nz_tap

    allocate(rmtap(nxt*nyt*nzt),rtmptap(nxt*nyt*nzt),rm(nx*ny*nz))
    call helicon_mod_init(sparseparam%zz)
    call helicon2_mod_init(sparseparam%xx)
    if (ny.ne.1) call helicon3_mod_init(sparseparam%yy)
    

    deallocate(rm,rmtap,rtmptap)
  end subroutine SparseRegularization_apply

  subroutine SparseRegularization_add_taper(sparseparam,arraytap,array,forw)
    type(SparseRegParam) ::                 sparseparam
    real    :: arraytap(-sparseparam%ntaperz+1:sparseparam%nz+sparseparam%ntaperz,-sparseparam%ntaperx+1:sparseparam%nx+sparseparam%ntaperx,-sparseparam%ntapery+1:sparseparam%ny+sparseparam%ntapery)
    real    :: array(sparseparam%nz,sparseparam%nx,sparseparam%ny)
    integer :: i,j,nx,ny,nz,ntaperx,ntapery,ntaperz
    logical :: forw
    
    nz=sparseparam%nz
    ny=sparseparam%ny
    nx=sparseparam%nx
    ntaperz=sparseparam%ntaperz
    ntapery=sparseparam%ntapery
    ntaperx=sparseparam%ntaperx

    if (forw) then
       arraytap=0.
       arraytap(1:nz,1:nx,1:ny)=array
       
       do j=-ntaperz+1,0
          arraytap(j,:,:)=arraytap(1,:,:)
       end do
       do j=nz+1,nz+ntaperz
          arraytap(j,:,:)=arraytap(nz,:,:)
       end do
       do j=-ntaperx+1,0
          arraytap(:,j,:)=arraytap(:,1,:)
       end do
       do j=nx+1,nx+ntaperx
          arraytap(:,j,:)=arraytap(:,nx,:)
       end do
       do j=-ntapery+1,0
          arraytap(:,:,j)=arraytap(:,:,1)
       end do
       do j=ny+1,ny+ntapery
          arraytap(:,:,j)=arraytap(:,:,ny)
       end do
    else
       array=arraytap(1:nz,1:nx,1:ny)
    end if
  end subroutine SparseRegularization_add_taper

end module Sparse_regularization_mod
