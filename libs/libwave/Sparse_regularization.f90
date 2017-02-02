module Sparse_regularization_mod

  use sep
  use helix
  use print
  use helicon_mod
  use helicon2_mod
  use helicon3_mod
  use createhelixmod
  use Norm_mod
  use ModelSpace_types
  
  implicit none

  type SparseRegParam
     logical          :: derivdx
     logical          :: derivdz
     logical          :: derivdy
     character(len=6) :: mod_nrm_type_char
     integer          :: mod_nrm_type 
     type(filter)     :: xx,zz,yy
     integer          :: ntaperz
     integer          :: ntaperx
     integer          :: ntapery
     real             :: eps
     real             :: ratio
     real             :: mod_thresh
     integer          :: nx,nz,ny
     integer          :: nx_tap,nz_tap,ny_tap
  end type SparseRegParam

contains

  subroutine Init_SparseRegularization(sparseparam,mod)
    type(SparseRegParam) ::            sparseparam
    type(ModelSpace)     ::                        mod

    call from_param('derivdx',sparseparam%derivdx,.false.)
    call from_param('derivdy',sparseparam%derivdy,.false.)
    call from_param('derivdz',sparseparam%derivdz,.false.)
    call from_param('model_nrm_type',sparseparam%mod_nrm_type_char,'Cauchy')
    call from_param('ratio',sparseparam%ratio,10.)
    call from_param('reg_threshold',sparseparam%mod_thresh,0.)
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

    if(sparseparam%mod_nrm_type_char(1:5).eq.'Cauch') sparseparam%mod_nrm_type=3
    if(sparseparam%mod_nrm_type_char(1:5).eq.'Huber') sparseparam%mod_nrm_type=12
    if(sparseparam%mod_nrm_type_char(1:5).eq.'L1nor') sparseparam%mod_nrm_type=1
    if(sparseparam%mod_nrm_type_char(1:5).eq.'L2nor') sparseparam%mod_nrm_type=2

    write(0,*) 'INFO:---------------------------'
    write(0,*) 'INFO: Regularization parameters '
    write(0,*) 'INFO:---------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO:  d/dx       = ',sparseparam%derivdx
    write(0,*) 'INFO:  d/dy       = ',sparseparam%derivdy
    write(0,*) 'INFO:  d/dz       = ',sparseparam%derivdz
    write(0,*) 'INFO:  nrm_type   = ',sparseparam%mod_nrm_type_char  
    if ((sparseparam%mod_nrm_type.eq.3).or.(sparseparam%mod_nrm_type.eq.12)) then
      write(0,*) 'INFO:  mod thresh = ',sparseparam%mod_thresh
    end if
    write(0,*) 'INFO:  ratio      = ',sparseparam%ratio
    write(0,*) 'INFO:  taperz     = ',sparseparam%ntaperz
    write(0,*) 'INFO:  taperx     = ',sparseparam%ntaperx
    write(0,*) 'INFO:  tapery     = ',sparseparam%ntapery
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
       if (sparseparam%derivdy) then  
          write(0,*) 'INFO: Filter y:'
          call printn(n0,center,y,sparseparam%yy) 
       end if
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

    if (sparseparam%derivdz) then
       write(0,*) 'INFO: Filter z:'
       call printn(n0,center,z,sparseparam%zz)
    end if
    if (sparseparam%derivdx) then
       write(0,*) 'INFO: Filter x:' 
       call printn(n0,center,x,sparseparam%xx) 
    end if
    write(0,*) 'INFO: ----------------------------'

    deallocate(x,y,z,center,gap,n0,npef)

  end subroutine SparseRegularization_filter_init

  subroutine SparseRegularization_filter_close(sparseparam)
    type(SparseRegParam) ::                    sparseparam
    call deallocatehelix(sparseparam%xx)
    call deallocatehelix(sparseparam%zz)
    if(sparseparam%ny_tap.ne.1) call deallocatehelix(sparseparam%yy)
  end subroutine SparseRegularization_filter_close

  subroutine SparseRegularization_apply(sparseparam,                  vel,g,f)
    type(SparseRegParam) ::             sparseparam
    real, dimension(sparseparam%nx*sparseparam%nz*sparseparam%ny) ::  vel
    real, dimension(sparseparam%nx*sparseparam%nz*sparseparam%ny) ::      g
    double precision     ::                                                 f
    
    call helicon_mod_init(sparseparam%zz)
    call helicon2_mod_init(sparseparam%xx)
    if (sparseparam%ny.ne.1) call helicon3_mod_init(sparseparam%yy)
    
    g=0.
    f=0.

    if (sparseparam%derivdx) call SparseRegularization_Compute_f_g(helicon2_mod_lop,sparseparam,vel,f,g)
    if (sparseparam%derivdz) call SparseRegularization_Compute_f_g(helicon_mod_lop,sparseparam,vel,f,g)

    if ((sparseparam%ny.ne.1).and.sparseparam%derivdy) call SparseRegularization_Compute_f_g(helicon2_mod_lop,sparseparam,vel,f,g)
    
  end subroutine SparseRegularization_apply

  subroutine SparseRegularization_Compute_f_g(op,sparseparam,vel,f,g)
    interface
       function op( adj, add, xx, yy) result(stat)
         integer            :: stat 
         logical,intent(in) :: adj,add 
         real,dimension(:)  :: xx,yy 
       end function op
    end interface
    type(SparseRegParam) ::                     sparseparam
    real, dimension(sparseparam%nx*sparseparam%nz*sparseparam%ny) :: vel
    real, dimension(sparseparam%nx*sparseparam%nz*sparseparam%ny) :: g
    double precision     ::                                          f

    real, dimension(:), allocatable :: rm   
    real, dimension(:), allocatable :: rmtap   
    real, dimension(:), allocatable :: rtmptap

    integer :: nxt,nyt,nzt
    integer :: nz,ny,nx,stat,nrm_type

    nz=sparseparam%nz
    ny=sparseparam%ny
    nx=sparseparam%nx
    nxt=sparseparam%nx_tap
    nyt=sparseparam%ny_tap
    nzt=sparseparam%nz_tap
    nrm_type=sparseparam%mod_nrm_type

    allocate(rmtap(nxt*nyt*nzt),rtmptap(nxt*nyt*nzt),rm(nx*ny*nz))
    
    rm=0.
    
    call SparseRegularization_add_taper(sparseparam,rtmptap,vel,forw=.true.)
    stat=op(.false.,.false.,rtmptap,rmtap)   
    call SparseRegularization_add_taper(sparseparam,rmtap,rm,forw=.false.)

    ! Add model fitting side of objective function
    ! --------------------------------------------  
    f=f+fct_compute(nrm_type,rm,nz*nx,sparseparam%mod_thresh)

    ! Now apply the adjoint to form the gradient of the model fitting
    ! ---------------------------------------------------------------
    stat=gdt_compute(nrm_type,rmtap,nzt*nxt,sparseparam%mod_thresh)  
    stat=op(.true.,.false.,rtmptap,rmtap)
    call SparseRegularization_add_taper(sparseparam,rtmptap,g,forw=.false.)

    deallocate(rmtap,rtmptap,rm)

  end subroutine SparseRegularization_Compute_f_g
  
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
       array=array+arraytap(1:nz,1:nx,1:ny)
    end if
  end subroutine SparseRegularization_add_taper

end module Sparse_regularization_mod
