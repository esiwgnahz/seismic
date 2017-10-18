program Compute_Norm

  use Norm_mod
  use sep

  implicit none

  real, dimension(:), allocatable :: res
  real, dimension(:), allocatable :: gdt
  real                            :: fct
  integer :: task ! compute of or gradient
  integer :: norm ! norm type

  integer, dimension(:), allocatable:: n
  real,    dimension(:), allocatable:: o
  real,    dimension(:), allocatable:: d
  character(len=1024)               :: label

  integer :: ndim,i,j,nd
  real    :: thresh

  call sep_init(SOURCE)
  ndim=sep_dimension()
  allocate(n(ndim),o(ndim),d(ndim))
  
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do

  nd=product(n)
  allocate(res(nd))
  call sreed('in',res,4*nd)

  call from_param('task',task,0)
  if (task.eq.0) then
     write(0,*) 'INFO: Compute objective function'
  else if (task.eq.1) then
     write(0,*) 'INFO: Compute gradient'
  else if (task.eq.2) then
     write(0,*) 'INFO: Compute gradient and function'
  end if

  thresh=0.
  call from_parm('norm',norm,2)
  if (norm.eq.2) then
     write(0,*) 'INFO: using L2 norm'
  else if (norm.eq.1) then
     write(0,*) 'INFO: using L1 norm'
  else if (norm.eq.12) then
     write(0,*) 'INFO: using Huber norm'
     call from_param('thresh',thresh)
     write(0,*) 'INFO: with threshold',thresh
  else if (norm.eq.3) then
     write(0,*) 'INFO: using Cauchy norm'
     call from_param('thresh',thresh)
     write(0,*) 'INFO: with hyperparameter',thresh
  end if

  if (task.eq.0) then
     fct=fct_compute(norm,res,nd,thresh)
     call srite('fct',fct,4)
     call to_history('n1',1,'fct')
     call to_history('n2',1,'fct')
     call to_history('n3',1,'fct')
  else if (task.eq.1) then
     allocate(gdt(nd)); gdt=0.
     gdt=gdt_compute(norm,res,nd,thresh)
     call srite('out',gdt,4*nd)
  else if (task.eq.2) then
     allocate(gdt(nd)); gdt=0.
     gdt=gdt_compute(norm,res,nd,thresh)
     fct=fct_compute(norm,res,nd,thresh)
     call srite('fct',fct,4)
     call to_history('n1',1,'fct')
     call to_history('n2',1,'fct')
     call to_history('n3',1,'fct')
     call srite('out',gdt,4*nd)
  end if
  
  if (allocated(gdt)) deallocate(gdt) 
  deallocate(res)
  write(0,*) 'INFO: Done with fct/gdt computation'
  
end program Compute_Norm
