program Compute_Norm

  use Norm_mod
  use sep

  implicit none

  real, dimension(:), allocatable :: res
  real                            :: fct
  integer :: task ! compute of or gradient
  integer :: norm ! norm type
  integer :: add  ! whether or not we cummulate fct

  integer, dimension(:), allocatable:: n
  real,    dimension(:), allocatable:: o
  real,    dimension(:), allocatable:: d
  character(len=1024)               :: label

  integer :: ndim,i,j,nd,stat
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
  call from_param('norm',norm,2)
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
  else if (norm.eq.4) then
     write(0,*) 'INFO: using Hyperbolic functional'
     call from_param('thresh',thresh)
     write(0,*) 'INFO: with hyperparameter',thresh
  end if

  call from_param('add',add,0)
  fct=0.

  if (add.eq.1) then
     write(0,*) 'INFO: Adding to existing of file'
     call sreed('fct',fct,4)
     call auxclose('fct')
     write(0,*) 'INFO: fct in=',fct
  end if

  if (task.eq.0) then
     fct=fct+fct_compute(norm,res,nd,thresh)
     call srite('fct',fct,4)
     call to_history('n1',1,'fct')
     call to_history('n2',1,'fct')
     call to_history('n3',1,'fct')
     write(0,*) 'INFO: fct out=',fct
  else if (task.eq.1) then
     stat=gdt_compute(norm,res,nd,thresh)
     call srite('out',res,4*nd)
  else if (task.eq.2) then
     fct=fct+fct_compute(norm,res,nd,thresh)
     stat=gdt_compute(norm,res,nd,thresh)
     write(0,*) 'INFO: fct out=',fct
     call srite('fct',fct,4)
     call to_history('n1',1,'fct')
     call to_history('n2',1,'fct')
     call to_history('n3',1,'fct')
     call srite('out',res,4*nd)
  end if
  
  deallocate(res)
  write(0,*) 'INFO: Done with fct/gdt computation'
  
end program Compute_Norm
