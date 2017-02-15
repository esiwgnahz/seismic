program Apply_Filter

  use Filter_types
  use GenParam_types_flt
  use DataSpace_types_flt
  use ReadData_mod
  use ReadParams_mod
  use BuildAdaptiveFilter_mod
  use sep
  use ncnhelicon

  implicit none

  type(cube)   :: input 
  type(cube)   :: output

  type(NSfilter)    :: nmatch
  character(len=4)  :: filttag
  character(len=7)  :: filtpchtag

  integer :: ndim,stat
  logical :: adj,add

  call sep_init(SOURCE)
  filttag='filt'
  filtpchtag='filtpch'
  
  ndim=sep_dimension()

  call from_param('adj',adj,.false.)
  call from_param('add',add,.false.)

  write(0,*) 'INFO: Reading input'
  call ReadData_dim('in',input,ndim)
  call ReadData_cube('in',input)
  write(0,*) 'INFO: Allocating output'
  allocate(output%dat(product(input%n)))
  allocate(output%n(3),output%o(3),output%d(3))
  output%n=input%n
  output%o=input%o
  output%d=input%d

  write(0,*) 'INFO: '
  write(0,*) 'INFO: ADJ/ADD parameters'
  write(0,*) 'INFO: ------------------'
  write(0,*) 'INFO: Adj=',adj
  write(0,*) 'INFO: Add=',add
  write(0,*) 'INFO: ------------------'
  write(0,*) 'INFO: '
  write(0,*) 'INFO: '
  if (.not.add) output%dat=0.
  write(0,*) 'INFO: Reading filter'
  call NSfilter_read_param_from_file(filttag,nmatch,ndim)
  call psize_init(input,ndim,nmatch)
  call pch_init(input,nmatch)
  call create_nsmatch_filter(input,ndim,nmatch)
  call NSfilter_read_from_file(filttag,filtpchtag,nmatch,input,ndim)
  
  write(0,*) 'INFO: Applying filter'
  call ncnhelicon_init(nmatch%nmatch)

  if (adj) then
     stat=ncnhelicon_lop(adj,add,output%dat,input%dat)
  else
     stat=ncnhelicon_lop(adj,add,input%dat,output%dat)
  end if
 
  write(0,*) 'INFO: Writing output'
  call WriteData_dim('out',output,ndim)
  call WriteData_cube('out',output)

  write(0,*) 'INFO: deallocating'
  call cube_deallocate(output)
  call cube_deallocate(input)
  call NSfilter_deallocate(nmatch)

end program Apply_Filter
