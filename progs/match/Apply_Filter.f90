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

  call sep_init(SOURCE)
  filttag='filt'
  filtpchtag='filtpch'
  
  ndim=sep_dimension()

  write(0,*) 'INFO: Reading input'
  call ReadData_dim('in',input,ndim)
  call ReadData_cube('in',input)
  write(0,*) 'INFO: Allocating output'
  allocate(output%dat(product(input%n)))
  output%dat=0.
  allocate(output%n(3),output%o(3),output%d(3))
  output%n=input%n
  output%o=input%o
  output%d=input%d

  write(0,*) 'INFO: Reading filter'
  call NSfilter_read_param_from_file(filttag,nmatch,ndim)
  call psize_init(input,ndim,nmatch)
  call pch_init(input,nmatch)
  call create_nsmatch_filter(input,ndim,nmatch)
  call NSfilter_read_from_file(filttag,filtpchtag,nmatch,input,ndim)
  
  write(0,*) 'INFO: Applying filter'
  call ncnhelicon_init(nmatch%nmatch)
  stat=ncnhelicon_lop(.false.,.false.,input%dat,output%dat)
 
  write(0,*) 'INFO: Writing output'
  call WriteData_dim('out',output,ndim)
  call WriteData_cube('out',output)

  write(0,*) 'INFO: deallocating'
  call cube_deallocate(output)
  call cube_deallocate(input)
  call NSfilter_deallocate(nmatch)

end program Apply_Filter
