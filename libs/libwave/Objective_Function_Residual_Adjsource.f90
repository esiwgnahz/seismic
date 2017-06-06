! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module OF_Res_AdjSrc_mod
  
  use sep
  use Norm_mod
  use Inversion_types
  use DataSpace_types
  
  implicit none

  contains

    subroutine Compute_OF_RES_ADJ(begi,invparam,dobsgath,dmodgath,mutegath,resigath,f)
      integer              ::     begi
      type(InversionParam) ::          invparam
      type(TraceSpace), dimension(:) ::         dobsgath,dmodgath,mutegath,resigath
      double precision ::                                                           f

      integer :: i,j,nt,nd,stat,nrm_type
      real    :: thresh

      thresh=invparam%dat_thresh
      nrm_type=invparam%dat_nrm_type

      nt=size(dmodgath(1)%trace(:,1))
      nd=nt*size(dmodgath)

      do j=1,size(dmodgath)
         dmodgath(j)%trace=mutegath(j)%trace*(dobsgath(j)%trace-dmodgath(j)%trace)
         f=f+fct_compute(nrm_type,dmodgath(j)%trace,nt,thresh)
         resigath(begi+j-1)%trace=dmodgath(j)%trace
         dmodgath(j)%trace=mutegath(j)%trace*dmodgath(j)%trace
         stat=gdt_compute(nrm_type,dmodgath(j)%trace,nt,thresh)
      end do
      
    end subroutine Compute_OF_RES_ADJ
end module OF_Res_AdjSrc_mod
