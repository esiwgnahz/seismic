module OF_Res_AdjSrc_mod
  
  use sep
  use DataSpace_types
  
  implicit none

  type OFParam
     integer :: dnorm ! data space norm
  end type OFParam

  contains

    subroutine Compute_OF_RES_ADJ(ofpar,shotgath,dmodgath,resigath,f)
      type(OFParam) ::            ofpar
      type(GatherSpace)              :: shotgath
      type(TraceSpace), dimension(:) ::          dmodgath,resigath
      double precision ::                                          f

      integer :: i,begi

      begi=shotgath%begi
      if (ofpar%dnorm.eq.2) then
         
      else if (ofpar%dnorm.eq.1) then
      end if

    end subroutine Compute_OF_RES_ADJ
end module OF_Res_AdjSrc_mod
