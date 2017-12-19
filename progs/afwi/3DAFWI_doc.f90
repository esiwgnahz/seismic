module AFWI_Doc

  use sep
  use Common_Doc

  implicit none

contains

  subroutine THREEDAFWI_doc()
    call sep_add_doc_line("NAME")
    call sep_add_doc_line("")
    call sep_add_doc_line("    3DAFWI.x - 3D/2D computation of gradients for acoustic FWI, one shot only - 2nd order Time/8th order Space domain FD scheme ")
    call sep_add_doc_line("")
    call sep_add_doc_line("SYNOPSIS")
    call sep_add_doc_line("")
    call sep_add_doc_line("    3DAFWI.x vel=vel.H rho=rho.H traces=traces.H coordfile=coord.H stdin=wavelet.H data=data.H gradient=grad.H function=f.H > /dev/null ")
    call sep_add_doc_line("")
    call sep_add_doc_line("INPUT PARAMETERS")
    call sep_add_doc_line("")
    call sep_add_doc_line("    data_norm_type - char")
    call sep_add_doc_line("                     L2norm= L2 norm, sum of squares ")
    call sep_add_doc_line("                     L1norm= L1 norm, sum of abs(r) ")
    call sep_add_doc_line("                     Hyper = hyperbolic function (see Claerbout chapter 6) ")
    call sep_add_doc_line("                     Cauch = Cauchy measure ")
    call sep_add_doc_line("                     Huber = Huber functional ")
    call sep_add_doc_line("")
    call sep_add_doc_line("    data_threshold - float")
    call sep_add_doc_line("                     To be used if data_norm_type is Hyper, Cauchy or Huber ")
    call sep_add_doc_line("")
    call sep_add_doc_line("    vprho_param - integer")
    call sep_add_doc_line("                  0 - gradient for Vp/Rho paramaterization")
    call sep_add_doc_line("                  1 - gradient for Vp/Imp paramaterization")
    call sep_add_doc_line("                  Works if withRho=1 only")
    call sep_add_doc_line("")
    call Com_doc()
    call sep_add_doc_line("OPTIONAL INPUT PARAMETERS")
    call sep_add_doc_line("")
    call sep_add_doc_line("")
    call sep_add_doc_line("OPTIONAL INPUT FILES")
    call sep_add_doc_line("")
    call sep_add_doc_line("   rho - sepfile")
    call sep_add_doc_line("         density.H (nz,nx,ny) = density file, same size as velocity model")
    call sep_add_doc_line("         Only if withRho=1")
    call sep_add_doc_line("")
    call sep_add_doc_line("   wave_fwd - sepfile")
    call sep_add_doc_line("         file for the downngoing wavefield if memory needed > maxsize")
    call sep_add_doc_line("")
    call sep_add_doc_line("OUTPUT FILES")
    call sep_add_doc_line("")
    call sep_add_doc_line("   gradient - sepfile")
    call sep_add_doc_line("         grad.H (nz,nx,ny,np) = gradient file, same size as velocity model")
    call sep_add_doc_line("                 np=2 if withRho=0")
    call sep_add_doc_line("                  n4=1 f4=0 Vp gradient")
    call sep_add_doc_line("                  n4=1 f4=1 Illum file ")    
    call sep_add_doc_line("                 np=3 if withRho=1")
    call sep_add_doc_line("                  n4=1 f4=0 Vp gradient")
    call sep_add_doc_line("                  n4=1 f4=1 Rho/Imp gradient, depending on vprho_param")
    call sep_add_doc_line("                  n4=1 f4=2 Illum file ") 
    call sep_add_doc_line("                  if withRho=1")
    call sep_add_doc_line("")
    call sep_add_doc_line("   function - sepfile")
    call sep_add_doc_line("         fct.H (1) = value of objective function")
    call sep_add_doc_line("")
    call sep_add_doc_line("DESCRIPTION")
    call sep_add_doc_line("")
    call sep_add_doc_line("      Computation of gradient and objective function for one shot only for AFWI")
    call sep_add_doc_line("")
    call sep_add_doc_line("COMMENTS")
    call sep_add_doc_line("")
    call sep_add_doc_line("     (1) This code works for one shot only. For many shots, loop in a script.")
    call sep_add_doc_line("     (2) The code first checks for stability and dispersion according to.")
    call sep_add_doc_line("     the velocity model, dx, dy, dz and dt.")
    call sep_add_doc_line("     If the paramaterization is deemed unsuitable, an error will occur")
    call sep_add_doc_line("     and you have to change dx, dy, dz, dt accordingly.")
    call sep_add_doc_line("     (3) This code works for any trace geometry, no need to regularize beforehand:")
    call sep_add_doc_line("     traces are injected with a sinc interpolation at their true spatial location.")
    call sep_add_doc_line("")
    call sep_add_doc_line("")
    
    
    
    

    call doc('SOURCE')
  end subroutine THREEDAFWI_doc

end module AFWI_Doc
