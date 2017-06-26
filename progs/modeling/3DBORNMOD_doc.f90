module THREED_Doc

  use sep
  use Common_Doc

  implicit none

contains

  subroutine BORNMOD_doc()
    call sep_add_doc_line("NAME")
    call sep_add_doc_line("    3DBORNMOD.x - 3D/2D Borm modelin - 2nd order Time/8th order Space domain FD scheme ")
    call sep_add_doc_line("")
    call sep_add_doc_line("SYNOPSIS")
    call sep_add_doc_line("    3DBORNMOD.x vel=vel.H ref=reflectivity.H traces=traces.H coordfile=coord.H stdin=wavelet.H data=data.H > /dev/null ")
    call sep_add_doc_line("")
    call sep_add_doc_line("INPUT PARAMETERS")
    call sep_add_doc_line("    ref - sepfile")
    call sep_add_doc_line("          reflectivity.H (nz,nx,ny) = reflectivity file, same size as velocity model")
    call sep_add_doc_line("")
    call Com_doc()
    call sep_add_doc_line("")
    call sep_add_doc_line("    Born -boolean")
    call sep_add_doc_line("               [1] Should always be one")
    call sep_add_doc_line("")
    call sep_add_doc_line("OUTPUT PARAMETERS")
    call sep_add_doc_line("    data - sepfile")
    call sep_add_doc_line("            data.H (nt,ntraces) = output")
    call sep_add_doc_line("")
    call sep_add_doc_line("DESCRIPTION")
    call sep_add_doc_line("      Acoustic Born modeling in 2D and 3D (no density)")
    call sep_add_doc_line("")
    call sep_add_doc_line("COMMENTS")
    call sep_add_doc_line("     (1) This code works for one shot only. For many shots, loop in a script.")
    call sep_add_doc_line("     (2) The code first checks for stability and dispersion according to.")
    call sep_add_doc_line("     the velocity model, dx, dy, dz and dt.")
    call sep_add_doc_line("     If the paramaterization is deemed unsuitable, an error will occur")
    call sep_add_doc_line("     and you have to change dx, dy, dz, dt accordingly.")
    call sep_add_doc_line("     (3) This code works for any trace geometry, no need to regularize beforehand:")
    call sep_add_doc_line("     traces are injected with a sinc interpolation at their true spatial location.")
    call sep_add_doc_line("")
    call doc('SOURCE')

  end subroutine BORNMOD_doc

end module THREED_Doc
