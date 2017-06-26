module THREED_Doc

  use sep
  use Common_Doc

  implicit none

contains

  subroutine RTM_doc()
    call sep_add_doc_line("NAME")
    call sep_add_doc_line("    3DRTM.x - 3D/2D RTM - 2nd order Time/8th order Space domain FD scheme ")
    call sep_add_doc_line("")
    call sep_add_doc_line("SYNOPSIS")
    call sep_add_doc_line("    3DRTM.x vel=vel.H traces=traces.H  coordfile=coord.H stdin=wavelet.H image=img.H > /dev/null ")
    call sep_add_doc_line("")
    call sep_add_doc_line("INPUT PARAMETERS")
    call Com_doc()
    call sep_add_doc_line("")
    call sep_add_doc_line("    LSRTM -boolean")
    call sep_add_doc_line("               [0] Whether or not we want the adjoint of Born modeling")
    call sep_add_doc_line("")
    call sep_add_doc_line("    CHEAPLSRTM -boolean")
    call sep_add_doc_line("               [0] Whether or not we want the fast version of the adjoint of Born modeling")
    call sep_add_doc_line("")
    call sep_add_doc_line("OUTPUT PARAMETERS")
    call sep_add_doc_line("    image - sepfile")
    call sep_add_doc_line("            img.H (nz,nx,ny,2) = output")
    call sep_add_doc_line("            n4=1 f1=0 RTM image, n4=1 f1=1 illumination")
    call sep_add_doc_line("")
    call sep_add_doc_line("DESCRIPTION")
    call sep_add_doc_line("      Acoustic Reverse Time Migration in 2D and 3D")
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

  end subroutine RTM_doc

end module THREED_Doc
