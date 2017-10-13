module LBFGS_doc_mod

  use sep
  
  implicit none

contains

  subroutine LBFGS_doc()
    call sep_add_doc_line("NAME")
    call sep_add_doc_line("    L_BFGS.x - Limited Memory BFGS ")
    call sep_add_doc_line("")
    call sep_add_doc_line("SYNOPSIS")
    call sep_add_doc_line("    L_BFGS.x Model=mod.H Gradient=grad.H Function=fct.H Initial_Model=initmod.H Gradient_Weight=gw.H xmin= xmax= lbfgs_type= bound_type= stp1_opt= > /dev/null ")
    call sep_add_doc_line("")
    call sep_add_doc_line("INPUT PARAMETERS")
    call sep_add_doc_line("")
    call sep_add_doc_line("")
    call sep_add_doc_line("OPTIONAL INPUT FILE")

    call sep_add_doc_line("")
    call sep_add_doc_line("OUTPUT PARAMETERS")

    call sep_add_doc_line("")
    call sep_add_doc_line("DESCRIPTION")
    call sep_add_doc_line("")
    call sep_add_doc_line("")
    call sep_add_doc_line("COMMENTS")
    call sep_add_doc_line("")

    call doc('SOURCE')

  end subroutine LBFGS_DOC

end module LBFGS_doc_mod
