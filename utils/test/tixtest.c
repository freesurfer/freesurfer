//
// tixtest.c
//

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <tix.h>

int main(int argc, char *argv[])
{
  int eTcl;

  Tcl_Interp * interp = 0;
  interp = Tcl_CreateInterp();
  
  /* read tcl/tk internal startup scripts */
  eTcl = Tcl_Init( interp );
  if( TCL_OK != eTcl ) {
    printf("Tcl_Init returned %d: %s\n", (int)eTcl, interp->result );
  }
  eTcl = Tk_Init( interp );
  if( TCL_OK != eTcl ) {
    printf("Tk_Init returned %d: %s\n", (int)eTcl, interp->result );
  }
  eTcl = Itcl_Init(interp);
  if (TCL_OK != eTcl)
  {
    printf("Itcl_Init reutrned %d: %s\n", (int)eTcl, interp->result);
  }
  eTcl = Tix_Init( interp );
  if( TCL_OK != eTcl ) {
    printf("Tix_Init returned %d: %s\n", (int)eTcl, interp->result );
  }
  return 0;
}
