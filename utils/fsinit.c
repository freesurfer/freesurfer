#include <locale.h>
#include "error.h"



int
FSinit(void)
{
  char *cp ;
  cp = setlocale(LC_NUMERIC, "en_US");
  return(NO_ERROR) ;
}

