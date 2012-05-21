#include <stdlib.h>
#include <locale.h>
#include "error.h"



int
FSinit(void)
{
  char *cp ;
  if (getenv("FS_DISABLE_LANG") == NULL)
    cp = setlocale(LC_NUMERIC, "en_US");
  return(NO_ERROR) ;
}

