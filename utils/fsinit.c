#include <stdlib.h>
#include <locale.h>
#include "error.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif



int
FSinit(void)
{
  char *cp ;
  if (getenv("FS_DISABLE_LANG") == NULL)
    cp = setlocale(LC_NUMERIC, "en_US");
  if (getenv("OMP_NUM_THREADS") == NULL)
    omp_set_num_threads(1);
  return(NO_ERROR) ;
}

