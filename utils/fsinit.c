#include <locale.h>
#include <stdlib.h>
#include "error.h"
#ifdef HAVE_OPENMP
#include "romp_support.h"
#endif

int FSinit(void)
{
  // char *cp;
  if (getenv("FS_DISABLE_LANG") == NULL)
    if (!setlocale(LC_NUMERIC, "en_US")) {
      ErrorPrintf(NO_ERROR, "Could not set locale");
    }
#ifdef HAVE_OPENMP
  if (getenv("OMP_NUM_THREADS") == NULL) omp_set_num_threads(1);
#endif
  return (NO_ERROR);
}
