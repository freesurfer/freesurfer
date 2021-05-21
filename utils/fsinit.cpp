#include <locale.h>
#include <stdlib.h>
#include "error.h"

#ifdef HAVE_OPENMP
#include "romp_support.h"
#endif


static void setLocale()
{
  if (getenv("FS_DISABLE_LANG") != NULL) return;
  if (setlocale(LC_NUMERIC, "en_US")) return;
  if (setlocale(LC_NUMERIC, "en_US.utf8")) return;
  if (setlocale(LC_NUMERIC, "en_US.UTF-8")) return;

  // use C.UTF-8 as a last resort
  if (setlocale(LC_NUMERIC, "C.UTF-8")) return;

  // but don't accept anything else
  fs::warning() << "Please install the en_US locale for numeric consistency";
}


int FSinit(void)
{
  setLocale();
#ifdef HAVE_OPENMP
  if (getenv("OMP_NUM_THREADS") == NULL) omp_set_num_threads(1);
#endif
  return (NO_ERROR);
}
