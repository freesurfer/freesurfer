#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "diag.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "label.h"
#include "mri.h"
#include "mrisegment.h"


const char *Progname ;

static int use_abs ;
static int size_thresh = 10 ;



int
main(int argc, char *argv[])
{
  double           thresh ;
  MRI              *mri, *mri_abs ;
  char             *out_stem, fname[STRLEN] ;
  MRI_SEGMENTATION *mriseg ;
  int              s ;
  LABEL            *area ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  mri = MRIread(argv[1]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load MRI from %s\n", Progname, argv[1]) ;

  if (use_abs)
    mri_abs = MRIabs(mri, NULL) ;
  else
    mri_abs = MRIcopy(mri, NULL) ;

  thresh = atof(argv[2]) ;
  out_stem = argv[3] ;

  mriseg = MRIsegment(mri, thresh, 1e10) ;
  MRIremoveSmallSegments(mriseg, size_thresh) ;
  printf("segmenting volume at threshold %2.1f yields %d segments\n", thresh, mriseg->nsegments) ;

  for (s = 0 ; s < mriseg->nsegments ; s++)
  {
    area = MRIsegmentToLabel(mriseg, mri_abs, s) ;
    sprintf(fname, "%s.%3.3d.label", out_stem, s) ;
    LabelWrite(area, fname) ;
    
  }
  return(0) ;
}
  

