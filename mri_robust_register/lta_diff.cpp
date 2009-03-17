#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include "Registration.h"

// all other software are all in "C"
#ifdef __cplusplus
extern "C" {
#endif
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"

#ifdef __cplusplus
}
#endif

using namespace std;

static char vcid[] = "$Id: lta_diff.cpp,v 1.2 2009/03/17 19:24:59 mreuter Exp $";
char *Progname = NULL;

double cornerdiff(LTA* lta)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta,LINEAR_VOX_TO_VOX);
  return -1;
}

double interpolationError(LTA* lta)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta,LINEAR_VOX_TO_VOX);

  // sample from dst back to src
  MATRIX *mAinv = MatrixInverse(lta->xforms[0].m_L, NULL) ;
  if (!mAinv)
    ErrorExit(ERROR_BADPARM, "interpolationError: xform is singular") ;
  int width  = lta->xforms[0].dst.width ;
  int height = lta->xforms[0].dst.height ;
  int depth  = lta->xforms[0].dst.depth ;
  VECTOR * v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  VECTOR * v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */
  int y3,y2,y1;
  Real  x1, x2, x3 ;
  
  double errorsum = 0;
  v_Y->rptr[4][1] = 1.0f ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    V3_Z(v_Y) = y3 ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      V3_Y(v_Y) = y2 ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        V3_X(v_Y) = y1 ;
        MatrixMultiply(mAinv, v_Y, v_X) ;

        x1 = V3_X(v_X) ;
        x2 = V3_Y(v_X) ;
        x3 = V3_Z(v_X) ;
	
	errorsum += 1; //wrong
      }
    }
  }

  return errorsum;

}

int main(int argc, char *argv[])
{

//   // Default initialization
//   int nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
//   if (nargs && argc - nargs == 1) exit (0);
//   argc -= nargs;
   Progname = argv[0] ;
//   argc --;
//   argv++;
   if (vcid) {};

  if (argc < 3)
  {
     cout << argv[0] << " file1.lta file2.lta " << endl;
     exit(1);
  }
  string lta1f = argv[1];
  string lta2f = argv[2];
  
  double d = 1.0;
  if (argc >3 )
     d = atof(argv[3]);
     

  LTA* lta1 = LTAreadEx(lta1f.c_str());
  LTA* lta2 = LTAreadEx(lta2f.c_str());
  
  Registration R;
  cout << sqrt(R.RigidTransDistSq(lta1->xforms[0].m_L, lta2->xforms[0].m_L))/d;
  

}
