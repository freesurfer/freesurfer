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

static char vcid[] = "$Id: lta_diff.cpp,v 1.3 2009/03/17 23:16:27 mreuter Exp $";
char *Progname = NULL;

double cornerdiff(LTA* lta1, LTA* lta2)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta1,LINEAR_VOX_TO_VOX);
  LTAchangeType(lta2,LINEAR_VOX_TO_VOX);

  VECTOR * v_X  = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  VECTOR * v_Y1 = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */
  VECTOR * v_Y2 = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */

  assert (lta1->xforms[0].src.depth  == lta2->xforms[0].src.depth);
  assert (lta1->xforms[0].src.height == lta2->xforms[0].src.height);
  assert (lta1->xforms[0].src.width  == lta2->xforms[0].src.width);

  int y3,y2,y1;
  double d = 0;
  for (y3 = 0 ; y3 < 2 ; y3++)
  {
    V3_Z(v_X) = y3 * (lta1->xforms[0].src.depth-1) ;
    for (y2 = 0 ; y2 < 2 ; y2++)
    {
      V3_Y(v_X) = y2 * (lta1->xforms[0].src.height-1);
      for (y1 = 0 ; y1 < 2 ; y1++)
      {
        V3_X(v_X) = y1* (lta1->xforms[0].src.width-1) ;  
        MatrixMultiply(lta1->xforms[0].m_L, v_X, v_Y1) ;
        MatrixMultiply(lta2->xforms[0].m_L, v_X, v_Y2) ;
	double d1 = V3_X(v_Y1) - V3_X(v_Y2);
	double d2 = V3_Y(v_Y1) - V3_Y(v_Y2);
	double d3 = V3_Z(v_Y1) - V3_Z(v_Y2);
        d += sqrt(d1*d1 + d2*d2 + d3*d3);
	//cout << " corner : " << V3_X(v_X) << " , " <<  V3_Y(v_X) << " , " <<  V3_Z(v_X) << endl;
	//cout << "   mapped to "<< V3_X(v_Y1) << " , " <<  V3_Y(v_Y1) << " , " <<  V3_Z(v_Y1) << endl;
	//cout << "   mapped to "<< V3_X(v_Y2) << " , " <<  V3_Y(v_Y2) << " , " <<  V3_Z(v_Y2) << endl;	
      }
    }
  }
  
  return d/8.0;
}

double interpolationError(LTA* lta)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta,LINEAR_VOX_TO_VOX);

  // sample from dst back to src
  MATRIX *mAinv = MatrixInverse(lta->xforms[0].m_L, NULL) ;
  if (!mAinv) ErrorExit(ERROR_BADPARM, "interpolationError: xform is singular") ;
  int width  = lta->xforms[0].dst.width ;
  int height = lta->xforms[0].dst.height ;
  int depth  = lta->xforms[0].dst.depth ;
  VECTOR * v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  VECTOR * v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */
  int y3,y2,y1;
  double  x, y, z ;
  int  xm, xp, ym, yp, zm, zp;
  double val, xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */
  
  MRI* mri_error = MRIalloc(width, height,depth,MRI_FLOAT);
  
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

        x = V3_X(v_X) ;
        y = V3_Y(v_X) ;
        z = V3_Z(v_X) ;
	
        xm = MAX((int)x, 0) ;
        xp = MIN(width-1, xm+1) ;
        ym = MAX((int)y, 0) ;
        yp = MIN(height-1, ym+1) ;
        zm = MAX((int)z, 0) ;
        zp = MIN(depth-1, zm+1) ;

        xmd = x - (float)xm ;
        ymd = y - (float)ym ;
        zmd = z - (float)zm ;
        xpd = (1.0f - xmd) ;
        ypd = (1.0f - ymd) ;
        zpd = (1.0f - zmd) ;
	
	val = 0; // sum of distance to each coordinate (use smallest)
	if (xmd < xpd) val += xmd;
	else val += xpd;
	if (ymd < ypd) val += ymd;
	else val += ypd;
	if (zmd < zpd) val += zmd;
	else val += zpd;
	MRIFvox(mri_error,y1,y2,y3) = (float)(val) ;
	errorsum += val;
      }
    }
  }
  MRIwrite(mri_error,"mri_error.mgz");
  MRIfree(&mri_error);
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
     cout << endl;
     cout << argv[0] << " file1.lta file2.lta [norm-div] [dist-type]" << endl;
     cout << endl;
     cout << "    norm-div  (=1)  divide final distance by this (e.g. step adjustment)" << endl;
     cout << "    dist-type " << endl;
     cout << "       1  (default) Rigid Transform Distance (||log(R)|| + ||T||)" << endl;
     cout << "       2            8-corners distance after Transform " << endl;
     cout << endl;
     exit(1);
  }
  string lta1f = argv[1];
  string lta2f = argv[2];
  
  double d = 1.0;
  int disttype = 1;
  if (argc >3 ) d = atof(argv[3]);
  if (argc >4 ) disttype = atoi(argv[4]);
     

  LTA* lta1 = LTAreadEx(lta1f.c_str());
  LTA* lta2 = LTAreadEx(lta2f.c_str());
  
  if (!lta1 || !lta2)
  {
     cerr << "Could not open one of the LTA input files" << endl;
     exit(1);
  }
  
  double dist = -1;
  Registration R;
  
  switch (disttype)
  {
  case 1 :
     dist = sqrt(R.RigidTransDistSq(lta1->xforms[0].m_L, lta2->xforms[0].m_L))/d;
     break;
  case 2 :
     dist =  cornerdiff(lta1,lta2)/d;
     break;
   assert(1==2);
  }
  cout << dist << endl;
}
