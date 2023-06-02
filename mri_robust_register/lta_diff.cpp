/**
 * @brief A programm to compute differences of lta files (transforms)
 *
 */

/*
 * Original Author: Martin Reuter
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include <limits>
#include <vcl_iostream.h>

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_matlab_print.h>
#undef export

#include "Registration.h"
#include "MyMatrix.h"

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"

using namespace std;

struct Parameters
{
  string progname;
  string t1name;
  string t2name;
  int disttype;
  double normdiv;
  bool invert1;
  bool invert2;
  bool vox2vox;
  double radius;
};

static struct Parameters P =
{ "","","",2,1.0,false,false,false, 100.0};



/*----------------------------------------------------------------------
 ----------------------------------------------------------------------*/
#include "lta_diff.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(lta_diff_help_xml, lta_diff_help_xml_len);
}


/*!
 \fn int parseNextCommand(int argc, char **argv)
 \brief Parses the command-line for next command
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       number of used arguments for this command
 */
static int parseNextCommand(int argc, char *argv[], Parameters & P)
{
  int nargs = 0;
  char *option;

  option = argv[0] + 1;                     // remove '-'
  if (option[0] == '-')
  {
    option = option + 1;  // remove second '-'
  }
  StrUpper(option);

  //cout << " option: " << option << endl;

  if (!strcmp(option, "DIST") || !strcmp(option, "D"))
  {
    P.disttype = atoi(argv[1]);
    nargs = 1;
    cout << "--dist: Computing distance type: " << P.disttype << " ." << endl;
  }
  else if (!strcmp(option, "NORMDIV") )
  {
    P.normdiv = atof(argv[1]);
    nargs = 1;
    cout << "--normdiv: divide final distance by " << P.normdiv << " ." << endl;
  }
  else if (!strcmp(option, "RADIUS") )
  {
    P.radius = atof(argv[1]);
    nargs = 1;
    cout << "--radius: use RMS radius " << P.radius << "mm." << endl;
  }
  else if (!strcmp(option, "INVERT1") )
  {
    P.invert1 = true;
    nargs = 0;
    cout << "--invert1: inverting first transform. " << endl;
  }
  else if (!strcmp(option, "INVERT2") )
  {
    P.invert2 = true;
    nargs = 0;
    cout << "--invert2: inverting second transform. " << endl;
  }
  else if (!strcmp(option, "VOX") )
  {
    P.vox2vox = true;
    nargs = 0;
    cout << "--vox: analysing VOX to VOX transform (after correcting for voxel size). " << endl;
  }
  else if (!strcmp(option, "HELP") || !strcmp(option, "H"))
  {
    printUsage();
    exit(1);
  }
  else
  {
    cerr << endl << endl << "ERROR: Option: " << argv[0] << " unknown !! "
        << endl << endl;
    exit(1);
  }

  fflush(stdout);

  return (nargs);
}

/*!
 \fn int parseCommandLine(int argc, char **argv)
 \brief Parses the command-line
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       if all necessary parameters were set
 */
static bool parseCommandLine(int argc, char *argv[], Parameters & P)
{
  // Check for version flag and exit if nothing else is passed:
  int nargs = handleVersionOption(argc, argv, "lta_diff");
  if (nargs && argc - nargs == 1)
  {
    exit(0);
  }
  argc -= nargs;
  
  // Read name of executable
  P.progname = argv[0];
  argc--;
  argv++;
  
  // Make sure we have at least 1 transform filename
  int inputargs = argc;
  if (inputargs == 0)
  {
    printUsage();
    exit(1);
  }
  
  cout << endl;
  // Read positional arguments
  if (ISOPTION(*argv[0]))
  {
    printUsage();
    cerr << endl << endl << "ERROR: Please specify a transform file first !  "
        << endl << endl;

    exit(1);
  }
  
  // read first filename
  P.t1name = argv[0];
  argc--;
  argv++;
  cout << "First transform is " << P.t1name << endl;
  
  // read second transform if passed
  if (argc > 0 && !ISOPTION(*argv[0]))
  {
    P.t2name = argv[0];
    if (P.t2name == "identity.nofile")
    {
      P.t2name = "";
      cout << "No need to pass identity.nofile as 2nd, omitting it will work the same." << endl;
    }
    else
      cout << "Second transform is " << P.t2name << endl;
      
    argc--;
    argv++;
  }

  // read disttype if passed as positional (to be backward compatible)
  if (argc > 0 && !ISOPTION(*argv[0]))
  {
    P.disttype= atoi(argv[0]);
    if (P.disttype == 0) 
    {
      cerr << "ERROR: --dist: \"" << P.disttype << "\" not valid, expecting <int> 1... !" << endl;
      exit(1);
    }
    cout << "Using distance " << P.disttype << " ..." << endl;
    argc--;
    argv++;
    cout << endl  <<"WARNING: passing 'distance type' as positional 3rd argument is deprecated and will be removed in future versions! Use --dist <int> ..." << endl << endl;
  }

  // read normdiv if passed as positional (to be backward compatible)
  if (argc > 0 && !ISOPTION(*argv[0]))
  {
    P.normdiv= atof(argv[0]);
    argc--;
    argv++;
    cout << "Using normdiv " << P.normdiv << " ..." << endl;
    cout << endl <<"WARNING: passing 'normdiv' as positional 4th argument is deprecated and will be removed in future versions! Use --normdiv <float> ..." << endl << endl;
  }

  // read invert if passed as positional (to be backward compatible)
  if (argc > 0 && !ISOPTION(*argv[0]))
  {
    int invert = atoi(argv[0]);
    if (invert == 1)
    {
      P.invert1 = true;
      cout << "Will invert first transform ..." << endl;
    }
    if (invert == 2)
    {
      P.invert2 = true;
      cout << "Will invert second transform ..." << endl;
    }
    argc--;
    argv++;
    cout << endl <<"WARNING: passing 'invert' as positional 5th argument is deprecated and will be removed in future versions! Use --invert1 or --invert2 ..." << endl << endl;
  }
  
  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv, P);
    argc -= nargs;
    argv += nargs;
  }


  bool test1 = (!P.invert2 || P.t2name != "");
  if (!test1)
  {
    printUsage();
    cerr << endl << endl << "ERROR: Please specify second transform with --invert2 !  "
        << endl << endl;
    exit(1);
  }

  if (P.t2name != "" && P.disttype == 6)
  {
    cerr << endl << endl << "ERROR: distance type 6 (interpolation) can only take a single transform." << endl << endl;
    exit(1);
  }

  if (P.vox2vox && P.disttype == 6)
  {
    cerr << endl << endl << "WARNING: distance type 6 (interpolation) is independent of VOX or RAS coords." << endl << endl;
  }

  return test1;
}

void writeVox2Vox(LTA * lta)
{
  cout << " convet to vox 2 vox" << endl;
  LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  cout << " writing" << endl;
  LTAwrite(lta, "test-vox2vox.lta");
}

double cornerdiff(LTA* lta1, LTA* lta2, bool vox2vox)
{
  // get vox2vox using lta geometry info
  if (vox2vox)
  {
    LTAchangeType(lta1, LINEAR_VOX_TO_VOX);
    LTAchangeType(lta2, LINEAR_VOX_TO_VOX);
  }
  else
  {
    LTAchangeType(lta1, LINEAR_RAS_TO_RAS);
    LTAchangeType(lta2, LINEAR_RAS_TO_RAS);
  }
  

  VECTOR * v_X  = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  VECTOR * v_XR = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  VECTOR * v_Y1 = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR * v_Y2 = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR_ELT(v_X,4) = 1;
  VECTOR_ELT(v_Y1,4)= 1;
  VECTOR_ELT(v_Y2,4)= 1;

  assert(lta1->xforms[0].src.depth  == lta2->xforms[0].src.depth);
  assert(lta1->xforms[0].src.height == lta2->xforms[0].src.height);
  assert(lta1->xforms[0].src.width  == lta2->xforms[0].src.width);

  int y3, y2, y1;
  double d = 0;
  for (y3 = 0; y3 < 2; y3++)
  {
    V3_Z(v_X)= y3 * (lta1->xforms[0].src.depth-1);
    for (y2 = 0; y2 < 2; y2++)
    {
      V3_Y(v_X) = y2 * (lta1->xforms[0].src.height-1);
      for (y1 = 0; y1 < 2; y1++)
      {
        V3_X(v_X) = y1* (lta1->xforms[0].src.width-1);
        if (vox2vox)
        {
          MatrixMultiply(lta1->xforms[0].m_L, v_X, v_Y1);
          MatrixMultiply(lta2->xforms[0].m_L, v_X, v_Y2);
        }
        else
        {
          MATRIX * mv2r;
          mv2r = vg_i_to_r(&lta1->xforms[0].src);
          MatrixMultiply(mv2r, v_X, v_XR);
          MatrixMultiply(lta1->xforms[0].m_L, v_XR, v_Y1);
          MatrixFree(&mv2r);
          
          mv2r = vg_i_to_r(&lta2->xforms[0].src);
          MatrixMultiply(mv2r, v_X, v_XR);
          MatrixMultiply(lta2->xforms[0].m_L, v_XR, v_Y2);
          MatrixFree(&mv2r);
        }
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

  VectorFree(&v_X);
  VectorFree(&v_XR);
  VectorFree(&v_Y1);
  VectorFree(&v_Y2);

  return d / 8.0;
}

double cornerdiff(LTA* lta1, bool vox2vox)
{
  // get vox2vox using lta geometry info
  
  if (vox2vox)
    LTAchangeType(lta1, LINEAR_VOX_TO_VOX);
  else
    LTAchangeType(lta1, LINEAR_RAS_TO_RAS);
  

  VECTOR * v_X  = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  VECTOR * v_Y1 = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR_ELT(v_X,4) = 1;
  VECTOR_ELT(v_Y1,4)= 1;

  int y3, y2, y1;
  double d = 0;
  double dmax = 0;
  for (y3 = 0; y3 < 2; y3++)
  {
    V3_Z(v_X)= y3 * (lta1->xforms[0].src.depth-1);
    for (y2 = 0; y2 < 2; y2++)
    {
      V3_Y(v_X) = y2 * (lta1->xforms[0].src.height-1);
      for (y1 = 0; y1 < 2; y1++)
      {
        V3_X(v_X) = y1* (lta1->xforms[0].src.width-1);
        
        if (vox2vox)
          MatrixMultiply(lta1->xforms[0].m_L, v_X, v_Y1);
        else
        {
          //map corner to ras, map it with Ras2ras and compute distance
          MATRIX * mv2r = vg_i_to_r(&lta1->xforms[0].src);
          MatrixMultiply(mv2r, v_X, v_X);
          MatrixMultiply(lta1->xforms[0].m_L, v_X, v_Y1);
          MatrixFree(&mv2r);
        }  
          
        double d1 = V3_X(v_Y1) - V3_X(v_X);
        double d2 = V3_Y(v_Y1) - V3_Y(v_X);
        double d3 = V3_Z(v_Y1) - V3_Z(v_X);
        double dd = sqrt(d1*d1 + d2*d2 + d3*d3);
        //cout << " dd: " << dd << endl;
        if ( dd > dmax) dmax = dd;
        d += dd;
        //cout << " corner : " << V3_X(v_X) << " , " <<  V3_Y(v_X) << " , " <<  V3_Z(v_X) << endl;
        //cout << "   mapped to "<< V3_X(v_Y1) << " , " <<  V3_Y(v_Y1) << " , " <<  V3_Z(v_Y1) << endl;
        //cout << "   mapped to "<< V3_X(v_Y2) << " , " <<  V3_Y(v_Y2) << " , " <<  V3_Z(v_Y2) << endl;
      }
    }
  }
  VectorFree(&v_X);
  VectorFree(&v_Y1);
  cout << " dmax: " << dmax << endl;
  return d / 8.0;
}

double determinant(MATRIX * M1, MATRIX* M2)
{

//   MATRIX* M = MatrixAlloc(4,4,MATRIX_REAL);
  MATRIX* M = MatrixCopy(M1, NULL);
  if (M2 != NULL)
  {
    //cout << " inverting" << endl;
    //M = MatrixInverse(M1,M);
    M = MatrixMultiply(M2, M, M);
  }

  double d = MatrixDeterminant(M);
  MatrixFree(&M);
  return d;
}

void testQuaternion(const vnl_matrix<double>& Rot)
{

  Quaternion Q;
  Q.importMatrix(Rot[0][0],Rot[0][1],Rot[0][2],Rot[1][0],Rot[1][1],Rot[1][2],Rot[2][0],Rot[2][1],Rot[2][2]);
  cout << "Quaternion: " << Q << endl;
  
  vnl_vector < double > v(3);
  v[0] = 1.2; v[1]=-0.5; v[2] = .7;
  cout << " v    :  " << v << endl;
  cout << " v rot:  " << Rot * v << endl;
  
  std::vector < double > v1 = Q.rotate(v[0],v[1],v[2]);
  cout << " v rotQ: " << v1[0] <<" " << v1[1] << " " << v1[2] << endl;
  
  std::vector < double > m = Q.getRotMatrix3d();
  cout << " M = [ " << m[0] << " " << m[1] << " " << m[2] << endl << m[3] << " " << m[4] << " " << m[5] << endl << m[6] << " " << m[7] << " " << m[8] << " ]" << endl; 

exit(1);
}

void decompose(MATRIX * M1, MATRIX* M2)
{
  vnl_matrix<double> m = MyMatrix::convertMATRIX2VNL(M1);

  if (M2 != NULL)
  {
    //cout << " inverting" << endl;
    //M = MatrixInverse(M1,M);
    vnl_matrix<double> m2 = MyMatrix::convertMATRIX2VNL(M2);
    m = m * m2;
  }

  cout << " Decompose M1*M2 into Rot * Shear * Scale + Trans: " << endl
      << endl;
      
  vnl_matrix<double> Rot, Shear;
  vnl_diag_matrix<double> Scale;
  MyMatrix::Polar2Decomposition(m.extract(3, 3), Rot, Shear, Scale);
  vnl_matlab_print(vcl_cout,Rot,"Rot",vnl_matlab_print_format_long);
  cout << endl;
  
  Quaternion Q;
  Q.importMatrix(Rot[0][0],Rot[0][1],Rot[0][2],Rot[1][0],Rot[1][1],Rot[1][2],Rot[2][0],Rot[2][1],Rot[2][2]);
  std::vector < double > v = Q.getRotVec();
  cout << "RotVec = [ " << v[0] << " " << v[1] << " " << v[2] << " ] " << endl;
  cout << endl;
  
  cout << "RotAngle = " << Q.getRotAngle() << endl;
  cout << endl;
  
  vnl_matlab_print(vcl_cout,Shear,"Shear",vnl_matlab_print_format_long);
  cout << endl;
  
  vnl_matlab_print(vcl_cout,Scale,"Scale",vnl_matlab_print_format_long);
  cout << endl;
  
  vnl_vector<double> t = m.extract(3, 1, 0, 3).get_column(0);
  vnl_matlab_print(vcl_cout,t,"Trans",vnl_matlab_print_format_long);
  cout << endl;
  
  cout << "AbsTrans = " << t.two_norm() << endl;
  cout << endl;
  
  cout << "Determinant = " << vnl_determinant(m) << endl << endl;

}

double sphereDiff(MATRIX * M1, MATRIX* M2, double r)
{

  MATRIX* M = MatrixAlloc(4, 4, MATRIX_REAL);
  if (M2 == NULL)
    M = MatrixCopy(M1, M);
  else
  {
    M = MatrixInverse(M1, M);
    M = MatrixMultiply(M2, M, M);
  }

  // MatrixPrintFmt(stdout,"% 2.8f",M);

  double dmax = 0;
  double dmin = std::numeric_limits<double>::infinity();
  double davg = 0;
  int counter = 0;

  VECTOR * v_A = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  VECTOR * v_B = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  VECTOR_ELT(v_A,4)= 1;
  VECTOR_ELT(v_B,4)= 1;

  // poles
  V3_X(v_A)= 0;
  V3_Y(v_A)= 0;
  V3_Z(v_A)= r;
  MatrixMultiply(M, v_A, v_B);
  // MatrixPrintFmt(stdout,"% 2.8f",v_A);
  // MatrixPrintFmt(stdout,"% 2.8f",v_B);

  double d1 = V3_X(v_B)- V3_X(v_A);
  double d2 = V3_Y(v_B)- V3_Y(v_A);
  double d3 = V3_Z(v_B)- V3_Z(v_A);
  double dd = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  dmax = dd;
  dmin = dd;
  davg = dd;
  counter++;
  V3_X(v_A)= 0;
  V3_Y(v_A)= 0;
  V3_Z(v_A)= -r;
  MatrixMultiply(M, v_A, v_B);
  d1 = V3_X(v_B)- V3_X(v_A);
  d2 = V3_Y(v_B)- V3_Y(v_A);
  d3 = V3_Z(v_B)- V3_Z(v_A);
  dd = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  if (dd > dmax)
    dmax = dd;
  if (dd < dmin)
    dmin = dd;
  davg += dd;
  counter++;

  // sample sphere with radius 1
  int max = 10;
  for (int i = -max + 1; i < max; i++)
  {
    double angle1 = (i * M_PI * 0.5) / max; // from -pi/2 to +pi/2
    // radius:
    double r1 = cos(angle1);
    double h = sin(angle1);
    // circumference is 2pi *r1, we want 4*max samples at aequator (where r=1 and cc 2pi)
    int max2 = int(4.0 * max * r1);
    for (int j = 0; j < max2; j++)
    {
      double angle2 = (2.0 * M_PI * j) / max2; // from 0 to 2pi
      V3_X(v_A)= r*r1*cos(angle2);
      V3_Y(v_A)= r*r1*sin(angle2);
      V3_Z(v_A)= r*h;
      MatrixMultiply(M, v_A, v_B);
      d1 = V3_X(v_B)- V3_X(v_A);
      d2 = V3_Y(v_B)- V3_Y(v_A);
      d3 = V3_Z(v_B)- V3_Z(v_A);
      dd = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
      if (dd > dmax)
        dmax = dd;
      if (dd < dmin)
        dmin = dd;
      davg += dd;
      counter++;
    }
  }
  davg = davg / counter;
//   cout << " max: " << dmax << " min: " << dmin << " avg: " << davg << endl;
  MatrixFree(&M);
  return dmax;
}

void testSphereDiff()
{

  MATRIX* T = MatrixIdentity(4, NULL);
  T->rptr[1][4] = 10.0;
  cout << sphereDiff(T, NULL, 100) << endl;

  // rotation 90degree around z axis
  MATRIX* R = MatrixIdentity(4, NULL);
  R->rptr[1][1] = 0;
  R->rptr[1][2] = -1;
  R->rptr[2][1] = 1;
  R->rptr[2][2] = 0;
  cout << sphereDiff(R, NULL, 100) << endl;

  exit(1);
}

double interpolationError2D(double angle)
{
  int side = 256;

  //MRI* mri_error = MRIalloc(side, side,1,MRI_FLOAT);
  MATRIX *a = MatrixAllocRotation(3, angle, Z_ROTATION);

  VECTOR * v_X = VectorAlloc(3, MATRIX_REAL);  // input (src) coordinates
  VECTOR * v_Y = VectorAlloc(3, MATRIX_REAL);  // transformed (dst) coordinates
  double errorsum = 0, x, y;
  int xm, ym;
  double val, xmd, ymd, xpd, ypd;  // d's are distances
  V3_Z(v_Y)= 0;
  for (int y1 = 0; y1 < side; y1++)
    for (int y2 = 0; y2 < side; y2++)
    {
      V3_X(v_Y)= y1 - 128;
      V3_Y(v_Y) = y2 - 128;
      MatrixMultiply(a, v_Y, v_X);

      x = V3_X(v_X) + 128;
      y = V3_Y(v_X) + 128;

//         xm = MAX((int)x, 0) ;
//         xp = MIN(side-1, xm+1) ;
//         ym = MAX((int)y, 0) ;
//         yp = MIN(side-1, ym+1) ;
      xm = (int)floor(x);
      ym = (int)floor(y);

      xmd = x - (float)xm;
      ymd = y - (float)ym;
      xpd = (1.0f - xmd);
      ypd = (1.0f - ymd);

      //cout << "x: " << x << " xm: " << xm << " xp: " << xp << " xmd: " << xmd << " xpd: " << xpd << endl;
      //cout << "y: " << y <<" ym: " << ym << " yp: " << yp << " ymd: " << ymd << " ypd: " << ypd << endl;
      //assert(x>=0);
      assert (xmd >= 0 && xpd >= 0);
      assert (ymd >= 0 && ypd >= 0);

      val = 0;// sum of distance to each coordinate (use smallest)
      if (xmd < xpd) val += xmd;
      else val += xpd;
      if (ymd < ypd) val += ymd;
      else val += ypd;
      //MRIFvox(mri_error,y1,y2,0) = (float)(val) ;
      errorsum += val;

    }
    //MRIwrite(mri_error,"mri_error.mgz");
    //MRIfree(&mri_error);
  return errorsum;
}

double interpolationError(LTA* lta)
{
  // get vox2vox using lta geometry info
  LTAchangeType(lta, LINEAR_VOX_TO_VOX);

  // sample from dst back to src
  MATRIX *mAinv = MatrixInverse(lta->xforms[0].m_L, NULL);
  if (!mAinv)
    ErrorExit(ERROR_BADPARM, "interpolationError: xform is singular");
  int width = lta->xforms[0].dst.width;
  int height = lta->xforms[0].dst.height;
  int depth = lta->xforms[0].dst.depth;
  VECTOR * v_X = VectorAlloc(4, MATRIX_REAL);  // input (src) coordinates
  VECTOR * v_Y = VectorAlloc(4, MATRIX_REAL);  // transformed (dst) coordinates
  int y3, y2, y1;
  double x, y, z;
  int xm, ym, zm;
  double val, xmd, ymd, zmd, xpd, ypd, zpd;  // d's are distances

  MRI* mri_error = MRIalloc(width, height, depth, MRI_FLOAT);

  double errorsum = 0;
  v_Y->rptr[4][1] = 1.0f;
  for (y3 = 0; y3 < depth; y3++)
  {
    V3_Z(v_Y)= y3;
    for (y2 = 0; y2 < height; y2++)
    {
      V3_Y(v_Y) = y2;
      for (y1 = 0; y1 < width; y1++)
      {
        V3_X(v_Y) = y1;
        MatrixMultiply(mAinv, v_Y, v_X);

        x = V3_X(v_X);
        y = V3_Y(v_X);
        z = V3_Z(v_X);

        xm = MAX((int)x, 0);
        ym = MAX((int)y, 0);
        zm = MAX((int)z, 0);

        xmd = x - (float)xm;
        ymd = y - (float)ym;
        zmd = z - (float)zm;
        xpd = (1.0f - xmd);
        ypd = (1.0f - ymd);
        zpd = (1.0f - zmd);

        val = 0; // sum of distance to each coordinate (use smallest)
        if (xmd < xpd) val += xmd;
        else val += xpd;
        if (ymd < ypd) val += ymd;
        else val += ypd;
        if (zmd < zpd) val += zmd;
        else val += zpd;
        MRIFvox(mri_error,y1,y2,y3) = (float)(val);
        errorsum += val;
      }
    }
  }
  MRIwrite(mri_error, "mri_error.mgz");
  MRIfree(&mri_error);
  return errorsum;

}

MATRIX* getIsoVOX(LTA *lta)
{
  LT *lt= &lta->xforms[0];
  
  if (lt->dst.valid == 0 || lt->src.valid == 0)
  {
    cerr << "ERROR:********************************************************\n";
    cerr << "ERROR: dst or src info invalid - cannot get voxel information.\n";
    cerr << "ERROR:********************************************************\n";
    exit(1);
  }
  
  MATRIX* Ms = MatrixAlloc(4, 4, MATRIX_REAL);
  Ms->rptr[1][1] = lt->src.xsize;
  Ms->rptr[2][2] = lt->src.ysize;
  Ms->rptr[3][3] = lt->src.zsize;
  Ms->rptr[4][4] = 1;
  
  MATRIX* Mt = MatrixAlloc(4, 4, MATRIX_REAL);
  Mt->rptr[1][1] = lt->dst.xsize;
  Mt->rptr[2][2] = lt->dst.ysize;
  Mt->rptr[3][3] = lt->dst.zsize;
  Mt->rptr[4][4] = 1;

  LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  MATRIX* VOX = MatrixCopy(lt->m_L,NULL);
  VOX = MatrixMultiply(Mt, VOX, VOX);
  VOX = MatrixMultiply(VOX, Ms, VOX);

  MatrixFree(&Ms);
  MatrixFree(&Mt);
  return VOX;
}

int main(int argc, char *argv[])
{
  if (!parseCommandLine(argc, argv, P)) exit(1);

  if (P.disttype == 100)
  {
    int steps = 200;
    double div = 16.0;
    vector<double> theta(steps);
    vector<double> err(steps);
    for (int i = 0; i < steps; i++)
    {
      // 0.. PI/div in 20 steps
      // -PI/div ..0 is symmetric
      theta[i] = M_PI * (i + 1) / ((steps) * div);
      err[i] = interpolationError2D(theta[i]);
    }
    ostringstream ss;
    ss << "interror-rot16";
    string fn = ss.str() + ".plot";
    ofstream f(fn.c_str(), ios::out);

    f << "set terminal postscript eps color" << endl;
    f << "set title \"Interpolation error when rotating \"" << endl;
    f << "set output \"" << ss.str() << ".eps\"" << endl;
    f << "plot  \"-\" notitle with lines 1" << endl;
    for (int i = 0; i < steps; i++)
    {
      cout << theta[i] << " " << err[i] << endl;
      f << theta[i] << " " << err[i] << endl;
    }
    f << "e" << endl;
    exit(0);
  }

  LTA* lta1 = LTAreadEx(P.t1name.c_str());
  if (!lta1)
  {
    cerr << "Could not open the first input file: " << P.t1name << endl;
    exit(1);
  }
  LTA* lta2 = NULL;
  if (P.t2name != "")
  {
    lta2 = LTAreadEx(P.t2name.c_str());
    if (!lta2)
    {
      cerr << "Could not open the second input file: " << P.t2name << endl;
      exit(1);
    }
  }
  //else
  //  lta2 = LTAreadEx("identity.nofile");


  if (P.invert1)
  {
    VOL_GEOM vgtmp;
    LT *lt;
    MATRIX *m_tmp = lta1->xforms[0].m_L;
    lta1->xforms[0].m_L = MatrixInverse(lta1->xforms[0].m_L, NULL);
    MatrixFree(&m_tmp);
    lt = &lta1->xforms[0];
    if (lt->dst.valid == 0 || lt->src.valid == 0)
    {
      cerr << "WARNING:********************************************************\n";
      cerr << "WARNING: dst or src volume is invalid.  Inverse likely wrong.\n";
      cerr << "WARNING:********************************************************\n";
    }
    //copyVolGeom(&lt->dst, &vgtmp);
    vgtmp = lt->dst;
    //copyVolGeom(&lt->src, &lt->dst);
    lt->dst = lt->src;
    //copyVolGeom(&vgtmp, &lt->src);
    lt->src = vgtmp;
  }
  if (P.invert2 )
  {
    if (!lta2)
    {
      cerr << "ERROR: cannot invert 2nd transform, as it is not given or identity!\n" ;
      exit(1);
    }
  
    VOL_GEOM vgtmp;
    LT *lt;
    MATRIX *m_tmp = lta2->xforms[0].m_L;
    lta2->xforms[0].m_L = MatrixInverse(lta2->xforms[0].m_L, NULL);
    MatrixFree(&m_tmp);
    lt = &lta2->xforms[0];
    if (lt->dst.valid == 0 || lt->src.valid == 0)
    {
      cerr << "WARNING:********************************************************\n";
      cerr << "WARNING:dst or src volume is invalid.  Inverse likely wrong.\n";
      cerr << "WARNING:********************************************************\n";
    }
    //copyVolGeom(&lt->dst, &vgtmp);
    vgtmp = lt->dst;
    //copyVolGeom(&lt->src, &lt->dst);
    lt->dst = lt->src;
    //copyVolGeom(&vgtmp, &lt->src);
    lt->src = vgtmp;
  }

  //LTAchangeType(lta1,LINEAR_VOX_TO_VOX);
  //MATRIX* VOX1 = MatrixCopy(lta1->xforms[0].m_L,NULL);
  MATRIX* VOX1 = getIsoVOX(lta1);
  LTAchangeType(lta1, LINEAR_RAS_TO_RAS);
  MATRIX* RAS1 = MatrixCopy(lta1->xforms[0].m_L,NULL);
  MATRIX* M1 = RAS1;
  if (P.vox2vox) M1 = VOX1;
  
  MATRIX* RAS2 = NULL;
  MATRIX* VOX2 = NULL;
  MATRIX* M2 = NULL;
  if (lta2)
  {
    //LTAchangeType(lta2,LINEAR_VOX_TO_VOX);
    //VOX2 = MatrixCopy(lta2->xforms[0].m_L,NULL);
    VOX2 = getIsoVOX(lta2);
    LTAchangeType(lta2, LINEAR_RAS_TO_RAS);
    RAS2 = MatrixCopy(lta2->xforms[0].m_L,NULL);
    M2 = RAS2;
    if (P.vox2vox) M2 = VOX2;
  }

  double dist = -1.0;
  switch (P.disttype)
  {
  case 1:
    dist = sqrt(MyMatrix::RigidTransDistSq(M1, M2));
    break;
  case 2:
    dist = sqrt(MyMatrix::AffineTransDistSq(M1, M2, P.radius));
    break;
  case 3:
    if (P.t2name == "")
      dist = cornerdiff(lta1, P.vox2vox);
    else
      dist = cornerdiff(lta1, lta2, P.vox2vox);
    break;
  case 4:
    dist = sphereDiff(M1, M2, 100);
    break;
  case 5:
    dist = determinant(M1, M2);
    break;
  case 6:
    dist = MyMatrix::getResampSmoothing(lta1); //independent of vox or ras
    break;
  case 7:
    decompose(M1, M2);
    exit(0);
    break;
  default:
    cerr << "ERROR: distance type " << P.disttype << " unknown!" << endl;
    exit(1);
    break;
  }
  if (P.disttype < 7)
    cout << dist/P.normdiv << endl;
    
  MatrixFree(&VOX1);
  MatrixFree(&RAS1);
  if (VOX2) MatrixFree(&VOX2);
  if (RAS2) MatrixFree(&RAS2);
}
