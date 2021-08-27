/**
 * @brief Creates a modified image with noise or transformed 
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

//
// mri_create_tests.cpp
//
// written by Martin Reuter
// Oct. 12th ,2009
//
////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include <ctime> 
#include <cstdlib>

#include "Quaternion.h"
#include "MyMRI.h"

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include <vnl/vnl_matrix.h>

using namespace std;

struct Parameters
{
  string in;
  string in_t;
  string outs;
  string outt;
  string ltain;
  string ltaout;
  string mask;
  MRI* mri_in;
  double noise;
  int outlier;
  int outlierbox;
  bool translation;
  bool rotation;
  double iscale;
  bool doiscale;
  string iscaleout;
  string ltaouts;
  string ltaoutt;
  double transdist;
  double maxdeg;
};

static struct Parameters P =
{ "", "", "", "", "", "", "", NULL, 0.0, 0, -1, false, false, 1.0, false, "",
    "", "" , 11, 25};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);

const char *Progname = NULL;

std::vector<int> get_random(int lowest, int highest, int num = 3)
{

  unsigned int ttt = time(0);
  //cout << " seed: " << ttt << endl;
  srand(ttt);
  vector<int> ret(num);
  //int lowest=1, highest=10; 
  int range = (highest - lowest) + 1;
  for (int index = 0; index < num; index++)
  {
    int r = rand();
    //cout << " rand: " << r/(RAND_MAX + 1.0) << endl;
    ret[index] = lowest + int(range * (r / (RAND_MAX + 1.0)));
  }
  return ret;
}

void testmalloc()
{

  double * a = (double*) malloc(sizeof(double) * 200 * 111322800);
  if (a == NULL)
    cout << " not enough mem" << endl;
  else
    cout << "OK" << endl;

  a[5053030 + 15 * 111322800] = 3.0;
  cout << "a: " << a[5053030 + 15 * 111322800] << endl;
  free(a);
}

void testvnl()
{

  vnl_matrix<double> A;
  bool OK = A.set_size(111322800, 7);
  if (OK)
    cout << " OK" << endl;
  else
    cout << "not OK" << endl;

  A.assert_size(111322800, 7);

  cout << "assert ok" << endl;

  A[0][0] = 3.0;
  cout << " success 0-0" << endl;
  A[0][5] = 3.0;
  cout << " success 0-5" << endl;
  A[1][0] = 3.0;
  cout << " success 1-0" << endl;
  A[10530300][0] = 3.0;
  cout << " success 10530300-0" << endl;
  A[10530300][5] = 3.0;
  cout << " success 10530300-5" << endl;

  A[50530303][0] = 3.0;
  cout << " success 505...0" << endl;
  A[50530303][1] = 3.0;
  cout << " success 1" << endl;
  A[50530303][5] = 3.0;
  cout << " success 5" << endl;
  exit(0);

}

int main(int argc, char *argv[])
{

  // testmalloc();
  // testvnl();

  { // for valgrind, so that everything is freed
    cout << getVersion() << endl;

    // Default initialization
    int nargs = handleVersionOption(argc, argv, "mri_create_tests");
    if (nargs && argc - nargs == 1)
      exit(0);
    argc -= nargs;
    Progname = argv[0];
    argc--;
    argv++;
    ErrorInit(NULL, NULL, NULL);
//  DiagInit(NULL, NULL, NULL) ;

    if (!parseCommandLine(argc, argv, P))
    {
      printUsage();
      exit(1);
    }

    // ================================================================================= 

    // read input
    MRI* mriS = MRIread(P.in.c_str());
    assert(mriS!= NULL);
    MRI* mriT = NULL;
    if (P.in_t == "")
      mriT = MRIcopy(mriS, NULL);
    else
      mriT = MRIread(P.in_t.c_str());
    assert(mriT!= NULL);

    // mask target
    if (P.mask != "")
    {
      cout << " Applying mask " << P.mask << " to SOURCE" << endl;
      MRI *mri_mask;
      mri_mask = MRIread(P.mask.c_str());
      if (!mri_mask)
        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
            Progname, P.mask.c_str());
      MRImask(mriS, mri_mask, mriS, 0, 0);
      MRIfree(&mri_mask);
    }

    // read/construct lta
    if (P.ltain != "" && (P.translation || P.rotation))
    {
      cerr << " Cannot specify lta-in AND (translation OR rotation)" << endl;
      exit(1);
    }
    LTA * lta = NULL;
    if (P.ltain != "")
    {
      // try to read other transform
      TRANSFORM * trans = TransformRead(P.ltain.c_str());
      lta = (LTA *) trans->xform;
      if (!lta)
        ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
            Progname, P.ltain.c_str());
      lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      if (lta->type != LINEAR_VOX_TO_VOX)
      {
        ErrorExit(ERROR_BADFILE, "%s: must be LINEAR_VOX_TO_VOX (=0), but %d",
            Progname, P.ltain.c_str(), lta->type);
      }
      //R.setMinit(lta->xforms[0].m_L);
      cout << " Read " << P.ltain << " transform." << endl;
    }
    else
    {
      lta = LTAalloc(1, mriS);
      lta->type = LINEAR_VOX_TO_VOX;
      //MatrixPrintFmt(stdout,"% 2.8f",lta->xforms[0].m_L); cout << endl <<endl;
    }

    if (P.translation)
    {
      vector<int> t = get_random(-100, 100, 3);
//    vector < int > t(3,4);
      //cout << " T: " << t[0] << " " << t[1] << " " << t[2] << endl;
      //cout << " length: " << sqrt(t[0]*t[0] + t[1] * t[1] + t[2] * t[2]) << endl;

      assert(t.size() == 3);
///    double transdist = 100 ; // large 100mm=10 cm
//    double transdist = 50 ;
//     double transdist = 11;
//      double transdist = 0.05;
      double ff = 0.5 * P.transdist
          / sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
      //cout << " ff: " << ff << endl;
      float t0 = (float) (ff * t[0]);
      float t1 = (float) (ff * t[1]);
      float t2 = (float) (ff * t[2]);
      //cout << " New length: " << sqrt(t[0]*t[0] + t[1] * t[1] + t[2] * t[2]) << endl;
      cout << " Random Translation: ( " << 2 * t[0] << " , " << 2 * t[1]
          << " , " << 2 * t[2] << " )  length: "
          << 2 * sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) << endl;
      *MATRIX_RELT(lta->xforms[0].m_L, 1, 4) = t0;
      *MATRIX_RELT(lta->xforms[0].m_L, 2, 4) = t1;
      *MATRIX_RELT(lta->xforms[0].m_L, 3, 4) = t2;
      //*MATRIX_RELT(lta->xforms[0].m_L, 1, 4) = transdist / 2.0;
      //*MATRIX_RELT(lta->xforms[0].m_L, 2, 4) = 0;
      //*MATRIX_RELT(lta->xforms[0].m_L, 3, 4) = 0;
    }

    cout<< endl << " Matrix so far: " << endl;
    MatrixPrintFmt(stdout, "% 2.8f", lta->xforms[0].m_L);
    cout << endl << endl;

    if (P.rotation) // create half the rotation 
    {
      vector<int> t = get_random(-100, 100, 4);
      double length = sqrt(t[1] * t[1] + t[2] * t[2] + t[3] * t[3]);
      Quaternion Q;
      //double maxdeg = 40.0; // large 40 degree
      //double maxdeg = 25.0;
      double maxrad = 2.0 * M_PI * P.maxdeg / 360;
      //double rot    = 0.5 * maxrad * t[0]/100.0;
      double rot = 0.5 * maxrad;
      Q.importRotVec(rot, t[1] / length, t[2] / length, t[3] / length);
      vector<double> R = Q.getRotMatrix3d();
      cout << " Random Rotation: " << endl;
      for (int r = 0; r < 3; r++)
      {
        for (int c = 0; c < 3; c++)
        {
          //cout << R[r*3+c] << " " << flush;
          *MATRIX_RELT(lta->xforms[0].m_L, r+1, c+1) = R[r * 3 + c];
        }
        //cout << endl;
      }
      //MatrixPrintFmt(stdout,"% 2.8f",lta->xforms[0].m_L); cout << endl <<endl;

      // should be around center of the image:
      MATRIX * T1 = MatrixIdentity(4, NULL);
      MATRIX * T2 = MatrixIdentity(4, NULL);
      *MATRIX_RELT(T1, 1, 4) = -mriS->width / 2;
      *MATRIX_RELT(T2, 1, 4) = mriS->width / 2;
      *MATRIX_RELT(T1, 2, 4) = -mriS->height / 2;
      *MATRIX_RELT(T2, 2, 4) = mriS->height / 2;
      *MATRIX_RELT(T1, 3, 4) = -mriS->depth / 2;
      *MATRIX_RELT(T2, 3, 4) = mriS->depth / 2;
      lta->xforms[0].m_L = MatrixMultiply(lta->xforms[0].m_L, T1,
          lta->xforms[0].m_L);
      lta->xforms[0].m_L = MatrixMultiply(T2, lta->xforms[0].m_L,
          lta->xforms[0].m_L);
      //MatrixPrintFmt(stdout,"% 2.8f",lta->xforms[0].m_L); cout << endl <<endl;

      MatrixFree(&T1);
      MatrixFree(&T2);

    }
    assert(lta != NULL);

    // apply lta to image:  // symmetric
    MATRIX* a = MatrixCopy(lta->xforms[0].m_L, NULL);
    MATRIX* ai = MatrixInverse(lta->xforms[0].m_L, NULL);
    lta->xforms[0].m_L = MatrixMultiply(a, a, lta->xforms[0].m_L);

    cout << " Final Transform Matrix: " << endl;
    MatrixPrintFmt(stdout, "% 2.8f", lta->xforms[0].m_L);
    cout << " Determinant: " << MatrixDeterminant(lta->xforms[0].m_L) << endl
        << endl;

    mriS = MRIlinearTransformInterp(mriS, NULL, ai, SAMPLE_TRILINEAR);
    mriT = MRIlinearTransformInterp(mriT, NULL, a, SAMPLE_TRILINEAR);

    // iscale random
    double iscale = P.iscale;
    if (P.doiscale)
    {
      vector<int> s = get_random(95, 105, 1);
      iscale = s[0] / 100.0;
    }

    if (iscale != 1.0) // symmetric
    {
      cout << " Scaling Intenisty: " << iscale << endl;
      mriS = MyMRI::MRIvalscale(mriS, mriS, 1.0 / sqrt(iscale));
      mriT = MyMRI::MRIvalscale(mriT, mriT, sqrt(iscale));
    }

    // noise
    if (P.noise > 0.0)
    {
      cout << " Applying noise to image: " << P.noise << endl;
      MRI * mri_noise = MRIrandn(mriT->width, mriT->height, mriT->depth, 1, 0,
          P.noise, NULL);
      MRImaskZero(mri_noise, mriT, mri_noise);
      MRIadd(mriT, mri_noise, mriT);
      MRIfree(&mri_noise);
    }

    // outlier
    if (P.outlier > 0)
    {
//     cout << " Setting " << P.outlier << " random voxels to [200...255]" << endl;
//     vector <int> p = get_random(0,255,P.outlier*3);
//     vector <int> t = get_random(0,255,P.outlier);
//     for (int i = 0;i<P.outlier;i++)
//        MRIvox(mriT,p[i*3],p[i*3+1],p[i*3+2]) = t[i];
//     cout << " Setting " << P.outlier << " random voxel boxes 20^3" << endl;
      int bsize = 30; // should be even number
      int bsizeh = bsize / 2;
      cout << " Creating " << P.outlier << " outlier boxes " << bsize
          << "^3 (image copies)" << endl;
      vector<int> p = get_random(bsizeh + 1, 254 - bsizeh, P.outlier * 3);
      vector<int> ps = get_random(bsizeh + 50, 200 - bsizeh, P.outlier * 3);
      //vector <int> t = get_random(0,255,P.outlier);
      for (int i = 0; i < P.outlier; i++)
        for (int x = 0; x < bsize; x++)
          for (int y = 0; y < bsize; y++)
            for (int z = 0; z < bsize; z++)
            {
              int xxs = ps[i * 3] + x - (bsizeh);
              int yys = ps[i * 3 + 1] + y - (bsizeh);
              int zzs = ps[i * 3 + 2] + z - (bsizeh);
              int xx = p[i * 3] + x - (bsizeh);
              int yy = p[i * 3 + 1] + y - (bsizeh);
              int zz = p[i * 3 + 2] + z - (bsizeh);
              float val;

              if (i < P.outlier / 2)
              {
                if (xx < 0 || yy < 0 || zz < 0 || xx >= mriT->width
                    || yy >= mriT->height || yy >= mriT->depth)
                  continue;
                if (xxs < 0 || yys < 0 || zzs < 0 || xxs >= mriT->width
                    || yys >= mriT->height || yys >= mriT->depth)
                  val = 0;
                else
                  val = MRIgetVoxVal(mriT, xxs, yys, zzs, 0);

                MRIsetVoxVal(mriT, xx, yy, zz, 0, val);
                //MRIvox(mriT,xx,yy,zz) = val;
              }
              else
              {
                if (xx < 0 || yy < 0 || zz < 0 || xx >= mriS->width
                    || yy >= mriS->height || yy >= mriS->depth)
                  continue;
                if (xxs < 0 || yys < 0 || zzs < 0 || xxs >= mriS->width
                    || yys >= mriS->height || yys >= mriS->depth)
                  val = 0;
                else
                  val = MRIgetVoxVal(mriS, xxs, yys, zzs, 0);

                MRIsetVoxVal(mriS, xx, yy, zz, 0, val);
//        MRIvox(mriS,xx,yy,zz) = val;
              }
            }
    }

    if (P.outlierbox > 0)
    {
      int offset = 128;
      cout << " Creating Single Outlier Box [ " << offset << " , "
          << P.outlierbox + offset << " ]^3 with values [200..255]" << endl;
      vector<int> t = get_random(200, 255,
          P.outlierbox * P.outlierbox * P.outlierbox);
      int count = 0;
      for (int x = 0; x < P.outlierbox; x++)
        for (int y = 0; y < P.outlierbox; y++)
          for (int z = 0; z < P.outlierbox; z++)
          {
            MRIvox(mriT, offset+x, offset+y, offset+z) = t[count];
            count++;
          }
        }

//====================== OUTPUT ==========================================

        //cout << " OUTPUT results ... " << endl;

        // output source and target
    cout << " OUTPUT source MRI : " << P.outs << endl;
    MRIwrite(mriS, P.outs.c_str());

    cout << " OUTPUT target MRI : " << P.outt << endl;
    MRIwrite(mriT, P.outt.c_str());

    // output lta:
    if (P.ltaout != "")
    {
      cout << " OUTPUT LTA : " << P.ltaout << endl;
      // add src and dst info
      getVolGeom(mriS, &lta->xforms[0].src);
      getVolGeom(mriT, &lta->xforms[0].dst);
      LTAwriteEx(lta, P.ltaout.c_str());
    }

    if (P.ltaouts != "")
    {
      cout << " OUTPUT LTA (input -> new source) : " << P.ltaouts << endl;
      lta->xforms[0].m_L = MatrixCopy(ai, lta->xforms[0].m_L);
      // add src and dst info
      getVolGeom(mriS, &lta->xforms[0].src);
      getVolGeom(mriS, &lta->xforms[0].dst);
      LTAwriteEx(lta, P.ltaouts.c_str());
    }

    if (P.ltaoutt != "")
    {
      cout << " OUTPUT LTA (input -> new target) : " << P.ltaoutt << endl;
      lta->xforms[0].m_L = MatrixCopy(a, lta->xforms[0].m_L);
      // add src and dst info
      getVolGeom(mriT, &lta->xforms[0].src);
      getVolGeom(mriT, &lta->xforms[0].dst);
      LTAwriteEx(lta, P.ltaoutt.c_str());
    }

    MRIfree(&mriS);
    MRIfree(&mriT);

    // output iscale:
    if (P.iscaleout != "")
    {
      if (iscale != 1.0)
      {
        cout << " OUTPUT iscale file : " << P.iscaleout << endl;
        ofstream f(P.iscaleout.c_str(), ios::out);
        f << iscale;
        f.close();
      }
      else
        cout << " No iscale used, therefore will not ouput (iscale = 1.0)"
            << endl;
    }

  }
}

/*----------------------------------------------------------------------
 ----------------------------------------------------------------------*/
static void printUsage(void)
{
  cout << endl << endl;
  cout << "Usage: mri_create_tests <required arguments>" << endl << endl;

  cout << "Creates test cases for the registraition by mapping" << endl;
  cout << "the input to a source (half way backward) and to a" << endl;
  cout << "target (half way forward)." << endl;
  cout << endl;

  cout << "Required arguments" << endl << endl;
  cout << "  --in   invol.mgz       input volume to be modified" << endl;
  cout << "  --outs outsrc.mgz      output source volume name" << endl;
  cout << "  --outt outtarget.mgz   output target volume name" << endl;
  cout << endl;

  cout << "Optional arguments" << endl << endl;
  cout << "  --int  intvol.mgz      input target volume to be modified" << endl;
  cout << "                           must be in same space as invol.mgz"
      << endl;
  cout << "                           default: use invol to create outtarget"
      << endl;
  cout << "  --lta-in lta           specify lta for mapping input to outtarget"
      << endl;
  cout
      << "                           (inverse will be used to create outsource)"
      << endl;
  cout
      << "                           cannot be used with --rotation or --translation"
      << endl;
  cout << "  --mask mask.mgz        mask src mri with mask.mgz" << endl;
  cout << "  --noise <float>        add global Gaussian noise" << endl;
  cout << "  --outlier <int>        add <int> outlier voxel randomly" << endl;
  cout << "  --outlier-box <int>    add box 0..<int> containing random voxels"
      << endl;
  cout << "  --translation          apply random translation" << endl;
  cout << "  --transdist            set maximal trans. distance in mm (default 11)" << endl;
  cout << "  --rotation             apply random rotation" << endl;
  cout << "  --maxdeg               maximal rotation in degree (default 25)" << endl;
  cout << "  --intensity            apply random intensity scaling" << endl;
  cout << "  --iscale <double>      use as fixed intensity scaling parameter"
      << endl;
  cout << "  --lta-out lta          write used random transform to lta" << endl;
  cout << "  --lta-outs lta         write half way lta for source" << endl;
  cout << "  --lta-outt lta         write half way lta for target" << endl;
  cout << "  --iscale-out <string>  write used intensity scaling parameter"
      << endl;

  cout << endl << endl;
  cout
      << " If --translation and/or --rotation is specified, a matrix A is generated "
      << endl;
  cout
      << " (for mapping the input to the outtarget), then the input is also mapped via "
      << endl;
  cout
      << " the inverse of A to outsource. Therefore, --lta-out is A*A (the map from"
      << endl;
  cout << " outsource to outtarget), --lta-outs is Inv(A) and --lta-outt is A."
      << endl;
  cout << " If the same transform is to be applied to a different input image, "
      << endl;
  cout
      << " you need to first output the --lta-outt (A) and then pass it for the "
      << endl;
  cout << " different input via --lta-in in a subsequent call." << endl;

  cout << endl << endl;

  cout << " Report bugs to: Freesurfer@nmr.mgh.harvard.edu" << endl;

  cout << endl;

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
    option = option + 1;  // remove second '-'
  StrUpper(option);

  //cout << " option: " << option << endl;

  if (!strcmp(option, "IN"))
  {
    P.in = string(argv[1]);
    nargs = 1;
    cout << "Using " << P.in << " as input volume." << endl;
  }
  else if (!strcmp(option, "INT"))
  {
    P.in_t = string(argv[1]);
    nargs = 1;
    cout << "Using " << P.in_t << " as input target volume." << endl;
  }
  else if (!strcmp(option, "OUTS"))
  {
    P.outs = string(argv[1]);
    nargs = 1;
    cout << "Using " << P.outs << " as output source volume." << endl;
  }
  else if (!strcmp(option, "OUTT"))
  {
    P.outt = string(argv[1]);
    nargs = 1;
    cout << "Using " << P.outt << " as output target volume." << endl;
  }
  else if (!strcmp(option, "LTA-IN"))
  {
    P.ltain = string(argv[1]);
    nargs = 1;
    cout << "Input transform as " << P.ltain << " . " << endl;
  }
  else if (!strcmp(option, "LTA-OUT"))
  {
    P.ltaout = string(argv[1]);
    nargs = 1;
    cout << "Storing transform as " << P.ltaout << " . " << endl;
  }
  else if (!strcmp(option, "LTA-OUTS"))
  {
    P.ltaouts = string(argv[1]);
    nargs = 1;
    cout << "Storing half way source lta (input -> out-source) " << P.ltaouts
        << " . " << endl;
  }
  else if (!strcmp(option, "LTA-OUTT"))
  {
    P.ltaoutt = string(argv[1]);
    nargs = 1;
    cout << "Storing half way target lta (input -> out-target) " << P.ltaoutt
        << " . " << endl;
  }
  else if (!strcmp(option, "MASK"))
  {
    P.mask = string(argv[1]);
    nargs = 1;
    cout << "Using mask " << P.mask << " . " << endl;
  }
  else if (!strcmp(option, "NOISE"))
  {
    P.noise = atof(argv[1]);
    nargs = 1;
    cout << "Using global Gaussian noise " << P.noise << " ." << endl;
  }
  else if (!strcmp(option, "OUTLIER"))
  {
    P.outlier = atoi(argv[1]);
    nargs = 1;
    cout << "Randomly inserting " << P.outlier << " outlier voxel." << endl;
  }
  else if (!strcmp(option, "OUTLIER-BOX"))
  {
    P.outlierbox = atoi(argv[1]);
    nargs = 1;
    cout << "Inserting outlier box at 128 .. " << P.outlierbox << " ." << endl;
  }
  else if (!strcmp(option, "TRANSLATION"))
  {
    P.translation = true;
    nargs = 0;
    cout << "Creating random translation." << endl;
  }
  else if (!strcmp(option, "TRANSDIST"))
  {
    P.transdist = atof(argv[1]);
    nargs = 1;
    cout << "Translation distance" << P.transdist <<" ." << endl;
  }
  else if (!strcmp(option, "ROTATION"))
  {
    P.rotation = true;
    nargs = 0;
    cout << "Creating random rotation." << endl;
  }
  else if (!strcmp(option, "MAXDEG"))
  {
    P.maxdeg = atof(argv[1]);
    nargs = 1;
    cout << "Max rotation degree" << P.maxdeg <<" ." << endl;
  }
  else if (!strcmp(option, "INTENSITY"))
  {
    P.doiscale = true;
    nargs = 0;
    cout << "Applying random intensity scaling." << endl;
  }
  else if (!strcmp(option, "ISCALE"))
  {
    P.iscale = atof(argv[1]);
    nargs = 1;
    cout << "Using intensity scale " << P.iscale << " . " << endl;
  }
  else if (!strcmp(option, "ISCALE-OUT"))
  {
    P.iscaleout = string(argv[1]);
    nargs = 1;
    cout << "Writing intensity scale as " << P.iscaleout << " . " << endl;
  }
  else
  {
    cerr << "Option: " << argv[0] << " unknown !! " << endl;
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
  int nargs;

  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv, P);
    argc -= nargs;
    argv += nargs;
  }

  return (P.in != "" && P.outs != "" && P.outt != "");
}

