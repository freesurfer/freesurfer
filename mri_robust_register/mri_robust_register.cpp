/**
 * @brief Linear registration of two volumes using robust statistics
 *
 * Based on ideas in "Robust Multiresolution Alignment of MRI Brain Volumes"
 * by Nestares and Heeger (2000)
 */

/*
 * Original Author: Martin Reuter, Nov. 4th ,2008
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

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#if ITK_VERSION_MAJOR >= 5
#include <vcl_compiler.h>
#else
#include <vcl_iostream.h>
#endif

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_matlab_print.h>
#undef export

#include "Registration.h"
#include "RegistrationStep.h"
#include "Regression.h"
#include "RegPowell.h"
#include "CostFunctions.h"
#include "MyMRI.h"
#include "MyMatrix.h"
#include "JointHisto.h"

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "mriBSpline.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"

using namespace std;

#define SAT -1 // leave blank, either passed by user or --satit
//#define SAT 4.685 // this is suggested for gaussian noise
//#define SAT 20
#define SSAMPLE -1
#define ERADIUS 5

struct Parameters
{
  string mov;
  string dst;
  string lta;
  string maskmov;
  string maskdst;
  string halfmov;
  string halfdst;
  string halfweights;
  string halfmovlta;
  string halfdstlta;
  string weightsout;
  bool satit;
  bool nomulti;
  bool conform;
  bool floattype;
  bool lta_vox2vox;
  bool affine;
  bool iscale;
  bool iscaleonly;
  bool transonly;
  bool isoscale;
  string transform;
  bool leastsquares;
  int iterate;
  double epsit;
  double sat;
  string warpout;
  string norlout;
  int subsamplesize;
  int debug;
  MRI* mri_mov;
  MRI* mri_dst;
  bool dosatest;
  bool initorient;
  bool inittrans;
  bool initscaling;
  int verbose;
  int highit;
  bool doubleprec;
  double wlimit;
  bool oneminusweights;
  bool symmetry;
  string iscaleout;
  string iscalein;
  int minsize;
  int maxsize;
  Registration::Cost cost;
//  int    bins;
  int finalsampletype;
  bool entropy;
  int entroradius;
  string entmov;
  string entdst;
  bool entball;
  bool entcorrection;
  bool sobel;
  double powelltol;
  bool whitebgmov;
  bool whitebgdst;
  bool uchartype;
};
static struct Parameters P =
{ "", "", "", "", "", "", "", "", "", "", "", false, false, false, false, false, false,
    false, false, false, false, "", false, 5, 0.01, SAT, "", "", SSAMPLE, 0,
    NULL, NULL, false, false, true, false, 1, -1, false, 0.16, true, true, "",
    "", -1, -1, Registration::ROB,
//  256,
    SAMPLE_CUBIC_BSPLINE, false, ERADIUS, "", "", false, false, false, 1e-5, false, false,false};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);
static void initRegistration(Registration & R, Parameters & P);

const char *Progname = NULL;

//static MORPH_PARMS  parms ;
//static FILE *diag_fp = NULL ;

void debug(Parameters &P)
{
  MRI * mriS = MRIread(P.mov.c_str());
  MRI * mriT = MRIread(P.dst.c_str());
  MRI * SmT = MRIallocSequence(mriS->width, mriS->height, mriS->depth, MRI_FLOAT, mriS->nframes);
  SmT = MRIsubtract(mriS, mriT, SmT);
  SmT = MyMRI::getBlur(SmT, SmT);

// MRI * mri = MRIread(P.mov.c_str());
// mri = MyMRI::getBlur(mri,mri);
// MRI * mri2 = MRIread(P.dst.c_str());
// MRI * mri = MRIallocHeader(256,256,256,MRI_UCHAR,1);
// MRI * mri2 = MRIallocHeader(122,256,256,MRI_UCHAR,1);
// MRIcopyFrame(mri,mri2,0,0);
  MRIfree(&mriS);
  MRIfree(&mriT);
  MRIfree(&SmT);
  //MRIfree(&mri2);
  exit(0);
}

void conv(MRI * i)
{
  cout << " adsf" << endl;
  RegRobust R;
  MRI * fmri = MRIalloc(i->width, i->height, i->depth, MRI_FLOAT);
  MRIcopyHeader(i, fmri);
  fmri->type = MRI_FLOAT;
  float f;
  for (int z = 0; z < i->depth; z++)
    for (int y = 0; y < i->height; y++)
      for (int x = 0; x < i->width; x++)
      {
        f = MRIgetVoxVal(i, x, y, z, 0);
        // cout << " f " << f << endl;
        MRIsetVoxVal(fmri, x, y, z, 0, f);
      }
  cout << "asdfasdf" << endl;
  MRIwrite(fmri, "float-1.mgz");
  MRI * sfmri;
  sfmri = MyMRI::MRIvalscale(fmri, NULL, 100);
  MRIwrite(sfmri, "float-100.mgz");
  sfmri = MyMRI::MRIvalscale(fmri, sfmri, 1000);
  MRIwrite(sfmri, "float-1000.mgz");
  exit(0);
}

void testRegression()
{

  int n = 200;
  vnl_matrix<double> A(n, 1);
  vnl_vector<double> b(n);

  for (int i = 0; i < n; i++)
  {
    A[i][0] = i;
    b[i] = 4 * i;
  }
  for (int i = 0; i < n; i += 5)
  {
    b[i] = 0;
  }

  Regression<double> R1(A, b);
  vnl_vector<double> M1 = R1.getLSEst();
  cout << M1 << endl;
  cout << endl << endl;

  Regression<double> R2(A, b);
  vnl_vector<double> M2 = R2.getRobustEst();
  cout << M1 << endl;
  cout << endl << endl;

  exit(0);

}

void testSubsamp(Parameters &P)
{

  MRI * mri1 = MRIread(P.mov.c_str());
//  MRI * mris = MRIdownsample2(mri1,NULL);
  //MRI * mris = MRIdownsample2BSpline(mri1,NULL);
  MRI * mris = MRIupsampleN(mri1, NULL, 2);
  MRIwrite(mris, "mri1up2.mgz");

  MRI * mril = MRIupsampleN(mri1, NULL, 3);
  MRIwrite(mril, "mri1up3.mgz");

  exit(0);
}

void entro(Parameters & P)
{

  //int sigma  = 7;
  int radius = 5;

//  cout << "Entropy sigma: " << sigma << "  radius: " << radius << endl;
  cout << "Entropy radius: " << radius << endl;

  cout << "Converting: " << P.mov.c_str() << endl;
  MRI * mri1 = MRIread(P.mov.c_str());
  //MRI * mri1e = MyMRI::entropyImage(mri1,radius,sigma);
  MRI * mri1e = MyMRI::entropyImage(mri1, radius);
  MRIwrite(mri1e, "mri1e.mgz");
  P.mov = "mri1e.mgz";

  cout << "Converting: " << P.dst.c_str() << endl;
  MRI * mri2 = MRIread(P.dst.c_str());
  //MRI * mri2e = MyMRI::entropyImage(mri2,radius,sigma);
  MRI * mri2e = MyMRI::entropyImage(mri2, radius);
  MRIwrite(mri2e, "mri2e.mgz");
  P.dst = "mri2e.mgz";

  exit(0);
}

void lnorm(Parameters & P)
{

  //int sigma  = 7;
  int radius = 2;

//  cout << "Entropy sigma: " << sigma << "  radius: " << radius << endl;
  cout << "Normalization Box radius: " << radius << endl;

  cout << "Converting: " << P.mov.c_str() << endl;
  MRI * mri1 = MRIread(P.mov.c_str());
  //MRI * mri1e = MyMRI::entropyImage(mri1,radius,sigma);
  MRI * mri1e = MyMRI::getNormalizedImage(mri1, radius);
  MRIwrite(mri1e, "mri1n.mgz");
  P.mov = "mri1n.mgz";

  cout << "Converting: " << P.dst.c_str() << endl;
  MRI * mri2 = MRIread(P.dst.c_str());
  //MRI * mri2e = MyMRI::entropyImage(mri2,radius,sigma);
  MRI * mri2e = MyMRI::getNormalizedImage(mri2, radius);
  MRIwrite(mri2e, "mri2n.mgz");
  P.dst = "mri2n.mgz";

  exit(0);
}


void jointhisto(Parameters & P)
{

  MRI* mriS = MRIread(P.mov.c_str());
  MyMRI::MRInorm255(mriS, mriS);
  MRI* mriT = MRIread(P.dst.c_str());
  MyMRI::MRInorm255(mriT, mriT);
  vnl_matrix_fixed<double, 4, 4> M;
  M.set_identity();
  JointHisto hm(mriS, mriT, M, M, 1, 1, 1);
  hm.smooth();
  hm.save("histo.m", "H");
  //hm.normalize();
  MRI * weights = hm.locate(mriS, mriT, M, M, 1, 1, 1);
  MRIwrite(weights, "weights.mgz");

//    int n = 10;
//    vnl_matrix < double > m(n,n);
//    int i,j;
//    for (i=0;i<n;i++)
//    for (j=0;j<n;j++)
//      m[i][j] = i+j+2;
//    m[2][3] = 15;
//    m[1][2] = 13;
//    m[6][9] = 4;
//    vnl_matlab_print(vcl_cerr,m,"m",vnl_matlab_print_format_long);std::cerr << std::endl;
//    
//    JointHisto h;
//    h.set(m); 
//    h.print("H");
//    cout << " NMI: " << h.computeNMI() << endl;
//    
//    h.smooth(1);
//    h.print("Hs");
//    cout << " NMI: " << h.computeNMI() << endl;
// //   
// //   cout << " MI : " << h.computeMI() << endl;
// //   cout << " ECC: " << h.computeECC() << endl;
// //   cout << " NMI: " << h.computeNMI() << endl;
// //   cout << " NCC: " << h.computeNCC() << endl;
// 
// exit(1);

//    MRI* mriS = MRIread(P.mov.c_str());
//    MyMRI::MRInorm255(mriS,mriS);
//    std::pair < float, float > mm = CostFunctions::minmax(mriS);
//    cout << " mriS   min: " << mm.first << "   max : " << mm.second << "   mean : " << CostFunctions::mean(mriS) << endl;
//    MRI* mriT = MRIread(P.dst.c_str());
//    MyMRI::MRInorm255(mriT,mriT);
//    mm = CostFunctions::minmax(mriT);
//    cout << " mriT   min: " << mm.first << "   max : " << mm.second << "   mean : " << CostFunctions::mean(mriT) << endl;
//    JointHisto hm(mriS,mriT);
//    //hm.print("H");
//    cout << " NMI : " << hm.computeNMI() << endl;
//    hm.smooth(7);
//    cout << " NMIs: " << hm.computeNMI() << endl;
// 

  exit(1);
}

void gradmag(Parameters & P)
{

  MRI* mri1 = MRIread(P.mov.c_str());
  MRI* mri_mag1 = MRIalloc(mri1->width, mri1->height, mri1->depth, MRI_FLOAT);
  MRI* mri_grad1 = MRIsobel(mri1, NULL, mri_mag1);
//  MRIwrite(mri_mag1,"mri1mag.mgz");
//  MRIwriteFrame(mri_grad1,"mri1sobel_grad1.mgz",0);
//  MRIwriteFrame(mri_grad1,"mri1sobel_grad2.mgz",1);
//  MRIwriteFrame(mri_grad1,"mri1sobel_grad3.mgz",2);

  int dd, hh, ww;
  float x, y, z, r, phi, psi;

//  MRI* mri_phi1  = MRIalloc(mri1->width, mri1->height, mri1->depth, MRI_FLOAT);
//  MRI* mri_psi1  = MRIalloc(mri1->width, mri1->height, mri1->depth, MRI_FLOAT);
  MRI* mri_pp1 = MRIalloc(mri1->width, mri1->height, mri1->depth, MRI_FLOAT);
  for (dd = 0; dd < mri_mag1->depth; dd++)
    for (hh = 0; hh < mri_mag1->height; hh++)
      for (ww = 0; ww < mri_mag1->width; ww++)
      {
        x = MRIgetVoxVal(mri_grad1, ww, hh, dd, 0);
        y = MRIgetVoxVal(mri_grad1, ww, hh, dd, 1);
        z = MRIgetVoxVal(mri_grad1, ww, hh, dd, 2);
        r = MRIgetVoxVal(mri_mag1, ww, hh, dd, 0);
        phi = 0;
        psi = 0;
        if (r > 10)
        {
          phi = acos(z / r);
          psi = atan2(y, x);
        }
//     MRIsetVoxVal(mri_phi1,ww,hh,dd,0,phi);
//     MRIsetVoxVal(mri_psi1,ww,hh,dd,0,psi);
        MRIsetVoxVal(mri_pp1, ww, hh, dd, 0, psi + phi);
      }
//  MRIwrite(mri_phi1,"mri1phi.mgz");
//  MRIwrite(mri_psi1,"mri1psi.mgz");
  MRIwrite(mri_pp1, "mri1phipsi.mgz");

  MRI* mri2 = MRIread(P.dst.c_str());
  MRI* mri_mag2 = MRIalloc(mri2->width, mri2->height, mri2->depth, MRI_FLOAT);
  MRI* mri_grad2 = MRIsobel(mri2, NULL, mri_mag2);
//  MRIwrite(mri_mag2,"mri2mag.mgz");
//  MRIwriteFrame(mri_grad2,"mri2sobel_grad1.mgz",0);
//  MRIwriteFrame(mri_grad2,"mri2sobel_grad2.mgz",1);
//  MRIwriteFrame(mri_grad2,"mri2sobel_grad3.mgz",2);
//  MRI* mri_phi2  = MRIalloc(mri2->width, mri2->height, mri2->depth, MRI_FLOAT);
//  MRI* mri_psi2  = MRIalloc(mri2->width, mri2->height, mri2->depth, MRI_FLOAT);
  MRI* mri_pp2 = MRIalloc(mri2->width, mri2->height, mri2->depth, MRI_FLOAT);
  for (dd = 0; dd < mri_mag2->depth; dd++)
    for (hh = 0; hh < mri_mag2->height; hh++)
      for (ww = 0; ww < mri_mag2->width; ww++)
      {
        x = MRIgetVoxVal(mri_grad2, ww, hh, dd, 0);
        y = MRIgetVoxVal(mri_grad2, ww, hh, dd, 1);
        z = MRIgetVoxVal(mri_grad2, ww, hh, dd, 2);
        r = MRIgetVoxVal(mri_mag2, ww, hh, dd, 0);
        phi = 0;
        psi = 0;
        if (r > 10)
        {
          phi = acos(z / r);
          psi = atan2(y, x);
        }
//     MRIsetVoxVal(mri_phi2,ww,hh,dd,0,phi);
//     MRIsetVoxVal(mri_psi2,ww,hh,dd,0,psi);
        MRIsetVoxVal(mri_pp2, ww, hh, dd, 0, psi + phi);
      }
//  MRIwrite(mri_phi2,"mri2phi.mgz");
//  MRIwrite(mri_psi2,"mri2psi.mgz");
  MRIwrite(mri_pp1, "mri1phipsi.mgz");

  MRIfree(&mri_grad1);
  MRIfree(&mri_grad2);

  exit(0);
}

int ReadByteImageRawData(float **Image, /* output image data */
long *Width, /* output image width */
long *Height /* output image height */
)

{ /* begin ReadByteImageRawData */

  //char  Filename[256];
  FILE *f = (FILE *) NULL;
  float *p;
  unsigned char *Line;
  long x, y;
  int Error;

  /* interactivity */
//   do {
//     printf("Give the input image file name (enter ? to cancel):\n");
//     printf("--> ");
//     scanf("%255s", Filename);
//     f = fopen(Filename, "rb");
//   } while (strcmp(Filename, "?") && (f == (FILE *)NULL));
//   if (!strcmp(Filename, "?")) {
//     if (f != (FILE *)NULL) {
//       fclose(f);
//       printf("Sorry: reserved file name\n");
//     }
//     printf("Cancel\n");
//     return(1);
//   }
//   do {
//     printf("Give the input image width:\n");
//     printf("--> ");
//     scanf("%ld", Width);
// 
//     printf("Give the input image height:\n");
//     printf("--> ");
//     scanf("%ld", Height);
//   } while ((*Width < 1L) || (*Height < 1L));
  f = fopen("lena.img", "rb");
  *Width = 256;
  *Height = 256;

  /* allocation of workspace */
  *Image = (float *) malloc((size_t) (*Width * *Height * (long) sizeof(float)));
  if (*Image == (float *) NULL)
  {
    fclose(f);
    printf("Allocation of input image failed\n");
    return (1);
  }
  Line = (unsigned char *) malloc(
      (size_t) (*Width * (long) sizeof(unsigned char)));
  if (Line == (unsigned char *) NULL)
  {
    free(*Image);
    *Image = (float *) NULL;
    fclose(f);
    printf("Allocation of buffer failed\n");
    return (1);
  }

  /* perform reading in raster fashion */
  p = *Image;
  for (y = 0L; y < *Height; y++)
  {
    Error = (*Width
        != (long) fread(Line, sizeof(unsigned char), (size_t) *Width, f));
    if (Error)
    {
      free(Line);
      free(*Image);
      *Image = (float *) NULL;
      fclose(f);
      printf("File access failed\n");
      return (1);
    }
    for (x = 0L; x < *Width; x++)
    {
      *p++ = (float) Line[x];
    }
  }

  free(Line);
  fclose(f);
  return (0);
} /* end ReadByteImageRawData */

void testcubic(Parameters & P)
{

  float *ImageRasterArray;
  long Width, Height;
  int Error;

  Error = ReadByteImageRawData(&ImageRasterArray, &Width, &Height);
  if (Error)
  {
    printf("Failure to import image data\n");
    exit(1);
  }

  MRI* img = MRIalloc(Width, Height, 1, MRI_FLOAT);
  int x, y;
  for (x = 0; x < Width; x++)
    for (y = 0; y < Height; y++)
      MRIsetVoxVal(img, x, y, 0, 0, ImageRasterArray[y * Height + x]);

  //MRIwrite(img,"coeff-orig.tif");
  //exit(1);
  MRI_BSPLINE * bspline = NULL;
  bspline = MRItoBSpline(img, NULL, 3);
  MRIwrite(bspline->coeff, "img-coeff.tif");

  double val;
  MRI* mri2 = MRIcopy(img, NULL);
  for (x = 0; x < Width; x++)
    for (y = 0; y < Height; y++)
    {
      MRIsampleBSpline(bspline, x, y, 0, 0, &val);
      MRIsetVoxVal(mri2, x, y, 0, 0, val);
    }

  MRIwrite(mri2, "img-id.tif");
  exit(1);

//   MRI* mri1 = MRIread(P.mov.c_str());
//   //LTA * lta = LTAread("identity.nofile");
//   
//   //MRI* mri2 = LTAtransformInterp(mri1,NULL,lta,SAMPLE_CUBIC_BSPLINE);
//  // MRIwrite(mri2,"idlta.mgz");
//  
//   MRI* mri2 = MRIcopy(mri1,NULL);
//   
//   // now manually do it:
//   MRI_BSPLINE * bspline = NULL;
//   bspline = MRItoBSpline(mri1,NULL,3);
//   MRIwrite(bspline->coeff,"coeff.mgz");
//   
// 
//   int r,c,s;
//   double val;
//   for (s=0;s<mri1->depth;s++)
//   for (c=0;c<mri1->height;c++)
//   for (r=0;r<mri1->width;r++)
//   {
//     MRIsampleBSpline(bspline, r, c, s, 0, &val);
//     MRIsetVoxVal(mri2, r, c,s, 0, val) ;
//   }
//   MRIwrite(mri2,"idmy.mgz");
//   exit(1);

}

int main(int argc, char *argv[])
{
  {
    // for valgrind, so that everything is freed
    cout << getVersion() << endl << endl;
//  setenv("SURFER_FRONTDOOR","",1) ;
    // set the environment variable
    // to store mri as chunk in memory:
//  setenv("FS_USE_MRI_CHUNK","",1) ;

    // Default initialization
    int nargs = handleVersionOption(argc, argv, "mri_robust_register");
    if (nargs && argc - nargs == 1)
    {
      exit(0);
    }
    argc -= nargs;
    Progname = argv[0];
    argc--;
    argv++;
    ErrorInit(NULL, NULL, NULL);
//  DiagInit(NULL, NULL, NULL) ;

    if (!parseCommandLine(argc, argv, P))
    {
      //printUsage();
      exit(1);
    }
// testSubsamp(P);
// testcubic(P);
//entro(P);
//gradmag(P);
//jointhisto(P);
//debug(P);
//lnorm(P);

// vnl_vector < double > p(6,0.0);
// p[3] = 0.1;
// p[4] = .2;
// p[5]=.3;
// vnl_matrix < double > M = (RegistrationStep<double>::convertP2Md(p,2)).first;
// vnl_matlab_print(vcl_cerr,M,"M",vnl_matlab_print_format_long);std::cerr << std::endl;
// exit(1);

    // Timer
    Timer start;
    int msec, minutes, seconds;
    start.reset();

    // init registration from Parameters
    Registration * Rp = NULL;
    if (P.cost == Registration::ROB) //|| P.cost == Registration::LS )
      Rp = new RegRobust;
    else
    {
      P.floattype = true; // bad way, I know: will allow type to switch inside Registration, and because of the constfun, will switch to uchar
      Rp = new RegPowell;
    }
    // keep as reference (in order not to modify everything below to pointer)
    Registration &R = *Rp;

//    Registration R;
    initRegistration(R, P);
//conv(P.mri_dst);

//  cout << " mean mov : " << CostFunctions::mean(P.mri_mov) << "  mean dst: " << CostFunctions::mean(P.mri_dst) << endl;
//  cout << " sdev mov : " << CostFunctions::sdev(P.mri_mov) << "  sdev dst: " << CostFunctions::sdev(P.mri_dst) << endl;
//  cout << " median mov : " << CostFunctions::median(P.mri_mov) << "  median dst: " << CostFunctions::median(P.mri_dst) << endl;
//  cout << " mad mov : " << CostFunctions::mad(P.mri_mov) << "  mad dst: " << CostFunctions::mad(P.mri_dst) << endl;
// does not work with different image dimensions:
//  cout << " LS difference before: " << CostFunctions::leastSquares(P.mri_mov,P.mri_dst) << endl;
//  cout << " NC difference before: " << CostFunctions::normalizedCorrelation(P.mri_mov,P.mri_dst) << endl;

    // compute Alignment
    //std::pair <MATRIX*, double> Md;
//  if (P.satest) R.computeSatEstimate(2,P.iterate,P.epsit);
////  else if (P.satit) Md = R.computeIterativeRegSat(P.iterate,P.epsit);
//  else
    if (P.satit && P.cost == Registration::ROB)
    {
      dynamic_cast<RegRobust*>(Rp)->findSaturation();
      R.computeMultiresRegistration(0, P.iterate, P.epsit);
    }
    else if (P.nomulti)
    {
      R.computeIterativeRegistration(P.iterate, P.epsit);
    }
    else
    {
      R.computeMultiresRegistration(0, P.iterate, P.epsit);
    }

//   if (P.satest) // old stuff, can be removed ?
//   {
//     cout << "run:" << endl;
//     cout << " gnuplot " << R.getName() << "-sat.plot ; \\ " << endl;
//     cout << " epstopdf " << R.getName() << "-sat.eps " << endl;
//     cout << " and view the pdf " << endl << endl;
//     msec = start.milliseconds() ;
//     seconds = nint((float)msec/1000.0f) ;
//     minutes = seconds / 60 ;
//     seconds = seconds % 60 ;
//     cout << "registration took "<<minutes<<" minutes and "<<seconds<<" seconds." << endl;
//     exit(0);
//   }

    //Md.first = MatrixReadTxt("xform.txt",NULL);
    //Md.second = 1;

    // Print results:
    std::pair<MATRIX*, double> Md;
    cout << endl << "Final Transform:" << endl;
    vnl_matrix<double> fMv2v(R.getFinalVox2Vox());
    vnl_matlab_print(vcl_cout,fMv2v,"M",vnl_matlab_print_format_long);
    cout << endl;
    Md.first = MyMatrix::convertVNL2MATRIX(fMv2v);
    Md.second = R.getFinalIscale();
//  MatrixPrintFmt(stdout,"% 2.8f",Md.first);

    cout << " Determinant : " << vnl_determinant(fMv2v) << endl << endl;
    if (P.affine || P.isoscale)
    {
      cout << " Decompose into Rot * Shear * Scale : " << endl << endl;
      vnl_matrix<double> Rot, Shear;
      vnl_diag_matrix<double> Scale;
      MyMatrix::Polar2Decomposition(fMv2v.extract(3, 3), Rot, Shear, Scale);
      vnl_matlab_print(vcl_cout,Rot,"Rot",vnl_matlab_print_format_long);
      cout << endl;
      vnl_matlab_print(vcl_cout,Shear,"Shear",vnl_matlab_print_format_long);
      cout << endl;
      vnl_matlab_print(vcl_cout,Scale,"Scale",vnl_matlab_print_format_long);
      cout << endl;

//    vnl_svd <double > svd(fMv2v.extract(3,3));
//    // svd: M = U * W * V'
//    //vnl_matlab_print(vcl_cout,svd.U(),"U");cout << endl;
//    //vnl_matlab_print(vcl_cout,svd.W(),"W");cout << endl;
//    //vnl_matlab_print(vcl_cout,svd.V(),"V");cout << endl;
//    // Polar: M = R*S = (U*V') * (V*W*V')
//    vnl_matrix < double > Rot(svd.U()*svd.V().transpose());
//    vnl_matrix < double > S(svd.V()*svd.W()*svd.V().transpose());
//    vnl_matlab_print(vcl_cout,Rot,"Rot");cout << endl;
//    //vnl_matlab_print(vcl_cout,S,"S");cout << endl;
//    // further decompose S into shear * diag(scales)
//    vnl_matrix < double > Shear(3,3);
//    vnl_diag_matrix < double > Scale(3);
//    for (unsigned int c=0;c<3;c++)
//    {
//      Scale[c] = S[c][c];
//      Shear.set_column(c,S.get_column(c) / Scale[c]);
//    }
//    vnl_matlab_print(vcl_cout,Shear,"Shear");cout << endl;
//    vnl_matlab_print(vcl_cout,Scale,"Scale");cout << endl;
//    //cout << " decompose error: " << (Rot*Shear*Scale - fMv2v.extract(3,3)).frobenius_norm() << endl;
    }

    if (R.isIscale())
    {
      cout << "Intenstiy Scale Factor: " << Md.second << endl;
    }
    cout << endl;

    // writing transform section here
    cout << "writing output transformation to " << P.lta << " ..." << endl;
    char reg[STRLEN];
    strcpy(reg, P.lta.c_str());
    LTA * lta = LTAalloc(1, P.mri_mov);
    if (!P.lta_vox2vox) // do ras to ras
    {
      cout << "converting VOX to RAS and saving RAS2RAS..." << endl;
      //cout << "VOX2VOX:" << endl ;
      //MatrixPrint(stdout, Md.first) ;
      lta->xforms[0].m_L = MRIvoxelXformToRasXform(P.mri_mov, P.mri_dst,
          Md.first, lta->xforms[0].m_L);
      //cout << "RAS2RAS:" << endl ;
      //MatrixPrint(stdout,lta->xforms[0].m_L) ;
      lta->type = LINEAR_RAS_TO_RAS;
    }
    else // vox to vox
    {
      cout << "saving VOX2VOX..." << endl;
      lta->xforms[0].m_L = MatrixCopy(Md.first, lta->xforms[0].m_L);
      lta->type = LINEAR_VOX_TO_VOX;
    }
    // add src and dst info
    getVolGeom(P.mri_mov, &lta->xforms[0].src);
    getVolGeom(P.mri_dst, &lta->xforms[0].dst);
    LTAwriteEx(lta, reg);

    if (R.isIscale() && Md.second > 0 && P.iscaleout != "")
//  if (R.isIscale() && Md.second >0)
    {
      //string fn;
      //if (P.iscaleout != "") fn = P.iscaleout;
      //else fn = R.getName() + "-intensity.txt";
      //ofstream f(fn.c_str(),ios::out);
      ofstream f(P.iscaleout.c_str(), ios::out);
      f << Md.second;
      f.close();
    }

    //  MatrixWriteTxt("xform.txt",Md.first);
    // end of writing transform

    // here do scaling of intensity values
    if (R.isIscale() && Md.second > 0 && !P.entropy)
    {
      cout << "Adjusting Intensity of MOV by " << Md.second << endl;
      P.mri_mov = MyMRI::MRIvalscale(P.mri_mov, P.mri_mov, Md.second);
    }

    // maybe map source to target (resample):
    if (P.warpout != "")
    {
      //cout << "using lta" << endl;
      cout << endl;
      cout << "mapmov: resampling MOV to DST ..." << endl;
      int nframes = P.mri_mov->nframes;
      if (P.mri_mov->nframes > 1)
      {
        cout
            << " WARNING: movable has more than one frame !!! Only map first ..."
            << endl;
      }
      P.mri_mov->nframes = 1; // only map frame 1

//       //switch type to dst type:
//       MRI *mri_aligned = MRIcloneDifferentType(P.mri_dst, MRI_FLOAT) ;
//       mri_aligned = LTAtransform(P.mri_mov,mri_aligned, lta);
//       if (mri_aligned->type != P.mri_dst->type)
//       {
//          int no_scale_flag = FALSE;
//          MRI* mri2 = MRISeqchangeType(mri_aligned, P.mri_dst->type, 0.0, 0.999, no_scale_flag);
//          if (mri2 == NULL)
//          {
//            printf("ERROR: MRISeqchangeType\n");
//            exit(1);
//          }
//          MRIfree(&mri_aligned);
//          mri_aligned = mri2;
//       }

      // keep mov type:
      cout << " copy ... " << endl;
      MRI *mri_aligned = MRIcloneDifferentType(P.mri_dst, P.mri_mov->type);
      cout << " transform ... " << endl;
      mri_aligned = LTAtransformInterp(P.mri_mov, mri_aligned, lta,
          P.finalsampletype);

      // reset mov n frames:
      P.mri_mov->nframes = nframes;

      // keep acquisition params:
      MRIcopyPulseParameters(P.mri_mov, mri_aligned);

//    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
//    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
//    sprintf(fname, "%s_target", parms.base_name) ;
//    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;

//      cout << " mean warp : " << CostFunctions::mean(mri_aligned) << "  mean dst: " << CostFunctions::mean(P.mri_dst) << endl;
//      cout << " sdev warp : " << CostFunctions::sdev(mri_aligned) << "  sdev dst: " << CostFunctions::sdev(P.mri_dst) << endl;
//      cout << " median warp : " << CostFunctions::median(mri_aligned) << "  median dst: " << CostFunctions::median(P.mri_dst) << endl;
//      cout << " mad warp : " << CostFunctions::mad(mri_aligned) << "  mad dst: " << CostFunctions::mad(P.mri_dst) << endl;
// does not work with different image dimensions:
//    cout << " LS difference after: " << CostFunctions::leastSquares(mri_aligned,P.mri_dst) << endl;
//    cout << " NC difference after: " << CostFunctions::normalizedCorrelation(mri_aligned,P.mri_dst) << endl;

      MRIwrite(mri_aligned, P.warpout.c_str());
      MRIfree(&mri_aligned);

      cout << endl;
      cout << "To check aligned result, run:" << endl;
      cout << "  freeview -v " << P.dst << " " << P.warpout << endl;
    }

    if (P.norlout != "") // map source to target (no resample only adjust header)
    {
      cout << endl;
      cout << "mapmovhdr: Changing vox2ras MOV header (to map to DST) ..."
          << endl;
      // Compute new vox2ras instead of resampling
      // vox2ras = Stemp * invTtemp * invR * Tin
      MATRIX * ras2ras = MRIvoxelXformToRasXform(P.mri_mov, P.mri_dst, Md.first,
          NULL);
      MATRIX * vox2ras = MRIgetVoxelToRasXform(P.mri_mov);

      // concat:
      vox2ras = MatrixMultiply(ras2ras, vox2ras, vox2ras);

      MRI *mri_aligned = MRIcopy(P.mri_mov, NULL);
      MRIsetVoxelToRasXform(mri_aligned, vox2ras);
      int err = MRIwrite(mri_aligned, P.norlout.c_str());
      MRIfree(&mri_aligned);
      if (err)
      {
        printf("ERROR: writing %s\n", P.norlout.c_str());
        exit(1);
      }

      cout << endl;
      cout << "To check aligned result, run:" << endl;
      cout << "  freeview -v " << P.dst << " " << P.norlout << endl;

    }

    // maybe write out weights in target space:
    if (P.weightsout != "")
    {

      MRI * mri_weights = R.getWeights(); // in target space
      if (mri_weights != NULL)
      {
        cout << endl;
        cout << "Writing out Weights ..." << endl;
        //MRIwrite(mri_weights,"temp.mgz") ;
        //cout << " mri_weights type: " << mri_weights->type << endl;
        if (P.oneminusweights)
        {
          MRI * mri_iweights = MRIlinearScale(mri_weights, NULL, -1, 1, 0);
          MRIwrite(mri_iweights, P.weightsout.c_str());
          MRIfree(&mri_iweights);
        }
        else
          MRIwrite(mri_weights, P.weightsout.c_str());

//       // map to target and use target geometry
//       std::pair < vnl_matrix_fixed < double, 4, 4>, vnl_matrix_fixed < double, 4, 4> > map2weights = R.getHalfWayMaps();
//       vnl_matrix_fixed < double, 4, 4> hinv = vnl_inverse(map2weights.second);
//       MRI * wtarg = MRIalloc(P.mri_dst->width,P.mri_dst->height,P.mri_dst->depth,MRI_FLOAT);
//       MRIcopyHeader(P.mri_dst,wtarg);
//       wtarg->type = MRI_FLOAT;
//       wtarg = MyMRI::MRIlinearTransform(mri_weights,wtarg, hinv);
//       MRIwrite(wtarg, P.weightsout.c_str()) ;
//       MRIfree(&wtarg);
//       //MatrixFree(&hinv);
        cout << "... overlay the weights:" << endl;
        if (P.warpout != "")
          cout << "  freeview -v " << P.dst << " " << P.warpout << " "
              << P.weightsout << ":colormap=heat" << endl;
        else if (P.norlout != "")
          cout << "  freeview -v " << P.dst << " " << P.norlout << " "
              << P.weightsout << ":colormap=heat" << endl;
//        cout << "  tkmedit -f "<< P.dst <<" -aux "<< P.warpout << " -overlay " << P.weightsout <<endl;
      }
      else
      {
        cout
            << "Warning: no weights could be created! Maybe you ran with --leastsquares??"
            << endl;
      }
    }

    // write out images in half way space
    if (P.halfmov != "" || P.halfdst != "" || P.halfweights != ""
        || P.halfdstlta != "" || P.halfmovlta != "")
    {
      cout << endl;

      if (!P.symmetry)
      {
        cout << "ERROR: no half way space created (symmetry was switched off)!"
            << endl;
      }
      else if (!R.getHalfWayGeom())
      {
        cout
            << "ERROR: no half way space created (not implemented for this cost function)!"
            << endl;
      }
      else
      {

        cout << "Creating half way data ..." << endl;
        std::pair<vnl_matrix_fixed<double, 4, 4>, vnl_matrix_fixed<double, 4, 4> > maps2weights =
            R.getHalfWayMaps();
        vnl_matlab_print(vcl_cerr,maps2weights.first,"movhw",vnl_matlab_print_format_long);
        std::cerr << std::endl;
        vnl_matlab_print(vcl_cerr,maps2weights.second,"dsthw",vnl_matlab_print_format_long);
        std::cerr << std::endl;

        MRI * mri_hwgeom = R.getHalfWayGeom();

        LTA * m2hwlta = LTAalloc(1, P.mri_mov);
        LTA * d2hwlta = LTAalloc(1, P.mri_dst);
        if (!P.lta_vox2vox) // do ras to ras
        {
          // cout << "converting VOX to RAS and saving RAS2RAS..." << endl ;
          // (use geometry of destination space for half-way)
          m2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform(P.mri_mov,
              mri_hwgeom, MyMatrix::convertVNL2MATRIX(maps2weights.first.as_matrix()),
              m2hwlta->xforms[0].m_L);
          m2hwlta->type = LINEAR_RAS_TO_RAS;
          d2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform(P.mri_dst,
              mri_hwgeom, MyMatrix::convertVNL2MATRIX(maps2weights.second.as_matrix()),
              d2hwlta->xforms[0].m_L);
          d2hwlta->type = LINEAR_RAS_TO_RAS;
        }
        else // vox to vox
        {
          // cout << "saving VOX2VOX..." << endl ;
          //m2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.first, m2hwlta->xforms[0].m_L) ;
          m2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(
              maps2weights.first.as_matrix(), m2hwlta->xforms[0].m_L);
          m2hwlta->type = LINEAR_VOX_TO_VOX;
          //d2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.second, d2hwlta->xforms[0].m_L) ;
          d2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(
              maps2weights.second.as_matrix(), d2hwlta->xforms[0].m_L);
          d2hwlta->type = LINEAR_VOX_TO_VOX;
        }
        // add src and dst info (use mri_weights as target geometry in both cases)
        getVolGeom(P.mri_mov, &m2hwlta->xforms[0].src);
        getVolGeom(mri_hwgeom, &m2hwlta->xforms[0].dst);
        getVolGeom(P.mri_dst, &d2hwlta->xforms[0].src);
        getVolGeom(mri_hwgeom, &d2hwlta->xforms[0].dst);

        // write lta to half way
        if (P.halfmovlta != "")
        {
          LTAwriteEx(m2hwlta, P.halfmovlta.c_str());
        }
        if (P.halfdstlta != "")
        {
          LTAwriteEx(d2hwlta, P.halfdstlta.c_str());
        }

        if (P.halfmov != "")
        {
          cout << " creating half-way movable ..." << endl;
          // take dst geometry info from lta:
          MRI* mri_Swarp = LTAtransformInterp(P.mri_mov, NULL, m2hwlta,
              P.finalsampletype);

          //cout << " MOV       RAS: " << P.mri_mov->c_r << " , " <<  P.mri_mov->c_a << " , " <<  P.mri_mov->c_s << endl;
          //cout << " DST       RAS: " << P.mri_dst->c_r << " , " <<  P.mri_dst->c_a << " , " <<  P.mri_dst->c_s << endl;
          //cout << " weights   RAS: " << mri_weights->c_r << " , " <<  mri_weights->c_a << " , " <<  mri_weights->c_s << endl;
          //cout << " Swarp_old RAS: " << mri_Swarp_old->c_r << " , " <<  mri_Swarp_old->c_a << " , " <<  mri_Swarp_old->c_s << endl;
          //MRI* mri_Swarp = MRIalloc(mri_weights->width, mri_weights->height, mri_weights->depth, P.mri_mov->type);
          //MRIcopyHeader(mri_weights,mri_Swarp);
          //mri_Swarp->type = P.mri_mov->type;
          //LTAtransform(P.mri_mov,mri_Swarp, m2hwlta);
          //cout << " Swarp     RAS: " << mri_Swarp->c_r << " , " <<  mri_Swarp->c_a << " , " <<  mri_Swarp->c_s << endl;
          MRIcopyPulseParameters(P.mri_mov, mri_Swarp);
          MRIwrite(mri_Swarp, P.halfmov.c_str());

          if (P.debug)
          {
            MRIiterator mw(R.getWeights());
            MRIiterator ms(mri_Swarp);
            double meanw1 = 0, meanw0 = 0, mean = 0, meanw = 0, countw = 0;
            int countw1 = 0, countw0 = 0, count = 0;
            for (ms.begin(); !ms.isEnd(); ms++)
            {
              if (fabs(*mw) > 0.0001)
              {
                meanw0 += (*ms);
                countw0++;
              }
              if (fabs(*mw - 1.0) < 0.0001)
              {
                meanw1 += *ms;
                countw1++;
              }

              mean += *ms;
              count++;

              meanw += *ms * *mw;
              countw += *mw;

              assert(! (mw.isEnd() && !ms.isEnd()));
              mw++;
            }

            cout << " mov int means: " << mean / count << " ( " << count
                << " )  w0: " << meanw0 / countw0 << " ( " << countw0
                << " ) w1: " << meanw1 / countw1 << " ( " << countw1
                << " )  weighted: " << meanw / countw << " ( " << countw << " )"
                << endl;
          }

          MRIfree(&mri_Swarp);

          //MRIwrite(P.mri_mov,"movable-original.mgz");
          //mri_Swarp = R.makeConform(P.mri_mov,NULL,false,true);
          //MRIwrite(mri_Swarp,"movable-uhar.mgz");
          //MRI * tttemp = MRIclone(mri_Swarp,NULL);
          //tttemp =  MRIlinearTransform(mri_Swarp,tttemp, mh);
          //MRIwrite(tttemp,"movable-uhar-half.mgz");
          //MRIfree(&mri_Swarp);
          //MRIfree(&tttemp);

        }
        if (P.halfdst != "")
        {
          cout << " creating half-way destination ..." << endl;
          MRI* mri_Twarp = LTAtransformInterp(P.mri_dst, NULL, d2hwlta,
              P.finalsampletype);
          MRIcopyPulseParameters(P.mri_dst, mri_Twarp);
          MRIwrite(mri_Twarp, P.halfdst.c_str());

          if (P.debug)
          {
            MRIiterator mw(R.getWeights());
            MRIiterator ms(mri_Twarp);
            double meanw1 = 0, meanw0 = 0, mean = 0, meanw = 0, countw = 0;
            int countw1 = 0, countw0 = 0, count = 0;
            for (ms.begin(); !ms.isEnd(); ms++)
            {
              if (fabs(*mw) > 0.0001)
              {
                meanw0 += (*ms);
                countw0++;
              }
              if (fabs(*mw - 1.0) < 0.0001)
              {
                meanw1 += *ms;
                countw1++;
              }

              mean += *ms;
              count++;

              meanw += *ms * *mw;
              countw += *mw;

              assert(! (mw.isEnd() && !ms.isEnd()));
              mw++;
            }
            cout << " mov int means: " << mean / count << " ( " << count
                << " )  w0: " << meanw0 / countw0 << " ( " << countw0
                << " ) w1: " << meanw1 / countw1 << " ( " << countw1
                << " )  weighted: " << meanw / countw << " ( " << countw << " )"
                << endl;
          }

          MRIfree(&mri_Twarp);
        }
        if (P.halfweights != "")
        {
          MRI * mri_weights = R.getWeights();
          if (mri_weights != NULL)
          {
            cout << " saving half-way weights ..." << endl;
            MRI* mri_wtemp = LTAtransformInterp(mri_weights, NULL, d2hwlta,
                P.finalsampletype);
            if (P.oneminusweights)
            {
              mri_wtemp = MRIlinearScale(mri_wtemp, mri_wtemp, -1, 1, 0);
            }
            MRIwrite(mri_wtemp, P.halfweights.c_str());
            MRIfree(&mri_wtemp);
            //MRIwrite(mri_weights,P.halfweights.c_str());
          }
          else
          {
            cout
                << "Warning: no weights have been computed! Maybe you ran with --leastsquares??"
                << endl;
          }
        }
      }
    }

    if (P.debug > 0)
    {
      cout << endl;
      cout << "To check debug output, run:" << endl;
      std::string name = R.getName();
      cout << "  tkmedit -f " << name << "-mriS-mapped.mgz -aux " << name
          << "-mriT-mapped.mgz -overlay " << name << "-mriS-weights.mgz"
          << endl;
    }

    cout << endl;
    //  cout << "To check transform, run:" << endl;
    //  cout << "  tkregister2 --mov "<< P.mov <<" --targ " << P.dst <<" --lta " << P.lta << " --reg " << R.getName() << ".reg" << endl;

    // cleanup
    if (Md.first)
    {
      MatrixFree(&Md.first);
    }
    if (P.mri_mov)
    {
      MRIfree(&P.mri_mov);
    }
    if (P.mri_dst)
    {
      MRIfree(&P.mri_dst);
    }
    if (Rp)
    {
      delete (Rp);
      Rp = NULL;
    }
    if (lta)
      LTAfree(&lta);

    ///////////////////////////////////////////////////////////////
    msec = start.milliseconds();
    seconds = nint((float) msec / 1000.0f);
    minutes = seconds / 60;
    seconds = seconds % 60;
    cout << endl << "Registration took " << minutes << " minutes and "
        << seconds << " seconds." << endl;

    cout << endl << " Thank you for using RobustRegister! " << endl;
    cout << " If you find it useful and use it for a publication, please cite: "
        << endl << endl;
    cout
        << " Highly Accurate Inverse Consistent Registration: A Robust Approach"
        << endl;
    cout
        << " M. Reuter, H.D. Rosas, B. Fischl.  NeuroImage 53(4):1181-1196, 2010."
        << endl;
    cout << " http://dx.doi.org/10.1016/j.neuroimage.2010.07.020" << endl;
    cout << " http://reuter.mit.edu/papers/reuter-robreg10.pdf" << endl << endl;
    
    if (P.cost == Registration::SB) 
    {
      cout << " And for the segmentation based cost function:" << endl;
      cout << " Multimodal image registration through simultaneous segmentation," << endl;
      cout << " I. Aganj and B. Fischl. IEEE Signal Processing Letters 24(11):1661-1665, 2017." << endl;
      cout << " https://doi.org/10.1109/LSP.2017.2754263" << endl << endl;
    }

    //if (diag_fp) fclose(diag_fp) ;
  } // for valgrind, so that everything is free
  exit(0);
  return (0);
}

// int main(int argc, char *argv[])
// {
//
//   char         *fname_src, *fname_dst, *fname_out, fname[STRLEN];
//   MRI          *mri_src, *mri_dst, *mri_tmp;
//
//   int          nargs,ninputs,i,msec,minutes,seconds;
//   Timer start ;
//
//   // defaults
//   Progname = argv[0] ;
//
//   Registration R; // sets its own default parameters
//
// //  int ac = argc ;
// //  char **av = argv ;
//   for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
//   {
//     nargs = get_option(argc, argv,R) ;
//     argc -= nargs ;
//     argv += nargs ;
//   }
//
//   if (argc < 4)
//   {
//     printUsage();
//     exit(1);
//   }
//
//   ninputs = argc-3 ;
//   //cout << "reading "<<ninputs<<" input volumes..."<< endl;
//   fname_dst = argv[ninputs+1] ;
//   fname_out = argv[ninputs+2] ;
//   FileNameOnly(fname_out, fname) ;
//   FileNameRemoveExtension(fname, fname) ;
//   strcpy(parms.base_name, fname) ;
//  // Gdiag |= DIAG_WRITE ;
//  // cout << "logging results to "<< parms.base_name <<".log" << endl;
//
//   start.reset() ;
//   ///////////  read MRI Target //////////////////////////////////////////////////
//   cout << endl << "reading target '"<<fname_dst<<"'..."<< endl;;
//   fflush(stdout) ;
//   mri_dst = MRIread(fname_dst) ;
//   if (mri_dst == NULL)
//     ErrorExit(ERROR_NOFILE, "%s: could not open MRI Target %s.\n",
//               Progname, fname_dst) ;
//
//    //////////////////////////////////////////////////////////////
//   // create a list of MRI volumes
//   cout << "reading "<<ninputs<<" source (movable) volumes..."<< endl;
//   for (i = 0 ; i < ninputs ; i++) {
//     fname_src = argv[i+1] ;
//     cout << "reading source '"<<fname_src<<"'..." << endl;
//     fflush(stdout) ;
//     mri_tmp = MRIread(fname_src) ;
//     if (!mri_tmp)
//       ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
//                 Progname, fname_src) ;
//
// //    TRs[i] = mri_tmp->tr ;
// //    fas[i] = mri_tmp->flip_angle ;
// //    TEs[i] = mri_tmp->te ;
//
// //    if (mask_fname) {
// //      MRI *mri_mask ;
// //
// //      mri_mask = MRIread(mask_fname) ;
// //      if (!mri_mask)
// //        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
// //                  Progname, mask_fname) ;
// //      MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
// //      MRIfree(&mri_mask) ;
// //    }
//     if (i == 0)
//     {
//       mri_src = MRIallocSequence(mri_tmp->width,
//                                 mri_tmp->height,
//                                 mri_tmp->depth,
//                                 mri_tmp->type,
//                                 ninputs) ;
//       MRIcopyHeader(mri_tmp, mri_src) ;
//     }
//     MRIcopyFrame(mri_tmp, mri_src, 0, i) ;
//     MRIfree(&mri_tmp) ;
//   }
//   cout << endl;
//   //////////////////////////////////////////////////////////////
//
// //  if (!FZERO(blur_sigma)) {
// //    MRI *mri_tmp, *mri_kernel ;
//
// //    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
// //    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
// //    MRIfree(&mri_in) ;
// //    mri_in = mri_tmp ;
// // }
//
//     ////////////////////////////////////////////////////
//     // now start working (remember this is vox-to-vox transform)
//   //  parms.lta->xforms[0].m_L = MatrixIdentity(4, NULL) ;
//
//   // real work done here
// //  MATRIX * Minit = MRIgetVoxelToVoxelXform(mri_src,mri_dst) ;
// //  cout << "initial transform:\n" ;
// //  MatrixPrintFmt(stdout,"% 2.8f",Minit);
// //  //std::pair <MATRIX*, double> Md = R.computeIterativeRegistration(10,mri_src,mri_dst,Minit);
// //  std::pair <MATRIX*, double> Md = R.computeMultiresRegistration(mri_src,mri_dst,Minit);
// //  std::pair <MATRIX*, double> Md = R.computeMultiresRegistration(mri_src,mri_dst);
//
//   std::pair <MATRIX*, double> Md;
//   if (nit > 0)  Md = R.computeIterativeRegistration(nit,mri_src,mri_dst);
//   else Md = R.computeMultiresRegistration(mri_src,mri_dst);
//
//
//   cout << "final transform:\n" ;
//   MatrixPrintFmt(stdout,"% 2.8f",Md.first);
//   if (R.isIscale()) cout << " iscale: " << Md.second << endl;
//   cout << endl ;
//
//   /////////////////////////diagnostics/////////////////////////////////
// //  if ((Gdiag & DIAG_WRITE) && (parms.write_iterations != 0)) {
//     MRI *mri_aligned ;
//     int nframes = mri_src->nframes;
//     mri_src->nframes = 1 ;
//     mri_aligned = MRIlinearTransform(mri_src, NULL,Md.first);
//     // here also do scaling of intensity values
//     mri_src->nframes = nframes ;
// //    sprintf(fname, "%s_after_final_alignment", parms.base_name) ;
// //    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
// //    sprintf(fname, "%s_target", parms.base_name) ;
// //    MRIwriteImageViews(mri_dst, fname, IMAGE_SIZE) ;
//     sprintf(fname, "%s.mgz", parms.base_name) ;
//     MRIwrite(mri_aligned, fname) ;
//     MRIfree(&mri_aligned) ;
//     cout << "To check results, run:" << endl;
//     cout << "  tkmedit -f "<< fname_dst <<" -aux " << fname << endl;
// //    cout << " or " << endl;
// //    cout << "  tkmedit -f "<< fname_src <<" -aux " << fname << endl;
//
// //  }
//
//    /////////////////////////////////////////////////////////////////////
//   cout << "writing output transformation to "<<fname_out<<"..." << endl;
//   // writing transform section here
//   // create gca volume for outputting dirction cosines and c_(ras)
//   //mri_dst = MRIallocHeader(gca->width, gca->height, gca->depth, mri_in->type);
//   //GCAcopyDCToMRI(gca, mri_dst);
//   //strcpy(mri_dst->fname,gca_fname); // copy gca name
//   LTA * lta = LTAalloc(1,mri_src);
//
//   if (!stricmp(fname_out+strlen(fname_out)-3, "XFM"))
//   {
//     cout << "converting xform to RAS..." << endl ;
//     cout << "initial:" << endl ;
//
//     MatrixPrint(stdout, Md.first) ;
//     lta->xforms[0].m_L = MRIvoxelXformToRasXform (mri_src, mri_dst, Md.first, lta->xforms[0].m_L) ;
//     cout << "final:" << endl ;
//     MatrixPrint(stdout,lta->xforms[0].m_L) ;
//     lta->type = LINEAR_RAS_TO_RAS ;
//   }
//   else
//   {
//     lta->xforms[0].m_L = MatrixCopy(Md.first, lta->xforms[0].m_L) ;
//     lta->type = LINEAR_VOX_TO_VOX ;
//   }
//
//   // add src and dst info
//     getVolGeom(mri_src, &lta->xforms[0].src);
//     getVolGeom(mri_dst, &lta->xforms[0].dst);
//     LTAwriteEx(lta, fname_out) ;
//
//   ///////////////////////////////////////////// end of writing transform
//
//    if (mri_src)
//     MRIfree(&mri_src) ;
//    if (mri_dst)
//     MRIfree(&mri_dst) ;
//    if (Md.first)
//     MatrixFree(&Md.first) ;
//
//
//   ///////////////////////////////////////////////////////////////
//   msec = start.milliseconds() ;
//   seconds = nint((float)msec/1000.0f) ;
//   minutes = seconds / 60 ;
//   seconds = seconds % 60 ;
//   cout << "registration took "<<minutes<<" minutes and "<<seconds<<" seconds." << endl;
//   if (diag_fp)
//     fclose(diag_fp) ;
//   exit(0) ;
//   return(0) ;
//
// }

/*----------------------------------------------------------------------
 ----------------------------------------------------------------------*/
#include "mri_robust_register.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(mri_robust_register_help_xml, mri_robust_register_help_xml_len);
}

/*!
 \fn void initRegistration(Registration & R, const Parameters & P)
 \brief Initializes a Registration with Parameters (affine, iscale, transonly, leastsquares, sat and trans)
 \param R  Registration to be initialized
 \param P  Paramters for the initialization
 */
static void initRegistration(Registration & R, Parameters & P)
{
  if (!P.affine)
    R.setRigid();
  else
    R.setAffine();
  if (P.transonly)
    R.setTransonly();
  if (P.isoscale)
    R.setIsoscale();
  R.setIscale(P.iscale);
  if (P.iscaleonly)
    R.setIscaleOnly();
    
  //R.setRobust(!P.leastsquares);
  //R.setSaturation(P.sat);
  R.setVerbose(P.verbose); // set before debug, as debug sets its own verbose level
  R.setDebug(P.debug);
  R.setHighit(P.highit);
  R.setInitTransform(P.inittrans);
  R.setInitOrient(P.initorient);
  R.setInitScaling(P.initscaling);
  R.setDoublePrec(P.doubleprec);
  //R.setWLimit(P.wlimit);
  R.setSymmetry(P.symmetry);
  R.setCost(P.cost);
  //R.setOutputWeights(P.weights,P.weightsout);
  // set only for robust registration
  //if (P.cost == Registration::ROB)
  if (R.getClassName() == "RegRobust")
  {
    dynamic_cast<RegRobust*>(&R)->setSaturation(P.sat);
    dynamic_cast<RegRobust*>(&R)->setWLimit(P.wlimit);
  }
  if (R.getClassName() == "RegPowell")
  {
    dynamic_cast<RegPowell*>(&R)->setTolerance(P.powelltol);
  }

  int pos = P.lta.rfind(".");
  if (pos > 0)
  {
    R.setName(P.lta.substr(0, pos));
  }
  else
  {
    R.setName(P.lta);
  }

  if (P.subsamplesize > 0)
  {
    R.setSubsampleSize(P.subsamplesize);
  }
  R.setMinSize(P.minsize);
  R.setMaxSize(P.maxsize);

//   //////////////////////////////////////////////////////////////
//   // create a list of MRI volumes
//   //cout << "reading "<<ninputs<<" source (movable) volumes..."<< endl;
//   int ninputs = 1; // later we might want to load several frames
//   MRI* mri_src = NULL;
//   MRI* mri_tmp = NULL;
//   if (P.mri_mov) MRIfree(&P.mri_mov);
//   for (int i = 0 ; i < ninputs ; i++)
//   {
//     cout << "reading source '"<<P.mov<<"'..." << endl;
//     fflush(stdout) ;
//
//     mri_tmp = MRIread(P.mov.c_str()) ;
//     if (!mri_tmp)
//     {
//       ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
//                 Progname, P.mov.c_str()) ;
//       //cerr << Progname << " could not open input volume " << P.mov << endl;
//       //exit(1);
//     }
//
//     if (i == 0)
//     {
//       mri_src = MRIallocSequence(mri_tmp->width,
//                                  mri_tmp->height,
//                                  mri_tmp->depth,
//                                  mri_tmp->type,
//                                  ninputs) ;
//       MRIcopyHeader(mri_tmp, mri_src) ;
//       P.mri_mov = MRIallocSequence(mri_tmp->width,
//                                    mri_tmp->height,
//                                    mri_tmp->depth,
//                                    mri_tmp->type,
//                                    ninputs) ;
//       MRIcopyHeader(mri_tmp, P.mri_mov) ;
//     }
//     MRIcopyFrame(mri_tmp, P.mri_mov, 0, i) ; // store input in P.mri_mov
//
//     if (P.maskmov != "") // work only on mri_src to init registration (not P.mri_mov)
//     {
//       MRI *mri_mask = MRIread(P.maskmov.c_str());
//       if (!mri_mask)
//         ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
//                   Progname, P.maskmov.c_str()) ;
//       MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
//       MRIfree(&mri_mask) ;
//     }
//
//     MRIcopyFrame(mri_tmp, mri_src, 0, i) ;
//     MRIfree(&mri_tmp) ;
//   }
////   R.setSource(mri_src,P.fixvoxel,P.fixtype);
////   MRIfree(&mri_src);

  ///////////  read MRI Source //////////////////////////////////////////////////
  cout << endl;
  cout << "reading source '" << P.mov << "'..." << endl;
  fflush(stdout);

  MRI* mri_mov = MRIread(P.mov.c_str());
  if (mri_mov == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open MRI source %s.\n", Progname,
        P.mov.c_str());
    //cerr << Progname << " could not open MRI Target " << P.mov << endl;
    //exit(1);
  }
  if (mri_mov->nframes != 1)
  {
    ErrorExit(ERROR_NOFILE, "%s: only pass single frame MRI source %s.\n",
        Progname, P.mov.c_str());
  }

  P.mri_mov = MRIcopy(mri_mov, P.mri_mov); // save dst mri
  
  // Load Mask
  MRI *mri_mask = NULL;
  if (P.maskmov != "")
  {
    cout << "reading source mask '" << P.maskmov << "'..." << endl;
    mri_mask = MRIread(P.maskmov.c_str());
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n", Progname,
          P.maskmov.c_str());
  }
  
  if (P.sobel)
  {
    cout << "Converting mov to sobel magnitute image ..." << endl;
    MRI * temp = mri_mov;
    mri_mov = MRIalloc(temp->width, temp->height, temp->depth, MRI_FLOAT);
    mri_mov = MRIcopyHeader(temp, mri_mov);
    mri_mov->type = MRI_FLOAT;
    MRIsobel(temp, NULL, mri_mov);
    //MRIwrite(mri_mov,"sobel-mov.mgz");
    MRIfree(&temp);
  }
  
  if (P.entropy)
  {
    MRI * temp = mri_mov;
    Timer start;
    int msec, seconds;
    start.reset();
    cout << "Converting mov to entropy image (radius " << P.entroradius
        << " ) ... (can take 1-2 min)" << endl;
    mri_mov = MyMRI::entropyImage(temp, P.entroradius, P.entball,
        P.entcorrection,mri_mask);
    if (P.entmov != "")
      MRIwrite(mri_mov, P.entmov.c_str());
    msec = start.milliseconds();
    seconds = nint((float) msec / 1000.0f);
//    minutes = seconds / 60;
    //seconds = seconds % 60 ;
    cout << " Entropy computation took " << seconds << " seconds." << endl;
    MRIfree(&temp);
  }
  // apply mask
  else if (mri_mask)
  {
    MRImask(mri_mov, mri_mask, mri_mov, 0, 0);
   }
   
  if (mri_mask)
  {
    MRIfree(&mri_mask);
  }
  ///////////  read MRI Target //////////////////////////////////////////////////
  cout << "reading target '" << P.dst << "'..." << endl;
  fflush(stdout);

  MRI* mri_dst = MRIread(P.dst.c_str());
  if (mri_dst == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open MRI target %s.\n", Progname,
        P.dst.c_str());
    //cerr << Progname << " could not open MRI Target " << P.dst << endl;
    //exit(1);
  }
  if (mri_dst->nframes != 1)
  {
    ErrorExit(ERROR_NOFILE, "%s: only pass single frame MRI target %s.\n",
        Progname, P.dst.c_str());
  }
  P.mri_dst = MRIcopy(mri_dst, P.mri_dst); // save dst mri
  
  if (P.maskdst != "")
  {
    cout << "reading target mask '" << P.maskdst << "'..." << endl;
    mri_mask = MRIread(P.maskdst.c_str());
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n", Progname,
          P.maskdst.c_str());
  }
  
  if (P.sobel)
  {
    cout << "Converting dst to sobel magnitute image ..." << endl;
    MRI * temp = mri_dst;
    mri_dst = MRIalloc(temp->width, temp->height, temp->depth, MRI_FLOAT);
    mri_dst = MRIcopyHeader(temp, mri_dst);
    mri_dst->type = MRI_FLOAT;
    MRIsobel(temp, NULL, mri_dst);
    //MRIwrite(mri_dst,"sobel-dst.mgz");
    MRIfree(&temp);
  }
  
  if (P.entropy)
  {
    MRI * temp = mri_dst;
    cout << "Converting dst to entropy image (radius " << P.entroradius
        << " ) ... (can take 1-2 min)" << endl;
    mri_dst = MyMRI::entropyImage(temp, P.entroradius, P.entball,
        P.entcorrection,mri_mask);
    if (P.entdst != "")
      MRIwrite(mri_dst, P.entdst.c_str());
    MRIfree(&temp);
  }
  else if (mri_mask)
  {
    MRImask(mri_dst, mri_mask, mri_dst, 0, 0);
  }
  if (mri_mask)
  {
    MRIfree(&mri_mask);
  }

  // Set initial transform //////////////////////////////////////////////////
  if (P.transform != "")
  {
    cout << endl << "reading initial transform '" << P.transform << "'..."
        << endl;

//     // try to read simple text
//     bool st = true;
//     MATRIX* mi = MatrixAlloc(4,4,MATRIX_REAL);
//     ifstream f;
//     while (1==1) // fake while loop (to be run once)
//     {
//       string sin (P.transform);
//       if (sin == "id" || sin == "identity.nofile")
//       {
//         mi = MatrixIdentity(4,mi);
//         break;
//       }
//
//       f.open(P.transform.c_str(),ios::in);
//       if (!f)
//       {
//         cerr <<" Read 4x4: could not open initial transform file " << P.transform <<endl;
//         st = false;
//         break;
//       }
//
//       int row, col;
//       for (row = 1 ; row <= 4 ; row++)
//       {
//         string s;
//         getline(f,s);
//         istringstream s1(s);
//         for (col = 1 ; col <= 4 ; col++)
//         {
//           s1>> mi->rptr[row][col];
//         }
//         if (!s1.good())
//         {
//           st = false;
//           break;
//         }
//       }
//       break; // quit fake while loop
//     }
//     f.close();
//
//
//     if (st)
//     {
//       R.setMinit(MyMatrix::convertMATRIX2VNL(mi));
//     }
//     MatrixFree(&mi);
//
//
//     if (!st)
//     {
    // try to read other transform
    TRANSFORM * trans = TransformRead(P.transform.c_str());
    if (!trans)
    {
      ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s", Progname,
          P.transform.c_str());
    }
    LTA* lta = (LTA *) trans->xform;
    if (!lta)
    {
      ErrorExit(ERROR_BADFILE, "%s: could not read lta transform file %s",
          Progname, P.transform.c_str());
    }
    if (!lta->xforms[0].src.valid)
    {
      cout
          << " WARNING: no source geometry (RAS) in transform, assuming movable !!!"
          << endl;
      getVolGeom(mri_mov, &lta->xforms[0].src);
    }
    if (!lta->xforms[0].dst.valid)
    {
      cout
          << " WARNING: no target geometry (RAS) in transform, assuming destination !!!"
          << endl;
      getVolGeom(mri_dst, &lta->xforms[0].dst);
    }

    // change to Ras2Ras, then swap geometries (this is important only, if the geometries in the lta
    // differ from the source and target passed on the command line):
    lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
    LTAmodifySrcDstGeom(lta, P.mri_mov, P.mri_dst);
    lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);
    if (lta->type != LINEAR_VOX_TO_VOX)
    {
      // should never happen:
      ErrorExit(ERROR_BADFILE, "%s: must be LINEAR_VOX_TO_VOX (=0), but %d",
          Progname, P.transform.c_str(), lta->type);
    }
    R.setMinitOrig(MyMatrix::convertMATRIX2VNL(lta->xforms[0].m_L));
    //if (P.debug) // apply init transform to input source image directly
    //{
    //  MRI * mri_tmp = LTAtransform(mri_mov,NULL, lta);
    //  string fn = R.getName() + "-source-init.mgz";
    //  MRIwrite(mri_tmp,fn.c_str());
    //  MRIfree(&mri_tmp);
    //}
//    }
  }

  // load initial iscale value
  if (P.iscalein != "")
  {
    double iscale = 1.0;
    ifstream f(P.iscalein.c_str(), ios::in);
    if (f.good())
    {
      f >> iscale;
      f.close();
    }
    else
    {
      ErrorExit(ERROR_BADFILE, "Load Iscale input: no such file ( %s )",
          P.iscalein.c_str());
    }
    R.setIscaleInit(iscale);
  }

  cout << endl;

  // set to uchar if passed
  if (P.uchartype)
  {
    MRI * mriuchar = MyMRI::setTypeUCHAR(mri_mov);
    MRIfree(&mri_mov);
    mri_mov = mriuchar;
    mriuchar = MyMRI::setTypeUCHAR(mri_dst);
    MRIfree(&mri_dst);
    mri_dst = mriuchar; 
  }

  // set oustide value to white if passed
  if (P.whitebgmov)
    MyMRI::setMaxOutsideVal(mri_mov);
  if (P.whitebgdst)
    MyMRI::setMaxOutsideVal(mri_dst);

  // now actually set source and target (and possibly reslice):
  // important that first everything else is set!
  R.setSourceAndTarget(mri_mov, mri_dst, !P.floattype);
  MRIfree(&mri_mov);
  MRIfree(&mri_dst);

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

  if (!strcmp(option, "MOV") || !strcmp(option, "M"))
  {
    P.mov = string(argv[1]);
    nargs = 1;
    cout << "--mov: Using " << P.mov << " as movable/source volume." << endl;
  }
  else if (!strcmp(option, "DST") || !strcmp(option, "D"))
  {
    P.dst = string(argv[1]);
    nargs = 1;
    cout << "--dst: Using " << P.dst << " as target volume." << endl;
  }
  else if (!strcmp(option, "LTA"))
  {
    P.lta = string(argv[1]);
    nargs = 1;
    cout << "--lta: Output transform as " << P.lta << " . " << endl;
  }
  else if (!strcmp(option, "VOX2VOX"))
  {
    P.lta_vox2vox = true;
    cout << "--vox2vox: Output transform as VOX2VOX. " << endl;
  }
  else if (!strcmp(option, "AFFINE") || !strcmp(option, "A"))
  {
    P.affine = true;
    cout << "--affine: Enabling affine transform!" << endl;
  }
  else if (!strcmp(option, "ISOSCALE"))
  {
    P.isoscale = true;
    cout << "--isoscale: Enabling isotropic scaling!" << endl;
  }
  else if (!strcmp(option, "INITSCALING"))
  {
    P.initscaling = true;
    cout
        << "--initscaling: Enabling initial scale adjustment based on image dimensions!"
        << endl;
  }
  else if (!strcmp(option, "ISCALE") || !strcmp(option, "I"))
  {
    P.iscale = true;
    cout << "--iscale: Enabling intensity scaling!" << endl;
  }
  else if (!strcmp(option, "ISCALEONLY") )
  {
    P.iscaleonly = true;
    P.iscale = true;
    P.affine = false;
    P.isoscale = false;
    P.initscaling = false;
    P.transonly = false;
    P.initorient = false;
    P.inittrans= false;
    cout << "--iscaleonly: No geometric transform, only iscale!" << endl;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
    P.transonly = true;
    cout << "--transonly: Using only translation!" << endl;
  }
  else if (!strcmp(option, "TRANSFORM") || !strcmp(option, "IXFORM"))
  {
    P.transform = string(argv[1]);
    nargs = 1;
    cout << "--ixform: Using previously computed initial transform: " << argv[1]
        << endl;
  }
  else if (!strcmp(option, "INITORIENT"))
  {
    P.initorient = true;
    cout << "--initorient: Using moments for initial orientation!" << endl;
  }
  else if (!strcmp(option, "NOINIT"))
  {
    P.inittrans = false;
    cout << "--noinit: Skipping init of transform !" << endl;
  }
  else if (!strcmp(option, "LEASTSQUARES"))
  {
    P.leastsquares = true;
    cout << "--leastsquares: Using standard least squares (non-robust)!"
        << endl;
  }
  else if (!strcmp(option, "COST"))
  {
    string cost(argv[1]);
    std::transform(cost.begin(), cost.end(), cost.begin(), ::toupper);

    nargs = 1;
    if (cost == "LS")
      P.cost = Registration::LS;
    else if (cost == "ROB")
      P.cost = Registration::ROB;
    else if (cost == "ROBENT")
    {
      P.cost = Registration::ROB;
      P.entropy = true;
    }
    else if (cost == "MI")
      P.cost = Registration::MI;
    else if (cost == "NMI")
      P.cost = Registration::NMI;
    else if (cost == "ECC")
      P.cost = Registration::ECC;
    else if (cost == "NCC")
      P.cost = Registration::NCC;
    else if (cost == "SCR")
      P.cost = Registration::SCR;
    else if (cost == "TB")
      P.cost = Registration::TB;
    else if (cost == "LNCC")
      P.cost = Registration::LNCC;
    else if (cost == "SAD")
      P.cost = Registration::SAD;
    else if (cost == "SB")
      P.cost = Registration::SB;
    else
    {
      cout << "ERROR: cost function " << cost << " unknown! " << endl;
      exit(1);
    }

    cout << "--cost: Using cost function: " << cost << " !" << endl;
  }
  else if (!strcmp(option, "ENTRADIUS")|| !strcmp(option, "RADIUS"))
  {
    P.entroradius = atoi(argv[1]);
    nargs = 1;
    cout << "--entradius: Using local boxes with radius " << P.entroradius << " = "
        << P.entroradius * 2 + 1 << " sides. " << endl;
  }
  else if (!strcmp(option, "ENTDST"))
  {
    P.entdst = string(argv[1]);
    nargs = 1;
    cout << "--entdst: Output entropy dst image as " << P.entdst << endl;
  }
  else if (!strcmp(option, "ENTBALL"))
  {
    P.entball = true;
    nargs = 0;
    cout << "--entball: Using ball instead of local boxes. " << endl;
  }
  else if (!strcmp(option, "ENTCORRECTION"))
  {
    P.entcorrection = true;
    nargs = 0;
    cout << "--entcorrection: Using 'correction' mode for entropy images. "
        << endl;
  }
  else if (!strcmp(option, "ENTMOV"))
  {
    P.entmov = string(argv[1]);
    nargs = 1;
    cout << "--entmov: Output entropy mov image as " << P.entmov << endl;
  }
  else if (!strcmp(option, "SOBEL"))
  {
    P.sobel = true;
    nargs = 0;
    cout << "--sobel: Using Sobel magnitude images."<< endl;
  }
  else if (!strcmp(option, "MAXIT"))
  {
    P.iterate = atoi(argv[1]);
    nargs = 1;
    cout << "--maxit: Performing maximal " << P.iterate
        << " iterations on each resolution" << endl;
  }
  else if (!strcmp(option, "HIGHIT"))
  {
    P.highit = atoi(argv[1]);
    nargs = 1;
    cout << "--highit: Performing maximal " << P.highit
        << " iterations on highest resolution" << endl;
  }
  else if (!strcmp(option, "EPSIT"))
  {
    P.epsit = atof(argv[1]);
    nargs = 1;
    cout << "--epsit: Stop iterations when change is less than " << P.epsit
        << " . " << endl;
  }
  else if (!strcmp(option, "NOMULTI"))
  {
    P.nomulti = true;
    nargs = 0;
    cout << "--nomulti: Will work on highest resolution only (nomulti)!"
        << endl;
  }
  else if (!strcmp(option, "MAXSIZE"))
  {
    P.maxsize = atoi(argv[1]);
    nargs = 1;
    cout << "--maxsize: Largest dimension < " << P.maxsize << " . " << endl;
  }
  else if (!strcmp(option, "MINSIZE"))
  {
    P.minsize = atoi(argv[1]);
    nargs = 1;
    cout << "--minsize: Smallest dimension > " << P.minsize << " . " << endl;
  }
  else if (!strcmp(option, "SAT"))
  {
    P.sat = atof(argv[1]);
    nargs = 1;
    cout << "--sat: Using saturation " << P.sat << " in M-estimator!" << endl;
  }
  else if (!strcmp(option, "WLIMIT"))
  {
    P.wlimit = atof(argv[1]);
    nargs = 1;
    cout << "--wlimit: Using wlimit in satit " << P.wlimit << endl;
  }
  else if (!strcmp(option, "SUBSAMPLE"))
  {
    if (argc == 1)
      ErrorExit(ERROR_BADPARM,
          "\nERROR: specify min side lenght to start subsampling, e.g. --subsample 200.\n",
          Progname);

    P.subsamplesize = atoi(argv[1]);
    nargs = 1;
    if (P.subsamplesize >= 0)
    {
      cout << "--subsample: Will subsample if size is larger than "
          << P.subsamplesize << " on all axes!" << endl;
    }
    else
    {
      cout << "--subsample -1: Will not subsample on any scale!" << endl;
    }
  }
  else if (!strcmp(option, "SATIT"))
  {
    P.satit = true;
    nargs = 0;
    cout
        << "--satit: Will iterate with different SAT to ensure outliers below wlimit!"
        << endl;
  }
//   else if (!strcmp(option, "SATEST") ) // old  remove
//   {
//     P.satest = true;
//     nargs = 0 ;
//     cout << "--satest: Will estimate SAT (never really tested, use --satit instead!)" << endl;
//   }
  else if (!strcmp(option, "SATEST")) // never reached???  - old remove
  {
    P.dosatest = true;
    nargs = 0;
    cout << "--satest: Trying to estimate SAT value!" << endl;
  }
  else if (!strcmp(option, "POWELLTOL"))
  {
    P.powelltol = atof(argv[1]);
    nargs = 1;
    cout << "--powelltol: Using tolerance " << P.powelltol << " in Powell optimizier!" << endl;
  }
  else if (!strcmp(option, "DOUBLEPREC"))
  {
    P.doubleprec = true;
    nargs = 0;
    cout
        << "--doubleprec: Will perform algorithm with double precision (higher mem usage)!"
        << endl;
  }
  else if (!strcmp(option, "DEBUG"))
  {
    P.debug = 1;
    nargs = 0;
    cout << "--debug: Will output debug info and files!" << endl;
  }
  else if (!strcmp(option, "VERBOSE"))
  {
    P.verbose = atoi(argv[1]);
    nargs = 1;
    cout << "--verbose: Will use verbose level : " << P.verbose << endl;
  }
  else if (!strcmp(option, "WEIGHTS"))
  {
    P.weightsout = string(argv[1]);
    nargs = 1;
    cout << "--weights: Will output weights (in target space) as "
        << P.weightsout << " !" << endl;
  }
  else if (!strcmp(option, "WARP") || !strcmp(option, "MAPMOV"))
  {
    P.warpout = string(argv[1]);
    nargs = 1;
    cout << "--mapmov: Will save resampled movable as " << P.warpout << " !"
        << endl;
  }
  else if (!strcmp(option, "MAPMOVHDR"))
  {
    P.norlout = string(argv[1]);
    nargs = 1;
    cout << "--mapmovhdr: Will save header adjusted movable as " << P.norlout
        << " !" << endl;
  }
  else if (!strcmp(option, "HALFMOV"))
  {
    P.halfmov = string(argv[1]);
    nargs = 1;
    cout << "--halfmov: Will output final half way MOV !" << endl;
  }
  else if (!strcmp(option, "HALFDST"))
  {
    P.halfdst = string(argv[1]);
    nargs = 1;
    cout << "--halfdst: Will output final half way DST !" << endl;
  }
  else if (!strcmp(option, "HALFWEIGHTS"))
  {
    P.halfweights = string(argv[1]);
    nargs = 1;
    cout << "--halfweights: Will output half way WEIGHTS from last step to "
        << P.halfweights << " !" << endl;
  }
  else if (!strcmp(option, "HALFMOVLTA"))
  {
    P.halfmovlta = string(argv[1]);
    nargs = 1;
    cout << "--halfmovlta: Will output half way transform (mov) "
        << P.halfmovlta << " !" << endl;
  }
  else if (!strcmp(option, "HALFDSTLTA"))
  {
    P.halfdstlta = string(argv[1]);
    nargs = 1;
    cout << "--halfdstlta: Will output half way transform (dst) "
        << P.halfdstlta << " !" << endl;
  }
  else if (!strcmp(option, "MASKMOV"))
  {
    P.maskmov = string(argv[1]);
    nargs = 1;
    cout << "--maskmov: Will apply " << P.maskmov << " to mask mov/src !"
        << endl;
  }
  else if (!strcmp(option, "MASKDST"))
  {
    P.maskdst = string(argv[1]);
    nargs = 1;
    cout << "--maskdst: Will apply " << P.maskdst << " to mask dst/target !"
        << endl;
  }
  else if (!strcmp(option, "TEST"))
  {
    cout << "--test: TEST-MODE " << endl;
    RegRobust R;
    R.testRobust(argv[2], atoi(argv[1]));
    nargs = 2;
    exit(0);
  }
  else if (!strcmp(option, "CONFORM"))
  {
    P.conform = true;
    nargs = 0;
    cout << "--conform: Will conform images to 256^3 and voxels to 1mm!"
        << endl;
  }
  else if (!strcmp(option, "FLOATTYPE"))
  {
    P.floattype = true;
    nargs = 0;
    cout << "--floattype: Use float images internally (independent of input)!"
        << endl;
  }
  else if (!strcmp(option, "ONEMINUSW"))
  {
    P.oneminusweights = false;
    nargs = 0;
    cout
        << "--oneminusw: Will output 1-weights (zero=outlier), as in earlier versions!"
        << endl;
  }
  else if (!strcmp(option, "NOSYM"))
  {
    P.symmetry = false;
    nargs = 0;
    cout << "--nosym: Will resample source to target (no half-way space)!"
        << endl;
  }
  else if (!strcmp(option, "WHITEBGMOV"))
  {
    P.whitebgmov = true;
    nargs = 0;
    cout << "--whitebgmov: Will assume white background in MOV!"
        << endl;
  }
  else if (!strcmp(option, "WHITEBGDST"))
  {
    P.whitebgdst = true;
    nargs = 0;
    cout << "--whitebgdst: Will assume white background in DST!"
        << endl;
  }
  else if (!strcmp(option, "UCHAR"))
  {
    P.uchartype = true;
    nargs = 0;
    cout << "--uchar: Will convert images to uchar (with re-scaling and histogram cropping)!"
        << endl;
  }
  else if (!strcmp(option, "ISCALEOUT"))
  {
    P.iscaleout = string(argv[1]);
    nargs = 1;
    P.iscale = true;
    cout << "--iscaleout: Will do --iscale and ouput intensity scale to "
        << P.iscaleout << endl;
  }
  else if (!strcmp(option, "ISCALEIN"))
  {
    nargs = 1;
    P.iscalein = string(argv[1]);
    cout << "--iscalein: Will use init intensity scale" << endl;
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
  int nargs;
  int inputargs = argc;
  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv, P);
    argc -= nargs;
    argv += nargs;
  }

  if (inputargs == 0)
  {
    printUsage();
    exit(1);
  }

  bool test1 = (P.mov != "" && P.dst != "");
  if (!test1)
  {
    printUsage();
    cerr << endl << endl << "ERROR: Please specify --mov and --dst !  "
        << endl << endl;
    exit(1);
  }
  bool test1b = (P.lta != "" || P.iscaleonly);
  if (!test1b)
  {
    printUsage();
    cerr << endl << endl << "ERROR: Please specify --lta !  "
        << endl << endl;
    exit(1);
  }
  bool test2 = (P.satit || P.sat > 0 || P.cost != Registration::ROB
      || P.leastsquares);
  if (!test2)
  {
    printUsage();
    cerr << endl << endl
        << "ERROR: Please specify either --satit or --sat <float> !  " << endl
        << endl;
    exit(1);
  }
  bool test3 = (P.iscaleout == "" || P.iscale);
  if (!test3)
  {
    printUsage();
    cerr << endl << endl
        << "ERROR: Please specify --iscale together with --iscaleout to compute and output global intensity scaling! "
        << endl << endl;
    exit(1);
  }
  bool test4 = (P.warpout == "" || (P.warpout != P.weightsout));
  if (!test4)
  {
    printUsage();
    cerr << endl << endl
        << "ERROR: Resampled input name (--mapmov) cannot be same as --weights output!"
        << endl << endl;
    exit(1);
  }
  bool test5 = (P.iterate > 0);
  if (!test5)
  {
    //printUsage();
    cerr << endl << endl
        << "ERROR: --maxit must be larger than 0 (max iterations on each resolution)"
        << endl << endl;
    exit(1);
  }

  return (test1 && test2 && test3 && test4 && test5);
}
