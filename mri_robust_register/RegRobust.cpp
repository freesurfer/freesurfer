/**
 * @brief A class to compute a registration using robust regression
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

#include "RegRobust.h"
#include "RegistrationStep.h"

#if ITK_VERSION_MAJOR >= 5
#include <iostream>
#include <vcl_compiler.h>
#else
#include <vcl_iostream.h>
#endif

RegRobust::~RegRobust()
{ // we cleanup our private variables
//std::cout << " Destroy Registration" << std::endl;
  if (mri_indexing)
    MRIfree(&mri_indexing);
  if (mri_weights)
    MRIfree(&mri_weights);
  if (mri_hweights)
    MRIfree(&mri_hweights);
  //std::cout << " Done " << std::endl;
}

void RegRobust::clear() // initialize registration (keep source and target and gauss pyramid)
// initialize registration (keep source and target and gauss pyramid)
// also keep Rsrc and Rtrg (resampling matrices, if exist).
{
  Registration::clear();
  sat = -1;
  if (mri_indexing)
    MRIfree(&mri_indexing);
  if (mri_weights)
    MRIfree(&mri_weights);
  if (mri_hweights)
    MRIfree(&mri_hweights);
}

void RegRobust::computeIterativeRegistrationFull(int nmax, double epsit, MRI * mriS,
    MRI* mriT, const vnl_matrix<double>& m, double scaleinit)
// private routine, as called from multiregistration (passing mriS and mriT...)
// computes iterative registration (recomputing A and b in each step)
// retruns 4x4 matrix Mfinal and iscalefinal (class member)
// The caller needs to retrieve any really final transform with getFinalVox2Vox
{

  // call helper to avoid code duplication:

  if (doubleprec)
    iterativeRegistrationHelper<double>(nmax, epsit, mriS, mriT, m, scaleinit);
  else
    iterativeRegistrationHelper<float>(nmax, epsit, mriS, mriT, m, scaleinit);

  return;

}
void RegRobust::findSatMultiRes(const vnl_matrix<double> &mi, double scaleinit)
// helper for findSaturation
// basically the code from multiresoltuion
// all kinds of stuff is initialized before (e.g. pyramid)
{
  // variables to store matrix m and scaling factor d:
  pair<vnl_matrix_fixed<double, 4, 4>, double> cmd;
  pair<vnl_matrix_fixed<double, 4, 4>, double> md(mi, scaleinit);

  // allow 2d case (depth == 1)
  if (gpS[0]->width < 16 || gpS[0]->height < 16
      || (gpS[0]->depth < 16 && gpS[0]->depth != 1))
  {
    ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n");
  }
  int resolution = gpS.size();
  assert(resolution >= 1);
  // otherwise we should have exited above
  int rstart = 1;  // at least 16^3

  // stop if we get larger than 64^3 or if we reach highest resolution:
  int stopres;
  for (stopres = resolution - rstart; stopres > 0; stopres--)
  {
    if (gpS[stopres]->width >= 64 || gpS[stopres]->height >= 64
        || gpS[stopres]->depth >= 64)
      break;
  }

//  bool iscaletmp = iscale;
//  iscale = false; //disable intensity scaling on low resolutions

  for (int r = resolution - rstart; r >= stopres; r--)
  {
    if (verbose > 1)
    {
      cout << endl << "Resolution: " << r << endl;
      cout << " gpS ( " << gpS[r]->width << " , " << gpS[r]->height << " , "
          << gpS[r]->depth << " )" << endl;
      cout << " gpT ( " << gpT[r]->width << " , " << gpT[r]->height << " , "
          << gpT[r]->depth << " )" << endl;
    }

//    if (r==2) iscale = iscaletmp; // set iscale if set by user

    // compute Registration
    if (verbose > 2)
      cout << "   - compute new iterative registration" << endl;

    int n = 3;
    if (r == stopres)
      n = 1;
    int vv = verbose;
    if (verbose == 1)
      verbose = 0;
    computeIterativeRegistrationFull(n, 0.05, gpS[r], gpT[r], md.first, md.second);
    cmd.first = Mfinal;
    cmd.second = iscalefinal;
    verbose = vv;

    if (verbose > 1)
    {
      cout << endl << " current Matrix: " << endl;
      vnl_matlab_print(vcl_cout,cmd.first,"Tc",vnl_matlab_print_format_long);
      cout << endl;
      cout << " intens: Ic = " << cmd.second << endl;

      // adjust to highest level for output only:
      double tx = cmd.first[0][3];
      double ty = cmd.first[1][3];
      double tz = cmd.first[2][3];
      for (int ll = r; ll > 0; ll--)
      {
        tx *= 2;
        ty *= 2;
        tz *= 2;
      }
      cout << " equiv trans on highres: " << tx << " " << ty << " " << tz
          << endl;
    }

//     if (r == stopres)
//     {
//        if (debug)
//        {
//          // write out wcheck
//          string fn = getName() + "-wcheck-est.txt";
//          ofstream f(fn.c_str(),ios::out);
//          f << sat << " " << wcheck << endl;
//          f.close();  
//          string fn2 = getName() + "-wchecksqrt-est.txt";
//          ofstream f2(fn.c_str(),ios::out);
//          f2 << sat << " " << wchecksqrt << endl;
//          f2.close();  
//        }
//        
// //        if (wcheck > wlimit)
// //        {
// //           sat = sat+0.5;
// //           if (verbose > 1) cout << "   - Weight check " << wcheck << " > "<< wlimit  << " increasing sat: " << sat << endl;
// //           md.first = firstbackup;
// //           md.second = scaleinit;
// //           r = resolution-rstart+1;
// //           continue;
// //         }
//     }

    if (r != 0) // adjust matrix to higher resolution level
    {
      for (int rr = 0; rr < 3; rr++)
      {
        cmd.first[rr][3] = 2.0 * cmd.first[rr][3];
      }
    }
    if (r == stopres) // passed the test, adjust to highest
    {
      for (; r > 1; r--) // adjusted once already in the lines above
      {
        for (int rr = 0; rr < 3; rr++)
        {
          cmd.first[rr][3] = 4.0 * cmd.first[rr][3];
        }
      }
    }
    md.first = cmd.first;
    md.second = cmd.second;
    if (verbose > 1)
    {
      cout << endl << " Matrix: " << endl;
      vnl_matlab_print(vcl_cout,md.first,"T",vnl_matlab_print_format_long);
      cout << endl;
      cout << " Intensity:  I = " << md.second << endl;
    }
  } // resolution loop

}

double RegRobust::findSaturation()
{
  if (verbose > 0)
    cout << endl << endl << " Registration::findSaturation " << endl;
//   if (!mriS) mriS = mri_source;
//   if (!mriT) mriT = mri_target;
  MRI * mriS = mri_source;
  MRI * mriT = mri_target;

  vnl_matrix_fixed<double, 4, 4> m;
  m.set_identity();

  // variables to store matrix m and scaling factor d:
  pair<vnl_matrix_fixed<double, 4, 4>, double> md(
      vnl_matrix_fixed<double, 4, 4>(), iscaleinit);

  if (!Minit.empty())
    md.first = getMinitResampled();
  else
    md.first = initializeTransform(mriS, mriT);

//   if (scaleinit != 1.0) md.second = scaleinit;
//   else md.second = iscaleinit;

  int MINS = 16;
  if (minsize > MINS)
    MINS = minsize; // use minsize, but at least 16
  pair<int, int> limits = getGPLimits(mriS, mriT, MINS, maxsize);
  if (gpS.size() == 0)
    gpS = buildGPLimits(mriS, limits);
  if (gpT.size() == 0)
    gpT = buildGPLimits(mriT, limits);
  assert(gpS.size() == gpT.size());
  if (gpS[0]->width < MINS || gpS[0]->height < MINS
      || (gpS[0]->depth < MINS && gpS[0]->depth != 1))
  {
    ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n");
  }
  if (gpT[0]->width < MINS || gpT[0]->height < MINS
      || (gpT[0]->depth < MINS && gpT[0]->depth != 1))
  {
    ErrorExit(ERROR_BADFILE, "Input images must be larger than 16^3.\n");
  }

  int resolution = gpS.size();

  assert(resolution >= 1);
  // otherwise we should have exited above
  int rstart = 1;  // at least 16^3, last and most coarse image

  // stop if we get larger than 64^3
  // should be as close as possible to 64
  // the problem is that on other scales the wcheck limit gets meaningless
  int stopres;
  for (stopres = resolution - rstart; stopres > 0; stopres--)
  {
    if (gpS[stopres]->width >= 64 || gpS[stopres]->height >= 64
        || gpS[stopres]->depth >= 64)
      break;
  }

  if (gpS[stopres]->width < 32 || gpS[stopres]->height < 32
      || gpS[stopres]->depth < 32)
  {
    cout << endl
        << " ========================================================================"
        << endl;
    cout
        << " WARNING: image might be too small (or ill shaped) for --satit to work."
        << endl;
    cout
        << "          Try to manually specify --sat # if not satisfied with result! "
        << endl;
    cout
        << " ========================================================================"
        << endl << endl;
    ;
  }

  cout << endl << "   - Max Resolution used: " << stopres << endl;
  cout << "     -- gpS ( " << gpS[stopres]->width << " , "
      << gpS[stopres]->height << " , " << gpS[stopres]->depth << " )" << endl;
  cout << "     -- gpT ( " << gpT[stopres]->width << " , "
      << gpT[stopres]->height << " , " << gpT[stopres]->depth << " )" << endl;

  if (verbose > 1)
  {
    cout << "   - initial transform:\n";
    vnl_matlab_print(vcl_cout,md.first,"Ti",vnl_matlab_print_format_long);
    cout << endl;
    cout << "   - initial iscale:   Ii = " << md.second << endl;
  }

  // adjust md.first to current (lowest) resolution:
  for (int r = 1; r <= resolution - rstart; r++)
    for (int rr = 0; rr < 3; rr++)
      md.first[rr][3] = 0.5 * md.first[rr][3];

  vnl_matrix_fixed<double, 4, 4> firstbackup = md.first;

  if (verbose > 1)
  {
    cout << "   - initial adjusted:\n";
    vnl_matlab_print(vcl_cout,md.first,"Tia",vnl_matlab_print_format_long);
    cout << endl;
  }

  // -------------------------------------------- RUN LOOP ----------------------------------
  // 
  cout << "   - running loop to estimate saturation parameter:\n";
  double satdiff = 0.5; // stop if we get closer than this
  double satmax = 0;
  double satmin = 0;
  double wmin = -1;
  double wmax = -1;
  int counter = 0;
  while (satmin == 0.0 || satmax == 0.0 || satmax - satmin > satdiff)
  {
    counter++;
    if (satmin == 0 && satmax == 0)
      sat = 16;
    else if (satmin == 0)
      sat = 0.5 * satmax;
    else if (satmax == 0)
      sat = 2 * satmin;
    else
      sat = 0.5 * (satmax + satmin);
    if (verbose > 0)
    {
      if (counter > 1)
        cout << "         min sat: " << satmin << " ( " << wmin
            << " ), max sat: " << satmax << " ( " << wmax << " ), sat diff: "
            << satmax - satmin << ", (wlimit=" << wlimit << ")"<< endl;
      cout << "     -- Iteration: " << counter << "  trying sat: " << sat
          << endl;
    }
    findSatMultiRes(md.first, md.second);
    if (wcheck > wlimit)
    {
      satmin = sat;
      wmin = wcheck;
    }
    else
    {
      satmax = sat;
      wmax = wcheck;
    }

    // if sat low (sensitive) and still not many outliers
    // quit, to prevent long (infinite) search
    // e.g. if source and target are same image
    if (sat < 6 && wcheck < 0.04)
    {
      satmax = sat;
      satmin = sat;
      wmax = wcheck;
      wmin = wcheck;
      break;
    }
  }

  // -------------------------------------------- SELECT FINAL ---------------------------------
  // 
  if (wmax <= wlimit)
  {
    sat = satmax;
    wcheck = wmax;
  }
  else
  {
    assert(wmin <= wlimit);
    sat = satmin;
    wcheck = wmin;
  }

  if (verbose > 0)
    cout << "   - final SAT: " << sat << " ( it: " << counter
        << " , weight check " << wcheck << " <= " << wlimit << " )" << endl;

  if (debug)
  {
    // write out wcheck
    string fn = getName() + "-wcheck-est.txt";
    ofstream f(fn.c_str(), ios::out);
    f << sat << " " << wcheck << endl;
    f.close();
  }

  return sat;
}

double RegRobust::estimateIScale(MRI *mriS, MRI *mriT)
{
  if (mriS->nframes > 1 || mriT->nframes > 1)
  {
    cerr << "RegRobust::estimateIScale multiple frames not supported!" << endl;
    exit(1);
    // see constructAB for modifications necessary to do multi frame
  }


  if (verbose > 1)
    cout << "   - estimateIScale: " << endl;

  assert(mriT != NULL);
  assert(mriS != NULL);
  assert(mriS->width == mriT->width);
  assert(mriS->height== mriT->height);
  assert(mriS->depth == mriT->depth);
  assert(mriS->type == mriT->type);
  //assert(mriS->width == mask->width);
  //assert(mriS->height== mask->height);
  //assert(mriS->depth == mask->depth);
  //assert(mask->type == MRI_INT);
  //MRIclear(mask);

  int z, y, x;
  long int ss = mriS->width * mriS->height * mriS->depth;
  if (mri_indexing)
    MRIfree(&mri_indexing);
  if (ss > std::numeric_limits<int>::max())
  {
    if (verbose > 1)
      cout << "     -- using LONG for indexing ... " << flush;
    mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth, MRI_LONG);
    if (mri_indexing == NULL)
      ErrorExit(ERROR_NO_MEMORY,
          "Registration::estimateIScale could not allocate memory for mri_indexing");
    if (verbose > 1)
      cout << " done!" << endl;
  }
  else
  {
    double mu = ((double) ss) * sizeof(int) / (1024.0 * 1024.0);
    if (verbose > 1)
      cout << "     -- allocating " << mu << "Mb mem for indexing ... "
          << flush;
    mri_indexing = MRIalloc(mriS->width, mriS->height, mriS->depth, MRI_INT);
    if (mri_indexing == NULL)
      ErrorExit(ERROR_NO_MEMORY,
          "Registration::estimateIScale could not allocate memory for mri_indexing");
    if (verbose > 1)
      cout << " done!" << endl;
  }

  for (z = 0; z < mriS->depth; z++)
    for (x = 0; x < mriS->width; x++)
      for (y = 0; y < mriS->height; y++)
        MRILvox(mri_indexing, x, y, z) = 0;

  bool dosubsample = false;
  if (subsamplesize > 0)
    dosubsample = (mriS->width > subsamplesize && mriS->height > subsamplesize
        && (mriS->depth > subsamplesize || mriS->depth == 1));
  dosubsample = false; // needs to be fixed below!!! indeces are now randomized

  // we will need the blurred images (as float):
  if (verbose > 1)
    cout << "     -- compute smoothie ... " << flush;
  MRI *Sbl = MRIalloc(mriS->width, mriS->height, mriS->depth, MRI_FLOAT);
  Sbl = MRIcopy(mriS, Sbl);
  Sbl = MyMRI::getBlur(Sbl, Sbl);
  MRI *Tbl = MRIalloc(mriT->width, mriT->height, mriT->depth, MRI_FLOAT);
  Tbl = MRIcopy(mriT, Tbl);
  Tbl = MyMRI::getBlur(Tbl, Tbl);

  if (verbose > 1)
    cout << " done!" << endl;

  if (dosubsample)
  {
    if (verbose > 1)
      cout << "     -- subsample ... " << flush;

    MRI * Sblt = Sbl;
    Sbl = MyMRI::subSample(Sblt);
    MRIfree(&Sblt);
    MRI * Tblt = Tbl;
    Tbl = MyMRI::subSample(Tblt);
    MRIfree(&Tblt);

    if (verbose > 1)
      cout << " done! " << endl;
  }

  // compute 'counti': the number of rows needed (zero elements need to be removed)
  int n = Sbl->width * Sbl->height * Sbl->depth;
  if (verbose > 1)
    cout << "     -- size " << Sbl->width << " x " << Sbl->height << " x "
        << Sbl->depth << " = " << n << flush;
  long int counti = 0;
  double eps = 0.00001;
  for (z = 0; z < Sbl->depth; z++)
    for (x = 0; x < Sbl->width; x++)
      for (y = 0; y < Sbl->height; y++)
      {
        if (isnan(MRIFvox(Sbl, x, y, z) )||isnan(MRIFvox(Tbl, x, y, z)))
        {
          if (verbose > 0) cout << " found a nan value!!!" << endl;
          continue;
        }
        if (fabs(MRIFvox(Sbl, x, y, z)) < eps && fabs(MRIFvox(Tbl, x, y, z)) <eps  )
        {
          //if (verbose > 0) cout << " found a zero element !!!" << endl;
          continue;
        }
        counti++; // start with 1
      }
  if (verbose > 1)
    cout << "  need only: " << counti << endl;
  if (counti == 0)
  {
    cerr << endl;
    cerr << " ERROR: All entries are zero! Images do not overlap (anymore?)."
        << endl;
    cerr
        << "    This can have several reasons (i.e. different modalities, different "
        << endl;
    cerr
        << "    intensity scales, large non-linearities, too diff. voxel sizes ...)"
        << endl;
    //cerr << "    Try calling with --noinit (if the original images are well aligned)" << endl;
    cerr << "    Maybe use --transform <init.lta> with an approx. alignment"
        << endl;
    cerr << "    obtained from tkregister or another registration program."
        << endl;
    cerr << "    Or do some prior intensity correction? " << endl;
    cerr << endl;
    exit(1);
  }

  // allocate the space for A and B
  double abmu = ((double) counti) * sizeof(double) / (1024.0 * 1024.0);
  if (verbose > 1)
    cout << "     -- allocating " << abmu << "Mb mem for A and b ... " << flush;
  pair<vnl_matrix<double>, vnl_vector<double> > Ab(
      vnl_matrix<double>(counti, 1), vnl_vector<double>(counti));
  if (verbose > 1)
    cout << " done! " << endl;
//      if (A == NULL || b == NULL) 
//         ErrorExit(ERROR_NO_MEMORY,"Registration::estimateIScale could not allocate memory for A or b") ;

  if (verbose > 1)
    cout << "     -- size " << Sbl->width << " " << Sbl->height << " "
        << Sbl->depth << flush;

  long int count = 0;
  int xp1, yp1, zp1;
  for (z = 0; z < Sbl->depth; z++)
    for (x = 0; x < Sbl->width; x++)
      for (y = 0; y < Sbl->height; y++)
      {
        if (isnan(MRIFvox(Sbl, x, y, z) )||isnan(MRIFvox(Tbl, x, y, z)) )
        {
          if (verbose > 0) cout << " found a nan value!!!" << endl;
          continue;
        }

        if (dosubsample)
        {
          xp1 = 2*x;
          yp1 = 2*y;
          zp1 = 2*z;
        }
        else
        {
          xp1 = x;
          yp1 = y;
          zp1 = z; // if not subsampled
          }
          assert(xp1 < mriS->width);
          assert(yp1 < mriS->height);
          assert(zp1 < mriS->depth);

          if (fabs(MRIFvox(Sbl, x, y, z)) < eps && fabs(MRIFvox(Tbl, x, y, z)) < eps )
          {
            //cout << " found a zero row!!!" << endl;
            MRILvox(mri_indexing, xp1, yp1, zp1) = -1;
            continue;
          }

          //count++; // start with 1

          if (xp1 >= mriS->width || yp1 >= mriS->height || zp1 >= mriS->depth)
          {

            cerr << " outside !!! " << xp1 << " " << yp1 << " " << zp1 << endl;
            assert(1==2);
          }

          MRILvox(mri_indexing, xp1, yp1, zp1) = count;

          //Ab.first[count][0]  = 0.5 / iscalefinal *( MRIFvox(Tbl,x,y,z) + MRIFvox(Sbl, x, y, z)); 
          //Ab.first[count][0]  = MRIFvox(Sbl, x, y, z);

          // intensity model: R(s,IS,IT) = exp(-0.5 s) IT - exp(0.5 s) IS
          //                  R'  = -0.5 ( exp(-0.5 s) IT + exp(0.5 s) IS)
          Ab.first[count][0] = 0.5 * (MRIFvox(Tbl,x,y,z) + MRIFvox(Sbl, x, y, z));

          Ab.second[count] = -(MRIFvox(Tbl, x, y, z) - MRIFvox(Sbl, x, y, z));

          count++;// start with 0

        }

        // free remaining MRI    
  if (Sbl)
    MRIfree(&Sbl);
  if (Tbl)
    MRIfree(&Tbl);

  assert(counti == count);

  Regression<double> R(Ab.first, Ab.second);
  R.setVerbose(verbose);
  R.setFloatSvd(true); // even for double, the svd can be float, better switch to float all toghether

  vnl_vector<double> p(R.getRobustEst());

  double is = p[0];
  double s = log(iscalefinal);
  s = s - is;
  iscalefinal = exp(s);
  cout << " ISCALE: " << iscalefinal << " returned: " << is << endl;

  return iscalefinal;
}
