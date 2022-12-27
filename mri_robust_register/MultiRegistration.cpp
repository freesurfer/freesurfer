/**
 * @brief A class to handle registration of multiple files
 *
 * MultiRegistration is a class to compute a robust registration
 *  of several images. It makes use routines from Registration 
 *
 * written by Martin Reuter
 *  Aug. 12th ,2009
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

#include "MultiRegistration.h"
#include "RegRobust.h"
#include "Regression.h"
#include "RobustGaussian.h"
#include "CostFunctions.h"
#include "MyMatrix.h"
#include "MyMRI.h"
#include "mriBSpline.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_determinant.h>

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"

using namespace std;

void MultiRegistration::clear()
{
  if (mri_mean)
    MRIfree(&mri_mean);

  for (unsigned int i = 0; i < mri_mov.size(); i++)
  {
    if (i < ltas.size() && ltas[i])
      LTAfree(&ltas[i]);
    if (i < mri_warps.size() && mri_warps[i])
      MRIfree(&mri_warps[i]);
    if (i < mri_weights.size() && mri_weights[i])
      MRIfree(&mri_weights[i]);
    if (i < mri_mov.size() && mri_mov[i])
      MRIfree(&mri_mov[i]);
    if (i < mri_bsplines.size() && mri_bsplines[i])
      MRIfreeBSpline(&mri_bsplines[i]);
  }

}

unsigned int MultiRegistration::getSeed()
// computes a seed for randomization based on the input images
// takes the middle slice 
{
  int n = (int) mri_mov.size();
  assert(n > 1);

  double dseed = 0.0;
  int x, y, z, i, p, xup, xdown, yup, ydown, zup, zdown;
  for (i = 0; i < n; i++)
  {
    // center voxel:
    x = mri_mov[i]->width / 2;
    y = mri_mov[i]->height / 2;
    z = mri_mov[i]->depth / 2;
    
    // sum up crossair
    for (p = 0; p < 20; p++)
    {
      xup = x + p;
      xdown = x - p;
      yup = y + p;
      ydown = y - p;
      zup = z + p;
      zdown = z - p;

      if (xdown >= 0)
        dseed += fabs(MRIgetVoxVal(mri_mov[i],xdown,y,z,0)) ;
      if (ydown >= 0)
        dseed += fabs(MRIgetVoxVal(mri_mov[i],x,ydown,z,0)) ;
      if (zdown >= 0)
        dseed += fabs(MRIgetVoxVal(mri_mov[i],x,y,zdown,0)) ;
      if (xup < mri_mov[i]->width)
        dseed += fabs(MRIgetVoxVal(mri_mov[i],xup,y,z,0)) ;
      if (yup < mri_mov[i]->height)
        dseed += fabs(MRIgetVoxVal(mri_mov[i],x,yup,z,0)) ;
      if (zup < mri_mov[i]->depth)
        dseed += fabs(MRIgetVoxVal(mri_mov[i],x,y,zup,0)) ;

      }
    }
  if (dseed == 0.0)
  {
    std::cout << "WARNING: image seed is zero, are images really empty at center?" << std::endl;
  } 
  else // bring into range from 10 ... 100  if dseed is small, as we convert to int
    while (dseed < 10) dseed = dseed * 10;
  
  //std::cout << "MultiRegistration::getSeed  = " << dseed << std::endl; 
  return (unsigned int)(dseed);
}

/*!
 \fn int loadMovables(const std::vector < std::string > pmov)
 \brief Loads the movable volumes as specified on command line
 \param P  Paramters for the initialization
 */
int MultiRegistration::loadMovables(const std::vector<std::string>& pmov,const std::vector<std::string>& pmasks)
{

  assert(mri_mov.size () == 0);
  int n = (int) pmov.size();
  int m = (int) pmasks.size();
  assert (n==m || m == 0);
  
  mov = pmov; // local copy of input filenames
  if (n == 1) // try read 4D volume
  {
    cout << "trying to read 4D source '" << mov[0] << "'..." << endl;
    MRI *mri4d = MRIread(mov[0].c_str());
    if (!mri4d)
    {
      ErrorExit(ERROR_NOFILE,
          "MultiRegistration::loadMovables: could not open input volume %s.\n",
          mov[0].c_str());
    }
    n = mri4d->nframes;
    cout << " Input has " << n << " frames." << endl;
    if (n == 1)
    {
      cerr << "ERROR: only single 3D volume passed!" << endl;
      cerr << "    Pass several inputs, or 4D input volume" << endl;
      exit(1);
    }
    
    MRI *mask4d = NULL;
    if (m == 1)
    {
      mask4d = MRIread(pmasks[0].c_str());
      if (!mri4d)
      {
        ErrorExit(ERROR_NOFILE,
          "MultiRegistration::loadMovables: could not open input volume %s.\n",
          pmasks[0].c_str());
      }
      m = mask4d->nframes;
      if (n != m)
      {
        cerr << "ERROR: frame count in mask does not agree with frames in mov!" << endl;
        exit(1);
      }
    } 
    mov.resize(n);
    mri_mov.resize(n);
    mri_bsplines.resize(n, NULL);
    for (int i = 0; i < n; i++)
    {
      mri_mov[i] = MRIcopyFrame(mri4d, NULL, i, 0);
      if (m>0)
      {
        cout << "masking source frame " << i << " ..." << endl;
        MRI *mask = MRIcopyFrame(mask4d, NULL, i, 0);
        mri_mov[i] = MRImaskZero(mri_mov[i], mask, mri_mov[i]);
        MRIfree(&mask);
      }
      if (sampletype == SAMPLE_CUBIC_BSPLINE)
      {
        cout << "converting source frame " << i << " to bspline ..." << endl;
        mri_bsplines[i] = MRItoBSpline(mri_mov[i], mri_bsplines[i], 3);
      }
      ostringstream oss;
      oss << "Frame_" << i;
      mov[i] = oss.str();
    }
    MRIfree(&mri4d);
    MRIfree(&mask4d);
  }
  else //several inputs;
  {
    assert(n>1);
    mri_mov.resize(n);
    mri_bsplines.resize(n, NULL);
    vector<double> msize(pmov.size());
    for (int i = 0; i < n; i++)
    {
      cout << "reading source '" << mov[i] << "'..." << endl;

      mri_mov[i] = MRIread(mov[i].c_str());
      if (!mri_mov[i])
      {
        ErrorExit(ERROR_NOFILE,
            "MultiRegistration::loadMovables: could not open input volume %s.\n",
            mov[i].c_str());
      }
      if (mri_mov[i]->nframes != 1)
      {
        ErrorExit(ERROR_NOFILE,
            "MultiRegistration::loadMovables: only pass single frame MRI %s.\n",
            mov[i].c_str());
      }
      msize[i] = mri_mov[i]->xsize;
      if (mri_mov[i]->ysize < msize[i])
        msize[i] = mri_mov[i]->ysize;
      if (mri_mov[i]->zsize < msize[i])
        msize[i] = mri_mov[i]->zsize;

      if (m>0)
      {
        cout << "masking source frame " << i << " ..." << endl;
        MRI *mask = MRIread(pmasks[i].c_str());
        if (!mask)
        {
          ErrorExit(ERROR_NOFILE,
            "MultiRegistration::loadMovables: could not open mask volume %s.\n",
            pmasks[i].c_str());
        }
        if (mask->nframes != 1)
        {
          ErrorExit(ERROR_NOFILE,
            "MultiRegistration::loadMovables: only pass single frame mask %s.\n",
            pmasks[i].c_str());
        }
        mri_mov[i] = MRImaskZero(mri_mov[i], mask, mri_mov[i]);
        MRIfree(&mask);
      }


      if (sampletype == SAMPLE_CUBIC_BSPLINE)
      {
        cout << "converting source '" << mov[i] << "' to bspline ..." << endl;
        mri_bsplines[i] = MRItoBSpline(mri_mov[i], mri_bsplines[i], 3);
      }

    }
    // user can use setResthresh() to specify EPS. the default is 0.001
    float EPS = resthresh;
    for (int i = 1; i < n; i++)
    {
      if (fabs(mri_mov[i]->xsize - mri_mov[0]->xsize) > EPS
          || fabs(mri_mov[i]->ysize - mri_mov[0]->ysize) > EPS
          || fabs(mri_mov[i]->zsize - mri_mov[0]->zsize) > EPS)
      {
        cerr
            << "ERROR: MultiRegistration::loadMovables: images have different voxel sizes.\n";
        cerr << "  Currently not supported, maybe first make conform?\n";
        cerr << "  Debug info: size(" << i << ") = " << mri_mov[i]->xsize
            << ", " << mri_mov[i]->ysize << ", " << mri_mov[i]->zsize
            << "   size(0) = " << mri_mov[0]->xsize << ", " << mri_mov[0]->ysize
            << ", " << mri_mov[0]->zsize << endl;
        ErrorExit(ERROR_BADFILE,
            "MultiRegistration::loadMovables: voxel size is different %s.\n",
            mov[i].c_str());
      }
    }
  }

  mri_warps.resize(n, NULL);
  intensities.resize(n, 1.0);
  ltas.resize(n, NULL);
  if (robust)
    mri_weights.resize(n, NULL);

  return mri_mov.size();
}

int MultiRegistration::loadLTAs(const std::vector<std::string> nltas)
{
  assert(nltas.size() == mov.size());
  assert(ltas.size() == mov.size());
  iltas = nltas; // copy of input filenames
  for (uint i = 0; i < nltas.size(); i++)
  {
    TRANSFORM * trans = TransformRead(nltas[i].c_str());
    LTA* lta = (LTA *) trans->xform;
    if (!lta)
      ErrorExit(ERROR_BADFILE,
          "MultiRegistration::loadLTAs could not read transform file %s",
          nltas[i].c_str());
    if (nltas[i] == "identity.nofile" && mri_mov[i] != NULL)
      LTAmodifySrcDstGeom(lta, mri_mov[i], mri_mov[i]);
    if (!lta->xforms[0].src.valid)
      ErrorExit(ERROR_BADFILE,
          "MultiRegistration::loadLTAs no source geometry, use lta with valid geometry ( %s )",
          nltas[i].c_str());
    if (!lta->xforms[0].dst.valid)
      ErrorExit(ERROR_BADFILE,
          "MultiRegistration::loadLTAs no target geometry, use lta with valid geometry ( %s )",
          nltas[i].c_str());

    lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);

    ltas[i] = lta;
  }
  return (int) ltas.size();
}

int MultiRegistration::loadIntensities(const std::vector<std::string> nintens)
{
  assert(nintens.size() == mov.size());
  assert(intensities.size() == mov.size());
  iintens = nintens; // copy of input filenames
  for (uint i = 0; i < nintens.size(); i++)
  {
    ifstream f(nintens[i].c_str(), ios::in);
    if (f.good())
    {
      f >> intensities[i];
      f.close();
    }
    else
    {
      ErrorExit(ERROR_BADFILE,
          "MultiRegistration::loadIntensities no such file ( %s )",
          nintens[i].c_str());
    };
  }
  return (int) intensities.size();
}

/*!
 \brief Initializes a Registration with Parameters (rigid, iscale, transonly, robust, sat and doubleprec)
 \param R  Registration to be initialized
 */
void MultiRegistration::initRegistration(RegRobust & R)
{
  // assert(n < (int) P.mov.size());

  if (rigid)
  {
    R.setRigid();
  }
  else
    R.setAffine();
    
  R.setIscale(iscale);
  if (transonly)
    R.setTransonly();
  if (iscaleonly)
    R.setIscaleOnly();
    
  R.setCost(Registration::ROB);
  R.setSaturation(sat);
  R.setDoublePrec(doubleprec);
  //R.setDebug(debug);

  if (subsamplesize > 0)
    R.setSubsampleSize(subsamplesize);
  if (highit >= 0)
    R.setHighit(highit);

//   int pos = P.mov[n].rfind(".");
//   if (pos > 0) R.setName(P.mov[n].substr(0,pos));
//   else  R.setName(P.mov[n]);
//
//  // if (P.subsamplesize > 0) R.setSubsampleSize(P.subsamplesize);
//
//
//   R.setSource(P.mri_mov[n],P.fixvoxel,P.keeptype);
//   R.setTarget(P.mri_mean,P.fixvoxel,P.keeptype);
}

/*!
 \fn void mapAndAverageMov(int itdebug)
 \brief  maps movables to template using lta's, adjusts intensities (if iscale) and creates average (mean,median)
 \param idebug  specify iteration number for debug output
 */
bool MultiRegistration::mapAndAverageMov(int itdebug)
// maps mov to template (using ltas)
// adjust intensities (if iscale)
// creates template acording to average:1 mean, 2 median..
{
  unsigned int nin = mri_mov.size();
  assert(nin > 1);
  assert(mri_warps.size() == nin);
  assert(ltas.size() == nin);
  assert(intensities.size() == nin);
  
  if (iscaleonly)
  {
    cout << "creating intitial template (iscale only)" << endl;
    for (unsigned int i = 0; i < nin; i++)
    {
      mri_warps[i] = MyMRI::MRIvalscale(mri_mov[i], mri_warps[i],
          intensities[i]);
      if (debug)
      {
        ostringstream oss;
        oss << outdir << "tp" << i + 1 << "_to_template-it" << itdebug << ".mgz";
        MRIwrite(mri_warps[i], oss.str().c_str());
      }
    }
  }
  else
  {
    cout << "mapping movs and creating initial template..." << endl;
    if (iscale)
      cout << " allow intensity scaling" << endl;
    for (unsigned int i = 0; i < nin; i++)
    {
      if (mri_warps[i])
        MRIfree(&mri_warps[i]);
      // use geometry from ltas
      // (if initXforms was called, this is the center of mass of all tps)
      if (sampletype == SAMPLE_CUBIC_BSPLINE)
        mri_warps[i] = LTAtransformBSpline(mri_bsplines[i], NULL, ltas[i]);
      else
        mri_warps[i] = LTAtransformInterp(mri_mov[i], NULL, ltas[i], sampletype);
      MRIcopyPulseParameters(mri_mov[i], mri_warps[i]);
      if (iscale)
      {
        mri_warps[i] = MyMRI::MRIvalscale(mri_warps[i], mri_warps[i],
            intensities[i]);
      }
      if (debug)
      {
        ostringstream oss;
        oss << outdir << "tp" << i + 1 << "_to_template-it" << itdebug << ".mgz";
        MRIwrite(mri_warps[i], oss.str().c_str());
      }
    }
  }
  mri_mean = averageSet(mri_warps, mri_mean, average, sat);
  return true;
}

MRI* MultiRegistration::averageConformSet(int itdebug)
// maps mov to template (using ltas)
// adjust intensities (if iscale)
// creates template acording to average:1 mean, 2 median..
{
  unsigned int nin = mri_mov.size();
  assert(nin > 1);
  assert(ltas.size() == nin);
  assert(intensities.size() == nin);

  std::vector<MRI*> mri_cwarps(nin, NULL);

  cout << "warping movs and creating conform template..." << endl;
  if (iscale)
    cout << " allow intensity scaling" << endl;
  for (unsigned int i = 0; i < nin; i++)
  {

    // use geometry from ltas
    // (if initXforms was called, this is the center of mass of all tps)
    mri_cwarps[i] = LTAtransform(mri_mov[i], NULL, ltas[i]);
    MRIcopyPulseParameters(mri_mov[i], mri_warps[i]);
    //????
    cerr << " MultiRegistration::averageConformSet not implemented, yet !!!"
        << endl;
    assert(1==2);

    if (iscale)
    {
      mri_cwarps[i] = MyMRI::MRIvalscale(mri_cwarps[i], mri_cwarps[i],
          intensities[i]);
    }
    if (debug)
    {
      ostringstream oss;
      oss << outdir << "tp" << i + 1 << "_to_template_conform-it" << itdebug
          << ".mgz";
      MRIwrite(mri_cwarps[i], oss.str().c_str());
    }
  }
  MRI* mri_cmean = averageSet(mri_cwarps, NULL, average, sat);
  return mri_cmean;
}

bool MultiRegistration::computeTemplate(int itmax, double eps, int iterate,
    double epsit)
// itmax : iterations for template creation
// iterate: iterations for registration
// epsit: 
// uses ltas (if set), mri_mov and parameters
// sets mri_mean, ltas, mri_warps (mri_weights and intensities if available)
{
  int nin = (int) mri_mov.size();
  assert(nin > 1);

  if (debug)
    cout << endl << endl << "MultiRegistration::computeTemplate(avitmax: "
        << itmax << ", aveps " << eps << ", regit: " << iterate << ", regeps: "
        << epsit << " )" << endl;

  cout << endl << "Computing first template" << endl;
  bool havexforms = (ltas[0] != NULL);
  if (!havexforms) // create simple initial average aligning centers (moments), blurry!
    initialAverageSet();
  else
    // we have initial transforms
    mapAndAverageMov(0);

  if (itmax == 0) // no iterations necessary just return with mean (and ltas and warps);
  {
    //strncpy(P.mri_mean->fname, P.mean.c_str(),STRLEN);
    return true;
  }

  strncpy(mri_mean->fname, (outdir + "template-it0.mgz").c_str(), STRLEN-1);
  if (debug)
  {
    cout << "debug: saving template-it0.mgz" << endl;
    MRIwrite(mri_mean, (outdir + "template-it0.mgz").c_str());
  }

  //cout << "template fname: " << mri_mean->fname << endl;

//  int itmax  = 10;
//  double eps = 0.025;

  int itcount = 0;
  double maxchange = 100;

//  vector < Registration > Rv(mri_mov.size()); //(removed to save mem)
//  for (int i = 0;i<nin;i++) Rv[i].setSource(mri_mov[i],fixvoxel,keeptype);

  vector<vnl_matrix_fixed<double, 4, 4> > transforms(mri_mov.size());
//  if (havexforms) //initial transforms
  {
    for (int i = 0; i < nin; i++)
    {
      // first change to Ras2Ras, then swap geometries
      // (this is important only, if ixforms were passed AND the geometries in the input lta
      //   differ from the source and target mean space specified on the command line):
      // if we constructed the initial xforms, geometries should already be accurate
      ltas[i] = LTAchangeType(ltas[i], LINEAR_RAS_TO_RAS);
      LTAmodifySrcDstGeom(ltas[i], mri_mov[i], mri_mean);
      ltas[i] = LTAchangeType(ltas[i], LINEAR_VOX_TO_VOX); //vox2vox for registration  

      // the transforms are maps from the original images
      // (not resampled to conform as inside the Registration)
      transforms[i] = MyMatrix::convertMATRIX2VNL(ltas[i]->xforms[0].m_L);
    }

  }

  // if we do not have good transforms, run special treatement
  // below on different resolutions, here we determine how often
  // these lowres registrations are run.
  // the methods are: maxit 3, maxit 2, maxit 1, subsample 180
  int noxformits[4] =
  { 3, 1, 0, 0 };
  while (itcount < itmax && maxchange > eps)
  {
    itcount++;

    // if we have prior ltas, allow termination
    // by initializing maxchange =0 so that it can be 
    // computed below
    if (ltas[0])
      maxchange = 0;

    cout << endl << "======================================================="
        << endl;
    cout << endl << "Working on global iteration : " << itcount << endl;
    cout << endl << "======================================================="
        << endl;
    //cout << "========================================" << endl;
    //printmemusage();
    //cout << "========================================" << endl << endl;

    if (itcount <= 4 && noxformits[itcount - 1] > 0)
      cout << "  noxformits = " << noxformits[itcount - 1] << endl;

    // register all inputs to mean
    vector<double> dists(nin, 1000); // should be larger than maxchange!
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static,1)
#endif
    for (int i = 0; i < nin; i++)
    {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
      cout << endl << "Working on TP " << i + 1 << endl << endl;
      RegRobust R; // create new registration each time to keep mem usage smaller
//      Rv[i].clear();
      R.setVerbose(0);
      initRegistration(R); //set parameters
//      R.setSource(mri_mov[i], fixvoxel, keeptype);
//      R.setTarget(mri_mean, fixvoxel, keeptype); // gaussian pyramid will be constructed for
//                                                 // each Rv[i], could be optimized
      R.setSourceAndTarget(mri_mov[i],mri_mean,keeptype);

      ostringstream oss;
      oss << outdir << "tp" << i + 1 << "_to_template-it" << itcount;
      R.setName(oss.str());

      // compute Alignment
      //std::pair <MATRIX*, double> Md;
      std::pair<vnl_matrix_fixed<double, 4, 4>, double> Md;
      //int iterate = P.iterate;
      //double epsit= P.epsit;
      int maxres = 0;
      // on higher iterations use subsamplesize as passed on commandline
      int subsamp = subsamplesize;
      // simplify first steps (only if we do not have good transforms):
      if (!iscaleonly && !havexforms && itcount <= 4 && noxformits[itcount - 1] > 0)
      {
        switch (itcount)
        {
        case 1:
          maxres = 3;
          break;
        case 2:
          maxres = 2;
          break;
        case 3:
          maxres = 1;
          break;
        case 4:
          subsamp = 180;
          break;
        }
      }

      R.setSubsampleSize(subsamp);
      R.setIscaleInit(intensities[i]);
      R.setMinitOrig(transforms[i]); // as the transforms are in the original space
      if (satit)
        R.findSaturation();

#ifdef HAVE_OPENMP
#pragma omp critical
#endif 
      if (nomulti || iscaleonly)
      {
        cout << " - running high-res registration on TP " << i + 1 << "..." << endl;
        R.computeIterativeRegistration(iterate, epsit); 
      }
      else
      {
        cout << " - running multi-resolutional registration on TP " << i + 1 << "..." << endl;
        R.computeMultiresRegistration(maxres, iterate, epsit);
      }

      Md.first = R.getFinalVox2Vox();
      Md.second = R.getFinalIscale();
      if (!R.getConverged())
      {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
        cout << "   *** WARNING: TP " << i + 1
            << " to template did not converge ***" << endl;
      }
       
      transforms[i] = Md.first;
      intensities[i] = Md.second;

      // convert Matrix to LTA ras to ras
      LTA * lastlta = NULL;
      if (ltas[i])
        lastlta = ltas[i];
      ltas[i] = MyMatrix::VOXmatrix2LTA(Md.first, mri_mov[i], mri_mean);
      //P.ltas[i] = LTAalloc(1,P.mri_mov[i]);
      //P.ltas[i]->xforms[0].m_L =
      //  MRIvoxelXformToRasXform (P.mri_mov[i],
      //                           P.mri_mean,
      //                           Md.first,
      //                           P.ltas[i]->xforms[0].m_L) ;
      //P.ltas[i]->type = LINEAR_RAS_TO_RAS ;
      //getVolGeom(P.mri_mov[i], &P.ltas[i]->xforms[0].src);
      //getVolGeom(P.mri_mean, &P.ltas[i]->xforms[0].dst);

      // compute maxchange
      if (lastlta)
      {
        LTAchangeType(lastlta, LINEAR_RAS_TO_RAS); //measure dist in RAS coords
        LTAchangeType(ltas[i], LINEAR_RAS_TO_RAS); //measure dist in RAS coords
        dists[i] = sqrt(
            MyMatrix::AffineTransDistSq(lastlta->xforms[0].m_L,
                ltas[i]->xforms[0].m_L));
        LTAfree(&lastlta);
        if (dists[i] > maxchange)
          maxchange = dists[i];
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
        cout << "   tp " << i + 1 << " distance: " << dists[i] << endl;
      }

      // create warps: warp mov to mean
      if (mri_warps[i])
        MRIfree(&mri_warps[i]);
      mri_warps[i] = MRIclone(mri_mean, mri_warps[i]);
      if (sampletype == SAMPLE_CUBIC_BSPLINE)
      {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
        cout << " - mapping tp " << i + 1 << " to template (cubic bspline) ..."
            << endl;
        mri_warps[i] = LTAtransformBSpline(mri_bsplines[i], mri_warps[i],
            ltas[i]);
      }
      else
      {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
        cout << " - mapping tp " << i + 1 << " to template..." << endl;
        mri_warps[i] = LTAtransformInterp(mri_mov[i], mri_warps[i], ltas[i],
            sampletype);
      }
      MRIcopyPulseParameters(mri_mov[i], mri_warps[i]);

      // here do scaling of intensity values
      if (R.isIscale() && Md.second > 0)
      {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
        cout << " - adjusting intensity of mapped tp " << i + 1 << " by "
            << Md.second << endl;
        mri_warps[i] = MyMRI::MRIvalscale(mri_warps[i], mri_warps[i],
            Md.second);
      }

      if (backupweights)
      {
        // copy weights (as RV will be cleared)
        //   (info: they are in original half way space)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
        cout << " - backup weights tp " << i + 1 << " ..." << endl;
        //if (mri_weights[i]) MRIfree(&mri_weights[i]); 
        mri_weights[i] = MRIcopy(R.getWeights(), mri_weights[i]);
      }

      //cout << " LS difference after: " <<
      //CF.leastSquares(mri_aligned,P.mri_dst) << endl;
      //cout << " NC difference after: " <<
      //CF.normalizedCorrelation(mri_aligned,P.mri_dst) << endl;

      if (debug)
      {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
        cout << " - debug tp " << i + 1 << " : writing transforms, warps, weights ..." << endl;

        LTAwriteEx(ltas[i], (oss.str() + ".lta").c_str());

        MRIwrite(mri_warps[i], (oss.str() + ".mgz").c_str());

        if (R.isIscale() && Md.second > 0)
        {
          string fn = oss.str() + "-intensity.txt";
          ofstream f(fn.c_str(), ios::out);
          f << Md.second;
          f.close();
        }

        // if we have weights:  
        if (mri_weights[i] != NULL)
        {
          std::pair<vnl_matrix_fixed<double, 4, 4>,
              vnl_matrix_fixed<double, 4, 4> > map2weights = R.getHalfWayMaps();
          vnl_matrix_fixed<double, 4, 4> hinv = vnl_inverse(map2weights.second);

#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
          {
            cout << endl;
            cout << map2weights.first << endl;
            cout << endl;
            cout << map2weights.second << endl;
            cout << endl;
            cout << hinv << endl;
            cout << endl;
          }
          MRI * wtarg = MRIallocSequence(mri_weights[i]->width, mri_weights[i]->height,
              mri_weights[i]->depth, MRI_FLOAT,mri_weights[i]->nframes);
          MRIcopyHeader(mri_weights[i], wtarg);
          MATRIX * v2r = MRIgetVoxelToRasXform(mri_mean);
          MRIsetVoxelToRasXform(wtarg, v2r);
          wtarg->type = MRI_FLOAT;
          wtarg->i_to_r__ = AffineMatrixCopy(mri_mean->i_to_r__,
              wtarg->i_to_r__);
          wtarg->r_to_i__ = MatrixCopy(mri_mean->r_to_i__, wtarg->r_to_i__);

          wtarg = MyMRI::MRIlinearTransform(mri_weights[i], wtarg, hinv);
          MRIwrite(wtarg, (oss.str() + "-weights.mgz").c_str());
          MRIwrite(mri_weights[i], (oss.str() + "-www.mgz").c_str());
          MRIfree(&wtarg);
          //MatrixFree(&hinv);
          MatrixFree(&v2r);
        }
      } // if debug end

      // clear to reduce memory usage:
      //Rv[i].clear();
      // Rv[i].freeGPT();

#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
      {
        cout << endl << "Finished TP : " << i + 1 << endl;
        cout << endl;
        cout << "=====================================================" << endl;
      }

    } // for loop end (all timepoints)

    // if we did not have initial transforms
    // allow for more iterations on different resolutions
    // based on noxformits vector defined above
    if (!havexforms && itcount <= 4 && noxformits[itcount - 1] > 1)
    {
      noxformits[itcount - 1]--;
      itcount--;
    }

    if (dists[0] <= maxchange) // it was computed, so print it:
    {
      cout << endl << "Global Iteration " << itcount << " Distances: " << endl;
      for (unsigned int i = 0; i < dists.size(); i++)
        cout << dists[i] << " ";
      cout << endl << "Maxchange: " << maxchange << endl << endl;
    }

    // create new mean
    cout << "Computing new template " << itcount << endl;
    mri_mean = averageSet(mri_warps, mri_mean, average, sat);
    if (debug)
    {
      ostringstream oss;
      oss << outdir << "template-it" << itcount << ".mgz";
      cout << "debug: writing template to " << oss.str() << endl;
      MRIwrite(mri_mean, oss.str().c_str());
      //strncpy(mri_mean->fname, oss.str().c_str(),STRLEN);
    }

  } // end while

  //strncpy(P.mri_mean->fname, P.mean.c_str(),STRLEN);

  cout << " DONE : computeTemplate " << endl;
  return true;
}

bool MultiRegistration::halfWayTemplate(int maxres, int iterate, double epsit,
    bool vox2vox)
// can be used only with two input images
{

  cerr << " Error MultiRegistration::halfWayTemplate is deprecated " << endl;
  exit(1);

  int nin = (int) mri_mov.size();
  if (nin != 2)
  {
    cerr << "Error, need 2 movs" << endl;
    exit(1);
  }

  // register 1 with 2

  RegRobust R;
  initRegistration(R); //set parameter
  R.setSourceAndTarget(mri_mov[0], mri_mov[1], keeptype);

  ostringstream oss;
  oss << outdir << "halfway_template.mgz";
  R.setName(oss.str().c_str());

  // use initial transform if given:
  //MATRIX *minit = NULL;
  vnl_matrix_fixed<double, 4, 4> minit;
  assert(ltas.size() ==2);
  if (ltas[0])
  {
    LTAchangeType(ltas[0], LINEAR_VOX_TO_VOX); //vox2vox for registration
    LTAchangeType(ltas[1], LINEAR_VOX_TO_VOX); //vox2vox for registration
    assert(ltas[0]->type == LINEAR_VOX_TO_VOX);
    assert(ltas[1]->type == LINEAR_VOX_TO_VOX);
    minit = vnl_inverse(MyMatrix::convertMATRIX2VNL(ltas[1]->xforms[0].m_L))
        * MyMatrix::convertMATRIX2VNL(ltas[0]->xforms[0].m_L);
    R.setDebug(1);
  }

  // compute Alignment
  //std::pair <MATRIX*, double> Md;
  std::pair<vnl_matrix_fixed<double, 4, 4>, double> Md;
  // adjust subsamplesize, if passed:
  if (subsamplesize > 0)
    R.setSubsampleSize(subsamplesize);
  if (satit)
    R.findSaturation();

  //!!!! what if iscale init was passed? needs fixing, if this is used at all?
  R.setMinitOrig(minit);
  if (nomulti)
  {
    R.computeIterativeRegistration(iterate, epsit); 
  }
  else
  {
    R.computeMultiresRegistration(maxres, iterate, epsit);
  }
  Md.first = R.getFinalVox2Vox();
  Md.second = R.getFinalIscale();
  //if (minit) MatrixFree(&minit);
  intensities[0] = Md.second;
  intensities[1] = 1;

  cout << "creating half way data ..." << endl;
  assert(MyMRI::isIsotropic(mri_mov[0]) && MyMRI::isIsotropic(mri_mov[1]));

  std::pair<vnl_matrix_fixed<double, 4, 4>, vnl_matrix_fixed<double, 4, 4> > maps2weights =
      R.getHalfWayMaps();
  LTA * m2hwlta = LTAalloc(1, mri_mov[0]);
  LTA * d2hwlta = LTAalloc(1, mri_mov[1]);
  if (!vox2vox) // do ras to ras
  {
    // cout << "converting VOX to RAS and saving RAS2RAS..." << endl ;
    // (use geometry of destination space for half-way) WIll NOT WORK FOR nonistorpic due to internal resampling
    m2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(maps2weights.first);
    m2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform(mri_mov[0], mri_mov[1],
        m2hwlta->xforms[0].m_L, m2hwlta->xforms[0].m_L);
    m2hwlta->type = LINEAR_RAS_TO_RAS;

    d2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(maps2weights.second);
    d2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform(mri_mov[1], mri_mov[1],
        d2hwlta->xforms[0].m_L, d2hwlta->xforms[0].m_L);
    d2hwlta->type = LINEAR_RAS_TO_RAS;
  }
  else // vox to vox
  {
    // cout << "saving VOX2VOX..." << endl ;
    m2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(maps2weights.first,
        m2hwlta->xforms[0].m_L);
    m2hwlta->type = LINEAR_VOX_TO_VOX;
    d2hwlta->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(maps2weights.second,
        m2hwlta->xforms[0].m_L);
    d2hwlta->type = LINEAR_VOX_TO_VOX;
  }
  // add src and dst info (use src as target geometry in both cases)
  getVolGeom(mri_mov[0], &m2hwlta->xforms[0].src);
  getVolGeom(mri_mov[0], &m2hwlta->xforms[0].dst);
  getVolGeom(mri_mov[1], &d2hwlta->xforms[0].src);
  getVolGeom(mri_mov[0], &d2hwlta->xforms[0].dst);

  ltas.resize(2, NULL);
  if (ltas[0])
    LTAfree(&ltas[0]);
  if (ltas[1])
    LTAfree(&ltas[1]);
  ltas[0] = m2hwlta;
  ltas[1] = d2hwlta;

  cout << " creating half-way template ..." << endl;
  // take dst info from lta:
  if (mri_warps[0])
    MRIfree(&mri_warps[0]);
  if (mri_warps[1])
    MRIfree(&mri_warps[1]);
  mri_warps[0] = LTAtransform(mri_mov[0], NULL, m2hwlta);
  MRIcopyPulseParameters(mri_mov[0], mri_warps[0]);
  mri_warps[1] = LTAtransform(mri_mov[1], NULL, d2hwlta);
  MRIcopyPulseParameters(mri_mov[1], mri_warps[1]);
  // here do scaling of intensity values
  if (R.isIscale() && Md.second > 0)
  {
    cout << "Adjusting Intensity of MOV by " << Md.second << endl;
    mri_warps[0] = MyMRI::MRIvalscale(mri_warps[0], mri_warps[0], Md.second);
  }
  mri_mean = averageSet(mri_warps, mri_mean, average, sat);

  mri_weights[0] = R.getWeights(); // in original half way space
  mri_weights[1] = mri_weights[0];

  // copy weights for both if we have them (as the R->weights will be destroyed):
  if (mri_weights[0])
    for (unsigned int i = 0; i < mri_weights.size(); i++)
      mri_weights[i] = MRIcopy(mri_weights[0], NULL);

  //MatrixFree(&Md.first);
  return true;
}

/** Flips/reorders images to be in LIA orientation (conform default)
   (future: maybe use the orientation of most of the inputs )
   then computes average cosine matrix. 
   mri_mov images need to be set.
 */
vnl_matrix_fixed<double, 3, 3> MultiRegistration::getAverageCosines()
{
  int nin = (int) mri_mov.size();

  std::vector <  vnl_matrix_fixed<double, 3, 3> > cosM(nin);
  vnl_matrix_fixed<double, 3, 3> meanr (0.0);

//  // conform slice orientation
//  // orientation: LIA
//  // primary slice direction: coronal
//  vnl_matrix_fixed<double, 3, 3> v2rconf (0.0);
//  v2rconf[0][0] = -1.0;
//  v2rconf[1][0] =  0.0;
//  v2rconf[2][0] =  0.0;
//  v2rconf[0][1] =  0.0;
//  v2rconf[1][1] =  0.0;
//  v2rconf[2][1] = -1.0;
//  v2rconf[0][2] =  0.0;
//  v2rconf[1][2] =  1.0;
//  v2rconf[2][2] =  0.0;
  
  // extract cosines
  bool same = true;
  double eps = 0.000000001;
  for (int i = 0; i<nin ; i++)
  {
    cosM[i][0][0] = mri_mov[i]->x_r;
    cosM[i][1][0] = mri_mov[i]->x_a;
    cosM[i][2][0] = mri_mov[i]->x_s;
    cosM[i][0][1] = mri_mov[i]->y_r;
    cosM[i][1][1] = mri_mov[i]->y_a;
    cosM[i][2][1] = mri_mov[i]->y_s;
    cosM[i][0][2] = mri_mov[i]->z_r;
    cosM[i][1][2] = mri_mov[i]->z_a;
    cosM[i][2][2] = mri_mov[i]->z_s;
    
    // for i>0 compare to first if all cos are the same
    if (i>0)
    {
      for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
        if (fabs(cosM[0][j][k] - cosM[i][j][k]) > eps)
        {
          same = false;
          std::cout << " Cosines of input transforms not identical !" << std::endl;
        }
    }
  }

  // if all are the same, return source cousine
  if (same)
    return cosM[0];

  // if all cos not the same, check reorder
  //bool reorder = false;
//  if (! same) for (int i = 0; i<nin ; i++)
  meanr = cosM[0];
  for (int i = 1; i<nin ; i++)
  {
//    // reorder axis to match conform ordering
    // we only swap dimensions around to make the cosine matrices comparable
    // (same main slice orientation). We will later average the cosines of all inputs
//    vnl_matrix_fixed<double, 3, 3> v2v = v2rconf.transpose() * cosM[i];
    vnl_matrix_fixed<double, 3, 3> v2v = cosM[0].transpose() * cosM[i];

    //cout << " v2v[" << i << "] = " << endl;
    //vnl_matlab_print(vcl_cout,v2v,"v2v",vnl_matlab_print_format_long);
    //cout << endl;
      
    // swap (and possibly invert) axis, so that rotation gets smaller
    int xd = 1;
    int yd = 1;
    int zd = 1;
      
    // determine max in each column:
    if (fabs(v2v[1][0]) > fabs(v2v[0][0]))
      xd = 2;
    if (fabs(v2v[2][0]) > fabs(v2v[0][0])
        && fabs(v2v[2][0]) > fabs(v2v[1][0]))
      xd = 3;
    if (fabs(v2v[1][1]) > fabs(v2v[0][1]))
      yd = 2;
    if (fabs(v2v[2][1]) > fabs(v2v[0][1])
        && fabs(v2v[2][1]) > fabs(v2v[1][1]))
      yd = 3;
    if (fabs(v2v[1][2]) > fabs(v2v[0][2]))
      zd = 2;
    if (fabs(v2v[2][2]) > fabs(v2v[0][2])
        && fabs(v2v[2][2]) > fabs(v2v[1][2]))
      zd = 3;
        
    // sign
    if (v2v[xd - 1][0] < 0.0)
      xd = -xd;
    if (v2v[yd - 1][1] < 0.0)
      yd = -yd;
    if (v2v[zd - 1][2] < 0.0)
      zd = -zd;
        
    // now reorder if necessary  
    if (xd != 1 || yd != 2 || zd != 3)
    {

      //reorder = true;
      if (abs(xd) + abs(yd) + abs(zd) != 6)
      {
        cout << "WARNING: reorder not clear ..." << endl;
        vnl_matlab_print(vcl_cout,v2v,"v2v",vnl_matlab_print_format_long);
        cout << endl;
        cout << " xd: " << xd << " yd: " << yd << " zd: " << zd << endl;
        if (vnl_determinant(v2v) < 0)
        {  // cannot run sqrt later if det < 0
          cout << "ERROR: vox2vox det: " << vnl_determinant(v2v) << " < 0"
              << endl;
          cout << "       Something might be wrong with RAS info in inputs."
              << endl;
          cout << "       Make sure volumes are in same voxel orientation." << endl;
          exit(1);
        }
        exit(1);
      }

      //cout << "   Reordering axes in mov "<< i << " to better fit first mov... ( " << xd << " " << yd
      //    << " " << zd << " )" << endl;

      vnl_matrix_fixed<double, 3, 3> Q(0.0);
      Q[0][abs(xd)-1] = (xd >= 0 ? 1 : -1);
      Q[1][abs(yd)-1] = (yd >= 0 ? 1 : -1);
      Q[2][abs(zd)-1] = (zd >= 0 ? 1 : -1);

      cosM[i] = cosM[i] * Q;

    }// end reorder
 
    //cout << " v2r[" << i << "] = " << endl;
    //vnl_matlab_print(vcl_cout,cosM[i],"cosM",vnl_matlab_print_format_long);
    //cout << endl;

    // once axes are all ordered the same, we can average cosines
    // by computing the mean and projecting it back to rotation space
    meanr += cosM[i];

  }

  // Project meanr back to SO(3) via polar decomposition:  
  meanr = (1.0 / nin) * meanr;
  //vnl_matlab_print(vcl_cout,meanr,"meanr",vnl_matlab_print_format_long);std::cout << std::endl;
  vnl_matrix<double> PolR(3, 3), PolS(3, 3);
  MyMatrix::PolarDecomposition(meanr, PolR, PolS);
  meanr = PolR;
  //vnl_matlab_print(vcl_cout,meanr,"meanrfinal",vnl_matlab_print_format_long);std::cout << std::endl;
  
  return meanr;
}


/** Creates the template image geometry.
   Sets template to isotropic voxels (select the largest over all images
   of the minimum side length). Then determines the widht height depth and
   sets all dimensions to the max (can treat 2D images with depth=1).
   Then sets average RAS cosines and sets the center to the average RAS center.
   Will return slices in COR (see getAverageCosines).
   mri_mov images need to be set.
 */
MRI * MultiRegistration::createTemplateGeo()
{
  int nin = (int) mri_mov.size();

  // determine isotropic voxel size for template
  // set to the largest (across images) of the
  // smallest voxel size
  double conform_size = 0.0;
  for (int i = 0; i<nin ; i++)
  {
    double xsize = fabs(mri_mov[i]->xsize);
    double ysize = fabs(mri_mov[i]->ysize);
    double zsize = fabs(mri_mov[i]->zsize);
    
    double minvox = xsize;
    if (ysize < minvox) minvox = ysize;
    if (zsize < minvox) minvox = zsize;
      
    if (minvox > conform_size) conform_size = minvox;
  }

  // determine width height depth
  // since inputs can be in different
  // slice orientations, set all to the max
  int maxdimw = 0;
  int maxdimh = 0;
  int maxdimd = 0;
  int count2d = 0;
  for (int i = 0; i<nin ; i++)
  {
    double xsize   = fabs(mri_mov[i]->xsize);
    double ysize   = fabs(mri_mov[i]->ysize);
    double zsize   = fabs(mri_mov[i]->zsize);
    double fwidth  = xsize * mri_mov[i]->width;
    double fheight = ysize * mri_mov[i]->height;
    double fdepth  = zsize * mri_mov[i]->depth;
    
    double eps =0.0001; // to prevent ceil(2.0*64 / 2.0) = ceil(64.000000000001) = 65
    int nw = (int) ceil((fwidth  / conform_size) - eps);
    int nh = (int) ceil((fheight / conform_size) - eps);
    int nd = (int) ceil((fdepth  / conform_size) - eps);
    if (mri_mov[i]->depth == 1)
      count2d++; // 2d image
    
    if (nw > maxdimw) maxdimw = nw;
    if (nh > maxdimh) maxdimh = nh;
    if (nd > maxdimd) maxdimd = nd;
  }

//  int depth = maxdimd;
//  if (count2d > 0)
//  {
//    if (count2d == nin)
//      depth = 1;
//    else
    if (count2d > 0 && count2d != nin)
    {
      cerr << " ERROR (createTemplateGeo) : mixing 2d and 3d input images not possible!" << endl;
      exit(1);
    }
//  }
//  MRI * template_geo   = MRIallocHeader(maxdimw, maxdimh, depth, MRI_FLOAT , 1);
  MRI * template_geo   = MRIallocHeader(maxdimw, maxdimh, maxdimd, MRI_FLOAT , 1);
  template_geo->xsize  = template_geo->ysize = template_geo->zsize = conform_size ;
  template_geo->imnr0  = 1;
  template_geo->imnr1  = maxdimh;
  template_geo->type   = MRI_UCHAR;
  template_geo->thick  = conform_size;
  template_geo->ps     = conform_size;
  template_geo->xstart = template_geo->ystart = - maxdimw/2;
  template_geo->zstart = -maxdimd/2;
  template_geo->xend   = template_geo->yend = maxdimw/2;
  template_geo->zend   = maxdimd/2;
  
  // set direction cosines to default directions:
  //template_geo->x_r = -1.0;
  //template_geo->x_a =  0.0;
  //template_geo->x_s =  0.0;
  //template_geo->y_r =  0.0;
  //template_geo->y_a =  0.0;
  //template_geo->y_s = -1.0;
  //template_geo->z_r =  0.0;
  //template_geo->z_a =  1.0;
  //template_geo->z_s =  0.0;
  
  vnl_matrix_fixed<double, 3, 3> avgcos = getAverageCosines();
  template_geo->x_r = avgcos[0][0];
  template_geo->x_a = avgcos[1][0];
  template_geo->x_s = avgcos[2][0];
  template_geo->y_r = avgcos[0][1];
  template_geo->y_a = avgcos[1][1];
  template_geo->y_s = avgcos[2][1];
  template_geo->z_r = avgcos[0][2];
  template_geo->z_a = avgcos[1][2];
  template_geo->z_s = avgcos[2][2];
  

  // set center: use the average cras coords of the input images
  double nr =0.0;
  double na =0.0;
  double ns =0.0;
  for (int i = 0; i < nin; i++)
  {
      nr += mri_mov[i]->c_r;
      na += mri_mov[i]->c_a;
      ns += mri_mov[i]->c_s;
  }
  nr /= (double)nin;
  na /= (double)nin;
  ns /= (double)nin;

  // set new center in geometry
  template_geo->c_r = nr;
  template_geo->c_a = na;
  template_geo->c_s = ns;
  template_geo->ras_good_flag = 1;
  
  MRIreInitCache(template_geo);

  return template_geo;
}

/** Will set initial LTA's from each tp to template space.
    First registers all inputs to tpi (1..N).
    maxres, iterate and epsit are passed to those pairwise registrations.
    If fixtp is true, it stops there and creates average in that tp space.
    Else it computes average rotation and translation,
    constructs average template geometry and creates
    LTA's to that mid space. 
    Uses mri_mov and other parameters (e.g. satit, iscale, crascenter)
*/
bool MultiRegistration::initialXforms(int tpi, bool fixtp, int maxres,
    int iterate, double epsit)
// will set ltas (as RAS to RAS ???)
// tpi 1....n
// uses outdir (really???)
{
  if (tpi <= 0)
    return false; // no init?
  cout << endl << "MultiRegistration::initializing Xforms (init " << tpi
      << " , maxres " << maxres << " , iterate " << iterate << " , epsit "
      << epsit << " ) : " << endl;
  assert(ltas.size() == mri_mov.size());
  assert(mri_mov.size() > 1);
  if (debug) printParams();
  
  int nin = (int) mri_mov.size();

  tpi--; // tpi 0....n-1
  assert(tpi>=0 && tpi <nin);

  // create index set for lookup (tpi should switch places with 0)
  vector<int> index(nin);
  for (int i = 0; i < nin; i++)
    index[i] = i;
  index[0] = tpi;
  index[tpi] = 0;
  vector<bool> converged(nin, true);

  // Register everything to tpi TP
//  vector < Registration > Rv(nin);
  vector<std::pair<vnl_matrix_fixed<double, 4, 4>, double> > Md(nin);
//  vector < double >  centroid;
  vnl_vector_fixed<double, 4> centroid(0.0);

  //Md[0].first = MatrixIdentity(4,NULL);
  Md[0].first.set_identity();
  Md[0].second = 1.0;
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static,1)
#endif
  for (int i = 1; i < nin; i++)
  {
    int j = index[i]; // use new index
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
    {
      cout << endl << "[init] ========================= TP " << j + 1
          << " to TP " << tpi + 1 << " ==============================" << endl;
      cout << "         Register TP " << j + 1 << " ( " << mov[j] << " )"
          << endl;
      cout << "          to      TP " << tpi + 1 << " ( " << mov[tpi] << " )"
          << endl << endl;
    }

    ostringstream oss;
    oss << outdir << "tp" << j + 1 << "_to_tp" << tpi;

    RegRobust R;
    initRegistration(R); //set parameter
    //R.setRigid(); // stay rigid for averaging initial template space, after that allow affine
    if (debug) R.setVerbose(1);
    else R.setVerbose(0);
    R.setSourceAndTarget(mri_mov[j], mri_mov[tpi], keeptype);
    R.setName(oss.str());

    // compute Alignment (maxres,iterate,epsit) are passed above
    if (satit)
      R.findSaturation();
    if (nomulti)
    {
      R.computeIterativeRegistration(iterate, epsit); 
    }
    else
    {
      R.computeMultiresRegistration(maxres, iterate, epsit);
    }
    Md[i].first = R.getFinalVox2Vox();
    Md[i].second = R.getFinalIscale();
    converged[i] = R.getConverged();

    // the centroid is the voxel coord where the moment based centroid is located
    vnl_vector_fixed<double, 4> centroid_temp(R.getCentroidSinT());
    // add centroid of tpi (target of the registration)
    // only do this once (when i==1) is enough
    if (i == 1)
    {
      centroid_temp += R.getCentroidT();
      if (debug)
      {
        vnl_matlab_print(vcl_cout,R.getCentroidT(),"CentroidT",vnl_matlab_print_format_long);
        std::cout << std::endl;
      }
    }
    if(debug)
    {
      vnl_matlab_print(vcl_cout,R.getCentroidSinT(),"CentroidSinT",vnl_matlab_print_format_long);
      std::cout << std::endl;
    }
#ifdef HAVE_OPENMP
#pragma omp critical
#endif  
    {
      centroid += centroid_temp;
    }
  } // end for loop (initial registration to inittp)

  centroid = (1.0 / nin) * centroid;
  if (debug)
  {
    vnl_matlab_print(vcl_cout,centroid,"Centroid",vnl_matlab_print_format_long);
    std::cout << std::endl;
  }
  
  for (int i = 1; i < nin; i++)
  {
    if (!converged[i])
      cout << "* WARNING: TP " << index[i] + 1 << " to " << tpi + 1
          << " did not converge !! " << endl;
  }

  // copy results in correct order back to global members
  // lta (transforms) and intensities:
  assert(ltas.size() == mri_mov.size());
  for (int i = 0; i < nin; i++)
  {
    int j = index[i]; // use new index

    assert(ltas[j] == NULL);
    ltas[j] = MyMatrix::VOXmatrix2LTA(Md[i].first, mri_mov[j], mri_mov[tpi]); // use geometry of tpi
    intensities[j] = Md[i].second;
    //LTAwrite(ltas[j],(mov[i]+"-temp.lta").c_str());

    if (debug && fixtp)
    {
      ostringstream oss;
      oss << outdir << "tp" << j + 1 << "_to_mcoord";
      LTAwriteEx(ltas[j], (oss.str() + ".lta").c_str());
      MRI * warped = NULL;
      if (sampletype == SAMPLE_CUBIC_BSPLINE)
        warped = LTAtransformBSpline(mri_bsplines[j], NULL, ltas[j]);
      else
        warped = LTAtransformInterp(mri_mov[j], NULL, ltas[j], sampletype);

      if (iscale)
      {
        string fn = oss.str() + "-intensity.txt";
        ofstream f(fn.c_str(), ios::out);
        f << Md[i].second;
        f.close();

        warped = MyMRI::MRIvalscale(warped, warped, Md[i].second);
      }
      MRIwrite(warped, (oss.str() + ".mgz").c_str());
      MRIfree(&warped);
    }

  }

  if (fixtp) // we are done, as we mapped to this (itp) TP
  {
    return true;
  }

  // find geometric mean intensity scale
  if (iscale)
  {
    // here normalize scales to geometric mean:
    normalizeIntensities();
  }
  
  // convert to ras2ras and invert (to map from tpi to each other image
  vector<vnl_matrix_fixed<double, 4, 4> > mras(nin);
  mras[0].set_identity();
  vector<vnl_matrix_fixed<double, 3, 3> > mras3(nin);
  mras3[0].set_identity();
  vnl_vector_fixed<double, 3> meant(0.0);
  vnl_matrix_fixed<double, 3, 3> lin;
  vnl_vector_fixed<double, 3> trans;
  vnl_matrix_fixed < double , 4 , 4 > mrasinv;

  for (int i = 1; i < nin; i++)
  {
    int j = index[i];
    //cout << "  computing coord of TP " << j + 1 << " ( " << mov[j] << " )" << endl;
    //cout << "   wrt to TP " << tpi + 1 << " ( " << mov[tpi] << " )" << endl;
    mras[i] = MyMRI::MRIvoxelXformToRasXform(mri_mov[j], mri_mov[tpi], Md[i].first);
    // mras is also needed later!
    
    ////vnl_matlab_print(vcl_cout,mras[i],"mras",vnl_matlab_print_format_long);std::cout << std::endl;
    
    // invert so that each transform points from tpi to other time points
    mrasinv = vnl_inverse(mras[i]);    
    
    // split translation and linear matrix (rotation, when rigid):
    MyMatrix::getRTfromM(mrasinv, lin, trans);
    //vnl_matlab_print(vcl_cout,trans,"trans",vnl_matlab_print_format_long);std::cout << std::endl;
    //vnl_matlab_print(vcl_cout,lin,"lin",vnl_matlab_print_format_long);std::cout << std::endl;

    // these point from TPI to each image as we inverted the RAS2RAS above
    // so no further processing is necessary as
    // M x = A x + t
    // M^-1 x = A^(-1) ( x - t ) = A^(-1) x - A^(-1) t (which is what we had summed up before)
    // This way the linear map A is done first and then things are moved
    // so we can simply average the linear (or rotation) part and the translation 
    // separately and put back together. Thus the linear average will be applied fist, then trans.

    // simply average translation part
    meant += trans;

    // store linear or rigid part for later:
    mras3[i] = lin;
    
  } 
  //average translation
  meant = (1.0 / nin) * meant;
  if (debug)
  {
    vnl_matlab_print(vcl_cout,meant,"meant",vnl_matlab_print_format_long);std::cout << std::endl;  
  }

  // find mean rotation or linear map (in RAS coordinates)
  vnl_matrix_fixed<double, 3, 3> meanr;
  if (rigid)
    meanr = MyMatrix::RotationMean(mras3, frobnormthresh);
  else
    meanr = MyMatrix::GeometricMean(mras3);
  
  if (debug)
  {
    vnl_matlab_print(vcl_cout,meanr,"meanr",vnl_matlab_print_format_long);std::cout << std::endl;  

    cout << " Determinant( meanr ) : " << vnl_determinant(meanr) << endl << endl;
    cout << " Decompose into Rot * Shear * Scale : " << endl << endl;
    vnl_matrix<double> Rot, Shear;
    vnl_diag_matrix<double> Scale;
    MyMatrix::Polar2Decomposition(meanr, Rot, Shear, Scale);
    vnl_matlab_print(vcl_cout,Rot,"Rot",vnl_matlab_print_format_long);
    cout << endl;
    vnl_matlab_print(vcl_cout,Shear,"Shear",vnl_matlab_print_format_long);
    cout << endl;
    vnl_matlab_print(vcl_cout,Scale,"Scale",vnl_matlab_print_format_long);
    cout << endl;
  }

  // put back together to matrix in homogeneous coords
  vnl_matrix_fixed<double, 4, 4> Mm(MyMatrix::getMfromRT(meanr, meant));


  // construct target geometry for the mean space 
  // (the mean space may be outside any of the input RAS geometries)
  // this averages direction cosines (after un-flipping/reordering axis)
  // it also sets the center to average of CRAS of inputs
  MRI * template_geo = createTemplateGeo(); 

  if (! crascenter) //default behavior: use average centroid
  {
    vnl_vector_fixed<double, 4> ncenter(0.0);
    
    // center at the mean of all tp centroids mapped to the mean space
    // map average centroid from TPI vox space to mean RAS space:
    MATRIX * mv2r_temp = MRIgetVoxelToRasXform(mri_mov[tpi]);
    vnl_matrix_fixed<double, 4, 4> tpi_v2r(MyMatrix::convertMATRIX2VNL(mv2r_temp));
    MatrixFree(&mv2r_temp);
    //vnl_matlab_print(vcl_cout,tpi_v2r,"tpiv2r",vnl_matlab_print_format_long);std::cout << std::endl;
  
    ncenter = centroid;
    //vnl_matlab_print(vcl_cout,centroid,"centroid",vnl_matlab_print_format_long);std::cout << std::endl;

    // map to RAS:
    ncenter = tpi_v2r * ncenter;
    //vnl_matlab_print(vcl_cout,ncenter,"centroidras",vnl_matlab_print_format_long);std::cout << std::endl;
  
    // map to mean space
    ncenter = Mm * ncenter;
    //vnl_matlab_print(vcl_cout,ncenter,"centroidmeanras",vnl_matlab_print_format_long);std::cout << std::endl;
  
    // set new center in geometry
    template_geo->c_r = ncenter[0];
    template_geo->c_a = ncenter[1];
    template_geo->c_s = ncenter[2];
    template_geo->ras_good_flag = 1;
    MRIreInitCache(template_geo);
  }
  
  // (old) set direction cosines back to tpi
  //template_geo->x_r = mri_mov[tpi]->x_r;
  //template_geo->x_a = mri_mov[tpi]->x_a;
  //template_geo->x_s = mri_mov[tpi]->x_s;
  //template_geo->y_r = mri_mov[tpi]->y_r;
  //template_geo->y_a = mri_mov[tpi]->y_a;
  //template_geo->y_s = mri_mov[tpi]->y_s;
  //template_geo->z_r = mri_mov[tpi]->z_r;
  //template_geo->z_a = mri_mov[tpi]->z_a;
  //template_geo->z_s = mri_mov[tpi]->z_s;
  
  
  

  // construct maps from each image to the mean space
  vnl_matrix_fixed<double, 4, 4> M;
  for (int i = 0; i < nin; i++)
  {
    int j = index[i];
    cout << "  computing mean coord of TP " << j + 1 << " ( " << mov[j] << " ) "
        << endl;
    // concat transforms: meant meanr mras
    // from right to left, first mras[j] (aligns to T1)
    // then do the mean rot and move to mean location
    M = Mm * mras[i];
//vnl_matlab_print(vcl_cout,M,"M",vnl_matlab_print_format_long);std::cout << std::endl;

    if (rigid)
    {
         vnl_matrix < double > R(3,3),S(3,3),A(3,3),I(3,3);
         I.set_identity();
         M.extract(A);
         MyMatrix::PolarDecomposition(A,R,S);
         if (S[0][0] < 0.0 || S[1][1] < 0.0 || S[2][2] < 0.0)
           ErrorExit(ERROR_OUT_OF_BOUNDS, "Internal Error:  produced reflection.\n") ;
         double eps = 0.0001;
         double fnorm1 = (S-I).frobenius_norm();
         if (fnorm1 > eps)
         {
           std::cerr << "Internal Error for tp " << j << " -> template" << std::endl;
           std::cerr << " Rotation should not scale ( "<< fnorm1 << " )" << std::endl;
           std::cerr << " Debug Info: " << std::endl;
           vnl_matlab_print(vcl_cerr,A,"A",vnl_matlab_print_format_long);std::cerr << std::endl;
           vnl_matlab_print(vcl_cerr,R,"R",vnl_matlab_print_format_long);std::cerr << std::endl;
           vnl_matlab_print(vcl_cerr,S,"S",vnl_matlab_print_format_long);std::cerr << std::endl;
           std::cerr << " Make sure input voxel sizes are identical for all images!" << std::endl;
           ErrorExit(ERROR_OUT_OF_BOUNDS, "Internal Error: Rotation should not scale.\n") ;
         }
     
        cout << "   mapping back to rot, err = " << fnorm1 << endl;
        M.update(R);
        M.set_row(3, 0.0);
        M[3][3] = 1.0;
         
         
    }

    // make lta from M (M is RAS to RAS)
    assert(ltas.size() == mri_mov.size());
    //assert(ltas[j] == NULL);
    if (ltas[j] != NULL)
      LTAfree(&ltas[j]);
    ltas[j] = MyMatrix::RASmatrix2LTA(M, mri_mov[j], template_geo); // use geometry of template_geo
    //LTAwrite(ltas[j],(mov[i]+"-temp2.lta").c_str());  
    LTAchangeType(ltas[j], LINEAR_VOX_TO_VOX);


    /* removed this now (as inputs and template may have different voxel sizes...
       not sure what effect that has on isotropic images , proably rounding errors are small???
       hope we switch to double in image header geometries and transform
    if (rigid) // map back to Rotation (RAS2RAS->VOX2VOX introduces scaling!)
    {
      vnl_matrix<double> MM(
          MyMatrix::convertMATRIX2VNL(ltas[j]->xforms[0].m_L));
//vnl_matlab_print(vcl_cout,MM,"MM",vnl_matlab_print_format_long);std::cout << std::endl;
      vnl_matrix<double> R(3, 3), S(3, 3), A(3, 3), I(3, 3);
      I.set_identity();
      MM.extract(A);
      MyMatrix::PolarDecomposition(A, R, S);
      if (S[0][0] < 0.0 || S[1][1] < 0.0 || S[2][2] < 0.0)
        ErrorExit(ERROR_OUT_OF_BOUNDS,
            "Internal InitialXforms Error:  produced reflection.\n");

//         double eps = 0.0000001;        
      double fnorm1 = (S - I).frobenius_norm();
      cout << "   mapping back to rot, err = " << fnorm1 << endl;
//         if (fnorm1 > eps)
//         {
//           std::cerr << "Internal Error: " << std::endl;
//           std::cerr << " Rotation should not scale ( "<< fnorm1 << " )" << std::endl;
//           std::cerr << " Debug Info: " << std::endl;
//           vnl_matlab_print(vcl_cerr,A,"A",vnl_matlab_print_format_long);std::cerr << std::endl;
//           vnl_matlab_print(vcl_cerr,R,"R",vnl_matlab_print_format_long);std::cerr << std::endl;
//           vnl_matlab_print(vcl_cerr,S,"S",vnl_matlab_print_format_long);std::cerr << std::endl;
//           
//         
//           ErrorExit(ERROR_OUT_OF_BOUNDS, "Internal Error: Sqrt of Rotation should not scale.\n") ;
//         }

      MM.update(R);
      MM.set_row(3, 0.0);
      MM[3][3] = 1.0;
//vnl_matlab_print(vcl_cout,MM,"MM2",vnl_matlab_print_format_long);std::cout << std::endl;
      ltas[j]->xforms[0].m_L = MyMatrix::convertVNL2MATRIX(MM,
          ltas[j]->xforms[0].m_L);
    //LTAwrite(ltas[j],(mov[i]+"-temp3.lta").c_str());

    }*/

    if (debug)
    {
      ostringstream oss;
      oss << outdir << "tp" << j + 1 << "_to_mcoord";
      LTAwriteEx(ltas[j], (oss.str() + ".lta").c_str());
      MRI * warped = NULL;
      if (sampletype == SAMPLE_CUBIC_BSPLINE)
        warped = LTAtransformBSpline(mri_bsplines[j], NULL, ltas[j]);
      else
        warped = LTAtransformInterp(mri_mov[j], NULL, ltas[j], sampletype);
      if (iscale)
      {
        string fn = oss.str() + "-intensity.txt";
        ofstream f(fn.c_str(), ios::out);
        f << Md[i].second;
        f.close();

        warped = MyMRI::MRIvalscale(warped, warped, Md[i].second);
      }
      MRIwrite(warped, (oss.str() + ".mgz").c_str());
      MRIfree(&warped);
    }
  }

  //cleanup
  MRIfree(&template_geo);

  return true;
}

void MultiRegistration::normalizeIntensities()
// normalizes class member intensities:
// computes the geometric mean intensity and sets the intensity scale factors so that 
// their geometric mean is equal to 1
{
  double mint = 1.0;
  int nin = (int) intensities.size();
  for (int i = 0; i < nin; i++)
    mint *= intensities[i];

  // geometric mean
  mint = pow(mint, 1.0 / nin);

  // set intenstiy factors so that geo-mean is 1:
  for (int i = 0; i < nin; i++)
    intensities[i] *= (1.0 / mint);

}

bool MultiRegistration::writeMean(const std::string& mean)
{
  if (!mri_mean)
  {
    cout << " ERROR: No average exists! Skipping output." << endl;
    return false;
  }
  strncpy(mri_mean->fname, mean.c_str(), STRLEN-1);
  return (MRIwrite(mri_mean, mean.c_str()) == 0);
}

bool MultiRegistration::writeConformMean(const std::string& mean)
{
  if (!mri_mean)
  {
    cout << " ERROR: No average exists! Skipping output." << endl;
    return false;
  }

  // create conform mean
  MRI * mri_cmean = averageConformSet(0);

  strncpy(mri_cmean->fname, mean.c_str(), STRLEN-1);
  int ok = MRIwrite(mri_cmean, mean.c_str());
  return (ok == 0);
}

bool MultiRegistration::writeLTAs(const std::vector<std::string> & nltas,
    bool vox2vox, const std::string & mean)
{
  assert(nltas.size() == ltas.size());
  int error = 0;
  for (unsigned int i = 0; i < nltas.size(); i++)
  {
    if (!ltas[i])
    {
      cout << " ERROR: No ltas exist! Skipping output." << endl;
      return false;
    }

    if (vox2vox)
    {
      error += (LTAchangeType(ltas[i], LINEAR_VOX_TO_VOX) == NULL);
    }
    else
    {
      error += (LTAchangeType(ltas[i], LINEAR_RAS_TO_RAS) == NULL);
    }
    strncpy(ltas[i]->xforms[0].dst.fname, mean.c_str(), STRLEN-1);
    strncpy(ltas[i]->xforms[0].src.fname, mov[i].c_str(), STRLEN-1);
    LTAwriteEx(ltas[i], nltas[i].c_str());

    vnl_matrix<double> fMv2v = MyMatrix::LTA2VOXmatrix(ltas[i]);
    cout << " Determinant( lta[ " << i << " ]) : " << vnl_determinant(fMv2v)
        << endl << endl;

    if (!rigid)
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
    }

  }
  return (error == 0);
}

bool MultiRegistration::writeMapMovHdr(
  const std::vector<std::string>& mapmovhdr)
{
  assert(mapmovhdr.size() == mri_mov.size());
  MATRIX* ras2ras = MatrixAlloc(4, 4, MATRIX_REAL);
  MATRIX* vox2ras;
  int error = 0;
  for (unsigned int i = 0; i < mapmovhdr.size(); i++)
  {
    if (!ltas[i])
    {
      std::cout << " ERROR: No LTAs exist! Skipping output.\n";
      error = 1;
      break;
    }
    vnl_matrix<double> fMr2r = MyMatrix::LTA2RASmatrix(ltas[i]);
    ras2ras = MyMatrix::convertVNL2MATRIX(fMr2r, ras2ras);
    vox2ras = MRIgetVoxelToRasXform(mri_mov[i]);
    vox2ras = MatrixMultiply(ras2ras, vox2ras, vox2ras);
    MRI *mri_aligned = MRIcopy(mri_mov[i], NULL);
    MRIsetVoxelToRasXform(mri_aligned, vox2ras);
    error = MRIwrite(mri_aligned, mapmovhdr[i].c_str());
    MRIfree(&mri_aligned);
    MatrixFree(&vox2ras);
    if (error)
    {
      std::cout << "ERROR: Can't write " << mapmovhdr[i].c_str() << '\n';
      break;
    }
  }
  MatrixFree(&ras2ras);
  return !error;
}

bool MultiRegistration::writeWarps(const std::vector<std::string>& nwarps)
{
  assert(nwarps.size() == mri_warps.size() || nwarps.size() == 1);
  int error = 0;
  //vnl_vector < double > ls (nwarps.size());
  //vnl_vector < double > lsa (nwarps.size(),0.0);
  //double lsavg = 0.0;
  //double lsavga = 0.0;
  //int count = 0;

  if (nwarps.size() == 1)
  {
    cout << " Writing single output (4D volume) ..." << endl;
    if (!mri_warps[0])
    {
      cout << " ERROR: Warps do not exist! Skipping output." << endl;
      return false;
    }

    MRI * mri4d = MRIallocSequence(mri_warps[0]->width, mri_warps[0]->height,
        mri_warps[0]->depth, mri_warps[0]->type, mri_warps.size());

    for (unsigned int i = 0; i < mri_warps.size(); i++)
    {
      if (!mri_warps[i])
      {
        cout << " ERROR: Warp " << i << " does not exist! Skipping output."
            << endl;
        return false;
      }
      MRIcopyFrame(mri_warps[i], mri4d, 0, i);
    }
    error = MRIwrite(mri4d, nwarps[0].c_str());
    MRIfree(&mri4d);
  }
  else // write individually
  {
    for (unsigned int i = 0; i < nwarps.size(); i++)
    {
      if (!mri_warps[i])
      {
        cout << " ERROR: Warp " << i << " does not exist! Skipping output."
            << endl;
        return false;
      }
      error += MRIwrite(mri_warps[i], nwarps[i].c_str());
      //ls[i] = CostFunctions::leastSquares(mri_warps[i],mri_mean);
      //for (unsigned int j = i+1;j<nwarps.size();j++)
      //{
      //  lsa[i] += CostFunctions::tukeyBiweight(mri_warps[i],mri_warps[j]);
      //  count++;
      //}
      //lsavga += lsa[i];
      //lsavg += ls[i];
    }
    //vnl_matlab_print(vcl_cout,ls,"tb",vnl_matlab_print_format_long);
    //cout << endl;
    //cout << " tbavg: " << lsavg/nwarps.size() << endl;
    //cout << " tbavga: " << lsavga/count << endl;
  }
  return (error == 0);
}

bool MultiRegistration::writeIntensities(
    const std::vector<std::string>& nintens)
{
  assert(nintens.size() == intensities.size());
  int error = 0;
  for (unsigned int i = 0; i < nintens.size(); i++)
  {
    assert(intensities[i] > 0);
    ofstream f(nintens[i].c_str(), ios::out);
    if (f.good())
    {
      f << intensities[i];
      f.close();
    }
    else
      error++;
  }
  return (error == 0);
}

bool MultiRegistration::writeWeights(const std::vector<std::string>& nweights,
    bool oneminusweights)
{
  assert(nweights.size() == mri_weights.size());
  int error = 0;
  for (unsigned int i = 0; i < nweights.size(); i++)
  {
    if (!mri_weights[i])
    {
      cout << " No weights constructed, skipping output of weights" << endl;
      return false;
    }
    MRI * wtmp = mri_weights[i];
    if (oneminusweights)
      wtmp = MRIlinearScale(wtmp, NULL, -1, 1, 0);
    error += MRIwrite(wtmp, nweights[i].c_str());
    if (oneminusweights)
      MRIfree(&wtmp);
  }
  return (error == 0);
}

/*!
 \fn MRI* averageSet(const vector < MRI * >& set, MRI* mean, int method, double sat)
 \brief Averages the passed and aligned volumes depending on method 
 \param set vector of movable volumes
 \param mean  output of mean volume
 \param method  0 = mean, 1 = median, 2 = tukey biweight (testing)
 \param sat     saturation for tukey biweight
 */
MRI* MultiRegistration::averageSet(const std::vector<MRI *>& set, MRI* mean,
    int method, double sat)
{
  assert(set.size() > 1);

//    for(uint i = 0;i<set.size();i++)
//    {
//       cout << " TP " << i+1 << endl;
//       cout << " mean   : " << CostFunctions::mean(set[i])   << endl;
//       cout << " sdev   : " << CostFunctions::sdev(set[i])   << endl;
//       cout << " median : " << CostFunctions::median(set[i]) << endl;
//       cout << " mad    : " << CostFunctions::mad(set[i])    << endl;
//     }

  if (method == 0)
  {
    // mean
    cout << "    using mean" << endl;
    for (unsigned int i = 0; i < set.size(); i++)
    {
      mean = MRIaverage(set[i], i, mean);
    }
  }
  else if (method == 1)
  {
    cout << " using median " << endl;
    // robust
    int x, y, z, i;
    assert(set.size() > 0);
    // float dd[set.size()];
    if (!mean)
      mean = MRIclone(set[0], NULL);
    //MRI * midx = MRIalloc(set[0]->width,set[0]->height,set[0]->depth, MRI_FLOAT);
    //MRIcopyHeader(mean,midx);
    //midx->type = MRI_FLOAT;
    //   pair < float, float > mm;
#ifdef HAVE_OPENMP
#pragma omp parallel for private(y,x,i) shared(set,mean) schedule(guided)
#endif
    for (z = 0; z < set[0]->depth; z++)
    {
      float dd[set.size()];
      for (y = 0; y < set[0]->height; y++)
        for (x = 0; x < set[0]->width; x++)
        {
          for (i = 0; i < (int) set.size(); i++)
            dd[i] = MRIgetVoxVal(set[i], x, y, z, 0);
          pair<float, float> mm = RobustGaussian<float>::medianI(dd,
              (int) set.size());
          MRIsetVoxVal(mean, x, y, z, 0, mm.first);
          //MRIsetVoxVal(midx,x,y,z,0,mm.second);
        }
      //MRIwrite(midx,"midx.mgz");
      //MRIwrite(mean,"m.mgz");
      //assert(1==2);
    }
  }
  else if (method == 2) // NEVER REALLY WORKED!! (not enough values and usually int)
  {
    cout << "    using tukey biweight" << endl;
    // robust tukey biweight
    int x, y, z, i;
//     MATRIX* b = MatrixAlloc(set.size(),1,MATRIX_REAL);
//     pair < MATRIX* , MATRIX* >  muw;
    vnl_vector<float> b(set.size());
    assert(set.size() > 0);
    if (!mean)
      mean = MRIclone(set[0], NULL);
    for (z = 0; z < set[0]->depth; z++)
      for (y = 0; y < set[0]->height; y++)
        for (x = 0; x < set[0]->width; x++)
        {
          // cout << " x: " << x << " y: " << y << " z: " <<z << "  size: " << set.size() <<endl;
          for (i = 0; i < (int) set.size(); i++)
          {
            //cout << "i: " << i << endl;
            //b->rptr[i+1][1] = MRIgetVoxVal(set[i],x,y,z,0);
            b[i] = MRIgetVoxVal(set[i], x, y, z, 0);
          }
          //cout << endl << "intgensities at this voxel:" << endl; ;
          //MatrixPrintFmt(stdout,"% 2.8f",b);

          Regression<float> R(b);
          vnl_vector<float> p = R.getRobustEst(sat);
          //muw = R.getRobustEstW(sat);
          //    cout << " tukey mean: " << muw.first->rptr[1][1] << endl;
          //MRIsetVoxVal(mean,x,y,z,0,muw.first->rptr[1][1]);
          MRIsetVoxVal(mean, x, y, z, 0, p[0]);
          //MatrixFree(&muw.first);
          //MatrixFree(&muw.second);
        }
  }
  else if (method == 3) // needs more development (sigma..)
  {
    cout << "    using tukey biweight (alltoghether)" << endl;
    // robust tukey biweight
    int x, y, z, i;
    int n = set[0]->depth * set[0]->height * set[0]->width;
//     MATRIX* b = MatrixAlloc(set.size()*n,1,MATRIX_REAL);
//     MATRIX* A = MatrixAlloc(set.size()*n,n,MATRIX_REAL);
//     A = MatrixZero(A->rows,A->cols,A);
//     pair < MATRIX* , MATRIX* >  muw;
    vnl_matrix<float> A(set.size() * n, n);
    vnl_vector<float> b(set.size() * n);
    A.fill(0.0);

    assert(set.size() > 0);
    if (!mean)
      mean = MRIclone(set[0], NULL);
    for (i = 0; i < (int) set.size(); i++)
    {
      assert(set[i]->width == set[0]->width);
      assert(set[i]->height == set[0]->height);
      assert(set[i]->depth == set[0]->depth);
      int pcount = 0;
      for (z = 0; z < set[0]->depth; z++)
        for (y = 0; y < set[0]->height; y++)
          for (x = 0; x < set[0]->width; x++)
          {
            // cout << " x: " << x << " y: " << y << " z: " <<z << "  size: " << set.size() <<endl;
            //cout << "i: " << i << endl;
//             b->rptr[pcount*set.size()+i+1][1] =
//               MRIgetVoxVal(set[i],x,y,z,0);
//             A->rptr[pcount*set.size()+i+1][pcount+1] = 1;
            b[pcount * set.size() + i] = MRIgetVoxVal(set[i], x, y, z, 0);
            A[pcount * set.size() + i][pcount] = 1.0;
            pcount++;
          }
    }
    //cout << endl << "intgensities at this voxel:" << endl; ;
    //MatrixPrintFmt(stdout,"% 2.8f",b);

    Regression<float> R(A, b);
    //muw = R.getRobustEstW(sat);
    vnl_vector<float> p = R.getRobustEst(sat);
    //    cout << " tukey mean: " << muw.first->rptr[1][1] << endl;
    int pcount = 0;
    for (z = 0; z < set[0]->depth; z++)
      for (y = 0; y < set[0]->height; y++)
        for (x = 0; x < set[0]->width; x++)
        {
          //pcount++;
          //MRIsetVoxVal(mean,x,y,z,0,muw.first->rptr[pcount][1]);
          MRIsetVoxVal(mean, x, y, z, 0, p[pcount]);
          pcount++;
        }
    //MatrixFree(&muw.first);
    //MatrixFree(&muw.second);

  }
  else
  {

    cerr << " averageSet  method " << method << " unknown" << endl;
    exit(1);
  }
//       cout << " Average Vol "  << endl;
//       cout << " mean   : " << CostFunctions::mean(mean) << endl;
//       cout << " sdev   : " << CostFunctions::sdev(mean) << endl;
//       cout << " median : " << CostFunctions::median(mean) << endl;
//       cout << " mad    : " << CostFunctions::mad(mean)<< endl;
  printf(
    "       -- Template : (%g, %g, %g)mm and (%d, %d, %d) voxels.\n",
      mean->xsize, mean->ysize, mean->zsize,
      mean->width, mean->height,mean->depth);


  return mean;
}

/*!
 \fn MRI* initialAverageSet()
 \brief Creates initial translation alignment based on centroids and averages the mapped movables
 */
bool MultiRegistration::initialAverageSet()
{
// the initial average can be any place as long as it is
// not identical to one of the
// initial images (to avoid introducing a bias).
// it can be the average of all images but that is extremely blurry, and if
// we have only few images, far apart, maybe the algorithm will not converge
// or converge very slowly (has not happened yet).
// therefore we find an initial coarse alignment by using moments
// initial ltas are also returned
  std::cout << " MultiRegistration::initialAverageSet " << std::endl;

  assert(mri_mov.size() > 1);
  int n = (int) mri_mov.size();
  if ((int) ltas.size() != n)
  {
    assert(ltas.size() == 0);
    ltas.resize(n, NULL);
  }

  if (iscaleonly)
  {
    std::cout << "    -- averaging without alignment (iscale only)" << std::endl;
  
    if (!mri_mean) // all inputs should be in same space for this!!!
    {
      mri_mean = MRIclone(mri_mov[0], NULL);   
    }
    
  
    MATRIX* Mtrans = MatrixIdentity(4, NULL);
    for (int i = 0; i < n; i++)
    {
      assert(ltas[i] == NULL);
      ltas[i] = MyMatrix::RASmatrix2LTA(Mtrans, mri_mov[i], mri_mean);
    }  
  
  }
  else
  {
    std::cout << "    -- aligning image centroids (translation)" << std::endl;
    // compute input centroids and common centroid
    vector<double> centroid(3, 0.0);
    vector<vector<double> > centroids(n);
    for (int i = 0; i < n; i++)
    {
      centroids[i] = CostFunctions::centroid(mri_mov[i]);
      MATRIX * mv2r_temp = MRIgetVoxelToRasXform(mri_mov[i]);
      vnl_matrix_fixed<double, 4, 4> tpi_v2r(
          MyMatrix::convertMATRIX2VNL(mv2r_temp));
      MatrixFree(&mv2r_temp);
      vnl_vector_fixed<double, 4> ncenter;
      for (uint ii = 0; ii < 3; ii++)
        ncenter[ii] = centroids[i][ii];
      ncenter[3] = 1.0;
      // map to RAS:
      ncenter = tpi_v2r * ncenter;
      for (uint ii = 0; ii < 3; ii++)
        centroids[i][ii] = ncenter[ii];
  
      //cout << " Centroid [ " << i << " ] = " << centroids[i][0] << " " << centroids[i][1] << " " << centroids[i][2] << endl;
      centroid[0] += centroids[i][0];
      centroid[1] += centroids[i][1];
      centroid[2] += centroids[i][2];
    }
    centroid[0] /= n;
    centroid[1] /= n;
    centroid[2] /= n;
    //cout << " Centroid : " << centroid[0] << " " << centroid[1] << " " << centroid[2]  << endl;
  
    if (!mri_mean) //default:
    {
      // take geometry from set[0] and set RAS center at joint center of mass
      mri_mean = MRIclone(mri_mov[0], NULL);
      // map average centroid from TPI vox space to mean RAS space:
      MATRIX * mv2r_temp = MRIgetVoxelToRasXform(mri_mov[0]);
      vnl_matrix_fixed<double, 4, 4> tpi_v2r(
          MyMatrix::convertMATRIX2VNL(mv2r_temp));
      MatrixFree(&mv2r_temp);
      vnl_vector_fixed<double, 4> ncenter;
      for (uint ii = 0; ii < 3; ii++)
        ncenter[ii] = centroid[ii];
      ncenter[3] = 1.0;
      // map to RAS:
      ncenter = tpi_v2r * ncenter;
      // set new center in geometry
      for (uint ii = 0; ii < 3; ii++)
        centroid[ii] = ncenter[ii];
      mri_mean->c_r = ncenter[0];
      mri_mean->c_a = ncenter[1];
      mri_mean->c_s = ncenter[2];
      mri_mean->ras_good_flag = 1;
      MRIreInitCache(mri_mean);
    }
  
    // create maps (lta) to centered image
    MATRIX* Mtrans = MatrixIdentity(4, NULL);
    for (int i = 0; i < n; i++)
    {
      Mtrans->rptr[1][4] = centroid[0] - centroids[i][0];
      Mtrans->rptr[2][4] = centroid[1] - centroids[i][1];
      Mtrans->rptr[3][4] = centroid[2] - centroids[i][2];
      // set initial ltas as maps to the blurry mean
      assert(ltas[i] == NULL);
      ltas[i] = MyMatrix::RASmatrix2LTA(Mtrans, mri_mov[i], mri_mean);
    }
  }

  // map source to average space and create mri_mean
  mapAndAverageMov(0);

  return true;
}
