/**
 * @file MultiRegistration.h
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
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/15 16:12:27 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2008-2009
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */



#include "MultiRegistration.h"
#include "Registration.h"
#include "Regression.h"
#include "RobustGaussian.h"
#include "CostFunctions.h"
#include "MyMatrix.h"
#include "MyMRI.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

// all other software are all in "C"
#ifdef __cplusplus
extern "C"
{
#endif
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#ifdef __cplusplus
}
#endif

using namespace std;

void MultiRegistration::clear()
{
  if (mri_mean) MRIfree(&mri_mean);

  for (unsigned int i = 0;i<mri_mov.size();i++)
  {
			if (i<ltas.size() && ltas[i]) LTAfree(&ltas[i]);
			if (i<mri_warps.size() && mri_warps[i]) MRIfree(&mri_warps[i]);
      if (i<mri_weights.size() && mri_weights[i]) MRIfree(&mri_weights[i]);
      if (i<mri_mov.size() && mri_mov[i]) MRIfree(&mri_mov[i]);

  }

}

/*!
  \fn int loadMovables(const std::vector < std::string > mov)
  \brief Loads the movable volumes as specified on command line
  \param P  Paramters for the initialization
*/
int MultiRegistration::loadMovables(const std::vector < std::string > pmov)
{

  assert (mri_mov.size () == 0);
  int n = (int) pmov.size();
  mri_mov.resize(n);
	mov = pmov; // copy of input filenames

  for (unsigned int i = 0;i<mov.size(); i++)
  {
    cout << "reading source '"<<mov[i]<<"'..." << endl;

    mri_mov[i] = MRIread(mov[i].c_str()) ;
    if (!mri_mov[i])
    {
      ErrorExit(ERROR_NOFILE, "MultiRegistration::loadMovables: could not open input volume %s.\n",
                mov[i].c_str()) ;
    }
  }

  mri_warps.resize(n,NULL);
  intensities.resize(n,1.0);
  ltas.resize(n,NULL);
  if (robust) mri_weights.resize(n,NULL);

  return mov.size();
}

int MultiRegistration::loadLTAs(const std::vector < std::string > nltas)
{
	 assert(nltas.size() == mov.size());
	 assert(ltas.size() == mov.size());
	 iltas = nltas; // copy of input filenames
	 for (uint i = 0;i<nltas.size();i++)
	 {
        TRANSFORM * trans = TransformRead(nltas[i].c_str());
        LTA* lta =  (LTA *)trans->xform ;
        if (!lta)
          ErrorExit(ERROR_BADFILE, "MultiRegistration::loadLTAs could not read transform file %s", nltas[i].c_str()) ;
        lta = LTAchangeType(lta,LINEAR_VOX_TO_VOX);

        ltas[i] = lta;
   }
   return (int)ltas.size();
}



/*!
  \fn void initRegistration(Registration & R, const Parameters & P)
  \brief Initializes a Registration with Parameters (rigid, iscale, transonly, robust, sat and trans)
  \param R  Registration to be initialized
*/
void MultiRegistration::initRegistration(Registration & R)
{
  // assert(n < (int) P.mov.size());

  R.setRigid(rigid);
  R.setIscale(iscale);
  R.setTransonly(transonly);
  R.setRobust(robust);
  R.setSaturation(sat);
  //R.setDebug(debug);
  
	
	//if (subsamplesize > 0) R.setSubsamplesize(subsamplesize);

//   int pos = P.mov[n].rfind(".");
//   if (pos > 0) R.setName(P.mov[n].substr(0,pos));
//   else  R.setName(P.mov[n]);
//
//  // if (P.subsamplesize > 0) R.setSubsamplesize(P.subsamplesize);
//
//
//   R.setSource(P.mri_mov[n],P.fixvoxel,P.fixtype);
//   R.setTarget(P.mri_mean,P.fixvoxel,P.fixtype);
}

bool MultiRegistration::averageSet(int itdebug)
// maps mov to template (using ltas)
// adjust intensities (if iscale)
// creates template acording to average:1 mean, 2 median..
{
  unsigned int nin = mri_mov.size();
	assert(nin > 1);
  assert (mri_warps.size() == nin);
	assert (ltas.size() == nin);
	assert (intensities.size() == nin);
	
    cout << "warping movs and creating initial template..." << endl;
		if (iscale) cout << " allow intensity scaling" << endl;
    for (unsigned int i = 0;i<nin;i++)
    {  
      if (mri_warps[i]) MRIfree(&mri_warps[i]); 
      mri_warps[i] = MRIclone(mri_mov[0],mri_warps[i]);
      mri_warps[i] = LTAtransform(mri_mov[i],mri_warps[i], ltas[i]);
			if (iscale)
			{
			   mri_warps[i] = MyMRI::MRIvalscale(mri_warps[i],mri_warps[i], intensities[i]);
			}
      if (debug)
      {
        ostringstream oss;
        oss << outdir << "tp" << i+1 << "_to_template-it"<<itdebug<<".mgz";
        MRIwrite(mri_warps[i], oss.str().c_str()) ;
      }
    }
    mri_mean = averageSet(mri_warps, mri_mean,average,sat);
    return true;
}

bool MultiRegistration::computeTemplate(int itmax, double eps , int iterate, double epsit)
// itmax : iterations for template creation
// iterate: iterations for registration
// epsit: 
// uses ltas (if set), mri_mov and parameters
// sets mri_mean, ltas, mri_warps (mri_weights and intensities if available)
{
  int nin = (int) mri_mov.size();
	assert(nin > 1);
  
	if (debug) cout << endl << endl << "MultiRegistration::computeTemplate(avitmax: "<< itmax << ", aveps " << eps << ", regit: " << iterate << ", regeps: " <<epsit<<" )" << endl;
	
	
  cout << "Computing first template" << endl;  
  bool havexforms = (ltas[0] != NULL);
  if (!havexforms) // create simple initial transforms (not recommended)
    mri_mean = initialAverageSet(mri_mov,NULL,average,sat);
  else // we have initial transforms
  {  

    averageSet(0);
//     cout << "warping movs and creating initial template..." << endl;
// 		if (iscale) cout << " allow intensity scaling" << endl;
//     for (int i = 0;i<nin;i++)
//     {  
//       if (mri_warps[i]) MRIfree(&mri_warps[i]); // should not happen
//       mri_warps[i] = MRIclone(mri_mov[0],mri_warps[i]);
//       mri_warps[i] = LTAtransform(mri_mov[i],mri_warps[i], ltas[i]);
// 			if (iscale)
// 			{
// 			   mri_warps[i] = MyMRI::MRIvalscale(mri_warps[i],mri_warps[i], intensities[i]);
// 			}
//       if (debug)
//       {
//         ostringstream oss;
//         oss << outdir << "tp" << i+1 << "_to_template-it0.mgz";
//         MRIwrite(mri_warps[i], oss.str().c_str()) ;
//       }
//     }
//     mri_mean = averageSet(mri_warps, mri_mean,average,sat);
  }
  
  if (itmax==0) // no iterations necessary just return with mean (and ltas and warps);
  {
    //strncpy(P.mri_mean->fname, P.mean.c_str(),STRLEN);
    return true;
  }

  
  strncpy(mri_mean->fname,(outdir+"template-it0.mgz").c_str(),STRLEN);
  if (debug)
  {
    cout << "debug: saving template-it0.mgz" << endl;
    MRIwrite(mri_mean,(outdir+"template-it0.mgz").c_str());
  }

   cout << "template fname: " << mri_mean->fname << endl;

//  int itmax  = 10;
//  double eps = 0.025;

  int itcount = 0;
  double maxchange = 100;

  vector < MATRIX * > transforms(mri_mov.size(),NULL);
  if (havexforms) //init transforms
  {
     for (int i=0;i<nin;i++)
     {
        LTAchangeType(ltas[i],LINEAR_VOX_TO_VOX); //vox2vox for registration
        transforms[i] = MatrixCopy(ltas[i]->xforms[0].m_L, NULL);
     }
	
  }
  
  vector < Registration > Rv(mri_mov.size());
  for (int i = 0;i<nin;i++) Rv[i].setSource(mri_mov[i],fixvoxel,fixtype);

  LTA * lastlta = NULL;
  while (itcount < itmax && maxchange > eps)
  {
    itcount++;
		
    // if we have prior ltas, allow termination
		// by initializing maxchange =0 so that it can be 
		// computed below
		if ( ltas[0] ) maxchange = 0;

    cout << endl << "=======================================================" << endl;
    cout << endl << "Working on global iteration : " << itcount << endl;
    cout << endl << "=======================================================" << endl;
    //cout << "========================================" << endl;
    //printmemusage();
    //cout << "========================================" << endl << endl;
		
    // register all inputs to mean
    vector < double > dists(nin,1000); // should be larger than maxchange!
    for (int i = 0;i<nin;i++)
    {
		  cout << endl << "Working on TP " << i+1 << endl;
      Rv[i].clear();
      initRegistration(Rv[i]); //set parameter
      Rv[i].setTarget(mri_mean,
                      fixvoxel,
                      fixtype); // gaussian pyramid will be constructed for
                                  // each Rv[i], could be optimized
      ostringstream oss;
      oss << outdir << "tp" << i+1 << "_to_template-it" << itcount;
      Rv[i].setName(oss.str());

      // compute Alignment
      std::pair <MATRIX*, double> Md;
      //int iterate = P.iterate;
      //double epsit= P.epsit;
      int maxres = 0;
			// on higher iterations use subsamplesize as passed on commandline
      int subsamp = subsamplesize;
			// simplify first steps (only if we do not have good transforms):
      switch (itcount)
      {
      case 1:
        maxres = 2;
        break;
      case 2:
        maxres = 2;
        break;
      case 3:
        maxres = 1;
        break;
      case 4:
        maxres = 1;
        break;
      case 5:
        subsamp = 180;
        break;
      }

      if (havexforms) // default case, as we use the new mean space init!
      {
        //// tried high res iteration only:
				// (does not make sense, we need 2 iter. on high res anyway)
        //Md= Rv[i].computeIterativeRegistration(1,epsit,NULL,NULL,transforms[i],P.intensities[i]);
        //dists[i] = sqrt(Rv[i].AffineTransDistSq(Md.first, transforms[i]));
        //if (itcount ==1) subsamp = 180;	// so use subsample for first step on high res
        maxres = 0; //go up to hig-res allways (skip first steps as above)
      }

      Rv[i].setSubsamplesize(subsamp);
      Md = Rv[i].computeMultiresRegistration(maxres,
                                               iterate,
                                               epsit,
                                               NULL,
                                               NULL,
                                               transforms[i],
                                               intensities[i]);
      if (transforms[i]) MatrixFree(&transforms[i]);
      transforms[i] = Md.first;
      intensities[i] = Md.second;

      // convert Matrix to LTA ras to ras
      if (lastlta)  LTAfree(&lastlta);
      if (ltas[i]) lastlta = ltas[i];
      ltas[i] = MyMatrix::VOXmatrix2LTA(Md.first,mri_mov[i],mri_mean);
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
        LTAchangeType(lastlta,LINEAR_RAS_TO_RAS); //measure dist in RAS coords
        dists[i] = sqrt(MyMatrix::AffineTransDistSq(lastlta->xforms[0].m_L,
                                                ltas[i]->xforms[0].m_L));
        LTAfree(&lastlta);
        if (dists[i] > maxchange) maxchange = dists[i];
        cout << endl << "tp " << i+1 << " distance: " << dists[i] << endl;
      }


      // create warps: warp mov to mean
      cout << "warping mov to template..." << endl;
      if (mri_warps[i]) MRIfree(&mri_warps[i]);
      mri_warps[i] = MRIclone(mri_mean,mri_warps[i]);
      mri_warps[i] = LTAtransform(mri_mov[i],mri_warps[i], ltas[i]);

      // here do scaling of intensity values
      if (Rv[i].isIscale() && Md.second > 0)
      {
        cout << "Adjusting Intensity of WARP by " << Md.second << endl;
        mri_warps[i] = MyMRI::MRIvalscale(mri_warps[i],
                                         mri_warps[i], Md.second);
      }
			
			// copy weights (as RV will be cleared)
      //   (info: they are in original half way space)
	     cout << "backup weights ..." << endl;
       mri_weights[i] = MRIcopy(Rv[i].getWeights(),mri_weights[i]);
			 
      //cout << " LS difference after: " <<
      //CF.leastSquares(mri_aligned,P.mri_dst) << endl;
      //cout << " NC difference after: " <<
      //CF.normalizedCorrelation(mri_aligned,P.mri_dst) << endl;


      if (debug)
      {
        cout << "debug: writing transforms, warps, weights ..." << endl;
        LTAwriteEx(ltas[i], (oss.str()+".lta").c_str()) ;

        MRIwrite(mri_warps[i], (oss.str()+".mgz").c_str()) ;

        if (Rv[i].isIscale() && Md.second >0)
        {
          string fn = oss.str() + "-intensity.txt";
          ofstream f(fn.c_str(),ios::out);
          f << Md.second;
          f.close();
        }

        // if we have weights:  
        if (mri_weights[i] != NULL)
        {
          std::pair <MATRIX*, MATRIX*> map2weights = Rv[i].getHalfWayMaps();
          MATRIX * hinv = MatrixInverse(map2weights.second,NULL);

          cout << endl;
          MatrixPrint(stdout,map2weights.first) ;
          cout << endl;
          MatrixPrint(stdout,map2weights.second) ;
          cout << endl;
          MatrixPrint(stdout,hinv) ;
          cout << endl;

          MRI * wtarg = MRIalloc(mri_weights[i]->width,
                                 mri_weights[i]->height,
                                 mri_weights[i]->depth,
                                 MRI_FLOAT);
          MRIcopyHeader(mri_weights[i],wtarg);
          MATRIX * v2r  = MRIgetVoxelToRasXform(mri_mean);
          MRIsetVoxelToRasXform(wtarg,v2r);
          wtarg->type = MRI_FLOAT;
          wtarg->i_to_r__ = MatrixCopy(mri_mean->i_to_r__, wtarg->i_to_r__);
          wtarg->r_to_i__ = MatrixCopy(mri_mean->r_to_i__, wtarg->r_to_i__);

          wtarg = MRIlinearTransform(mri_weights[i],wtarg, hinv);
          MRIwrite(wtarg, (oss.str()+"-weights.mgz").c_str()) ;
          MRIwrite(mri_weights[i], (oss.str()+"-www.mgz").c_str());
          MRIfree(&wtarg);
          MatrixFree(&hinv);
          MatrixFree(&v2r);
        }
      } // if debug end
    
		  // clear to reduce memory usage:
      Rv[i].clear();
      Rv[i].freeGPT();
			
      cout << endl << "Finished TP : " << i+1 << endl;
      cout << endl;
      cout << "=====================================================" << endl;
      //printmemusage();
      //cout << "========================================" << endl << endl;

    } // for loop end (all timepoints)

    if (dists[0] <= maxchange) // it was computed, so print it:
    {
      cout << endl << "Global Iteration " << itcount <<" Distances: " << endl;
      for (unsigned int i=0;i<dists.size();i++) cout << dists[i] <<" ";
      cout << endl << "Maxchange: " << maxchange << endl << endl;
    }

    // create new mean
    cout << "Computing new template " <<itcount << endl;
    mri_mean = averageSet(mri_warps, mri_mean,average,sat);
    if (debug)
    {
      ostringstream oss;
      oss << outdir << "template-it" << itcount << ".mgz";
      cout << "debug: writing template to " << oss.str() << endl;
      MRIwrite(mri_mean,oss.str().c_str());
      //strncpy(mri_mean->fname, oss.str().c_str(),STRLEN);
    }

  } // end while

  //strncpy(P.mri_mean->fname, P.mean.c_str(),STRLEN);
  
	cout << " DONE : computeTemplate " << endl;
  return true;
}


bool MultiRegistration::halfWayTemplate(int maxres, int iterate, double epsit, bool vox2vox)
// can be used only with two input images
{
  int nin = (int) mri_mov.size();
  assert (nin == 2);

  // register 1 with 2

  Registration R;
  initRegistration(R); //set parameter
  R.setSource(mri_mov[0],fixvoxel,fixtype);
  R.setTarget(mri_mov[1],fixvoxel,fixtype);

  ostringstream oss;
  oss << outdir << "halfway_template.mgz";
  R.setName(oss.str().c_str());

  // use initial transform if given:
  MATRIX *minit = NULL;
	assert(ltas.size() ==2);
  if (ltas[0])
  {
    LTAchangeType(ltas[0],LINEAR_VOX_TO_VOX); //vox2vox for registration
    LTAchangeType(ltas[1],LINEAR_VOX_TO_VOX); //vox2vox for registration
    assert (ltas[0]->type == LINEAR_VOX_TO_VOX);
    assert (ltas[1]->type == LINEAR_VOX_TO_VOX);
    minit = MatrixInverse(ltas[1]->xforms[0].m_L,NULL);
    minit = MatrixMultiply(minit,ltas[0]->xforms[0].m_L,minit);
    R.setDebug(1);
  }	
		
  // compute Alignment
  std::pair <MATRIX*, double> Md;
	// adjust subsamplesize, if passed:
	if (subsamplesize > 0 ) R.setSubsamplesize(subsamplesize);
  Md = R.computeMultiresRegistration(maxres,iterate,epsit,NULL,NULL,minit);
	if (minit) MatrixFree(&minit);
  intensities[0] = Md.second;
  intensities[1] = 1;

  cout << "creating half way data ..." << endl;
  std::pair < MATRIX*, MATRIX*> maps2weights = R.getHalfWayMaps();
  LTA * m2hwlta = LTAalloc(1,mri_mov[0]);
  LTA * d2hwlta = LTAalloc(1,mri_mov[1]);
  if (!vox2vox) // do ras to ras
  {
    // cout << "converting VOX to RAS and saving RAS2RAS..." << endl ;
    // (use geometry of destination space for half-way)
    m2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform (mri_mov[0],
                                                      mri_mov[1],
                                                      maps2weights.first,
                                                      m2hwlta->xforms[0].m_L) ;
    m2hwlta->type = LINEAR_RAS_TO_RAS ;
    d2hwlta->xforms[0].m_L = MRIvoxelXformToRasXform (mri_mov[1],
                                                      mri_mov[1],
                                                      maps2weights.second,
                                                      d2hwlta->xforms[0].m_L) ;
    d2hwlta->type = LINEAR_RAS_TO_RAS ;
  }
  else // vox to vox
  {
    // cout << "saving VOX2VOX..." << endl ;
    m2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.first,
                                        m2hwlta->xforms[0].m_L) ;
    m2hwlta->type = LINEAR_VOX_TO_VOX ;
    d2hwlta->xforms[0].m_L = MatrixCopy(maps2weights.second,
                                        m2hwlta->xforms[0].m_L) ;
    d2hwlta->type = LINEAR_VOX_TO_VOX ;
  }
  // add src and dst info (use src as target geometry in both cases)
  getVolGeom(mri_mov[0], &m2hwlta->xforms[0].src);
  getVolGeom(mri_mov[0], &m2hwlta->xforms[0].dst);
  getVolGeom(mri_mov[1], &d2hwlta->xforms[0].src);
  getVolGeom(mri_mov[0], &d2hwlta->xforms[0].dst);

  ltas.resize(2,NULL);
  if (ltas[0]) LTAfree(&ltas[0]);
  if (ltas[1]) LTAfree(&ltas[1]);
  ltas[0] = m2hwlta;
  ltas[1] = d2hwlta;

  cout << " creating half-way template ..." << endl;
  // take dst info from lta:
  if (mri_warps[0]) MRIfree(&mri_warps[0]);
  if (mri_warps[1]) MRIfree(&mri_warps[1]);
  mri_warps[0] = LTAtransform(mri_mov[0],NULL, m2hwlta);
  mri_warps[1] = LTAtransform(mri_mov[1],NULL, d2hwlta);
  // here do scaling of intensity values
  if (R.isIscale() && Md.second > 0)
  {
    cout << "Adjusting Intensity of MOV by " << Md.second << endl;
    mri_warps[0] = MyMRI::MRIvalscale(mri_warps[0], mri_warps[0], Md.second);
  }
  mri_mean = averageSet(mri_warps, mri_mean,average,sat);


  mri_weights[0] = R.getWeights(); // in original half way space
  mri_weights[1] = mri_weights[0];


  // copy weights for both if we have them (as the R->weights will be destroyed):
  if (mri_weights[0])
    for (unsigned int i = 0;i<mri_weights.size();i++)
      mri_weights[i] = MRIcopy(mri_weights[0],NULL);
			
  MatrixFree(&Md.first);
	return true;
}


bool MultiRegistration::initialXforms(int tpi, int maxres, int iterate, double epsit)
// will set ltas (as RAS to RAS )
// tpi 1....n
// uses outdir
{
  if (tpi <= 0) return false; // no init?
  cout << endl << "MultiRegistration::initializing Xforms (init " << tpi << " , maxres "<< maxres << " , iterate " << iterate << " , epsit " << epsit<< " ) : " << endl;
  assert(ltas.size() == mri_mov.size());
  assert(mri_mov.size() > 1);
  
  int nin = (int) mri_mov.size();
	
  tpi--; // tpi 0....n-1
  assert (tpi>=0 && tpi <nin);
	
  // create index set for lookup (tpi should switch places with 0)
  vector < int > index(nin);
  for (int i = 0;i<nin;i++) index[i] = i;
  index[0] = tpi;
  index[tpi] = 0;


  // Register everything to first TP (coarse)
  vector < Registration > Rv(nin);
  vector < std::pair <MATRIX*, double> > Md(nin);
  Md[0].first = MatrixIdentity(4,NULL);
  Md[0].second= 1.0;
  for (int i = 1;i<nin;i++) 
  {
	  int j = index[i]; // use new index
    cout << endl << "  computing initial registration of TP "<< j+1 <<" ( "<<mov[j]<<" )" << endl;
		cout << "   wrt to TP "<<tpi+1<<" ( "<<mov[tpi]<<" )" << endl;

    ostringstream oss;
    oss << outdir << "tp" << j+1 << "_to_tp" << tpi;
    
    Registration R;
    R.setSource(mri_mov[j],fixvoxel,fixtype);
    initRegistration(R); //set parameter
    R.setTarget(mri_mov[tpi], fixvoxel, fixtype);
    R.setName(oss.str());
		
    // compute Alignment (maxres,iterate,epsit) are passed above
    Md[i] = R.computeMultiresRegistration(maxres,iterate,epsit);

//    ostringstream oss2;
//    oss2 << "tp" << j+1;    
//    MatrixWrite(Md[i].first,(oss2.str()+".fsmat").c_str(),oss2.str().c_str());
//    Md[i].first = MatrixRead((oss2.str()+".fsmat").c_str());
//    Md[i].second = 1.0;

//     if (P.debug)
//     {
//       LTA * nlta = VOXmatrix2LTA(Md[i].first,P.mri_mov[j],P.mri_mov[tpi]);
//       LTAwriteEx(nlta, (oss.str()+".lta").c_str()) ;
// 
//       MRI* warped = MRIclone(P.mri_mov[tpi],NULL);
//       warped = LTAtransform(P.mri_mov[j],warped, nlta);
// 
//       if (R.isIscale() && Md[i].second >0)
//       {
//         string fn = oss.str() + "-intensity.txt";
//         ofstream f(fn.c_str(),ios::out);
//         f << Md[i].second;
//         f.close();
// 	
//         warped = R.MRIvalscale(warped,warped, Md[i].second);
//       }
//       MRIwrite(warped, (oss.str()+".mgz").c_str()) ;
//       // todo: maybe output weights (not necessary so far) 
//       LTAfree(&nlta);
//       MRIfree(&warped);     
//     }
  }
  
  // find mean center
  vector <MATRIX* > mras(nin);
  mras[0] = MatrixIdentity(4,NULL);
  MATRIX* rot = MatrixIdentity(3,NULL);
  MATRIX* roti = NULL;
  MATRIX* trans = VectorAlloc(3,MATRIX_REAL);
  MATRIX* meant = MatrixZero(3,1,NULL);
  MATRIX* meanr = MatrixIdentity(3,NULL);
  for (int i = 1;i<nin;i++) 
  {
	  int j = index[i];
    cout << "  computing coord of TP "<< j+1 <<" ( "<<mov[j]<<" )" << endl;
		cout << "   wrt to TP "<<tpi+1<<" ( "<<mov[tpi]<<" )" << endl;
    mras[i] = MRIvoxelXformToRasXform (mri_mov[j],mri_mov[tpi],Md[i].first,mras[i]);
    //MatrixPrintFmt(stdout,"% 2.8f",mras);cout << endl;
    // split into rotation translation
    MyMatrix::getRTfromM(mras[i],rot,trans);
   //MatrixPrintFmt(stdout,"% 2.8f",trans); cout << endl;
   //cout << " Mras: " << endl;
   //MatrixPrintFmt(stdout,"% 2.8f",mras); cout << endl;

    // reverse order (first translate, then rotate)
    // rot stays, trans changes: Rx+t = R (x + R^{-1} t)
    roti  = MatrixTranspose(rot,roti); // inverse
    //cout << " roti: " << endl;MatrixPrintFmt(stdout,"% 2.8f",roti); cout << endl;
    
    trans = MatrixMultiply(roti,trans,trans);
    //cout << "transi: - " << endl; MatrixPrintFmt(stdout,"% 2.8f",trans[i]); cout << endl;
    
//     if (P.debug) // output transonly
//     {
//       ostringstream oss;
//       oss << outdir << "tp" << j+1 << "_to_tp"<<tpi<<"-transonly";
//       MATRIX *Mtemp = getMfromRT(NULL,trans,NULL);
//       LTA * nlta = RASmatrix2LTA(Mtemp,P.mri_mov[j],P.mri_mov[tpi]);
//       MRI* warped = LTAtransform(P.mri_mov[j],NULL, nlta);
//       MRIwrite(warped, (oss.str()+".mgz").c_str()) ;
//       LTAfree(&nlta);
//       MRIfree(&warped);
//       MatrixFree(&Mtemp);         
//     }
    meant = MatrixSubtract(meant,trans,meant);
    meanr = MatrixAdd(meanr,roti,meanr);
  }
  
  MatrixFree(&roti);
  MatrixFree(&trans);
  MatrixFree(&rot);

  //average
  meant = MatrixScalarMul(meant,1.0/nin,meant);
  //cout << "meant: " << endl; MatrixPrintFmt(stdout,"% 2.8f",meant); cout << endl;
  meanr = MatrixScalarMul(meanr,1.0/nin,meanr);
  //cout << "meanr: " << endl;MatrixPrintFmt(stdout,"% 2.8f",meanr); cout << endl;
  
  // project meanr back to SO(3) (using polar decomposition)
  VECTOR * vz = VectorAlloc(3,MATRIX_REAL);
  MATRIX * mv = MatrixSVD(meanr,vz,NULL);
  //MatrixPrintFmt(stdout,"% 2.8f",meanr); cout << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",vz); cout << endl;
  //MatrixPrintFmt(stdout,"% 2.8f",mv); cout << endl;
  mv = MatrixTranspose(mv,mv);
  //MatrixPrintFmt(stdout,"% 2.8f",mv); cout << endl;
  meanr = MatrixMultiply(meanr,mv,meanr);
  //cout << " meanr.proj: " << endl;MatrixPrintFmt(stdout,"% 2.8f",meanr); cout << endl;
  MatrixFree(&mv);
  MatrixFree(&vz);

  // construct maps from each image to the mean space
  MATRIX * Mm = MyMatrix::getMfromRT(meanr,meant,NULL);
  MATRIX * M = NULL;
  //cout << " Mm: " << endl; MatrixPrintFmt(stdout,"% 2.8f",Mm); cout << endl;
  for (int i = 0;i<nin;i++) 
  {
	  int j = index[i];
    cout << "  computing mean coord of TP "<< j+1 <<" ( "<<mov[j]<<" ) " << endl;
    // concat transforms: meant meanr mras
    // from right to left: mras[j] (aligns to T1)
    //                     then do the mean rot and move to mean location

    //MatrixPrintFmt(stdout,"% 2.8f",mras[i]); cout << endl;
    M = MatrixMultiply(Mm,mras[i],M);
    //MatrixPrintFmt(stdout,"% 2.8f",M); cout << endl;
    
    // make lta from M (M is RAS to RAS)
    assert(ltas.size() == mri_mov.size());
    assert(ltas[j] == NULL);
    ltas[j] = MyMatrix::RASmatrix2LTA(M,mri_mov[j],mri_mov[tpi]); // use geometry of tpi

    // store intensities
		if (iscale)
		{
      if (i==0) intensities[j] = 1.0;
  		else
  		{
  		   assert (Md[i].second > 0);
  		   intensities[j] = Md[i].second;
  		}
	  }
			
    if (debug)
    {
      ostringstream oss;
      oss << outdir << "tp" << j+1 << "_to_mcoord";
      LTAwriteEx(ltas[j], (oss.str()+".lta").c_str()) ;
      MRI * warped = LTAtransform(mri_mov[j],NULL, ltas[j]);
      if (iscale)
      {
        string fn = oss.str() + "-intensity.txt";
        ofstream f(fn.c_str(),ios::out);
        f << Md[i].second;
        f.close();
	
        warped = MyMRI::MRIvalscale(warped,warped, Md[i].second);
      }
      MRIwrite(warped, (oss.str()+".mgz").c_str()) ;
      MRIfree(&warped);     
    }
    
    // cleanup
    MatrixFree(&mras[i]);
    MatrixFree(&Md[i].first);
  }
  //cleanup
  MatrixFree(&M);
  MatrixFree(&Mm);
  return true;
}

bool MultiRegistration::writeMean(const std::string& mean)
{
	    if (! mri_mean ) 
			{
			   cout << " ERROR: No average exists! Skipping output." << endl;
				 return false;
			}
  strncpy(mri_mean->fname, mean.c_str(),STRLEN);
  return (MRIwrite(mri_mean,mean.c_str()) == 0);
}

bool MultiRegistration::writeLTAs(const std::vector < std::string > & nltas, bool vox2vox,const string & mean)
{
   assert (nltas.size() == ltas.size());
	 int error = 0;
   for (unsigned int i = 0;i<nltas.size();i++)
	 {
	    if (! ltas[i] ) 
			{
			   cout << " ERROR: No ltas exist! Skipping output." << endl;
				 return false;
			}

      if (vox2vox)
      {
        error += (LTAchangeType(ltas[i], LINEAR_VOX_TO_VOX)==NULL);
      }
      else assert(ltas[i]->type == LINEAR_RAS_TO_RAS);
      strncpy(ltas[i]->xforms[0].dst.fname, mean.c_str(),STRLEN);
      LTAwriteEx(ltas[i], nltas[i].c_str()) ;	 
	 }
	 return (error == 0);
}

bool MultiRegistration::writeWarps(const std::vector <  std::string >& nwarps)
{
   assert (nwarps.size() == mri_warps.size());
	 int error = 0;
   for (unsigned int i = 0;i<nwarps.size();i++)
	 {
	    if (! mri_warps[i] ) 
			{
			   cout << " ERROR: No warps exist! Skipping output." << endl;
				 return false;
			}
		  error += MRIwrite(mri_warps[i],nwarps[i].c_str()) ;
	 }
	 return (error == 0);
}

bool MultiRegistration::writeIntensities(const std::vector < std::string >& nintens)
{
   assert (nintens.size() == intensities.size());
	 int error = 0;
   for (unsigned int i = 0;i<nintens.size();i++)
	 {
	     assert (intensities[i] > 0);
          ofstream f(nintens[i].c_str(),ios::out);
					if (f.good())
					{
            f << intensities[i];
            f.close();
					} else error++;
	 }
	 return (error == 0);
}

bool MultiRegistration::writeWeights(const std::vector < std::string >& nweights)
{
   assert (nweights.size() == mri_weights.size());
	 int error = 0;
   for (unsigned int i = 0;i<nweights.size();i++)
	 {
	    if (! mri_weights[i])
			{
			   cout << " No weights constructed, skipping output of weights" << endl;
				 return false;
			}
		  error += MRIwrite(mri_weights[i], nweights[i].c_str()) ;
	 }
	 return (error == 0);
}


/*!
  \fn MRI* averageSet(const vector < MRI * >& set, MRI* mean, int method, double sat)
  \brief Averages the movable volumes depending on method
  \param set vector of movable volumes
  \param mean  output of mean volume
  \param method  0 = mean, 1 = median, 2 = tukey biweight (testing)
  \param sat     saturation for tukey biweight
*/
MRI* MultiRegistration::averageSet(const vector < MRI * >& set,
                       MRI* mean, int method, double sat)
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
    for (unsigned int i = 0;i<set.size();i++)
    {
      mean = MRIaverage(set[i],i,mean);
    }
  }
  else if (method ==1)
  {
    cout << " using median " << endl;
    // robust
    int x,y,z,i;
    assert(set.size() > 0);
    double dd[set.size()];
    if (!mean) mean = MRIclone(set[0],NULL);
    //MRI * midx = MRIalloc(set[0]->width,set[0]->height,set[0]->depth, MRI_FLOAT);
    //MRIcopyHeader(mean,midx);
    //midx->type = MRI_FLOAT;
    pair < double, double > mm;
    for (z = 0 ; z < set[0]->depth ; z++)
      for (y = 0 ; y < set[0]->height ; y++)
        for (x = 0 ; x < set[0]->width ; x++)
        {
          for (i=0; i<(int) set.size();i++)
            dd[i] = MRIgetVoxVal(set[i],x,y,z,0);
          mm = RobustGaussian::medianI(dd,(int)set.size());
          MRIsetVoxVal(mean,x,y,z,0,mm.first);
	  //MRIsetVoxVal(midx,x,y,z,0,mm.second);
        }
    //MRIwrite(midx,"midx.mgz");
    //MRIwrite(mean,"m.mgz");
    //assert(1==2);
  }
  else if (method ==2)
  {
    cout << "    using tukey biweight" << endl;
    // robust tukey biweight
    int x,y,z,i;
    MATRIX* b = MatrixAlloc(set.size(),1,MATRIX_REAL);
    pair < MATRIX* , MATRIX* >  muw;
    assert(set.size() > 0);
    if (!mean) mean = MRIclone(set[0],NULL);
    for (z = 0 ; z < set[0]->depth ; z++)
      for (y = 0 ; y < set[0]->height ; y++)
        for (x = 0 ; x < set[0]->width ; x++)
        {
          // cout << " x: " << x << " y: " << y << " z: " <<z << "  size: " << set.size() <<endl;
          for (i=0; i<(int) set.size();i++)
          {
            //cout << "i: " << i << endl;
            b->rptr[i+1][1] = MRIgetVoxVal(set[i],x,y,z,0);
          }
          //cout << endl << "intgensities at this voxel:" << endl; ;
          //MatrixPrintFmt(stdout,"% 2.8f",b);

          Regression R(b);
          muw = R.getRobustEstW(sat);
          //    cout << " tukey mean: " << muw.first->rptr[1][1] << endl;
          MRIsetVoxVal(mean,x,y,z,0,muw.first->rptr[1][1]);
          MatrixFree(&muw.first);
          MatrixFree(&muw.second);
        }
  }
  else if (method ==3) // needs more development (sigma..)
  {
    cout << "    using tukey biweight (alltoghether)" << endl;
    // robust tukey biweight
    int x,y,z,i;
    int n = set[0]->depth * set[0]->height *  set[0]->width ;
    MATRIX* b = MatrixAlloc(set.size()*n,1,MATRIX_REAL);
    MATRIX* A = MatrixAlloc(set.size()*n,n,MATRIX_REAL);
    A = MatrixZero(A->rows,A->cols,A);
    pair < MATRIX* , MATRIX* >  muw;
    assert(set.size() > 0);
    if (!mean) mean = MRIclone(set[0],NULL);
    for (i=0; i<(int) set.size();i++)
    {
      assert(set[i]->width == set[0]->width);
      assert(set[i]->height == set[0]->height);
      assert(set[i]->depth == set[0]->depth);
      int pcount = 0;
      for (z = 0 ; z < set[0]->depth ; z++)
        for (y = 0 ; y < set[0]->height ; y++)
          for (x = 0 ; x < set[0]->width ; x++)
          {
            // cout << " x: " << x << " y: " << y << " z: " <<z << "  size: " << set.size() <<endl;
            //cout << "i: " << i << endl;
            b->rptr[pcount*set.size()+i+1][1] =
              MRIgetVoxVal(set[i],x,y,z,0);
            A->rptr[pcount*set.size()+i+1][pcount+1] = 1;
            pcount++;
          }
    }
    //cout << endl << "intgensities at this voxel:" << endl; ;
    //MatrixPrintFmt(stdout,"% 2.8f",b);

    Regression R(A,b);
    muw = R.getRobustEstW(sat);
    //    cout << " tukey mean: " << muw.first->rptr[1][1] << endl;
    int pcount = 0;
    for (z = 0 ; z < set[0]->depth ; z++)
      for (y = 0 ; y < set[0]->height ; y++)
        for (x = 0 ; x < set[0]->width ; x++)
        {
          pcount++;
          MRIsetVoxVal(mean,x,y,z,0,muw.first->rptr[pcount][1]);
        }
    MatrixFree(&muw.first);
    MatrixFree(&muw.second);

  }
  else
  {

    cerr <<  " averageSet  method " << method << " unknown" << endl;
    exit(1);
  }
//       cout << " Average Vol "  << endl;
//       cout << " mean   : " << CostFunctions::mean(mean) << endl;
//       cout << " sdev   : " << CostFunctions::sdev(mean) << endl;
//       cout << " median : " << CostFunctions::median(mean) << endl;
//       cout << " mad    : " << CostFunctions::mad(mean)<< endl;

  return mean;
}

/*!
  \fn MRI* initialAverageSet(const vector < MRI * >& set, MRI* mean, int method, double sat)
  \brief Averages the movable volumes depending on method
  \param set vector of movable volumes
  \param mean  output of mean volume
  \param method  0 = mean, 1 = median, 2 = tukey biweight (testing)
  \param sat     saturation for tukey biweight
*/
MRI* MultiRegistration::initialAverageSet(const vector < MRI * >& set,
                              MRI* mean, int method, double sat)
{
// the initial average can be any place as long as it is
// not identical to one of the
// initial images (to avoid introducing a bias).
// it can be the average of all images but that is extremely blurry, and if
// we have only few images, far apart, maybe the algorithm will not converge
// or converge very slowly (has not happened yet).
// therefore we find an initial coarse alignment by using moments

  assert(set.size() > 1);
  int n = (int) set.size();
  vector < double >  centroid(3,0);
  vector < vector < double > >centroids(n);
  for (int i = 0;i<n;i++)
  {
    centroids[i] = CostFunctions::centroid(set[i]);
    centroid[0] += centroids[i][0];
    centroid[1] += centroids[i][1];
    centroid[2] += centroids[i][2];
  }
  centroid[0] /= n;
  centroid[1] /= n;
  centroid[2] /= n;

  // map set to mean of centroids
  MATRIX* Mtrans = MatrixIdentity(4,NULL);
  vector < MRI* > newset(n);
  for (int i = 0;i<n;i++)
  {
    Mtrans->rptr[1][4] = centroid[0]-centroids[i][0];
    Mtrans->rptr[2][4] = centroid[1]-centroids[i][1];
    Mtrans->rptr[3][4] = centroid[2]-centroids[i][2];
    newset[i] = MRIlinearTransform(set[i],NULL,Mtrans);
  }

  return averageSet(newset,mean,method,sat);
}
