/**
 * @file  gcamregisterlevel_cuda.cu
 * @brief Implementation of GCAMregisterLevel for the GPU
 *
 * Reference:
  * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
  * Structures in the Human Brain", Fischl et al.
  * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.7 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifdef GCAMORPH_ON_GPU

#include "macros.h"
#include "error.h"

#include "gcamorph.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphenergy.hpp"

#include "gcamcomputegradient_cuda.hpp"
#include "gcamfots_cuda.hpp"

#include "gcamremovenegativenodes_cuda.hpp"

#include "gcamregisterlevel_cuda.hpp"

// ========================================================================

template<typename T, typename U>
int RegisterLevel( GPU::Classes::GCAmorphGPU& gcam,
                   const GPU::Classes::MRIframeGPU<T>& mri,
                   const GPU::Classes::MRIframeGPU<U>& mri_smooth,
                   GCA_MORPH_PARMS *parms )
{
  GPU::Algorithms::GCAmorphEnergy gcamEnergy;

  int which = GCAM_INTEGRATE_OPTIMAL;
  int done = 0;

  int n, nsmall, increasing;
  int good_step;
  int max_small;
  double pct_change, last_pct_change, tol, rms, last_rms;
  double orig_j;
  double orig_dt, min_dt;

  int myGinvalid = 0;

  if( parms->integration_type == GCAM_INTEGRATE_FIXED )
  {
    which = GCAM_INTEGRATE_FIXED;
  }

  max_small = parms->nsmall;

  gcam.ClearMomentum();


  orig_dt = parms->dt ;
  nsmall = 0 ;
  pct_change = 0.0 ;

  if( parms->integration_type == GCAM_INTEGRATE_FIXED )
  {
    increasing = 0;  /* will be set to 0 if new step is
                              smaller than last */
  }
  else
  {
    increasing = 0;  /* don't use it for optimal time-step stuff */
  }



  if( parms->uncompress )
  {
    // gcamRemoveCompressedNodes(gcam, mri, parms, parms->ratio_thresh);
    std::cerr << __FUNCTION__
              << ": gcamRemoveCompressedNodes not implemented on GPU"
              << std::endl;
    abort();
  }

  last_rms = gcamEnergy.ComputeRMS( gcam, mri, parms );


  {
    printf("%04d: dt=%2.3f, rms=%2.3f, neg=%d, invalid=%d\n",
           0, 0.0f, last_rms, gcam.neg, myGinvalid) ;

    fflush(stdout);
  }

  orig_j = parms->l_jacobian ;
  tol = parms->tol ;
  good_step = 0;


  for( n = parms->start_t; n<parms->start_t+parms->niterations; n++)
  {


    gcam.RemoveStatus( GCAM_LABEL_NODE );
    gcam.RemoveStatus( GCAM_IGNORE_LIKELIHOOD );

    ComputeGradient( gcam, mri, mri_smooth, parms );

    parms->l_jacobian = orig_j;

    switch( parms->integration_type )
    {
    case GCAM_INTEGRATE_OPTIMAL:
      parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ; /* will search around
                                                         this value */
      min_dt = FindOptimalTimestep(gcam, parms, mri) ;
      parms->dt = min_dt ;
      break ;

    case GCAM_INTEGRATE_FIXED:
      min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
      break ;

    case GCAM_INTEGRATE_BOTH:
      if (which == GCAM_INTEGRATE_OPTIMAL)
      {
        parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ; /* will search around
                                                           this value */
        min_dt = FindOptimalTimestep( gcam, parms, mri);
        parms->dt = min_dt ;
        max_small = parms->nsmall ;
        tol = parms->tol ;
      }
      else
      {
        /* take some momentum steps */
        max_small = 2*parms->nsmall ;
        tol = parms->tol/2 ;
      }
      break ;

    default:
      min_dt = parms->dt ;
      std::cerr << __FUNCTION__
                << ": Unknown integration type = "
                << parms->integration_type
                << std::endl;
      abort();
    }

    gcam.CopyNodePositions( CURRENT_POSITIONS, SAVED2_POSITIONS );
    gcam.ApplyGradient( parms );
    gcam.ComputeMetricProperties( myGinvalid );

    if( parms->constrain_jacobian )
    {
      //gcamConstrainJacobian(gcam, mri, parms) ;
      std::cerr << __FUNCTION__
                << ": gcamConstrainJacobian not implemented on GPU"
                << std::endl;
      abort();
    }

    if (gcam.neg > 0 && parms->noneg == True)
    {
      int i = 0;

      gcam.CopyNodePositions( SAVED2_POSITIONS, CURRENT_POSITIONS );
      gcam.ClearMomentum();
      gcam.ComputeMetricProperties( myGinvalid );
      gcam.ApplyGradient( parms );
      gcam.ComputeMetricProperties( myGinvalid );
      while( gcam.neg > 0 )
      {

        parms->dt *= 0.5 ;


        gcam.UndoGradient();
        gcam.ComputeMetricProperties( myGinvalid );
        gcam.ApplyGradient( parms );
        gcam.ComputeMetricProperties( myGinvalid );

        if( ++i > 50 )
        {
          if( gcam.neg > 0 )
          {
            /* couldn't find a  step without folds */
            gcam.CopyNodePositions( SAVED2_POSITIONS, CURRENT_POSITIONS );
            gcam.ComputeMetricProperties( myGinvalid );
          }
          break ;
        }
      }
    }
    min_dt = parms->dt;

    RemoveNegativeNodes( gcam, mri, parms );


    if( parms->uncompress )
    {
      //gcamRemoveCompressedNodes(gcam, mri, parms, parms->ratio_thresh);
      std::cerr << __FUNCTION__
                << ": gcamRemoveCompressedNodes not implemented on GPU"
                << std::endl;
    }

    if (gcam.neg > 0)
    {
      std::cout << "---------- unfolding failed"
                << " - restoring original position --------------------"
                << std::endl;

      gcam.CopyNodePositions( SAVED2_POSITIONS, CURRENT_POSITIONS );
      gcam.ComputeMetricProperties( myGinvalid );
    }



    rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
    last_pct_change = pct_change;

    if( FZERO(last_rms) )
    {
      pct_change = 0.0;
    }
    else
    {
      pct_change = 100.0*(last_rms-rms)/last_rms;
    }

    if( (pct_change < last_pct_change) || FZERO(pct_change) )
    {
      increasing = 0 ;
    }
    else
    {
      /* could check for last_pct_change == 0 here */
      increasing = 1 ;
    }


    if( pct_change <= 0 )
    {

      increasing = 0 ;
      if (parms->constrain_jacobian == 0)
      {
        if (parms->uncompress == False)
        {
          // otherwise it could be the
          // uncompressing that caused the sse to dec.
          done = 1 ;
        }

        gcam.CopyNodePositions( SAVED2_POSITIONS, CURRENT_POSITIONS );
        gcam.ComputeMetricProperties( myGinvalid );
        rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
      }
    }



    //  if (Gdiag & DIAG_SHOW)
    {
      printf("%04d: dt=%2.6f, rms=%2.3f (%2.3f%%), neg=%d, invalid=%d\n",
             n+1, min_dt, rms, pct_change, gcam.neg, myGinvalid ) ;
      fflush(stdout);
    }


    if ((pct_change < tol) && !increasing)
    {

      // if we ever took a good step since the last regridding,
      // and we tried to uncomress and failed, regrid

      if( (++nsmall >= max_small) || (pct_change <= 0) )
      {
        if (parms->integration_type == GCAM_INTEGRATE_BOTH)
        {
          if (!good_step)
          {
            done++;
          }

          if( done >= 2 )
          {
            /* couldn't take a step with either technique */
            n++ ;
            break ;

          }
          else
          {
            /* switch integration types */

            if( which == GCAM_INTEGRATE_FIXED )
            {
              increasing = 0 ;
              which = GCAM_INTEGRATE_OPTIMAL ;
              parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
              if (parms->uncompress)
              {
                std::cerr << __FUNCTION__
                          << ": gcamRemoveCompressedNodes unimplemented on GPU"
                          << std::endl;
                abort();
                //gcamRemoveCompressedNodes(gcam, mri,parms,parms->ratio_thresh);
                /* will search around this value */
              }
            }
            else
            {


              increasing = /*1 */0;
              pct_change = 0.0 ;
              which = GCAM_INTEGRATE_FIXED ;

              if( !DZERO(min_dt) )
              {
                parms->dt = min_dt;
              }
              else
              {
                min_dt = parms->dt = (sqrt(parms->navgs)+1.0f)*orig_dt ;
              }

            }
            good_step = nsmall = 0 ;
            gcam.ClearMomentum();
          }
        }
        else
        {
          n++;
          break ;
        }
      }
    }
    else if ( pct_change >= tol )
    {
      /* took at least one good step */

      good_step = 1 ;
      done = 0 ;    /* for integration type == BOTH, apply both types again */
      nsmall = 0 ;  /* start counting small steps again */
    }

    last_rms = rms ;
    /*          parms->label_dist -= 0.5 ;*/
    if (parms->label_dist < 0 )
    {
      parms->label_dist = 0 ;
    }

  }

  parms->start_t = n;
  parms->dt = orig_dt;

  SetGinvalid( myGinvalid );

  return( NO_ERROR );
}


// ==============





template<typename T, typename U>
void
gcamRLfinalDispatch( GCA_MORPH *gcam,
                     const MRI *mri,
                     const MRI *mri_smooth,
                     GCA_MORPH_PARMS *parms )
{

  GPU::Classes::GCAmorphGPU myGCAM;
  GPU::Classes::MRIframeGPU<T> myMRI;
  GPU::Classes::MRIframeGPU<U> myMRIsmooth;

  // Handle the MRIs
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );

  myMRIsmooth.Allocate( mri_smooth );
  myMRIsmooth.Send( mri_smooth, 0 );

  // Put the GCAM on the GPU
  myGCAM.CheckIntegrity(); // Shouldn't be necessary....
  myGCAM.SendAll( gcam );

  // Run the computation
  RegisterLevel( myGCAM, myMRI, myMRIsmooth, parms );

  // Retrieve results
  myGCAM.RecvAll( gcam );
}


template<typename T>
void
gcamRLsmoothDispatch( GCA_MORPH *gcam,
                      const MRI *mri,
                      const MRI *mri_smooth,
                      GCA_MORPH_PARMS *parms )
{

  switch( mri_smooth->type )
  {

  case MRI_UCHAR:
    gcamRLfinalDispatch<T,unsigned char>( gcam, mri, mri_smooth, parms );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri_smooth "
              << mri_smooth->type << std::endl;
    abort();
  }

}





void gcamRegisterLevelGPU( GCA_MORPH *gcam,
                           const MRI *mri,
                           const MRI *mri_smooth,
                           GCA_MORPH_PARMS *parms )
{

  switch( mri->type )
  {

  case MRI_UCHAR:
    gcamRLsmoothDispatch<unsigned char>( gcam, mri, mri_smooth, parms );
    break;

  case MRI_FLOAT:
    gcamRLsmoothDispatch<float>( gcam, mri, mri_smooth, parms );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised type for mri "
              << mri->type << std::endl;
    abort();
  }

}




#endif
