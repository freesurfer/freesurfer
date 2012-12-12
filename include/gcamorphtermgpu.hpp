/**
 * @file  gcamorphtermgpu.hpp
 * @brief Holds routines to compute GCAmorph terms on the GPU (header)
 *
 * This file holds a variety of routines which compute terms for
 * mri_ca_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.24 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef GCA_MORPH_TERM_GPU_HPP
#define GCA_MORPH_TERM_GPU_HPP

#include "macros.h"
#include "error.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "ctfactory.hpp"

#include "gcamorphgpu.hpp"
#include "gcamorphcpu.hpp"

namespace GPU
{
namespace Algorithms
{

//! Class to hold GCAMorph 'term' computations
class GCAmorphTerm
{

public:

  //! Constructor (stream defaults to zero)
  GCAmorphTerm( const cudaStream_t s = 0 ) : stream(s) {};

  //! Destructor
  ~GCAmorphTerm( void ) {};

  //! Routine to print timing information to std.out
  static void ShowTimings( void );


  //! Computes the smoothness term
  void Smoothness( GPU::Classes::GCAmorphGPU& gcam,
                   const float l_smoothness ) const;

  //! Computes the Jacobian term
  void Jacobian( GPU::Classes::GCAmorphGPU& gcam,
                 const float l_jacobian,
                 const float jac_scale ) const;

  //! Computes the Log Likelihood term
  template<typename T, typename U>
  void LogLikelihood( GPU::Classes::GCAmorphGPU& gcam,
                      const GPU::Classes::MRIframeGPU<T>& mri,
                      const GPU::Classes::MRIframeGPU<U>& mri_smooth,
                      double l_log_likelihood ) const;

  //! Dispatch routine for the Log Likelihood term
  template<typename T, typename U>
  void LLtermDispatch( GCA_MORPH *gcam,
                       const MRI*  mri,
                       const MRI* mri_smooth,
                       double l_log_likelihood ) const;

  //! Basic Dispatch routine for LLT
  void LLTDispatch( GCA_MORPH *gcam,
                    const MRI*  mri,
                    const MRI* mri_smooth,
                    double l_log_likelihood ) const;


  template<typename T>
  GPU::Classes::CTfactory* BindMRI( const GPU::Classes::MRIframeGPU<T>& mri ) const;

  template<typename T>
  GPU::Classes::CTfactory* BindMRIsmooth( const GPU::Classes::MRIframeGPU<T>& mri ) const;



  //! Posterior/Anterior consistency check for Label term
  int LabelPostAntConsistency( GPU::Classes::GCAmorphGPU& gcam,
                               GPU::Classes::MRIframeGPU<float>& mri_dist ) const;

  //! Final update stage for Label term
  int LabelFinalUpdate( GPU::Classes::GCAmorphGPU& gcam,
                        const GPU::Classes::MRIframeGPU<float>& mri_dist,
                        const float l_label ) const;


  //! Copy deltas for Label term
  void LabelCopyDeltas( GPU::Classes::GCAmorphGPU& gcam,
                        const GPU::Classes::MRIframeGPU<float>& mri_dist,
                        const float l_label ) const;


  //! Wrapper for Remove Label Outliers
  int RemoveLabelOutliersDispatch( GPU::Classes::GCAmorphGPU& gcam,
                                   MRI *mri_dist,
                                   const int whalf,
                                   const double thresh ) const;

  //! Wrapper for Label Main Loop
  void LabelMainLoopDispatch( GPU::Classes::GCAmorphGPU& gcam,
                              const MRI *mri,
                              MRI *mri_dist,
                              const double l_label,
                              const double label_dist ) const;


  //! Label Term for the GPU
  template<typename T>
  int LabelTerm( GPU::Classes::GCAmorphGPU& gcam,
                 const GPU::Classes::MRIframeGPU<T>& mri,
                 double l_label, double label_dist ) const
  {

    int num;
    MRI *mri_dist, *mriCPU;
    int nremoved;

    if( DZERO(l_label) )
    {
      return( NO_ERROR );
    }

    GCAmorphTerm::tLabelTot.Start();

    gcam.CheckIntegrity();

    const dim3 gcamDims = gcam.d_rx.GetDims();

    mri_dist = MRIalloc( gcamDims.x, gcamDims.y, gcamDims.z, MRI_FLOAT );
    //MRIsetResolution( mri_dist, gcam.spacing, gcam.spacing, gcam.spacing );
    std::cerr << __FUNCTION__ << ": Did not call MRIsetResolution" << std::endl;

    // GCAMresetLabelNodeStatus
    gcam.RemoveStatus( GCAM_LABEL_NODE );
    gcam.RemoveStatus( GCAM_IGNORE_LIKELIHOOD );



    // Set up the CPU copies
    static Freesurfer::GCAmorphCPU gcamCPU;
    gcamCPU.AllocateFromTemplate( gcam );
    gcamCPU.GetFromGPU( gcam );

    const dim3 mriDims = mri.GetDims();
    mriCPU = MRIalloc( mriDims.x, mriDims.y, mriDims.z, mri.MRItype() );
    mri.Recv( mriCPU, 0 );

    // Do the main loop
    this->LabelMainLoop( gcamCPU, mriCPU, mri_dist, l_label, label_dist );



    // Remove the outliers
    nremoved = RemoveLabelOutliers( gcamCPU, mri_dist, 2, 3 );




    SetInconsistentLabelNodes( nremoved );

    // Put data back on the GPU
    gcamCPU.PutOnGPU( gcam );

    // Put mri_dist on the GPU
    GPU::Classes::MRIframeGPU<float> mriDistGPU;
    mriDistGPU.Allocate( mri_dist );
    mriDistGPU.Send( mri_dist, 0 );

    // Copy the deltas
    this->LabelCopyDeltas( gcam, mriDistGPU, l_label );

    /*
      This call causes differences between the CPU and GPU.
      See notes in corresponding kernel.
    */
    nremoved += this->LabelPostAntConsistency( gcam, mriDistGPU );



    num = this->LabelFinalUpdate( gcam, mriDistGPU, l_label );

    MRIfree( &mri_dist );
    MRIfree( &mriCPU );

    GCAmorphTerm::tLabelTot.Stop();

    // Shut up compiler warning
    if( false )
    {
      std::cout << __FUNCTION__ << num << std::endl;
    }

    return( NO_ERROR );
  }


  //! Dispatch wrapper from CPU types
  template<typename T>
  int LabelTermDispatch( GCA_MORPH *gcam, const MRI *mri,
                         double l_label, double label_dist ) const
  {

    int retVal;

    GPU::Classes::MRIframeGPU<T> mriGPU;

    mriGPU.Allocate( mri );
    mriGPU.Send( mri, 0 );

    GPU::Classes::GCAmorphGPU gcamGPU;
    gcamGPU.SendAll( gcam );

    retVal = this->LabelTerm( gcamGPU, mriGPU, l_label, label_dist );

    gcamGPU.RecvAll( gcam );


    return( retVal );
  }


  // ######################################################
private:

  //! Stream to use for operations
  cudaStream_t stream;

  //! Timer for smoothness term
  static SciGPU::Utilities::Chronometer tSmoothTot;
  //! Timer for smoothness term subtraction
  static SciGPU::Utilities::Chronometer tSmoothSubtract;
  //! Timer for smoothness computation itself
  static SciGPU::Utilities::Chronometer tSmoothCompute;

  //! Timer for the Jacobian term
  static SciGPU::Utilities::Chronometer tJacobTot;
  //! Timer for norm calculation of Jacobian term
  static SciGPU::Utilities::Chronometer tJacobMaxNorm;
  //! Timer for jacobian computation itself
  static SciGPU::Utilities::Chronometer tJacobCompute;

  //! Timer for Log likelihood term
  static SciGPU::Utilities::Chronometer tLogLikelihoodTot;
  //! Timer for Log likelihood term computation
  static SciGPU::Utilities::Chronometer tLogLikelihoodCompute;

  //! Timer for the main loop of LabelTerm
  static SciGPU::Utilities::Chronometer tLabelMainLoop;
  //! Timer for the computation portion of RemoveOutliers
  static SciGPU::Utilities::Chronometer tRemoveOutliers;
  //! Timer for copy deltas portion of LabelTerm
  static SciGPU::Utilities::Chronometer tLabelCopyDeltas;
  //! Timer for post/ant consistency check of LabelTerm
  static SciGPU::Utilities::Chronometer tLabelPostAntConsistency;
  //! Timer for final update of Label term
  static SciGPU::Utilities::Chronometer tLabelFinal;
  //! Timer for complete Label term
  static SciGPU::Utilities::Chronometer tLabelTot;

  // ---------------------

  template<typename T>
  void LLTmrismoothDispatch( GCA_MORPH *gcam,
                             const MRI*  mri,
                             const MRI* mri_smooth,
                             double l_log_likelihood ) const;

  //! Remove Label Outliers for Label term
  int RemoveLabelOutliers( Freesurfer::GCAmorphCPU& gcam,
                           MRI *mri_dist,
                           const int whalf,
                           const double thresh ) const;

  //! Main loop of Label term
  void LabelMainLoop( Freesurfer::GCAmorphCPU& gcam,
                      const MRI *mri,
                      MRI *mri_dist,
                      const double l_label,
                      const double label_dist ) const;

  //! Variation on is_temporal_wm
  static int IsTemporalWM( const MRI *mri,
                           const GCA_NODE *gcan,
                           float xf, float yf, float zf );
};

}
}

#endif

#endif
