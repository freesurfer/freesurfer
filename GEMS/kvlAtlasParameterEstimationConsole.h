/**
 * @file  kvlAtlasParameterEstimationConsole.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.3 $
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
#ifndef __kvlAtlasParameterEstimationConsole_h
#define __kvlAtlasParameterEstimationConsole_h

#include "kvlAtlasParameterEstimationConsoleGUI.h"
#include <vector>
#include <string>
#include "kvlAtlasParameterEstimator.h"



namespace kvl
{

class AtlasParameterEstimationConsole : public kvlAtlasParameterEstimationConsoleGUI
{

public:

  //
  AtlasParameterEstimationConsole();

  //
  virtual ~AtlasParameterEstimationConsole();

  //
  void SetLabelImages( const std::vector< std::string >&  fileNames, bool useGaussians = false, bool ignoreLastImage = false );

  //
  void Show();

  //
  void Estimate();

  //
  void HandleEstimatorEvent( itk::Object* object, const itk::EventObject & event );

protected:

  //
  void DisplayLabelImage( unsigned int labelImageNumber );

#if 0
  //
  void SelectTriangleContainingPoint( float x, float y );
#endif

  //
  void InitializeMesh();

  //
  void Step();

  //
  void Interrupt();

  //
  void Continue();

  //
  void SetPositionEstimationResolution( unsigned int positionEstimationResolution );

private:

  AtlasParameterEstimator::Pointer  m_Estimator;

  std::string  m_TotalProgressLabel;
  std::string  m_SubProgressLabel;

  bool  m_Interrupted;
  bool  m_Stepping;


};



} // end namespace kvl

#endif
