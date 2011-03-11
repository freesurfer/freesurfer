/**
 * @file  SurfaceOverlayProperties.h
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:42 $
 *    $Revision: 1.1 $
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

#ifndef SurfaceOverlayProperties_h
#define SurfaceOverlayProperties_h

#include "vtkSmartPointer.h"
#include "Broadcaster.h"
#include <string>

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
class SurfaceOverlay;

class SurfaceOverlayProperties  : public Broadcaster
{
public:
  SurfaceOverlayProperties ( SurfaceOverlay* overlay);
  ~SurfaceOverlayProperties ();

  enum COLOR_SCALE  { CS_GreenRed = 0, CS_Heat, CS_BlueRed, CS_ColorWheel, CS_RYGBWheel, CS_TwoCondGR };
  enum COLOR_METHOD { CM_Linear = 0, CM_LinearOpaque, CM_Piecewise };
  
  double GetOpacity() const;
  void SetOpacity( double opacity );
  
  int GetColorScale() const;
  void SetColorScale( int nScale );

  void SetSurfaceOverlay( SurfaceOverlay* overlay );
  
  void SetMinPoint( double dValue );
  double GetMinPoint();
  
  void SetMidPoint( double dValue );
  double GetMidPoint();
  
  void SetMaxPoint( double dValue );
  double GetMaxPoint();
  
  int GetColorMethod();
  void SetColorMethod( int n );
  
  bool GetColorInverse();
  void SetColorInverse( bool bInverse );
  
  bool GetColorTruncate();
  void SetColorTruncate( bool bTruncate );
  
  void MapOverlayColor( float* data, unsigned char* colordata, int nPoints );
  void MapOverlayColor( unsigned char* colordata, int nPoints );
  
  vtkRGBAColorTransferFunction* GetLookupTable()
  {
    return m_lut;
  }
  
private:
  void ColorMapChanged();
  
  double      m_dOpacity;
  int         m_nColorScale;
  int         m_nColorMethod;
  double      m_dMinPoint;
  double      m_dMidPoint;
  double      m_dMaxPoint;
  int         m_colorMin[3];
  int         m_colorMid[3];
  int         m_colorMax[3];
  bool        m_bColorInverse;
  bool        m_bColorTruncate;
  
  vtkRGBAColorTransferFunction* m_lut;
  
  SurfaceOverlay* m_overlay;
};

#endif
