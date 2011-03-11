/**
 * @file  SurfaceOverlay.h
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

#ifndef SurfaceOverlay_h
#define SurfaceOverlay_h

#include "vtkSmartPointer.h"
#include "Broadcaster.h"
#include "Listener.h"
#include <string>

extern "C"
{
#include "mri.h"
}

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
class LayerSurface;
class SurfaceOverlayProperties;

class SurfaceOverlay  : public Broadcaster, public Listener
{
  friend class SurfaceOverlayProperties;
public:
  SurfaceOverlay ( LayerSurface* surf );
  ~SurfaceOverlay ();

  void SetSurface( LayerSurface* surf );
  
  SurfaceOverlayProperties* GetProperties()
  {
    return m_properties;
  }
  
  const char* GetName();
  
  void SetName( const char* name );
  
  void InitializeData();
  
  void MapOverlay( unsigned char* colordata );
  
  float* GetData()
  { 
    return m_fData;
  }
  
  int GetDataSize()
  {
    return m_nDataSize;
  }
  
  double GetDataAtVertex( int nVertex );
  
  void GetRange( double* range )
  {
    range[0] = m_dMinValue;
    range[1] = m_dMaxValue;
  }
  
  bool LoadCorrelationData( const char* filename );
  
  bool HasCorrelationData()
  {
    return m_bCorrelationData;
  }
  
  void UpdateCorrelationAtVertex( int nVertex );
  
protected:
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
private:
  float*        m_fData;
  int           m_nDataSize;
  double        m_dMaxValue;
  double        m_dMinValue;
  
  std::string   m_strName;
  SurfaceOverlayProperties* m_properties;
  LayerSurface* m_surface;
  
  bool        m_bCorrelationData;
  bool        m_bCorrelationDataReady;
  MRI*        m_mriCorrelation;
};

#endif
