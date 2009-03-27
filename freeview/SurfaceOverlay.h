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
 *    $Author: rpwang $
 *    $Date: 2009/03/27 21:25:11 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2007-2009,
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

#ifndef SurfaceOverlay_h
#define SurfaceOverlay_h

#include "vtkSmartPointer.h"
#include "Broadcaster.h"
#include "Listener.h"
#include <string>

class vtkLookupTable;
class vtkRGBATransferFunction;
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
};

#endif
