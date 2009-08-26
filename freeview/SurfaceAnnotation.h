/**
 * @file  SurfaceAnnotation.h
 * @brief Data handler for surface annotation
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/08/26 19:59:03 $
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

#ifndef SurfaceAnnotation_h
#define SurfaceAnnotation_h

#include "vtkSmartPointer.h"
#include "Broadcaster.h"
#include "Listener.h"
#include <string>

extern "C"
{
#include "colortab.h"
}

class vtkLookupTable;
class vtkRGBATransferFunction;
class LayerSurface;
class SurfaceAnnotationProperties;

class SurfaceAnnotation  : public Broadcaster, public Listener
{
public:
  SurfaceAnnotation ( LayerSurface* surf );
  ~SurfaceAnnotation ();

  void SetSurface( LayerSurface* surf );
  
  const char* GetName();
  
  void SetName( const char* name );
  
  bool LoadAnnotation( const char* fn );
    
  int* GetIndices()
  { 
    return m_nIndices;
  }
  
  int GetIndexSize()
  {
    return m_nIndexSize;
  }
  
  COLOR_TABLE* GetColorTable()
  {
    return m_lut;
  }
  
  int GetIndexAtVertex( int nVertex );
  
  std::string GetAnnotationNameAtVertex( int nVertex );
   
protected:
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
private:
  int*          m_nIndices;
  int           m_nIndexSize;
  
  std::string   m_strName;
  COLOR_TABLE*  m_lut;
  LayerSurface* m_surface;
};

#endif
