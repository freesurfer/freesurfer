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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:03 $
 *    $Revision: 1.9 $
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

#ifndef SurfaceAnnotation_h
#define SurfaceAnnotation_h

#include "Broadcaster.h"
#include "Listener.h"
#include <string>

extern "C"
{
#include "colortab.h"
}

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
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
  
  int GetNumberOfAnnotations()
  {
    return m_nAnnotations;
  }
  
  COLOR_TABLE* GetColorTable()
  {
    return m_lut;
  }
  
  int GetIndexAtVertex( int nVertex );
  
  void GetAnnotationPoint( int nIndex, double* pt_out );
  
  std::string GetAnnotationNameAtIndex( int nIndex );
  
  std::string GetAnnotationNameAtVertex( int nVertex );
  
  void GetAnnotationColorAtIndex( int nIndex, int* rgb );
   
protected:
  void Reset();
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
private:
  int*          m_nIndices;
  int           m_nIndexSize;
  int*          m_nCenterVertices;  // center vertex of each annotation
  int           m_nAnnotations;     // number of valid annotations
  
  std::string   m_strName;
  COLOR_TABLE*  m_lut;
  LayerSurface* m_surface;
};

#endif
