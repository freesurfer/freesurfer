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
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:53 $
 *    $Revision: 1.10 $
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

#include <QObject>

extern "C"
{
#include "colortab.h"
}

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
class LayerSurface;;

class SurfaceAnnotation  : public QObject
{
public:
  SurfaceAnnotation ( LayerSurface* surf );
  ~SurfaceAnnotation ();

  void SetSurface( LayerSurface* surf );
  
  QString GetName();
  
  void SetName( const QString& name );
  
  bool LoadAnnotation( const QString& fn );
    
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
  
  QString GetAnnotationNameAtIndex( int nIndex );
  
  QString GetAnnotationNameAtVertex( int nVertex );
  
  void GetAnnotationColorAtIndex( int nIndex, int* rgb );
   
protected:
  void Reset();
  
private:
  int*          m_nIndices;
  int           m_nIndexSize;
  int*          m_nCenterVertices;  // center vertex of each annotation
  int           m_nAnnotations;     // number of valid annotations
  
  QString       m_strName;
  COLOR_TABLE*  m_lut;
  LayerSurface* m_surface;
};

#endif
