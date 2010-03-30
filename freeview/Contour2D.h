/**
 * @file  Contour2D.h
 * @brief Contour2D data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/03/30 18:31:03 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2008-2009,
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

#ifndef Contour2D_h
#define Contour2D_h

#include <wx/wx.h>
#include "Broadcaster.h"
#include "Listener.h"
#include <vtkSmartPointer.h>

class RenderView2D;
class vtkRenderer;
class vtkImageData;
class vtkActor;
class vtkImageActor;
class vtkImageResample;
class vtkImageThreshold;
class vtkSimpleLabelEdgeFilter;
class vtkImageMapToColors;
class vtkImageMask;
class vtkImageGaussianSmooth;
class vtkImageLogic;

class Contour2D : public Broadcaster, public Listener
{
public:
  Contour2D( RenderView2D* view );
  virtual ~Contour2D();

  vtkImageData* GetInputImage()
  {
    return m_imageInput;
  }
  
  void SetInput( vtkImageData* imagedata, double dContourValue, double dSliceLocation );
  
  void SetContourValue( double dvalue );
  
  double GetContourValue()
  {
    return m_dContourValue;
  }
  
  int GetPlane()
  {
    return m_nPlane;
  }
  
  void UpdateSliceLocation( double slice_location );
  
  void Reset();

  vtkImageActor* GetActor();
  
  vtkImageData* GetThresholdedImage();
  
  void SetVisible ( bool visible );
  
  bool IsVisible();
  
  void AddLine( double* ras1, double* ras2 );
  
  void RemoveLine( double* ras1, double* ras2 );
  
  bool GetSmooth()
  {
    return m_bSmooth;
  }
  
  void SetSmooth( bool bSmooth );
  
  double GetSmoothSD();
  
  void SetSmoothSD( double sd );
  
  void SetContourColor( double r, double g, double b );
  
  double* GetContourColor()
  {
    return m_dContourColor;
  }
  
protected:
  void DrawPatchLineOnMask( vtkImageData* image, double* ras1, double* ras2, int nDrawValue );
  
  RenderView2D*   m_view;
  int             m_nPlane;
  double          m_dSliceLocation;
  double          m_dContourValue;
  bool            m_bSmooth;
  vtkImageData*   m_imageInput;
  double          m_dContourColor[3];
  
  vtkSmartPointer<vtkImageActor>      m_actorContour;
  vtkSmartPointer<vtkImageThreshold>  m_filterThreshold;
  vtkSmartPointer<vtkImageResample>   m_filterResample;
  vtkSmartPointer<vtkSimpleLabelEdgeFilter> m_filterEdge;
  vtkSmartPointer<vtkImageMapToColors>      m_colormap;
  vtkSmartPointer<vtkImageGaussianSmooth>   m_filterSmooth;
  vtkSmartPointer<vtkImageMask>             m_filterMask;
  vtkSmartPointer<vtkImageLogic>            m_filterLogic;
  vtkSmartPointer<vtkImageData>             m_imageMaskAdd;     // to add pixels
  vtkSmartPointer<vtkImageData>             m_imageMaskRemove;  // to remove pixels
};

#endif


