/**
 * @brief Contour2D data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#ifndef Contour2D_h
#define Contour2D_h

#include <QObject>
#include <QColor>
#include <vtkSmartPointer.h>

class RenderView2D;
class vtkRenderer;
class vtkImageData;
class vtkActor;
class vtkImageActor;
class vtkImageReslice;
class vtkImageThreshold;
class vtkSimpleLabelEdgeFilter;
class vtkImageMapToColors;
class vtkImageMask;
class vtkImageGaussianSmooth;
class vtkImageLogic;

class Contour2D : public QObject
{
  Q_OBJECT
public:
  Contour2D( RenderView2D* view );
  virtual ~Contour2D();

  vtkImageData* GetInputImage();

  void SetInput( vtkImageData* imagedata, double dContourValue, double dSliceLocation, int active_frame = 0 );

  double GetContourValue()
  {
    return m_dContourValue;
  }

  int GetPlane()
  {
    return m_nPlane;
  }

  void SetPlane( int nPlane )
  {
    m_nPlane = nPlane;
  }

  void UpdateSliceLocation( double slice_location, bool bForced = false );

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

  double GetSmoothSD();

  void SetColor( double r, double g, double b );

  double* GetColor()
  {
    return m_dContourColor;
  }

Q_SIGNALS:
  void ValueChanged();
  void ColorChanged();

public slots:
  void SetContourValue( double dvalue );
  void SetSmooth( bool bSmooth );
  void SetSmoothSD( double sd );
  void SetColor(const QColor& color)
  {
    SetColor(color.redF(), color.greenF(), color.blueF());
  }

protected:
  void DrawPatchLineOnMask( vtkImageData* image, double* ras1, double* ras2, int nDrawValue );

  RenderView2D*   m_view;
  int             m_nPlane;
  double          m_dSliceLocation;
  double          m_dContourValue;
  bool            m_bSmooth;
  vtkSmartPointer<vtkImageData>   m_imageInput;
  double          m_dContourColor[3];

  vtkSmartPointer<vtkImageActor>      m_actorContour;
  vtkSmartPointer<vtkImageThreshold>  m_filterThreshold;
  vtkSmartPointer<vtkImageReslice>   m_filterResample;
  vtkSmartPointer<vtkSimpleLabelEdgeFilter> m_filterEdge;
  vtkSmartPointer<vtkImageMapToColors>      m_colormap;
  vtkSmartPointer<vtkImageGaussianSmooth>   m_filterSmooth;
  vtkSmartPointer<vtkImageMask>             m_filterMask;
  vtkSmartPointer<vtkImageLogic>            m_filterLogic;
  vtkSmartPointer<vtkImageData>             m_imageMaskAdd;     // to add pixels
  vtkSmartPointer<vtkImageData>             m_imageMaskRemove;  // to remove pixels
};

#endif


