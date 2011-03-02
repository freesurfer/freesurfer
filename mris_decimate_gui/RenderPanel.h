/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:31 $
 *    $Revision: 1.4 $
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
#ifndef __RenderPanel__
#define __RenderPanel__

#include <wx/panel.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCommand.h>
#include <vtkWidgetEvent.h>
#include <vtkCallbackCommand.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSliderWidget.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>
#include <vtkTextActor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkDecimatePro.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXYPlotActor.h>
#include <vtkImageAccumulate.h>
#include <vtkScalarBarActor.h>

#include "wxVTKRenderWindowInteractor.h"

extern "C"
{
#include "mrisurf.h"
}

class RenderPanel : public wxPanel
{
public:

    /// Render mode
    typedef enum
    {
        e_Surface = 0,
        e_SurfaceAndWireframe,
        e_Wireframe,

        e_NumRenderModes
    } RenderMode;

    RenderPanel( wxWindow* parent );
    virtual ~RenderPanel();

    void SetSurface( MRI_SURFACE *mris);

    void OnSize(wxSizeEvent& event);

    void Render();

    void SetRenderMode(RenderMode renderMode);

    void ResetCamera();

    bool SaveScreenshot(const wxString& fileName);

    void ShowCurvature(bool enable);

    void SetMinMaxCurvatureRange(float minVal, float maxVal);

    void ShowHistogram(bool enable);

    void ShowColorBar(bool enable);

    ///
    /// Get the camera
    ///
    vtkCamera* GetCamera() const { return m_camera; } 

    ///
    /// Add an observer that gets notified when there are changes
    /// to the camera.
    ///        
    void AddCameraObserver(vtkCommand *observer);

    ///
    /// Automatically determine the min/max values to use for the
    ///	histogram based on the mean/std. dev. of the positive and
    /// negative lobes of the curvature histograms.
    ///
    void GetHistogramAutoRange(float& minVal, float &maxVal) const;

    ///
    /// Rotate the camera about the focal point by an angle.  The
    /// default position has the left side of the left hemisphere facing
    /// the camera.
    ///
    void AzimuthCamera(double angle);

    ///
    /// Record the current camera coordinates for resetting them
    ///
    void RecordCameraCoordinates();

    ///
    /// Reset the camera to the default view
    ///
    void ResetCameraToDefaultView();

    ///
    /// Reset the camera clipping range (required after adjusting
    /// the view transform)
    ///
    void ResetCameraClippingRange();


protected:

    ///
    ///	Given a buffer of scalar values, compute a histogram
    /// using VTK
    ///
    vtkSmartPointer<vtkImageAccumulate> ComputeHistogram(vtkSmartPointer<vtkFloatArray> scalars,
            int numScalars,
            double range[2]);

    vtkPolyData *m_mesh;
    vtkPolyData *m_wireframe;
    vtkActor *m_polyActor;
    vtkActor *m_wireframeActor;
    vtkPolyDataMapper *m_polyMapper;
    vtkPolyDataMapper *m_wireframeMapper;
    vtkRenderer *m_renderer;
    vtkCamera *m_camera;
    wxVTKRenderWindowInteractor *m_renWin;
    vtkXYPlotActor *m_histogramPlot;
    vtkPolyData *m_histogramMinLine;
    vtkPolyData *m_histogramMaxLine;
    vtkScalarBarActor *m_scalarBarActor;
    float m_histogramMaxY;

    float m_histogramPosMean;
    float m_histogramPosStdDev;
    float m_histogramNegMean;
    float m_histogramNegStdDev;

    double m_origCamFocalPoint[3];
    double m_origCamPosition[3];
    double m_origCamUpVector[3];

    MRI_SURFACE *m_surface;

    DECLARE_EVENT_TABLE()
};

#endif // __RenderPanel__
