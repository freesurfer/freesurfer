/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: ginsburg $
 *    $Date: 2010/08/04 20:38:52 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2010,
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

#include "wxVTKRenderWindowInteractor.h"

extern "C"
{
#include "mrisurf.h"
}

class RenderPanel : public wxPanel
{
public:
    RenderPanel( wxWindow* parent );
    virtual ~RenderPanel();

    void SetSurface( MRI_SURFACE *mris);

    void OnSize(wxSizeEvent& event);

    void Render();

    void SetWireframe(bool enable);

    void ResetCamera();

    bool SaveScreenshot(const wxString& fileName);

    void ShowCurvature(bool enable);

    void SetMinMaxCurvatureRange(float minVal, float maxVal);

    void SetHistogramMinMaxLines(float minVal, float maxVal);

    void ShowHistogram(bool enable);


protected:
    vtkPolyData *m_mesh;
    vtkActor *m_polyActor;
    vtkPolyDataMapper *m_polyMapper;
    vtkRenderer *m_renderer;
    vtkCamera *m_camera;
    wxVTKRenderWindowInteractor *m_renWin;
    vtkXYPlotActor *m_histogramPlot;
    vtkPolyData *m_histogramMinLine;
    vtkPolyData *m_histogramMaxLine;
    float m_histogramMaxY;


    MRI_SURFACE *m_surface;

    DECLARE_EVENT_TABLE()
};

#endif // __RenderPanel__
