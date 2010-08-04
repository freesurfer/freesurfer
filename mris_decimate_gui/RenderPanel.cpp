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
#include "RenderPanel.h"

#include "wxVTKRenderWindowInteractor.h"
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkVRMLExporter.h>
#include <vtkRenderLargeImage.h>
#include <vtkBMPWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkImageImport.h>
#include <vtkImageData.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageAccumulate.h>
#include <vtkXYPlotActor.h>

BEGIN_EVENT_TABLE( RenderPanel, wxPanel )
    EVT_SIZE        ( RenderPanel::OnSize )
END_EVENT_TABLE()

///////////////////////////////////////////////////////////////////////////////////
//  
//  Constructor/Destructor
//
//

///
/// Constructor
///
RenderPanel::RenderPanel( wxWindow* parent ) :
        wxPanel( parent ),
        m_histogramMaxY(0.0),
        m_surface( NULL )
{

    wxBoxSizer *sizer = new wxBoxSizer(wxHORIZONTAL);

    // Create variables to hold mesh data
    m_mesh = vtkPolyData::New();

    m_polyMapper = vtkPolyDataMapper::New();
    m_polyMapper->SetInput(m_mesh);
    m_polyMapper->SetScalarRange(-1,1);
    m_polyActor = vtkActor::New();
    m_polyActor->SetMapper(m_polyMapper);

    // The usual rendering stuff.
    m_camera = vtkCamera::New();
    m_camera->SetPosition(-1,0,0);
    m_camera->SetFocalPoint(0,0,0);
    m_camera->Roll(90.0);


    m_renderer = vtkRenderer::New();
    m_renWin = new wxVTKRenderWindowInteractor(this,-1, this->GetPosition(), this->GetSize());

    m_renWin->GetRenderWindow()->AddRenderer(m_renderer);
    sizer->Add(m_renWin, 1, wxEXPAND);

    m_renderer->AddActor(m_polyActor);
    m_renderer->SetActiveCamera(m_camera);
    m_renderer->ResetCamera();
    m_renderer->SetBackground(0,0,0);

    this->SetSizer(sizer);


    m_histogramPlot = vtkXYPlotActor::New();
    m_histogramPlot->ExchangeAxesOff();
    m_histogramPlot->SetLabelFormat( "%g" );
    m_histogramPlot->SetXTitle( "" );
    m_histogramPlot->SetYTitle( "" );
    m_histogramPlot->SetXValuesToValue();
    m_renderer->AddActor(m_histogramPlot);
    m_histogramPlot->VisibilityOff();

    m_histogramMinLine = vtkPolyData::New();
    m_histogramMaxLine = vtkPolyData::New();

}


///
/// Destructor
///
RenderPanel::~RenderPanel()
{
    m_renWin->Delete();
    m_mesh->Delete();
    m_polyMapper->Delete();
    m_polyActor->Delete();
    m_renderer->Delete();
    m_histogramPlot->Delete();
    m_histogramMinLine->Delete();
    m_histogramMaxLine->Delete();
}

///////////////////////////////////////////////////////////////////////////////////
//  
//  Public Methods
//
//


///
/// Set the MRIS and create the VTK structure for rendering it
///
void RenderPanel::SetSurface( MRI_SURFACE *mris )
{
    m_surface = mris;

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

    for ( int v = 0; v < mris->nvertices; v++)
    {
        points->InsertPoint(v, mris->vertices[v].x,
                            mris->vertices[v].y,
                            mris->vertices[v].z);
    }

    for (int f = 0; f < mris->nfaces; f++)
    {
        vtkIdType pts[3];
        pts[0] =  mris->faces[f].v[0];
        pts[1] =  mris->faces[f].v[1];
        pts[2] =  mris->faces[f].v[2];

        polys->InsertNextCell(3, pts);
    }

    m_mesh->SetPoints(points);
    m_mesh->SetPolys(polys);

    Render();
}

///////////////////////////////////////////////////////////////////////////////////
//  
//  wxWidgets Handlers
//
//

void RenderPanel::OnSize(wxSizeEvent& event)
{
    event.Skip();
    m_renWin->OnSize(event);
}

///////////////////////////////////////////////////////////////////////////////////
//  
//  Public Methods
//
//

///
/// Re-render the scene
///
void RenderPanel::Render()
{
    m_renWin->Render();
}

///
/// Turn wireframe rendering on
///
void RenderPanel::SetWireframe(bool enable)
{
    if (m_polyActor != NULL)
    {
        if (enable)
        {
            m_polyActor->GetProperty()->SetRepresentationToWireframe();
            m_polyActor->GetProperty()->ShadingOff();
        }
        else
        {
            m_polyActor->GetProperty()->SetRepresentationToSurface();
            m_polyActor->GetProperty()->ShadingOn();
        }
        Render();
    }
}

///
/// Reset the camera
///
void RenderPanel::ResetCamera()
{
    m_renderer->ResetCamera();
    Render();
}


///
/// Helper routine to compute file extension of images
///
bool HasExtension( const wxString& filename, const wxString& ext )
{
    return ( filename.Lower().Right( ext.Len() ) == ext.Lower() );
}

///
/// Save a screenshot, code adapted from freeview
///
bool RenderPanel::SaveScreenshot(const wxString& filename)
{
    wxString fn = filename;
    vtkImageWriter* writer = 0;
    if ( HasExtension( fn, _("wrl") ) )
    {
        vtkVRMLExporter* exporter = vtkVRMLExporter::New();
        exporter->SetFileName( fn );
        exporter->SetRenderWindow( m_renWin->GetRenderWindow() );
        exporter->Write();
        exporter->Delete();
    }
    else if ( HasExtension( fn, _("jpg") ) ||
              HasExtension( fn, _("jpeg") ) )
        writer = vtkJPEGWriter::New();
    else if ( HasExtension( fn, _("bmp") ) )
        writer = vtkBMPWriter::New();
    else if ( HasExtension( fn, _("ps") ) )
        writer = vtkPostScriptWriter::New();
    else if ( HasExtension( fn, _("tif") ) ||
              HasExtension( fn, _("tiff") ) )
        writer = vtkTIFFWriter::New();
    else
    {
        writer = vtkPNGWriter::New();
        if ( !HasExtension( fn, _("png")) )
            fn += _(".png");
    }

    bool ret = true;
    if (writer)
    {
        vtkRenderLargeImage* image = vtkRenderLargeImage::New();
        image->SetInput( m_renderer );
        writer->SetInput( image->GetOutput() );
        writer->SetFileName( fn.mb_str() );
        writer->Write();
        if ( writer->GetErrorCode() != 0 )
            ret = false;
        image->Delete();
        writer->Delete();

    }
    return ret;
}

///
/// If enabled, generate a histogram for the curvature values and set the scalar
/// range to map it to colors.  The min/max display range of the histogram can
/// be set via the GUI.
///
void RenderPanel::ShowCurvature(bool enable)
{
    if (enable && m_surface != NULL)
    {
        vtkSmartPointer<vtkFloatArray> curvs = vtkSmartPointer<vtkFloatArray>::New();
        curvs->Allocate( m_surface->nvertices );
        curvs->SetNumberOfComponents( 1 );
        for ( int vno = 0; vno < m_surface->nvertices; vno++ )
        {
            curvs->InsertNextValue( m_surface->vertices[vno].curv );
        }
        m_mesh->GetPointData()->SetScalars( curvs );

        // Compute a histogram for the curvs values
        vtkSmartPointer<vtkImageData> curvData = vtkSmartPointer<vtkImageData>::New();
        vtkSmartPointer<vtkImageImport> curvImport = vtkSmartPointer<vtkImageImport>::New();

        curvImport->SetDataOrigin(0, 0, 0);
        curvImport->SetWholeExtent(0, m_surface->nvertices - 1, 0, 0, 0, 0);
        curvImport->SetDataExtentToWholeExtent();
        curvImport->SetDataScalarTypeToFloat();
        curvImport->SetNumberOfScalarComponents(1);
        curvImport->SetImportVoidPointer(curvs->GetPointer(0));

        vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
        extract->SetInputConnection( curvImport->GetOutputPort() );
        extract->SetComponents( 0 );
        extract->Update();

        double range[2];
        extract->GetOutput()->GetScalarRange(range);

        vtkSmartPointer<vtkImageAccumulate> histogram = vtkSmartPointer<vtkImageAccumulate>::New();
        histogram->SetInputConnection( extract->GetOutputPort() );
        histogram->SetComponentExtent( 0, 255,0,0,0,0 );
        histogram->SetComponentOrigin( range[0],0,0 );
        histogram->SetComponentSpacing( (range[1] - range[0]) / 255.0f, 0, 0 );
        histogram->Update();

        m_histogramMaxY = histogram->GetOutput()->GetScalarRange()[1];

        SetHistogramMinMaxLines(range[0], range[1]);
        m_histogramPlot->RemoveAllInputs();
        m_histogramPlot->AddInput(m_histogramMinLine);
        m_histogramPlot->AddInput(m_histogramMaxLine);
        m_histogramPlot->AddInput(histogram->GetOutput());
        m_histogramPlot->SetPlotColor(0, 1, 0, 0);
        m_histogramPlot->SetPlotColor(1, 1, 0, 0);
        m_histogramPlot->SetPlotColor(2, 0, 1, 0);
        m_histogramPlot->SetXRange(range[0], range[1]);
        m_histogramPlot->SetYRange(0, m_histogramMaxY);

        m_histogramPlot->SetPosition( 0.5, 0.70 );
        m_histogramPlot->SetPosition2( 0.5, 0.25 );

    }
    else
    {
        m_mesh->GetPointData()->SetScalars(NULL);
    }

    Render();
}

///
/// Set min/max of histogram curvature range to display
///
void RenderPanel::SetMinMaxCurvatureRange(float minVal, float maxVal)
{
    m_polyMapper->SetScalarRange(minVal, maxVal);
    SetHistogramMinMaxLines(minVal, maxVal);
    Render();
}

///
/// Compute lines to draw on the 2D histogram plot for the min/max histogram
///
void RenderPanel::SetHistogramMinMaxLines(float minVal, float maxVal)
{
    // Create vertical lines for min/max histogram
    vtkSmartPointer<vtkPoints> vertPointsMin = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> vertScalarsMin = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkPoints> vertPointsMax = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> vertScalarsMax = vtkSmartPointer<vtkFloatArray>::New();


    float x[3] = {0, 0, 0};
    float y0 = 0.0;
    float y1 = m_histogramMaxY-1;
    x[0] = minVal;
    vertPointsMin->InsertPoint(0, x);
    vertScalarsMin->InsertTuple1(0, y0);
    vertPointsMin->InsertPoint(1, x);
    vertScalarsMin->InsertTuple1(1, y1);

    m_histogramMinLine->SetPoints(vertPointsMin);
    m_histogramMinLine->GetPointData()->SetScalars(vertScalarsMin);

    x[0] = maxVal;
    vertPointsMax->InsertPoint(0, x);
    vertScalarsMax->InsertTuple1(0, y0);
    vertPointsMax->InsertPoint(1, x);
    vertScalarsMax->InsertTuple1(1, y1);

    m_histogramMaxLine->SetPoints(vertPointsMax);
    m_histogramMaxLine->GetPointData()->SetScalars(vertScalarsMax);
}

void RenderPanel::ShowHistogram(bool enabled)
{
    if (enabled)
    {
        m_histogramPlot->VisibilityOn();
    }
    else
    {
        m_histogramPlot->VisibilityOff();
    }

    Render();
}

