/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/24 16:25:54 $
 *    $Revision: 1.6 $
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

#include "DecimatePanel.h"
#include "RenderPanel.h"
#include <wx/valtext.h>
#include <wx/progdlg.h>
#include <wx/filedlg.h>
#include <wx/filename.h>
#include <wx/tokenzr.h>


// Simple functions of principal curvatures (from mris_curvature_stats.c)
float f_sharpnessCurvature(float af_k1, float af_k2)
{
    return ((af_k1 - af_k2)*(af_k1 - af_k2));
}
float f_bendingEnergyCurvature(float af_k1, float af_k2)
{
    return (af_k1*af_k1 + af_k2*af_k2);
}
float f_curvednessCurvature(float af_k1, float af_k2)
{
    return (sqrt(0.5*(af_k1*af_k1 + af_k2*af_k2)));
}
float f_shapeIndexCurvature(float af_k1, float af_k2)
{
    return (af_k1 == af_k2 ? 0 : atan((af_k1+af_k2)/(af_k2 - af_k1)));
}
float f_foldingIndexCurvature(float af_k1, float af_k2)
{
    return (fabs(af_k1)*(fabs(af_k1) - fabs(af_k2)));
}

///
/// Callback function to display progress bar when decimating surface
/// 
void DecimateProgressCallback(float percent, const char *msg, void *userData)
{
    wxProgressDialog *progress = (wxProgressDialog*)userData;
    if (progress)
    {
        progress->Update((int)(percent * 100.0f), wxString::FromAscii(msg));
    }
}

///////////////////////////////////////////////////////////////////////////////////
//  
//  Constructor/Destructor
//
//


///
/// Constructor
///
DecimatePanel::DecimatePanel( wxWindow* parent, RenderPanel *renderPanel) :
    DecimatePanelBase( parent ),
    m_renderPanel(renderPanel),
    m_origSurface(NULL),
    m_decimatedSurface(NULL),
    m_lastSaveDir(NULL)
{
    m_curvatureTypes[e_None] = wxT("None");
    m_curvatureTypes[e_Gaussian] = wxT("Gaussian");
    m_curvatureTypes[e_Mean] = wxT("Mean");
    m_curvatureTypes[e_K1] = wxT("K1");
    m_curvatureTypes[e_K2] = wxT("K2");
    m_curvatureTypes[e_S] = wxT("S");
    m_curvatureTypes[e_C] = wxT("C");
    m_curvatureTypes[e_SI] = wxT("SI");
    m_curvatureTypes[e_BE] = wxT("BE");
    m_curvatureTypes[e_FI] = wxT("FI");

    ResetToDefault();

    for (int i = 0; i < (int)e_numCurvTypes; i++)
    {
        m_curvatureChoice->Append(m_curvatureTypes[i]);
    }
    m_curvatureChoice->Select(0);

    m_minimumValue = -1.0;
    m_maximumValue = 1.0;
    SetMinValue(0.0);
    SetMaxValue(0.0);
    m_minValueSlider->Disable();
    m_maxValueSlider->Disable();
    m_histogramCheckBox->Disable();
    m_colorBarCheckBox->Disable();
    m_upVectorText->Disable();
    m_cameraPositionText->Disable();
    m_focalPointText->Disable();
    m_setCameraButton->Disable();
    m_resetCameraButton->Disable();

    // Add validator
    wxTextValidator floatValidator(wxFILTER_NUMERIC);
    m_decimationLevelTextCtrl->SetValidator(floatValidator);
}

///////////////////////////////////////////////////////////////////////////////////
//  
//  Public Methods
//
//

///
/// Reset Everything to default state
///
void DecimatePanel::ResetToDefault()
{
    m_decimationLevelSlider->SetValue(100);
    SetDecimationLevel(1.0);
    m_minimumAngleSlider->SetValue(100);
    SetMinimumAngle(1.0);
    m_defaultButton->Disable();
}

///
/// Return the current decimation options set in the panel
///
DECIMATION_OPTIONS DecimatePanel::GetDecimationOptions() const
{
    DECIMATION_OPTIONS options;

    memset(&options, 0, sizeof(options));

    options.decimationLevel = m_decimationLevel;
    options.setMinimumAngle = true;
    options.minimumAngle = m_minimumAngle;
    return options;
}

///
/// Perform decimation on the original surface
///
void DecimatePanel::DoDecimate()
{
    wxProgressDialog *progress = new wxProgressDialog(wxT("Decimation Progress"),
        wxT("Loading Surface..."), 100, this);
    progress->SetSize(300, 100);

    if (m_decimatedSurface != NULL)
    {
        MRISfree(&m_decimatedSurface);
        m_decimatedSurface = NULL;
    }

    m_decimatedSurface = MRISclone(m_origSurface);
    DECIMATION_OPTIONS options = GetDecimationOptions();
    decimateSurface(&m_decimatedSurface, options, DecimateProgressCallback, progress);
    m_renderPanel->SetSurface(m_decimatedSurface);

	m_decTrianglesText->SetLabel(wxString::Format(_T("%d"), m_decimatedSurface->nfaces));
	m_decVerticesText->SetLabel(wxString::Format(_T("%d"), m_decimatedSurface->nvertices));
	m_decAvgVertexAreaText->SetLabel(wxString::Format(_T("%.3f mm^2"), m_decimatedSurface->avg_vertex_area));
	m_decAvgVertexDistText->SetLabel(wxString::Format(_T("%.3f mm"), m_decimatedSurface->avg_vertex_dist));
	    
	MRIScomputeMetricProperties(m_decimatedSurface);
	MRISsetNeighborhoodSize(m_decimatedSurface, 2);
	MRIScomputeSecondFundamentalFormDiscrete(m_decimatedSurface, 0);
	
    SetCurvatureType(m_curvatureChoice->GetSelection());
}

///
/// Set the original surface
///
void DecimatePanel::SetOrigSurface(MRI_SURFACE *origSurface)
{
    m_origSurface = origSurface;
    if (m_origSurface != NULL)
    {
		m_origTrianglesText->SetLabel(wxString::Format(wxT("%d"), m_origSurface->nfaces));
		m_origVerticesText->SetLabel(wxString::Format(wxT("%d"), m_origSurface->nvertices));
		m_origAvgVertexAreaText->SetLabel(wxString::Format(wxT("%.3f mm^2"), m_origSurface->avg_vertex_area));
		m_origAvgVertexDistText->SetLabel(wxString::Format(wxT("%.3f mm"), m_origSurface->avg_vertex_dist));	

		MRIScomputeMetricProperties(m_origSurface);
		MRISsetNeighborhoodSize(m_origSurface, 2);
	    MRIScomputeSecondFundamentalFormDiscrete(m_origSurface, 0);
	}
    else
    {
		m_origTrianglesText->SetLabel(wxT(""));
		m_origVerticesText->SetLabel(wxT(""));
		m_origAvgVertexAreaText->SetLabel(wxT(""));
		m_origAvgVertexDistText->SetLabel(wxT(""));		
    }
}

///
/// Save decimated surface to a file
///
int DecimatePanel::SaveDecimatedSurface(const char* filePath)
{
    char outFpath[STRLEN];

    FileNameAbsolute(filePath, outFpath);
    return MRISwrite(m_decimatedSurface, outFpath);
}

///
///	Update the decimation level internally along with on the GUI
///
void DecimatePanel::UpdateDecimationLevel(float val)
{
    int maxLevel = m_decimationLevelSlider->GetMax();
    int minLevel = m_decimationLevelSlider->GetMin();
    int valToSet = ((int)(val * (float) (maxLevel - minLevel))) + minLevel;
    m_decimationLevelSlider->SetValue( valToSet );
    SetDecimationLevel( val );
    m_defaultButton->Enable();
}

///
///	Update the curvature type by name.  Returns true if it is a known
///	curvature type or false if it is not found
///
bool DecimatePanel::UpdateCurvature( const wxString& curvature )
{
    for (int i = 0; i < (int)e_numCurvTypes; i++)
    {
        if (m_curvatureTypes[i].Lower() == curvature.Lower())
        {
            m_curvatureChoice->Select(i);
            SetCurvatureType(i);
            return true;
        }
    }

    return false;
}
///////////////////////////////////////////////////////////////////////////////////
//  
//  Protected Methods
//
//

///
/// Set the level of decimation
///
void DecimatePanel::SetDecimationLevel(float val)
{
    wxString valStr = wxString::Format(wxT("%.2f"), val);
    m_decimationLevelTextCtrl->Clear();
    m_decimationLevelTextCtrl->AppendText(valStr);
    m_decimationLevel = val;
}

///
/// Set the minimum angle between faces for decimation
///
void DecimatePanel::SetMinimumAngle(float val)
{
    wxString valStr = wxString::Format(wxT("%2.2f"), val);
    m_minimumAngleTextCtrl->Clear();
    m_minimumAngleTextCtrl->AppendText(valStr);
    m_minimumAngle = val;
}

///
/// Set the minimum value from the curvature histogram to color
///
void DecimatePanel::SetMinValue(float val)
{
    wxString valStr = wxString::Format(wxT("%.2f"), val);
    m_minValueTextCtrl->Clear();
    m_minValueTextCtrl->AppendText(valStr);
    m_currentMinValue = val;
}

///
/// Set the maximum value from the curvature histogram to color
///
void DecimatePanel::SetMaxValue(float val)
{
    wxString valStr = wxString::Format(wxT("%.2f"), val);
    m_maxValueTextCtrl->Clear();
    m_maxValueTextCtrl->AppendText(valStr);
    m_currentMaxValue = val;
}

///
/// Set a pointer to the variable used to store the last directory a file was saved in
///
void DecimatePanel::SetLastSaveDir(wxString *lastSaveDir)
{
    m_lastSaveDir = lastSaveDir;
}

///
///	Compute curvature stats (min, max, mean, sigma)
///
void DecimatePanel::ComputeStats(MRI_SURFACE *mris, float& minVal, float& maxVal,
			      			     float& mean, float& sigma)
{
	minVal = 100000.0f;
	maxVal = -100000.0f;
	mean = 0.0;
	sigma = 0.0;

	float total = 0.0;
	float totalSq = 0.0;
	float n = 0.0;
	

	for (int vno = 0; vno < mris->nvertices; vno++)
	{
		total += mris->vertices[vno].curv;
		totalSq += (mris->vertices[vno].curv * mris->vertices[vno].curv);
		n += 1.0;
		if (mris->vertices[vno].curv < minVal)
		{
			minVal = mris->vertices[vno].curv;
		}

		if (mris->vertices[vno].curv > maxVal)
		{
			maxVal = mris->vertices[vno].curv;
		}
	}

	if (n > 0.0)
	{
		mean = total / n;
		sigma = sqrt(totalSq / n - mean * mean);		
	}
}


///
/// Set the type of curvature to display and compute the
/// curvature using MRIS library functions as is done in mris_curvature_stats
///
void DecimatePanel::SetCurvatureType(int type)
{
	MRI_SURFACE *mris[2] = { m_origSurface, m_decimatedSurface };

	for (int i = 0; i < 2; i++)
	{		
    	if (mris[i] != NULL)
	    {
    	    short ret = -1;

			switch (type)
		    {
		        case e_Gaussian:
		            ret = MRISuseGaussianCurvature(mris[i]);
		            break;
		        case e_Mean:
		            ret = MRISuseMeanCurvature(mris[i]);
		            break;
		        case e_K1:
		            ret = MRISuseK1Curvature(mris[i]);
		            break;
		        case e_K2:
		            ret = MRISuseK2Curvature(mris[i]);
		            break;
		        case e_S:
		            ret = MRISusePrincipalCurvatureFunction(mris[i],
		                f_sharpnessCurvature);
		            break;
		        case e_C:
		            ret = MRISusePrincipalCurvatureFunction(mris[i],
		                f_curvednessCurvature);
		            break;
		        case e_BE:
		            ret = MRISusePrincipalCurvatureFunction(mris[i],
		                f_bendingEnergyCurvature);
		            break;
		        case e_SI:
		            ret = MRISusePrincipalCurvatureFunction(mris[i],
		                f_shapeIndexCurvature);
		            break;
		        case e_FI:
		            ret = MRISusePrincipalCurvatureFunction(mris[i],
		                f_foldingIndexCurvature);
		            break;
		        default:
		            m_renderPanel->ShowCurvature(false);
		            m_saveCurvatureButton->Disable();
		            m_minValueSlider->Disable();
		            m_maxValueSlider->Disable();
		            m_histogramCheckBox->Disable();
		            m_colorBarCheckBox->Disable();
		            m_renderPanel->ShowHistogram(false);
		            m_renderPanel->ShowColorBar(false);
					m_decCurvatureRangeText->SetLabel(wxT(""));
					m_decCurvatureMeanText->SetLabel(wxT(""));
					m_decCurvatureStdDevText->SetLabel(wxT(""));
					m_origCurvatureRangeText->SetLabel(wxT(""));
					m_origCurvatureMeanText->SetLabel(wxT(""));
					m_origCurvatureStdDevText->SetLabel(wxT(""));
					return;
		    }
		}		

        float minVal = 100000.0f;
        float maxVal = -100000.0f;
		float mean = 0.0;
		float sigma = 0.0;

		// Compute original surface stats and show values
		ComputeStats(m_origSurface, minVal, maxVal, mean, sigma);
		m_origCurvatureRangeText->SetLabel(wxString::Format(wxT("(%.2f, %.2f)"), minVal, maxVal));
		m_origCurvatureMeanText->SetLabel(wxString::Format(wxT("%.3f"), mean));
		m_origCurvatureStdDevText->SetLabel(wxString::Format(wxT("%.3f"), sigma));

		
		// Compute decimated surface stats and show values
		ComputeStats(m_decimatedSurface, minVal, maxVal, mean, sigma);
		m_decCurvatureRangeText->SetLabel(wxString::Format(wxT("(%.2f, %.2f)"), minVal, maxVal));
		m_decCurvatureMeanText->SetLabel(wxString::Format(wxT("%.3f"), mean));
		m_decCurvatureStdDevText->SetLabel(wxString::Format(wxT("%.3f"), sigma));
	
        m_saveCurvatureButton->Enable();
        m_renderPanel->ShowCurvature(true);
        m_minimumValue = minVal;
        m_maximumValue = maxVal;
		
		// Set the min/max based on the automatic analysis of the histogram
        float autoMinVal,
        autoMaxVal;
        m_renderPanel->GetHistogramAutoRange( autoMinVal, autoMaxVal );
        SetMinValue(autoMinVal);
        SetMaxValue(autoMaxVal);
        m_minValueSlider->Enable();
        m_maxValueSlider->Enable();


        // Adjust the slider to the automatic value
        if ((m_maximumValue - m_minimumValue) != 0.0)
        {
            float sliderMinVal = (autoMinVal - m_minimumValue) / (m_maximumValue - m_minimumValue);
            float sliderMaxVal = (autoMaxVal - m_minimumValue) / (m_maximumValue - m_minimumValue);

            float sliderMinRange = m_minValueSlider->GetMax() - m_minValueSlider->GetMin();
            float sliderMaxRange = m_maxValueSlider->GetMax() - m_maxValueSlider->GetMin();

            m_minValueSlider->SetValue( (int)(sliderMinVal * sliderMinRange) + m_minValueSlider->GetMin() );
            m_maxValueSlider->SetValue( (int)(sliderMaxVal * sliderMaxRange) + m_maxValueSlider->GetMin() );

        }

        m_renderPanel->SetMinMaxCurvatureRange(autoMinVal, autoMaxVal);
        m_histogramCheckBox->Enable();
        m_colorBarCheckBox->Enable();
        m_renderPanel->ShowHistogram(m_histogramCheckBox->GetValue() == 0 ? false : true);
        m_renderPanel->ShowColorBar(m_colorBarCheckBox->GetValue() == 0 ? false : true);

    }
}


///
/// Update the camera parameters from the current camera location
///
void DecimatePanel::UpdateCameraInfo()
{
    vtkCamera *camera = m_renderPanel->GetCamera();
    if (camera)
    {
        double vx, vy, vz;
        camera->GetViewUp( vx, vy, vz );
        m_upVectorText->SetValue(wxString::Format(wxT("%.2f %.2f %.2f"), vx, vy, vz));
        m_upVectorText->Enable();

        double px, py, pz;
        camera->GetPosition( px, py, pz );
        m_cameraPositionText->SetValue(wxString::Format(wxT("%.2f %.2f %.2f"), px, py, pz));
        m_cameraPositionText->Enable();

        double fx, fy, fz;
        camera->GetFocalPoint( fx, fy, fz );
        m_focalPointText->SetValue(wxString::Format(wxT("%.2f %.2f %.2f"), fx, fy, fz));
        m_focalPointText->Enable();

        m_setCameraButton->Enable();
        m_resetCameraButton->Enable();
    }
}

///
/// Given an input text string such as (1.0 2.0 3.0) convert
/// the strings to double values
/// \return True if succesfully converted, false otherwise
///
bool DecimatePanel::ConvertStringVecToDoubles( const wxString& str, double *vec, int components)
{
    wxStringTokenizer tokenizer(str, wxT(" ,"));
    int count = 0;

    while ( tokenizer.HasMoreTokens() && count < components )
    {
        wxString token = tokenizer.GetNextToken();
        double val;
        if (token.ToDouble(&val))
        {
            vec[count++] = val;
        }
        else
        {
            return false;
        }
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////////
//  
//  wxWidgets Handlers
//
//

void DecimatePanel::OnDecimationLevelChanged( wxScrollEvent& event )
{
    int maxLevel = m_decimationLevelSlider->GetMax();
    int minLevel = m_decimationLevelSlider->GetMin();
    float valToSetf = (float)(event.GetPosition() - minLevel) / (float)(maxLevel - minLevel);

    SetDecimationLevel(valToSetf);
    m_defaultButton->Enable();
}

void DecimatePanel::OnDecimationText( wxCommandEvent& event )
{
    double val;
    if (event.GetString().ToDouble(&val))
    {
        if (val < 0.0)
        {
            SetDecimationLevel(0.0);
            return;
        }
        else if( val > 1.0)
        {
            SetDecimationLevel(1.0);
            return;
        }

        int maxLevel = m_decimationLevelSlider->GetMax();
        int minLevel = m_decimationLevelSlider->GetMin();
        int valToSet = ((int)(val * (float) (maxLevel - minLevel))) + minLevel;

        m_decimationLevelSlider->SetValue(valToSet);    
        m_decimationLevel = val;   
    }
}


void DecimatePanel::OnMinimumAngleChanged( wxScrollEvent& event )
{
    int maxLevel = m_minimumAngleSlider->GetMax();
    int minLevel = m_minimumAngleSlider->GetMin();
    float valToSetf = (float)(event.GetPosition() - minLevel) / (float)(maxLevel - minLevel);

    SetMinimumAngle(valToSetf * 90.0f);
    m_defaultButton->Enable();
}

void DecimatePanel::OnMinimumAngleText( wxCommandEvent& event )
{
    double val;
    if (event.GetString().ToDouble(&val))
    {
        if (val < 0.0)
        {
            SetMinimumAngle(0.0);
            return;
        }
        else if(val > 90.0)
        {
            SetMinimumAngle(90.0);
            return;
        }

        int maxLevel = m_minimumAngleSlider->GetMax();
        int minLevel = m_minimumAngleSlider->GetMin();
        int valToSet = ((int)(val / 90.0 * (float) (maxLevel - minLevel))) + minLevel;

        m_minimumAngleSlider->SetValue(valToSet);  
        m_minimumAngle = val;
             
    }
}

void DecimatePanel::OnMinimumValueChanged( wxScrollEvent& event )
{
    int maxLevel = m_minValueSlider->GetMax();
    int minLevel = m_minValueSlider->GetMin();
    float valToSetf = (float)(event.GetPosition() - minLevel) / (float)(maxLevel - minLevel);

    SetMinValue((valToSetf * (m_maximumValue - m_minimumValue)) + m_minimumValue);

    if (m_currentMinValue > m_currentMaxValue)
    {
        m_maxValueSlider->SetValue(event.GetPosition());
        SetMaxValue(m_currentMinValue);
    }
    m_renderPanel->SetMinMaxCurvatureRange(m_currentMinValue, m_currentMaxValue);
}


void DecimatePanel::OnMinimumValueText( wxCommandEvent& event )
{
    double val;
    if (event.GetString().ToDouble(&val))
    {
        if (val < m_minimumValue)
        {
            SetMinValue(m_minimumValue);
            return;
        }
        else if(val > m_maximumValue)
        {
            SetMinValue(m_maximumValue);
            return;
        }

        int maxLevel = m_minValueSlider->GetMax();
        int minLevel = m_minValueSlider->GetMin();
        int valToSet = ((int)((val - m_minimumValue) / (m_maximumValue - m_minimumValue) * 
                            (float) (maxLevel - minLevel))) + minLevel;

        m_minValueSlider->SetValue(valToSet);  
        m_currentMinValue = val;
        m_renderPanel->SetMinMaxCurvatureRange(m_currentMinValue, m_currentMaxValue);     
    }
}

void DecimatePanel::OnMaximumValueChanged( wxScrollEvent& event )
{
    int maxLevel = m_maxValueSlider->GetMax();
    int minLevel = m_maxValueSlider->GetMin();
    float valToSetf = (float)(event.GetPosition() - minLevel) / (float)(maxLevel - minLevel);

    SetMaxValue((valToSetf * (m_maximumValue - m_minimumValue)) + m_minimumValue);
    if (m_currentMaxValue < m_currentMinValue)
    {
        m_minValueSlider->SetValue(event.GetPosition());
        SetMinValue(m_currentMaxValue);
    }

    m_renderPanel->SetMinMaxCurvatureRange(m_currentMinValue, m_currentMaxValue);
}

void DecimatePanel::OnMaximumValueText( wxCommandEvent& event )
{
    double val;
    if (event.GetString().ToDouble(&val))
    {
        if (val < m_minimumValue)
        {
            SetMaxValue(m_minimumValue);
            return;
        }
        else if(val > m_maximumValue)
        {
            SetMaxValue(m_maximumValue);
            return;
        }

        int maxLevel = m_maxValueSlider->GetMax();
        int minLevel = m_maxValueSlider->GetMin();
        int valToSet = ((int)((val - m_minimumValue) / (m_maximumValue - m_minimumValue) * 
                                (float) (maxLevel - minLevel))) + minLevel;

        m_maxValueSlider->SetValue(valToSet);     
        m_currentMaxValue = val;
        m_renderPanel->SetMinMaxCurvatureRange(m_currentMinValue, m_currentMaxValue);       
    }
}

void DecimatePanel::OnApplyButtonClick( wxCommandEvent& event )
{
    DoDecimate();
}

void DecimatePanel::OnDefaultButtonClick( wxCommandEvent& event )
{
    ResetToDefault();
}

void DecimatePanel::OnRenderModeChoice( wxCommandEvent& event )
{
    m_renderPanel->SetRenderMode((RenderPanel::RenderMode)event.GetSelection());
}

void DecimatePanel::OnSaveCurvatureClick( wxCommandEvent& event )
{
    wxFileDialog *saveDialog = new wxFileDialog(
        this, _("Save Curvature File As..."), *m_lastSaveDir, wxEmptyString,
        _("Curvature files(*.crv)|*.crv"),
        wxFD_SAVE | wxFD_OVERWRITE_PROMPT, wxDefaultPosition);

    if (saveDialog->ShowModal() == wxID_OK)
    {
        wxString currentFilePath = saveDialog->GetPath();
        wxFileName fileName(currentFilePath);
        *m_lastSaveDir = fileName.GetPath();

        if (MRISwriteCurvature(m_decimatedSurface, currentFilePath.mb_str(wxConvUTF8)) != 0)
        {
            wxMessageBox(wxString::Format(wxT("ERROR: Saving file '%s'"), currentFilePath.c_str()),
                wxT("ERROR"));
        }
    }
}

void DecimatePanel::OnHistogramCheckBox( wxCommandEvent& event )
{
    m_renderPanel->ShowHistogram(m_histogramCheckBox->GetValue() == 0 ? false : true);
}


void DecimatePanel::OnCurvatureChoice( wxCommandEvent& event )
{
    SetCurvatureType(event.GetSelection());
}

void DecimatePanel::OnColorBarCheckBox( wxCommandEvent& event )
{
    m_renderPanel->ShowColorBar(m_colorBarCheckBox->GetValue() == 0 ? false : true);
}

void DecimatePanel::OnUpVectorText( wxCommandEvent& event )
{
}

void DecimatePanel::OnCameraPositionText( wxCommandEvent& event )
{
}

void DecimatePanel::OnFocalPointText( wxCommandEvent& event )
{
}


void DecimatePanel::OnSetCameraClick( wxCommandEvent& event )
{
    double pos[3], focalPoint[3], upVector[3];
    
    if( !ConvertStringVecToDoubles( m_cameraPositionText->GetValue(), pos, 3) )
    {
        wxMessageBox(wxString::Format(wxT("ERROR: Could not convert to vector '%s'"), m_cameraPositionText->GetValue().c_str()),
                    wxT("ERROR"));
        return;
    }

    if( !ConvertStringVecToDoubles( m_focalPointText->GetValue(), focalPoint, 3) )
    {
        wxMessageBox(wxString::Format(wxT("ERROR: Could not convert to vector '%s'"), m_focalPointText->GetValue().c_str()),
                    wxT("ERROR"));
        return;
    }

    if( !ConvertStringVecToDoubles( m_upVectorText->GetValue(), upVector, 3) )
    {
        wxMessageBox(wxString::Format(wxT("ERROR: Could not convert to vector '%s'"), m_upVectorText->GetValue().c_str()),
                    wxT("ERROR"));
        return;
    }

    vtkCamera *camera = m_renderPanel->GetCamera();
    camera->SetFocalPoint(focalPoint);
    camera->SetPosition(pos);
    camera->SetViewUp(upVector);    
    UpdateCameraInfo();
    m_renderPanel->ResetCameraClippingRange();
}

void DecimatePanel::OnResetCameraClick( wxCommandEvent& event )
{
    m_renderPanel->ResetCameraToDefaultView();
}

void DecimatePanel::OnSaveScreenshotClick( wxCommandEvent& event )
{
    wxFileDialog *saveDialog = new wxFileDialog(
        this, _("Save Screenshot As..."), *m_lastSaveDir, wxEmptyString,
        _("Image files (*.png;*.gif;*.jpg;*.bmp:*.ps:*.tiff)|*.png;*.gif;*.jpg;*.bmp:*.ps:*.tiff"),
        wxFD_SAVE | wxFD_OVERWRITE_PROMPT, wxDefaultPosition);

    if (saveDialog->ShowModal() == wxID_OK)
    {
        wxString currentFilePath = saveDialog->GetPath();
        wxFileName fileName(currentFilePath);
        *m_lastSaveDir = fileName.GetPath();

        if ( !m_renderPanel->SaveScreenshot(currentFilePath) )
        {
            wxMessageBox(wxString::Format(wxT("ERROR: Failed to write image '%s'"), currentFilePath.c_str()),
                    wxT("ERROR"));
            return;
        }
    }
}

