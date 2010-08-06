/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: ginsburg $
 *    $Date: 2010/08/06 14:23:57 $
 *    $Revision: 1.3 $
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

#include "DecimatePanel.h"
#include "RenderPanel.h"
#include <wx/valtext.h>
#include <wx/progdlg.h>
#include <wx/filedlg.h>
#include <wx/filename.h>


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
        progress->Update((int)(percent * 100.0f), msg);
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
    m_curvatureTypes[e_None] = "None";
    m_curvatureTypes[e_Gaussian] = "Gaussian";
    m_curvatureTypes[e_Mean] = "Mean";
    m_curvatureTypes[e_K1] = "K1";
    m_curvatureTypes[e_K2] = "K2";
    m_curvatureTypes[e_S] = "S";
    m_curvatureTypes[e_C] = "C";
    m_curvatureTypes[e_SI] = "SI";
    m_curvatureTypes[e_BE] = "BE";
    m_curvatureTypes[e_FI] = "FI";

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
    m_decimationLevelSlider->SetValue(50);
    SetDecimationLevel(0.5);
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
    wxProgressDialog *progress = new wxProgressDialog("Decimation Progress",
        "Loading Surface...", 100, this);
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

    wxString stats = wxString::Format("Decimated Surface Statistics:\n"
        "Triangles: %d\n"
        "Vertices: %d\n", m_decimatedSurface->nfaces,
        m_decimatedSurface->nvertices);
    m_decimatedStatsText->SetLabel(stats);

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
        wxString stats = wxString::Format("Original Surface Statistics:\n"
            "Triangles: %d\n"
            "Vertices: %d\n", m_origSurface->nfaces,
            m_origSurface->nvertices);
        m_origStatsText->SetLabel(stats);
    }
    else
    {
        m_origStatsText->SetLabel("");
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
    wxString valStr = wxString::Format("%.2f", val);
    m_decimationLevelText->SetLabel(valStr);
    m_decimationLevel = val;
}

///
/// Set the minimum angle between faces for decimation
///
void DecimatePanel::SetMinimumAngle(float val)
{
    wxString valStr = wxString::Format("%2.2f", val);
    m_minimumAngleText->SetLabel(valStr);
    m_minimumAngle = val;
}

///
/// Set the minimum value from the curvature histogram to color
///
void DecimatePanel::SetMinValue(float val)
{
    wxString valStr = wxString::Format("%6.1f", val);
    m_minValueText->SetLabel(valStr);
    m_currentMinValue = val;
}

///
/// Set the maximum value from the curvature histogram to color
///
void DecimatePanel::SetMaxValue(float val)
{
    wxString valStr = wxString::Format("%6.1f", val);
    m_maxValueText->SetLabel(valStr);
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
/// Set the type of curvature to display and compute the
/// curvature using MRIS library functions as is done in mris_curvature_stats
///
void DecimatePanel::SetCurvatureType(int type)
{
    if (m_decimatedSurface != NULL)
    {
        short ret = -1;

        if (type != e_None)
        {
            MRIScomputeMetricProperties(m_decimatedSurface);
            MRISsetNeighborhoodSize(m_decimatedSurface, 2);
            MRIScomputeSecondFundamentalFormDiscrete(m_decimatedSurface, 0);
        }

        switch (type)
        {
            case e_Gaussian:
                ret = MRISuseGaussianCurvature(m_decimatedSurface);
                break;
            case e_Mean:
                ret = MRISuseMeanCurvature(m_decimatedSurface);
                break;
            case e_K1:
                ret = MRISuseK1Curvature(m_decimatedSurface);
                break;
            case e_K2:
                ret = MRISuseK2Curvature(m_decimatedSurface);
                break;
            case e_S:
                ret = MRISusePrincipalCurvatureFunction(m_decimatedSurface,
                    f_sharpnessCurvature);
                break;
            case e_C:
                ret = MRISusePrincipalCurvatureFunction(m_decimatedSurface,
                    f_curvednessCurvature);
                break;
            case e_BE:
                ret = MRISusePrincipalCurvatureFunction(m_decimatedSurface,
                    f_bendingEnergyCurvature);
                break;
            case e_SI:
                ret = MRISusePrincipalCurvatureFunction(m_decimatedSurface,
                    f_shapeIndexCurvature);
                break;
            case e_FI:
                ret = MRISusePrincipalCurvatureFunction(m_decimatedSurface,
                    f_foldingIndexCurvature);
                break;
            default:
                m_renderPanel->ShowCurvature(false);
                m_saveCurvatureButton->Disable();
                m_minValueSlider->Disable();
                m_maxValueSlider->Disable();
                m_histogramCheckBox->Disable();
                m_renderPanel->ShowHistogram(false);
                return;
        }

        float minVal = 100000.0f;
        float maxVal = -100000.0f;
        for (int vno = 0; vno < m_decimatedSurface->nvertices; vno++)
        {
            if (m_decimatedSurface->vertices[vno].curv < minVal)
            {
                minVal = m_decimatedSurface->vertices[vno].curv;
            }

            if (m_decimatedSurface->vertices[vno].curv > maxVal)
            {
                maxVal = m_decimatedSurface->vertices[vno].curv;
            }
        }

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
        m_renderPanel->ShowHistogram(m_histogramCheckBox->GetValue() == 0 ? false : true);

    }
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



void DecimatePanel::OnMinimumAngleChanged( wxScrollEvent& event )
{
    int maxLevel = m_minimumAngleSlider->GetMax();
    int minLevel = m_minimumAngleSlider->GetMin();
    float valToSetf = (float)(event.GetPosition() - minLevel) / (float)(maxLevel - minLevel);

    SetMinimumAngle(valToSetf * 90.0f);
    m_defaultButton->Enable();
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

void DecimatePanel::OnApplyButtonClick( wxCommandEvent& event )
{
    DoDecimate();
}

void DecimatePanel::OnDefaultButtonClick( wxCommandEvent& event )
{
    ResetToDefault();
}

void DecimatePanel::OnWireframeCheck( wxCommandEvent& event )
{
    m_renderPanel->SetWireframe(event.IsChecked());
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

        if (MRISwriteCurvature(m_decimatedSurface, currentFilePath) != 0)
        {
            wxMessageBox(wxString::Format("ERROR: Saving file '%s'", currentFilePath.c_str()),
                "ERROR");
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


