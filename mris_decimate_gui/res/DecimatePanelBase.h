///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Dec 29 2008)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#ifndef __DecimatePanelBase__
#define __DecimatePanelBase__

#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/slider.h>
#include <wx/textctrl.h>
#include <wx/sizer.h>
#include <wx/button.h>
#include <wx/choice.h>
#include <wx/checkbox.h>
#include <wx/statbox.h>
#include <wx/panel.h>
#include <wx/statline.h>

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class DecimatePanelBase
///////////////////////////////////////////////////////////////////////////////
class DecimatePanelBase : public wxPanel 
{
	private:
	
	protected:
		wxStaticText* m_staticText3;
		wxSlider* m_decimationLevelSlider;
		wxTextCtrl* m_decimationLevelTextCtrl;
		wxStaticText* m_staticText5;
		wxSlider* m_minimumAngleSlider;
		wxTextCtrl* m_minimumAngleTextCtrl;
		
		wxButton* m_applyButton;
		wxButton* m_defaultButton;
		wxStaticText* m_staticText7;
		wxChoice* m_curvatureChoice;
		wxButton* m_saveCurvatureButton;
		wxStaticText* m_staticText8;
		wxSlider* m_minValueSlider;
		wxTextCtrl* m_minValueTextCtrl;
		wxStaticText* m_staticText10;
		wxSlider* m_maxValueSlider;
		wxTextCtrl* m_maxValueTextCtrl;
		
		wxCheckBox* m_histogramCheckBox;
		
		wxCheckBox* m_colorBarCheckBox;
		wxStaticText* m_staticText15;
		wxChoice* m_choice2;
		
		wxButton* m_saveScreenshotButton;
		wxPanel* m_panel1;
		wxStaticText* m_staticText29;
		wxTextCtrl* m_upVectorText;
		wxStaticText* m_staticText30;
		wxTextCtrl* m_cameraPositionText;
		wxStaticText* m_staticText14;
		wxTextCtrl* m_focalPointText;
		wxButton* m_setCameraButton;
		wxButton* m_resetCameraButton;
		wxPanel* m_panel2;
		
		wxStaticText* m_staticText16;
		wxStaticText* m_staticText17;
		
		wxStaticLine* m_staticline2;
		wxStaticLine* m_staticline3;
		wxStaticText* m_staticText18;
		wxStaticText* m_origTrianglesText;
		wxStaticText* m_decTrianglesText;
		wxStaticText* m_staticText21;
		wxStaticText* m_origVerticesText;
		wxStaticText* m_decVerticesText;
		wxStaticText* m_staticText25;
		wxStaticText* m_origAvgVertexAreaText;
		wxStaticText* m_decAvgVertexAreaText;
		wxStaticText* m_staticText28;
		wxStaticText* m_origAvgVertexDistText;
		wxStaticText* m_decAvgVertexDistText;
		wxStaticText* m_staticText31;
		wxStaticText* m_origCurvatureRangeText;
		wxStaticText* m_decCurvatureRangeText;
		wxStaticText* m_staticText34;
		wxStaticText* m_origCurvatureMeanText;
		wxStaticText* m_decCurvatureMeanText;
		wxStaticText* m_staticText37;
		wxStaticText* m_origCurvatureStdDevText;
		wxStaticText* m_decCurvatureStdDevText;
		
		// Virtual event handlers, overide them in your derived class
		virtual void OnDecimationLevelChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnDecimationText( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnMinimumAngleChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnMinimumAngleText( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnApplyButtonClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnDefaultButtonClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnCurvatureChoice( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnSaveCurvatureClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnMinimumValueChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnMinimumValueText( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnMaximumValueChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnMaximumValueText( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnHistogramCheckBox( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnColorBarCheckBox( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnRenderModeChoice( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnSaveScreenshotClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnUpVectorText( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnCameraPositionText( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnFocalPointText( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnSetCameraClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnResetCameraClick( wxCommandEvent& event ){ event.Skip(); }
		
	
	public:
		DecimatePanelBase( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 374,821 ), long style = wxTAB_TRAVERSAL );
		~DecimatePanelBase();
	
};

#endif //__DecimatePanelBase__
