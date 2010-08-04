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
#include <wx/sizer.h>
#include <wx/button.h>
#include <wx/choice.h>
#include <wx/checkbox.h>
#include <wx/statline.h>
#include <wx/panel.h>

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
		wxStaticText* m_decimationLevelText;
		wxStaticText* m_staticText5;
		wxSlider* m_minimumAngleSlider;
		wxStaticText* m_minimumAngleText;
		
		wxButton* m_applyButton;
		wxButton* m_defaultButton;
		wxStaticText* m_staticText7;
		wxChoice* m_curvatureChoice;
		wxButton* m_saveCurvatureButton;
		wxStaticText* m_staticText8;
		wxSlider* m_minValueSlider;
		wxStaticText* m_minValueText;
		wxStaticText* m_staticText10;
		wxSlider* m_maxValueSlider;
		wxStaticText* m_maxValueText;
		
		wxCheckBox* m_histogramCheckBox;
		wxPanel* m_panel1;
		wxStaticLine* m_staticline1;
		wxCheckBox* m_wireframeCheckBox;
		wxStaticText* m_origStatsText;
		wxPanel* m_panel2;
		wxStaticText* m_decimatedStatsText;
		
		// Virtual event handlers, overide them in your derived class
		virtual void OnDecimationLevelChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnMinimumAngleChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnApplyButtonClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnDefaultButtonClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnCurvatureChoice( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnSaveCurvatureClick( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnMinimumValueChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnMaximumValueChanged( wxScrollEvent& event ){ event.Skip(); }
		virtual void OnHistogramCheckBox( wxCommandEvent& event ){ event.Skip(); }
		virtual void OnWireframeCheck( wxCommandEvent& event ){ event.Skip(); }
		
	
	public:
		DecimatePanelBase( wxWindow* parent, wxWindowID id = wxID_ANY, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 389,457 ), long style = wxTAB_TRAVERSAL );
		~DecimatePanelBase();
	
};

#endif //__DecimatePanelBase__
