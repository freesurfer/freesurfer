///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Dec 29 2008)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include "DecimatePanelBase.h"

///////////////////////////////////////////////////////////////////////////

DecimatePanelBase::DecimatePanelBase( wxWindow* parent, wxWindowID id, const wxPoint& pos, const wxSize& size, long style ) : wxPanel( parent, id, pos, size, style )
{
	wxBoxSizer* bSizer17;
	bSizer17 = new wxBoxSizer( wxVERTICAL );
	
	wxFlexGridSizer* fgSizer3;
	fgSizer3 = new wxFlexGridSizer( 7, 2, 0, 0 );
	fgSizer3->AddGrowableCol( 1 );
	fgSizer3->SetFlexibleDirection( wxBOTH );
	fgSizer3->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	
	m_staticText3 = new wxStaticText( this, wxID_ANY, wxT("Decimation Level"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText3->Wrap( -1 );
	m_staticText3->SetToolTip( wxT("The level of decimation applied to the mesh.  The resulting number of edges will be the Decimation Level times the original number of edges.") );
	
	fgSizer3->Add( m_staticText3, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );
	
	wxBoxSizer* bSizer5;
	bSizer5 = new wxBoxSizer( wxHORIZONTAL );
	
	m_decimationLevelSlider = new wxSlider( this, wxID_ANY, 50, 0, 100, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL );
	m_decimationLevelSlider->SetToolTip( wxT("The level of decimation applied to the mesh.  The resulting number of edges will be the Decimation Level times the original number of edges.") );
	
	bSizer5->Add( m_decimationLevelSlider, 2, wxALL, 5 );
	
	m_decimationLevelText = new wxStaticText( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decimationLevelText->Wrap( -1 );
	bSizer5->Add( m_decimationLevelText, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer5, 1, wxEXPAND, 5 );
	
	m_staticText5 = new wxStaticText( this, wxID_ANY, wxT("Minimum Angle"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText5->Wrap( -1 );
	m_staticText5->SetToolTip( wxT("The minimum angle in degrees between two neighboring triangles.") );
	
	fgSizer3->Add( m_staticText5, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );
	
	wxBoxSizer* bSizer6;
	bSizer6 = new wxBoxSizer( wxHORIZONTAL );
	
	m_minimumAngleSlider = new wxSlider( this, wxID_ANY, 100, 0, 900, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL );
	m_minimumAngleSlider->SetToolTip( wxT("The minimum angle in degrees between two neighboring triangles.") );
	
	bSizer6->Add( m_minimumAngleSlider, 2, wxALL, 5 );
	
	m_minimumAngleText = new wxStaticText( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_minimumAngleText->Wrap( -1 );
	bSizer6->Add( m_minimumAngleText, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer6, 1, wxEXPAND, 5 );
	
	
	fgSizer3->Add( 0, 0, 1, wxEXPAND, 5 );
	
	wxBoxSizer* bSizer14;
	bSizer14 = new wxBoxSizer( wxHORIZONTAL );
	
	m_applyButton = new wxButton( this, wxID_ANY, wxT("Decimate"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer14->Add( m_applyButton, 0, wxALL, 5 );
	
	m_defaultButton = new wxButton( this, wxID_ANY, wxT("Defaults"), wxDefaultPosition, wxDefaultSize, 0 );
	bSizer14->Add( m_defaultButton, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer14, 1, wxEXPAND, 5 );
	
	m_staticText7 = new wxStaticText( this, wxID_ANY, wxT("Display Curvature"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText7->Wrap( -1 );
	fgSizer3->Add( m_staticText7, 0, wxALL, 5 );
	
	wxBoxSizer* bSizer11;
	bSizer11 = new wxBoxSizer( wxHORIZONTAL );
	
	wxArrayString m_curvatureChoiceChoices;
	m_curvatureChoice = new wxChoice( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, m_curvatureChoiceChoices, 0 );
	m_curvatureChoice->SetSelection( 0 );
	bSizer11->Add( m_curvatureChoice, 0, wxALL, 5 );
	
	m_saveCurvatureButton = new wxButton( this, wxID_ANY, wxT("Save..."), wxDefaultPosition, wxDefaultSize, 0 );
	m_saveCurvatureButton->Enable( false );
	
	bSizer11->Add( m_saveCurvatureButton, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer11, 1, wxEXPAND, 5 );
	
	m_staticText8 = new wxStaticText( this, wxID_ANY, wxT("Minimum Value"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText8->Wrap( -1 );
	fgSizer3->Add( m_staticText8, 0, wxALL, 5 );
	
	wxBoxSizer* bSizer8;
	bSizer8 = new wxBoxSizer( wxHORIZONTAL );
	
	m_minValueSlider = new wxSlider( this, wxID_ANY, 500, 0, 1000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL );
	bSizer8->Add( m_minValueSlider, 2, wxALL, 5 );
	
	m_minValueText = new wxStaticText( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_minValueText->Wrap( -1 );
	bSizer8->Add( m_minValueText, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer8, 1, wxEXPAND, 5 );
	
	m_staticText10 = new wxStaticText( this, wxID_ANY, wxT("Maximum Value"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText10->Wrap( -1 );
	fgSizer3->Add( m_staticText10, 0, wxALL, 5 );
	
	wxBoxSizer* bSizer9;
	bSizer9 = new wxBoxSizer( wxHORIZONTAL );
	
	m_maxValueSlider = new wxSlider( this, wxID_ANY, 500, 0, 1000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL );
	bSizer9->Add( m_maxValueSlider, 2, wxALL, 5 );
	
	m_maxValueText = new wxStaticText( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_maxValueText->Wrap( -1 );
	bSizer9->Add( m_maxValueText, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer9, 1, wxEXPAND, 5 );
	
	
	fgSizer3->Add( 0, 0, 1, wxEXPAND, 5 );
	
	m_histogramCheckBox = new wxCheckBox( this, wxID_ANY, wxT("Show Histogram"), wxDefaultPosition, wxDefaultSize, 0 );
	m_histogramCheckBox->SetValue(true);
	
	fgSizer3->Add( m_histogramCheckBox, 0, wxALL, 5 );
	
	bSizer17->Add( fgSizer3, 1, wxEXPAND, 5 );
	
	wxBoxSizer* bSizer18;
	bSizer18 = new wxBoxSizer( wxVERTICAL );
	
	m_panel1 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	wxBoxSizer* bSizer19;
	bSizer19 = new wxBoxSizer( wxVERTICAL );
	
	m_staticline1 = new wxStaticLine( m_panel1, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer19->Add( m_staticline1, 0, wxEXPAND | wxALL, 5 );
	
	m_wireframeCheckBox = new wxCheckBox( m_panel1, wxID_ANY, wxT("Render Wireframe"), wxDefaultPosition, wxDefaultSize, 0 );
	
	bSizer19->Add( m_wireframeCheckBox, 0, wxALL, 5 );
	
	m_origStatsText = new wxStaticText( m_panel1, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origStatsText->Wrap( -1 );
	bSizer19->Add( m_origStatsText, 0, wxALL, 5 );
	
	m_panel1->SetSizer( bSizer19 );
	m_panel1->Layout();
	bSizer19->Fit( m_panel1 );
	bSizer18->Add( m_panel1, 1, wxEXPAND | wxALL, 5 );
	
	m_panel2 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	wxBoxSizer* bSizer20;
	bSizer20 = new wxBoxSizer( wxVERTICAL );
	
	m_decimatedStatsText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decimatedStatsText->Wrap( -1 );
	bSizer20->Add( m_decimatedStatsText, 0, wxALL, 5 );
	
	m_panel2->SetSizer( bSizer20 );
	m_panel2->Layout();
	bSizer20->Fit( m_panel2 );
	bSizer18->Add( m_panel2, 1, wxEXPAND | wxALL, 5 );
	
	bSizer17->Add( bSizer18, 1, wxEXPAND, 5 );
	
	this->SetSizer( bSizer17 );
	this->Layout();
	
	// Connect Events
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Connect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_applyButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnApplyButtonClick ), NULL, this );
	m_defaultButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnDefaultButtonClick ), NULL, this );
	m_curvatureChoice->Connect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( DecimatePanelBase::OnCurvatureChoice ), NULL, this );
	m_saveCurvatureButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnSaveCurvatureClick ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Connect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_histogramCheckBox->Connect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnHistogramCheckBox ), NULL, this );
	m_wireframeCheckBox->Connect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnWireframeCheck ), NULL, this );
}

DecimatePanelBase::~DecimatePanelBase()
{
	// Disconnect Events
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_decimationLevelSlider->Disconnect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnDecimationLevelChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_applyButton->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnApplyButtonClick ), NULL, this );
	m_defaultButton->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnDefaultButtonClick ), NULL, this );
	m_curvatureChoice->Disconnect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( DecimatePanelBase::OnCurvatureChoice ), NULL, this );
	m_saveCurvatureButton->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnSaveCurvatureClick ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_minValueSlider->Disconnect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMinimumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_histogramCheckBox->Disconnect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnHistogramCheckBox ), NULL, this );
	m_wireframeCheckBox->Disconnect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnWireframeCheck ), NULL, this );
}
