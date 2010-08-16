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
	fgSizer3 = new wxFlexGridSizer( 10, 2, 0, 0 );
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
	
	m_decimationLevelTextCtrl = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	bSizer5->Add( m_decimationLevelTextCtrl, 0, wxALL, 5 );
	
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
	
	m_minimumAngleTextCtrl = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	bSizer6->Add( m_minimumAngleTextCtrl, 0, wxALL, 5 );
	
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
	
	m_minValueSlider = new wxSlider( this, wxID_ANY, 5000, 0, 10000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL );
	bSizer8->Add( m_minValueSlider, 2, wxALL, 5 );
	
	m_minValueTextCtrl = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	bSizer8->Add( m_minValueTextCtrl, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer8, 1, wxEXPAND, 5 );
	
	m_staticText10 = new wxStaticText( this, wxID_ANY, wxT("Maximum Value"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText10->Wrap( -1 );
	fgSizer3->Add( m_staticText10, 0, wxALL, 5 );
	
	wxBoxSizer* bSizer9;
	bSizer9 = new wxBoxSizer( wxHORIZONTAL );
	
	m_maxValueSlider = new wxSlider( this, wxID_ANY, 5000, 0, 10000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL );
	bSizer9->Add( m_maxValueSlider, 2, wxALL, 5 );
	
	m_maxValueTextCtrl = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	bSizer9->Add( m_maxValueTextCtrl, 0, wxALL, 5 );
	
	fgSizer3->Add( bSizer9, 1, wxEXPAND, 5 );
	
	
	fgSizer3->Add( 0, 0, 1, wxEXPAND, 5 );
	
	m_histogramCheckBox = new wxCheckBox( this, wxID_ANY, wxT("Show Histogram"), wxDefaultPosition, wxDefaultSize, 0 );
	m_histogramCheckBox->SetValue(true);
	
	fgSizer3->Add( m_histogramCheckBox, 0, wxALL, 5 );
	
	
	fgSizer3->Add( 0, 0, 1, wxEXPAND, 5 );
	
	m_colorBarCheckBox = new wxCheckBox( this, wxID_ANY, wxT("Show Color Bar"), wxDefaultPosition, wxDefaultSize, 0 );
	m_colorBarCheckBox->SetValue(true);
	
	fgSizer3->Add( m_colorBarCheckBox, 0, wxALL, 5 );
	
	m_staticText15 = new wxStaticText( this, wxID_ANY, wxT("Render Mode"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText15->Wrap( -1 );
	fgSizer3->Add( m_staticText15, 0, wxALL, 5 );
	
	wxString m_choice2Choices[] = { wxT("Surface"), wxT("Surface and Wireframe"), wxT("Wireframe") };
	int m_choice2NChoices = sizeof( m_choice2Choices ) / sizeof( wxString );
	m_choice2 = new wxChoice( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, m_choice2NChoices, m_choice2Choices, 0 );
	m_choice2->SetSelection( 0 );
	fgSizer3->Add( m_choice2, 0, wxALL|wxEXPAND, 5 );
	
	
	fgSizer3->Add( 0, 0, 1, wxEXPAND, 5 );
	
	m_saveScreenshotButton = new wxButton( this, wxID_ANY, wxT("Save Screenshot..."), wxDefaultPosition, wxDefaultSize, 0 );
	fgSizer3->Add( m_saveScreenshotButton, 0, wxALL, 5 );
	
	bSizer17->Add( fgSizer3, 0, wxEXPAND, 5 );
	
	wxBoxSizer* bSizer18;
	bSizer18 = new wxBoxSizer( wxVERTICAL );
	
	m_panel1 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	wxBoxSizer* bSizer19;
	bSizer19 = new wxBoxSizer( wxVERTICAL );
	
	wxStaticBoxSizer* sbSizer1;
	sbSizer1 = new wxStaticBoxSizer( new wxStaticBox( m_panel1, wxID_ANY, wxT("Camera") ), wxVERTICAL );
	
	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer( 4, 2, 0, 0 );
	fgSizer2->AddGrowableCol( 1 );
	fgSizer2->SetFlexibleDirection( wxBOTH );
	fgSizer2->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	
	m_staticText29 = new wxStaticText( m_panel1, wxID_ANY, wxT("Up Vector"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText29->Wrap( -1 );
	fgSizer2->Add( m_staticText29, 0, wxALL, 5 );
	
	m_upVectorText = new wxTextCtrl( m_panel1, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	fgSizer2->Add( m_upVectorText, 2, wxALL|wxEXPAND, 5 );
	
	m_staticText30 = new wxStaticText( m_panel1, wxID_ANY, wxT("Position"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText30->Wrap( -1 );
	fgSizer2->Add( m_staticText30, 0, wxALL, 5 );
	
	m_cameraPositionText = new wxTextCtrl( m_panel1, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	fgSizer2->Add( m_cameraPositionText, 2, wxALL|wxEXPAND, 5 );
	
	m_staticText14 = new wxStaticText( m_panel1, wxID_ANY, wxT("Focal Point"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText14->Wrap( -1 );
	fgSizer2->Add( m_staticText14, 0, wxALL, 5 );
	
	m_focalPointText = new wxTextCtrl( m_panel1, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	fgSizer2->Add( m_focalPointText, 0, wxALL|wxEXPAND, 5 );
	
	m_setCameraButton = new wxButton( m_panel1, wxID_ANY, wxT("Set Camera"), wxDefaultPosition, wxDefaultSize, 0 );
	fgSizer2->Add( m_setCameraButton, 0, wxALL, 5 );
	
	m_resetCameraButton = new wxButton( m_panel1, wxID_ANY, wxT("Reset Camera"), wxDefaultPosition, wxDefaultSize, 0 );
	fgSizer2->Add( m_resetCameraButton, 1, wxALL, 5 );
	
	sbSizer1->Add( fgSizer2, 1, wxEXPAND, 5 );
	
	bSizer19->Add( sbSizer1, 1, wxEXPAND, 5 );
	
	m_panel1->SetSizer( bSizer19 );
	m_panel1->Layout();
	bSizer19->Fit( m_panel1 );
	bSizer18->Add( m_panel1, 0, wxEXPAND | wxALL, 5 );
	
	m_panel2 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	wxStaticBoxSizer* sbSizer2;
	sbSizer2 = new wxStaticBoxSizer( new wxStaticBox( m_panel2, wxID_ANY, wxT("Surface Statistics") ), wxVERTICAL );
	
	wxFlexGridSizer* fgSizer31;
	fgSizer31 = new wxFlexGridSizer( 9, 3, 0, 0 );
	fgSizer31->AddGrowableCol( 1 );
	fgSizer31->AddGrowableCol( 2 );
	fgSizer31->SetFlexibleDirection( wxBOTH );
	fgSizer31->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	
	
	fgSizer31->Add( 0, 0, 1, wxEXPAND, 5 );
	
	m_staticText16 = new wxStaticText( m_panel2, wxID_ANY, wxT("Original"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText16->Wrap( -1 );
	m_staticText16->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText16, 0, wxALL, 5 );
	
	m_staticText17 = new wxStaticText( m_panel2, wxID_ANY, wxT("Decimated"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText17->Wrap( -1 );
	m_staticText17->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText17, 0, wxALL, 5 );
	
	
	fgSizer31->Add( 0, 0, 1, wxEXPAND, 5 );
	
	m_staticline2 = new wxStaticLine( m_panel2, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	fgSizer31->Add( m_staticline2, 0, wxEXPAND | wxALL, 5 );
	
	m_staticline3 = new wxStaticLine( m_panel2, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	fgSizer31->Add( m_staticline3, 0, wxEXPAND | wxALL, 5 );
	
	m_staticText18 = new wxStaticText( m_panel2, wxID_ANY, wxT("Triangles"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText18->Wrap( -1 );
	m_staticText18->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText18, 0, wxALL, 5 );
	
	m_origTrianglesText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origTrianglesText->Wrap( -1 );
	m_origTrianglesText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_origTrianglesText, 0, wxALL, 5 );
	
	m_decTrianglesText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decTrianglesText->Wrap( -1 );
	m_decTrianglesText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_decTrianglesText, 0, wxALL, 5 );
	
	m_staticText21 = new wxStaticText( m_panel2, wxID_ANY, wxT("Vertices"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText21->Wrap( -1 );
	m_staticText21->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText21, 0, wxALL, 5 );
	
	m_origVerticesText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origVerticesText->Wrap( -1 );
	m_origVerticesText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_origVerticesText, 0, wxALL, 5 );
	
	m_decVerticesText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decVerticesText->Wrap( -1 );
	m_decVerticesText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_decVerticesText, 0, wxALL, 5 );
	
	m_staticText25 = new wxStaticText( m_panel2, wxID_ANY, wxT("Avg. Vertex Area"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText25->Wrap( -1 );
	m_staticText25->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText25, 0, wxALL, 5 );
	
	m_origAvgVertexAreaText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origAvgVertexAreaText->Wrap( -1 );
	m_origAvgVertexAreaText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_origAvgVertexAreaText, 0, wxALL, 5 );
	
	m_decAvgVertexAreaText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decAvgVertexAreaText->Wrap( -1 );
	m_decAvgVertexAreaText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_decAvgVertexAreaText, 0, wxALL, 5 );
	
	m_staticText28 = new wxStaticText( m_panel2, wxID_ANY, wxT("Avg. Vertex Dist"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText28->Wrap( -1 );
	m_staticText28->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText28, 0, wxALL, 5 );
	
	m_origAvgVertexDistText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origAvgVertexDistText->Wrap( -1 );
	m_origAvgVertexDistText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_origAvgVertexDistText, 0, wxALL, 5 );
	
	m_decAvgVertexDistText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decAvgVertexDistText->Wrap( -1 );
	m_decAvgVertexDistText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_decAvgVertexDistText, 0, wxALL, 5 );
	
	m_staticText31 = new wxStaticText( m_panel2, wxID_ANY, wxT("Curv. Range"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText31->Wrap( -1 );
	m_staticText31->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText31, 0, wxALL, 5 );
	
	m_origCurvatureRangeText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origCurvatureRangeText->Wrap( -1 );
	m_origCurvatureRangeText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_origCurvatureRangeText, 0, wxALL, 5 );
	
	m_decCurvatureRangeText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decCurvatureRangeText->Wrap( -1 );
	m_decCurvatureRangeText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_decCurvatureRangeText, 0, wxALL, 5 );
	
	m_staticText34 = new wxStaticText( m_panel2, wxID_ANY, wxT("Curv. Mean"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText34->Wrap( -1 );
	m_staticText34->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText34, 0, wxALL, 5 );
	
	m_origCurvatureMeanText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origCurvatureMeanText->Wrap( -1 );
	m_origCurvatureMeanText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_origCurvatureMeanText, 0, wxALL, 5 );
	
	m_decCurvatureMeanText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decCurvatureMeanText->Wrap( -1 );
	m_decCurvatureMeanText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_decCurvatureMeanText, 0, wxALL, 5 );
	
	m_staticText37 = new wxStaticText( m_panel2, wxID_ANY, wxT("Curv. StdDev"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText37->Wrap( -1 );
	m_staticText37->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_staticText37, 0, wxALL, 5 );
	
	m_origCurvatureStdDevText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_origCurvatureStdDevText->Wrap( -1 );
	m_origCurvatureStdDevText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_origCurvatureStdDevText, 0, wxALL, 5 );
	
	m_decCurvatureStdDevText = new wxStaticText( m_panel2, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	m_decCurvatureStdDevText->Wrap( -1 );
	m_decCurvatureStdDevText->SetFont( wxFont( 8, 70, 90, 90, false, wxEmptyString ) );
	
	fgSizer31->Add( m_decCurvatureStdDevText, 0, wxALL, 5 );
	
	sbSizer2->Add( fgSizer31, 1, wxEXPAND, 5 );
	
	m_panel2->SetSizer( sbSizer2 );
	m_panel2->Layout();
	sbSizer2->Fit( m_panel2 );
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
	m_decimationLevelTextCtrl->Connect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnDecimationText ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Connect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleTextCtrl->Connect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnMinimumAngleText ), NULL, this );
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
	m_minValueTextCtrl->Connect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnMinimumValueText ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Connect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueTextCtrl->Connect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnMaximumValueText ), NULL, this );
	m_histogramCheckBox->Connect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnHistogramCheckBox ), NULL, this );
	m_colorBarCheckBox->Connect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnColorBarCheckBox ), NULL, this );
	m_choice2->Connect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( DecimatePanelBase::OnRenderModeChoice ), NULL, this );
	m_saveScreenshotButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnSaveScreenshotClick ), NULL, this );
	m_upVectorText->Connect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnUpVectorText ), NULL, this );
	m_cameraPositionText->Connect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnCameraPositionText ), NULL, this );
	m_focalPointText->Connect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnFocalPointText ), NULL, this );
	m_setCameraButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnSetCameraClick ), NULL, this );
	m_resetCameraButton->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnResetCameraClick ), NULL, this );
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
	m_decimationLevelTextCtrl->Disconnect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnDecimationText ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleSlider->Disconnect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMinimumAngleChanged ), NULL, this );
	m_minimumAngleTextCtrl->Disconnect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnMinimumAngleText ), NULL, this );
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
	m_minValueTextCtrl->Disconnect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnMinimumValueText ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_TOP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_BOTTOM, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_LINEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_LINEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_PAGEUP, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_PAGEDOWN, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_THUMBTRACK, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_THUMBRELEASE, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueSlider->Disconnect( wxEVT_SCROLL_CHANGED, wxScrollEventHandler( DecimatePanelBase::OnMaximumValueChanged ), NULL, this );
	m_maxValueTextCtrl->Disconnect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnMaximumValueText ), NULL, this );
	m_histogramCheckBox->Disconnect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnHistogramCheckBox ), NULL, this );
	m_colorBarCheckBox->Disconnect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnColorBarCheckBox ), NULL, this );
	m_choice2->Disconnect( wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler( DecimatePanelBase::OnRenderModeChoice ), NULL, this );
	m_saveScreenshotButton->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnSaveScreenshotClick ), NULL, this );
	m_upVectorText->Disconnect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnUpVectorText ), NULL, this );
	m_cameraPositionText->Disconnect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnCameraPositionText ), NULL, this );
	m_focalPointText->Disconnect( wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler( DecimatePanelBase::OnFocalPointText ), NULL, this );
	m_setCameraButton->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnSetCameraClick ), NULL, this );
	m_resetCameraButton->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( DecimatePanelBase::OnResetCameraClick ), NULL, this );
}
