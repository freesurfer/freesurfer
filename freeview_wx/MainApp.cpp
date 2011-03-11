/**
 * @file  MainApp.cpp
 * @brief The application that starts all.
 *
 * Manages the 'current' tool, view, and layer. Holds one or more
 * views and lays them out. Top level document loading and creation of
 * layers. Runs the info area and gets updates from mouseovered
 * objects. Runs the toolbars. Manages the UI panels for settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
 *    $Revision: 1.1 $
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

#ifdef __WXMAC__
#include <ApplicationServices/ApplicationServices.h>
#endif

#include "wx/wxprec.h"
#include "MainApp.h"
#include "MainWindow.h"
#include "LayerPropertiesMRI.h"
#include "LayerMRI.h"
#include "res/resource.h"
#include "wxColorIndicatorXmlHandler.h"
#include "wxHistogramWidgetXmlHandler.h"
#include "MyCmdLineParser.h"
#include "CursorFactory.h"
#include "RenderView2D.h"
#include <wx/wx.h>
#include <wx/xrc/xmlres.h>
#include <wx/config.h>
#include <wx/tooltip.h>
#include <wx/datetime.h>
#include <wx/fs_mem.h>
#include <stdlib.h>
extern "C"
{
#include "error.h"
#include "mri.h"
}

void MyErrorExit( int nCode )
{
  throw nCode;
}

wxString Progname = _("freeview");

IMPLEMENT_APP( MainApp )

BEGIN_EVENT_TABLE( MainApp, wxApp )
// EVT_KEY_UP( MainApp::OnKeyDown )
END_EVENT_TABLE()

MainApp::MainApp() :
  m_wndMain( NULL )
{
#ifdef __WXMAC__
  // In order to run executable without a '.app' file, this code is required
  // as described here:
  //  http://wiki.wxwidgets.org/WxMac-specific_topics
  // This I find quite useful for debugging so that you don't need to always
  // build the .app and can just run "./freeview"
  ProcessSerialNumber PSN;
  GetCurrentProcess(&PSN);
  TransformProcessType(&PSN,kProcessTransformToForegroundApplication);

  // Set about menu info
  wxApp::s_macHelpMenuTitleName = "Help";
#endif
}

MainApp::~MainApp()
{}

bool MainApp::OnInit()
{
  // set error exit
  ErrorSetExitFunc( &MyErrorExit );
  
  setenv("SURFER_FRONTDOOR","",1);
  
  CmdLineEntry cmdLineDesc[] =
  {
    CmdLineEntry( CMD_LINE_OPTION, "v", "volume", "<FILE>...", "Load one or multiple volume files. Available sub-options are: \n\n':colormap=name' Set colormap for display. Valid names are grayscale/lut/heat/jet/gecolor/nih. \n\n':grayscale=min,max' Set grayscale window values.\n\n':heatscale=min,mid,max' Set heat scale values.\n\n':heatscaleoptions=option1[,option2]' Set heat scale options. Options can be 'truncate','invert', or both.\n\n':colorscale=min,max' Set generic colorscale values for jet/gecolor/nih.\n\n':lut=name' Set lookup table to the given name. Name can be the name of a stock color table or the filename of a color table file.\n\n':vector=flag' Display 3 frame volume as vectors. flag can be 'yes', 'true' or '1'.\n\n':tensor=flag' Display 9 frame volume as tensors. flag can be 'yes', 'true' or '1'.\n\n':render=flag' When displaying as vectors or tensors, render the glyph in the given form. For vector, flag can be 'line' as simple line or 'bar' as 3D bar (might be slow). For tensor, flag can be 'boxoid' or 'ellipsoid' (slow!).\n\n':inversion=flag' When displaying as vectors or tensors, invert the given component of the vectors. Valid flags are 'x', 'y' and 'z'.\n\n':reg=reg_filename' Set registration file for the volume. reg_filename can contain relative path to the volume file.\n\n':sample=method' Set the sample method when resampling is necessary. method can be 'nearest' (default) or 'trilinear'.\n\n':opacity=value' Set the opacity of the volume layer. value ranges from 0 to 1.\n\n':isosurface=low_threshold,high_threshold' Set 3D display as isosurface. High_threshold is optional. If no threshold or simply 'on' is given, threshold will be either automatically determined or retrieved from the save previously settings.\n\n':color=name' Set color of the isosurface. Name can be a generic color name such as 'red' or 'lightgreen', or three integer values as RGB values ranging from 0 to 255. For example '255,0,0' is the same as 'red'.\n\n':surface_region=file' Load isosurface region(s) from the given file. isosurface display will automatically be turned on.\n\n':name=display_name' Set the display name of the volume.\n\n':lock=lock_status' Lock the volume layer so it will not be moved in the layer stack. Status can be '1' or 'true'.\n\n':visible=visibility' Set the initial visibility of the volume. Visibility can be '1' or '0' or 'true' or 'false'.\n\nExample:\nfreeview -v T1.mgz:colormap=heatscale:heatscale=10,100,200\n", 1, 100 ),
    CmdLineEntry( CMD_LINE_SWITCH, "r", "resample", "", "Resample oblique data to standard RAS." ),
    CmdLineEntry( CMD_LINE_SWITCH, "conform", "conform", "", "Conform the volume to the first loaded volume." ),
    CmdLineEntry( CMD_LINE_SWITCH, "trilinear", "trilinear", "", "Use trilinear as the default resample method." ),
    CmdLineEntry( CMD_LINE_OPTION, "dti", "dti", "<VECTOR> <FA>...", "Load one or more dti volumes. Need two files for each dti volume. First one is vector file. Second one is FA (brightness) file.", 2, 100 ),
    CmdLineEntry( CMD_LINE_OPTION, "f", "surface", "<FILE>...", "Load one or multiple surface files. Available sub-options are:\n\n':curvature=curvature_filename' Load curvature data from the given curvature file. By default .curv file will be loaded if available.\n\n':overlay=overlay_filename' Load overlay data from file.\n\n':color=colorname' Set the base color of the surface. Color can be a color name such as 'red' or 3 values as RGB components of the color, e.g., '255,0,0'.\n\n':edgecolor=colorname' Set the color of the slice intersection outline on the surface. \n\n':edgethickness=thickness' Set the thickness of the slice intersection outline on the surface. set 0 to hide it.\n\n':annot=filenames' Set annotation files to load.\n\n':name=display_name' Set the display name of the surface.\n\n':offset=x,y,z' Set the position offset of the surface. Useful for connectivity display.\n\n':visible=visibility' Set the initial visibility of the surface. Visibility can be '1' or '0' or 'true' or 'false'.\n\n':vector=filename' Load a vector file for display.\n\n':target_surf=filename' Load a target surface file for vectors to project on for 2D display.\n\n':label=filename' Load a surface label file.\n", 1, 100 ),
    CmdLineEntry( CMD_LINE_OPTION, "l", "label", "<FILE>...", "Load one or multiple label(ROI) files. Available sub-options are:\n\n':ref=ref_volume' Enter the name of the reference volume for this label file. The volume is one of the volumes given by -v option. \n", 1, 100 ),
    CmdLineEntry( CMD_LINE_OPTION, "w", "way-points", "<FILE>...", "Load one or multiple way points files. Available sub-options are:\n\n':color=name' Set color of the way points. Name can be a generic color name such as 'red' or 'lightgreen', or three integer values as RGB values ranging from 0 to 255. For example '255,0,0' is the same as 'red'.\n\n':splinecolor=name' Set color of the spline.\n\n':radius=value' Set radius of the way points.\n\n':splineradius=value' Set radius of the spline tube.\n\n':name=display_name' Set the display name of the way points.\n\n':visible=visibility' Set the initial visibility of the way points. Visibility can be '1' or '0' or 'true' or 'false'.\n", 1, 100 ),
    CmdLineEntry( CMD_LINE_OPTION, "c", "control-points", "<FILE>...", "Load one or multiple control points files. Available sub-options are:\n\n':color=name' Set color of the control points. Name can be a generic color name such as 'red' or 'lightgreen', or three integer values as RGB values ranging from 0 to 255. For example '255,0,0' is the same as 'red'.\n\n':radius=value' Set radius of the control points.\n\n':name=display_name' Set the display name of the control points.\n\n':visible=visibility' Set the initial visibility of the control points. Visibility can be '1' or '0' or 'true' or 'false'.\n", 1, 100 ),
    CmdLineEntry( CMD_LINE_OPTION, "p-labels", "p-labels", "<FILES>...", "Load multiple p-label volume files.\n", 1, 100 ),
    CmdLineEntry( CMD_LINE_OPTION, "p-prefix", "p-prefix", "<PREFIX>...", "Set the file name prefix for p-label volume. program will use this to figure out label name from file name.\n", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "p-lut", "p-lut", "<NAME>...", "Set the look up table name to use for p-label display. name can be the name of a stock lookup table or the file name of a lookup table file. default is the default freesurfer look up table.\n", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "conn", "connectivity", "<DATA> <COLORMAP>", "Load connectivity data files.\n", 2, 2 ),
    CmdLineEntry( CMD_LINE_OPTION, "ss", "screenshot", "<FILE>", "Take a screen shot of the main viewport and then quit the program.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "viewport", "viewport", "<NAME>", "Set the main viewport as given. Accepted names are 'sagittal' or 'x', 'coronal' or 'y', 'axial' or 'z' and '3d'.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "viewsize", "viewsize", "<width> <height>", "Set the size of the main viewport. The size of the whole window will be changed accordingly.", 2, 2 ),
    CmdLineEntry( CMD_LINE_OPTION, "zoom", "zoom", "<FACTOR>", "Set zoom factor of the main viewport.", 1, 1 ),     
    CmdLineEntry( CMD_LINE_OPTION, "ras", "ras", "<X> <Y> <Z>", "Set cursor location at the given RAS coordinate.", 3, 3 ),    
    CmdLineEntry( CMD_LINE_OPTION, "slice", "slice", "<X> <Y> <Z>", "Set cursor location at the given slice numbers of the first loaded volume.", 3, 3 ),
 //   CmdLineEntry( CMD_LINE_SWITCH, "disable3d", "disable3d", "", "Disable rendering in 3D view (to save graphics resource)." ),
    CmdLineEntry( CMD_LINE_NONE )
  };
  
  char progDesc[] = "Volume/Surface viewer for freesurfer.";

  MyCmdLineParser cmd( (const char*)"freeview", (CmdLineEntry*)cmdLineDesc );
  cmd.SetProgramDescription( progDesc );
  if ( !cmd.Parse( argc, (char**)argv ) )
  {
    return false;
  }

  InitializeGUI();
  
  string_array sa;
  string_array floatingArgs = cmd.GetFloatingArguments(); 
  
  if ( cmd.Found( "disable3d" ) )
  {
    m_wndMain->GetRenderView(3)->SetRenderDisabled(true);
  }
  if ( cmd.Found( "trilinear" ) )
  {
    m_wndMain->SetDefaultSampleMethod( SAMPLE_TRILINEAR );
  }
  if ( cmd.Found( "conform" ) )
  {
    m_wndMain->SetDefaultConform( true );
  }
  if ( cmd.Found( "viewport", &sa ) )
  {
    wxString strg = wxString::FromAscii( sa[0].c_str() ).Lower();
    wxArrayString script;
    script.Add( _( "setviewport" ) );
    if ( strg == _( "sagittal" ) || strg == _( "x" ) )
      script.Add( _( "x" ) );
    else if ( strg == _( "coronal" ) || strg == _( "y" ) )
      script.Add( _( "y" ) );
    else if ( strg == _( "axial" ) || strg == _( "z" ) )
      script.Add( _( "z" ) );
    else if ( strg == _( "3d" ) )
      script.Add( _( "3d" ) );
    else
    {
      cerr << "Unrecognized viewport name '" << sa[0].c_str() << "'." << endl;
      return false;
    }
    m_wndMain->AddScript( script );
  }
  
  if ( cmd.Found( "viewsize", &sa ) )
  {
    wxArrayString script;
    script.Add( _( "setviewsize" ) );
    script.Add( wxString::FromAscii(sa[0].c_str()) );
    script.Add( wxString::FromAscii(sa[1].c_str()) );
    m_wndMain->AddScript( script );
  }
  
  
  if ( floatingArgs.size() > 0 )
  {
    for ( size_t i = 0; i < floatingArgs.size(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadvolume") );
      script.Add( wxString::FromAscii(floatingArgs[i].c_str()) );

      if ( cmd.Found( "r" ) )
        script.Add( _("r") );

      m_wndMain->AddScript( script );
    }
  }  
  
  int nRepeats = cmd.GetNumberOfRepeats( "v" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    cmd.Found( "v", &sa, n );
    for ( size_t i = 0; i < sa.size(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadvolume") );
      script.Add( wxString::FromAscii(sa[i].c_str()) );

      if ( cmd.Found( "r" ) )
        script.Add( _("r") );

      m_wndMain->AddScript( script );
    }
  }
  
  nRepeats = cmd.GetNumberOfRepeats( "dti" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    cmd.Found( "dti", &sa, n );
    for ( size_t i = 0; i < sa.size()/2; i++ )
    {
      wxArrayString script;
      script.Add( _("loaddti") );
      script.Add( wxString::FromAscii(sa[i*2].c_str()) );
      script.Add( wxString::FromAscii(sa[i*2+1].c_str()) );
      if ( cmd.Found( "r" ) )
        script.Add( _("r") );
      m_wndMain->AddScript( script );
    }
  }
  
  nRepeats = cmd.GetNumberOfRepeats( "l" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    cmd.Found( "l", &sa, n );
    for (size_t i = 0; i < sa.size(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadroi") );
      script.Add( wxString::FromAscii(sa[i].c_str()) );
      m_wndMain->AddScript( script );
    }
  }
  
  
  nRepeats = cmd.GetNumberOfRepeats( "f" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    cmd.Found( "f", &sa, n );
    for (size_t i = 0; i < sa.size(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadsurface") );
      script.Add( wxString::FromAscii(sa[i].c_str()) );
      if ( cmd.Found( "r" ) )
        script.Add( _("r") );
      m_wndMain->AddScript( script );
    }
  }
  
  nRepeats = cmd.GetNumberOfRepeats( "w" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    cmd.Found( "w", &sa, n );
    for ( size_t i = 0; i < sa.size(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadwaypoints") );
      script.Add( wxString::FromAscii(sa[i].c_str()) );
      if ( cmd.Found( "r" ) )
        script.Add( _("r") );
      m_wndMain->AddScript( script );
    }
  }
  
  nRepeats = cmd.GetNumberOfRepeats( "c" );
  for ( int n = 0; n < nRepeats; n++ )
  {
    cmd.Found( "c", &sa, n );
    for ( size_t i = 0; i < sa.size(); i++ )
    {
      wxArrayString script;
      script.Add( _("loadcontrolpoints") );
      script.Add( wxString::FromAscii(sa[i].c_str()) );
      if ( cmd.Found( "r" ) )
        script.Add( _("r") );
      m_wndMain->AddScript( script );
    }
  }
  
  if ( cmd.Found( "p-labels", &sa ) )
  {
    wxString filenames = sa[0].c_str();

    for ( size_t i = 1; i < sa.size(); i++ )
    {
      filenames += ";" + wxString(sa[i].c_str());
    }
    wxArrayString script;
    script.Add( _( "loadpvolumes" ) );
    script.Add( filenames );
    sa.clear();
    if ( cmd.Found( "p-prefix", &sa ) )
      script.Add( sa[0].c_str() );
    else
      script.Add( _("") );
    sa.clear();
    if ( cmd.Found( "p-lut", &sa ) )
      script.Add( sa[0].c_str() );
    m_wndMain->AddScript( script );
  }
  
  if ( cmd.Found( "conn", &sa ) )
  {
    wxArrayString script;
    script.Add( _( "loadconnectivity" ) );
    script.Add( sa[0].c_str() );
    script.Add( sa[1].c_str() );
    m_wndMain->AddScript( script );
  }
  
  if ( cmd.Found( "ras", &sa ) )
  {
    double ras[3];
    if ( !wxString::FromAscii( sa[0].c_str() ).ToDouble( &(ras[0]) ) || 
          !wxString::FromAscii( sa[1].c_str() ).ToDouble( &(ras[1]) ) || 
          !wxString::FromAscii( sa[2].c_str() ).ToDouble( &(ras[2]) ) )
    {
      cerr << "Invalid argument for 'ras'. Arguments must be valid float values." << endl;
      return false;
    }
    wxArrayString script;
    script.Add( _( "ras" ) );
    script.Add( wxString::FromAscii(sa[0].c_str()) );
    script.Add( wxString::FromAscii(sa[1].c_str()) );
    script.Add( wxString::FromAscii(sa[2].c_str()) );
    m_wndMain->AddScript( script );
  }  
  
  if ( cmd.Found( "slice", &sa ) )
  {
    long slice[3];
    if ( !wxString::FromAscii( sa[0].c_str() ).ToLong( &(slice[0]) ) || 
          !wxString::FromAscii( sa[1].c_str() ).ToLong( &(slice[1]) ) || 
          !wxString::FromAscii( sa[2].c_str() ).ToLong( &(slice[2]) ) )
    {
      cerr << "Invalid argument for 'slice'. Arguments must be valid integers." << endl;
      return false;
    }
    wxArrayString script;
    script.Add( _( "slice" ) );
    script.Add( wxString::FromAscii(sa[0].c_str()) );
    script.Add( wxString::FromAscii(sa[1].c_str()) );
    script.Add( wxString::FromAscii(sa[2].c_str()) );
    m_wndMain->AddScript( script );
  }
  
  if ( cmd.Found( "zoom", &sa ) )
  {
    double dValue;
    if ( !wxString::FromAscii( sa[0].c_str() ).ToDouble( &dValue ) || dValue == 0 )
    {
      cerr << "Invalid argument for 'zoom'. Argument must be a valid float value." << endl;
      return false;
    }
    wxArrayString script;
    script.Add( _( "zoom" ) );
    script.Add( wxString::FromAscii(sa[0].c_str()) );
    m_wndMain->AddScript( script );
  }
  
  if ( cmd.Found( "ss", &sa ) )
  {
    wxArrayString script;
    script.Add( _( "screencapture" ) );
    script.Add( wxString::FromAscii(sa[0].c_str()) );
    m_wndMain->AddScript( script );
    script.Clear();
    script.Add( _( "quit" ) );
    m_wndMain->AddScript( script );
  }

  m_wndMain->RunScript();
  return true;
}

void MainApp::InitializeGUI()
{
  wxXmlResource::Get()->InitAllHandlers();
  wxXmlResource::Get()->AddHandler(new wxColorIndicatorXmlHandler);
  wxXmlResource::Get()->AddHandler(new wxHistogramWidgetXmlHandler);

#if 0 // for development

  wxXmlResource::Get()->Load( "res/MainWindow.xrc" );
  wxXmlResource::Get()->Load( "res/PanelVolume.xrc" );
  wxXmlResource::Get()->Load( "res/PanelROI.xrc" );
  wxXmlResource::Get()->Load( "res/PanelSurface.xrc" );
  wxXmlResource::Get()->Load( "res/PanelWayPoints.xrc" );
  wxXmlResource::Get()->Load( "res/PanelSceneSetting.xrc" );
  wxXmlResource::Get()->Load( "res/DialogEditLookupTable.xrc" );
  wxXmlResource::Get()->Load( "res/DialogLoadVolume.xrc" );
  wxXmlResource::Get()->Load( "res/DialogLoadDTI.xrc" );
  wxXmlResource::Get()->Load( "res/DialogNewVolume.xrc" );
  wxXmlResource::Get()->Load( "res/DialogNewROI.xrc" );
  wxXmlResource::Get()->Load( "res/DialogNewWayPoints.xrc" );
  wxXmlResource::Get()->Load( "res/DialogOptimalVolume.xrc" );
  wxXmlResource::Get()->Load( "res/DialogPreferences.xrc" );
  wxXmlResource::Get()->Load( "res/WindowQuickReference.xrc" );
  wxXmlResource::Get()->Load( "res/ToolWindowEdit.xrc" );
  wxXmlResource::Get()->Load( "res/ToolWindowMeasure.xrc" );
  wxXmlResource::Get()->Load( "res/WindowHistogram.xrc" );
  wxXmlResource::Get()->Load( "res/WindowOverlayConfiguration.xrc" );
  wxXmlResource::Get()->Load( "res/WindowConnectivityConfiguration.xrc" );
  wxXmlResource::Get()->Load( "res/DialogWriteMovieFrames.xrc" );
  wxXmlResource::Get()->Load( "res/DialogRepositionSurface.xrc" );
  wxXmlResource::Get()->Load( "res/DialogCropVolume.xrc" );
     
#else  // for release
  InitXmlResource();

#ifdef __WXMAC__
  // Set the about item for the Mac menu
  wxApp::s_macAboutMenuItemId = XRCID("ID_HELP_ABOUT");
#endif


#endif
  SetVendorName( _("Massachusetts General Hospital") );
  SetAppName( Progname );

  wxImage::AddHandler( new wxPNGHandler );
  wxFileSystem::AddHandler( new wxMemoryFSHandler );
  CursorFactory::Initialize();

  CreateMainWindow();
}

void MainApp::CreateMainWindow()
{
  // create the main application window
  m_wndMain = new MainWindow();
  m_wndMain->SetTitle( Progname );

  // and show it (the frames, unlike simple controls, are not shown when
  // created initially)
  SetTopWindow( m_wndMain );
  m_wndMain->Show( true );

// wxToolTip::SetDelay(0);
}

int MainApp::OnExit()
{
  // clean up: Set() returns the active config object as Get() does, but unlike
  // Get() it doesn't try to create one if there is none (definitely not what
  // we want here!)
  delete wxConfigBase::Set( (wxConfigBase *) NULL );

  return 0;
}

int MainApp::FilterEvent( wxEvent& event )
{
  if ( !m_wndMain )
    return -1;
  
  if ( event.GetEventType() == wxEVT_KEY_DOWN )
  {
    wxKeyEvent& e = (wxKeyEvent&)event;
    // catch PageUp/Down key for slice scrolling
    if ( e.GetKeyCode() == WXK_PAGEUP || e.GetKeyCode() == WXK_PAGEDOWN )
    {
      RenderView* view = m_wndMain->GetActiveView();
   //   if ( !view )
   //     view = m_wndMain->GetPreviousActiveView();
        
      // view is not 3D view 
      if ( view && view != m_wndMain->GetRenderView( 3 ) )
      {
        ( (RenderView2D*)view )->MoveSlice( e.GetKeyCode() == WXK_PAGEUP ? 1 : -1 );
        return 1;
      }
    }
    // catch transparency setting hot key
    else if ( e.GetModifiers() == wxMOD_ALT &&
              ( e.GetKeyCode() == 'A' || e.GetKeyCode() == 'S' ) )
    {
      LayerMRI* mri = (LayerMRI*)m_wndMain->GetActiveLayer( "MRI" );
      if ( mri )
      {
        double dOpacity = mri->GetProperties()->GetOpacity();
        dOpacity += ( e.GetKeyCode() == 'A' ? -0.1 : 0.1 );
        dOpacity = ( dOpacity > 1. ? 1. : (dOpacity < 0. ? 0. : dOpacity ) );
        mri->GetProperties()->SetOpacity( dOpacity ); 
        return 1;
      }
    }
    else if ( e.GetModifiers() == wxMOD_ALT &&
              e.GetKeyCode() == 'L' )
    {
      LayerMRI* mri = (LayerMRI*)m_wndMain->GetActiveLayer( "MRI" );
      if ( mri /*&& mri->GetProperties()->GetColorMap() == LayerPropertiesMRI::LUT*/ )
      {
        mri->GetProperties()->SetShowLabelOutline( !mri->GetProperties()->GetShowLabelOutline() ); 
        return 1;
      }
    }
  }
  
  return -1;
}
