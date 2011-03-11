/**
 * @file  CursorFactory.cpp
 * @brief Cursor creator/loader for 2D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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

#include <wx/xrc/xmlres.h>
#include <wx/image.h>
#include <wx/mstream.h>
#include "CursorFactory.h"

#include "res/CursorFill_png.h"
#include "res/CursorPencil_png.h"
#include "res/CursorPolyline_png.h"
#include "res/CursorPan_png.h"
#include "res/CursorZoom_png.h"
#include "res/CursorColorPicker_png.h"
#include "res/CursorMeasureLine_png.h"
#include "res/CursorMeasureRectangle_png.h"
#include "res/CursorMeasurePolyline_png.h"
#include "res/CursorGrab_png.h"
#include "res/CursorContour_png.h"
#include "res/CursorAdd_png.h"
#include "res/CursorRemove_png.h"

wxCursor CursorFactory::CursorPencil    = wxCursor();
wxCursor CursorFactory::CursorFill      = wxCursor();
wxCursor CursorFactory::CursorPolyline  = wxCursor();
wxCursor CursorFactory::CursorPan       = wxCursor();
wxCursor CursorFactory::CursorZoom      = wxCursor();
wxCursor CursorFactory::CursorGrab      = wxCursor();
wxCursor CursorFactory::CursorColorPicker  = wxCursor();
wxCursor CursorFactory::CursorMeasureLine  = wxCursor();
wxCursor CursorFactory::CursorMeasureRectangle = wxCursor();
wxCursor CursorFactory::CursorMeasurePolyline  = wxCursor();
wxCursor CursorFactory::CursorContour   = wxCursor();
wxCursor CursorFactory::CursorAdd       = wxCursor();
wxCursor CursorFactory::CursorRemove    = wxCursor();

CursorFactory::CursorFactory()
{
  Initialize();
}

void CursorFactory::Initialize()
{
  wxMemoryInputStream s1( CursorPencil_png, CursorPencil_png_LEN );
  wxImage img( s1, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorPencil = wxCursor( img );

  wxMemoryInputStream s2( CursorFill_png, CursorFill_png_LEN );
  img = wxImage( s2, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 2 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 23 );
  CursorFill = wxCursor( img );

  wxMemoryInputStream s3( CursorPolyline_png, CursorPolyline_png_LEN );
  img = wxImage( s3, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorPolyline = wxCursor( img );

  wxMemoryInputStream s4( CursorPan_png, CursorPan_png_LEN );
  img = wxImage( s4, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 11 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 11 );
  CursorPan = wxCursor( img );

  wxMemoryInputStream s5( CursorZoom_png, CursorZoom_png_LEN );
  img = wxImage( s5, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 11 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 11 );
  CursorZoom = wxCursor( img );
  
  wxMemoryInputStream s6( CursorColorPicker_png, CursorColorPicker_png_LEN );
  img = wxImage( s6, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 23 );
  CursorColorPicker = wxCursor( img );
  
  wxMemoryInputStream s7( CursorMeasureLine_png, CursorMeasureLine_png_LEN );
  img = wxImage( s7, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorMeasureLine = wxCursor( img );
  
  wxMemoryInputStream s8( CursorGrab_png, CursorGrab_png_LEN );
  img = wxImage( s8, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 11 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 11 );
  CursorGrab = wxCursor( img );
    
  wxMemoryInputStream s9( CursorMeasureRectangle_png, CursorMeasureRectangle_png_LEN );
  img = wxImage( s9, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorMeasureRectangle = wxCursor( img );
  
  wxMemoryInputStream s10( CursorContour_png, CursorContour_png_LEN );
  img = wxImage( s10, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorContour = wxCursor( img );
  
  wxMemoryInputStream s11( CursorMeasurePolyline_png, CursorMeasurePolyline_png_LEN );
  img = wxImage( s11, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorMeasurePolyline = wxCursor( img );
  
  wxMemoryInputStream s12( CursorAdd_png, CursorAdd_png_LEN );
  img = wxImage( s12, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorAdd = wxCursor( img );
  
  wxMemoryInputStream s13( CursorRemove_png, CursorRemove_png_LEN );
  img = wxImage( s13, wxBITMAP_TYPE_PNG );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_X, 0 );
  img.SetOption( wxIMAGE_OPTION_CUR_HOTSPOT_Y, 0 );
  CursorRemove = wxCursor( img );
}
