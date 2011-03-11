/**
 * @file  wxColorIndicator.h
 * @brief wxColorIndicator class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:43 $
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

#ifndef wxColorIndicator_h
#define wxColorIndicator_h

#include <wx/panel.h>
#include <wx/colour.h>

class WXDLLEXPORT wxColorIndicator : public wxPanel
{
public:
  wxColorIndicator();

  wxColorIndicator( wxWindow *parent,
                    wxWindowID id = wxID_ANY,
                    const wxPoint& pos = wxDefaultPosition,
                    const wxSize& size = wxDefaultSize,
                    long style = wxTAB_TRAVERSAL,
                    const wxString& name = wxT( "ColorIndicator" ) );

  virtual ~wxColorIndicator();

  bool Create( wxWindow *parent,
               wxWindowID id = wxID_ANY,
               const wxPoint& pos = wxDefaultPosition,
               const wxSize& size = wxDefaultSize,
               long style = wxTAB_TRAVERSAL,
               const wxString& name = wxT( "ColorIndicator" ) );

  wxColour GetColor();
  wxColour GetBorderColor();
  void SetColor( const wxColour& color );
  void SetBorderColor( const wxColour& color );

protected:
  void OnPaint( wxPaintEvent &event );

private:
  wxColour m_colorFill;
  wxColour m_colorBorder;

  DECLARE_DYNAMIC_CLASS( wxColorIndicator )
  DECLARE_EVENT_TABLE()
};

#endif


