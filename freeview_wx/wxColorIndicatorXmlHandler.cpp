/**
 * @file  wxColorIndicatorXmlHandler.cpp
 * @brief wxColorIndicatorXmlHandler class.
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


#include "wxColorIndicator.h"
#include "wxColorIndicatorXmlHandler.h"


IMPLEMENT_DYNAMIC_CLASS( wxColorIndicatorXmlHandler, wxXmlResourceHandler )

wxColorIndicatorXmlHandler::wxColorIndicatorXmlHandler() : wxXmlResourceHandler()
{
  AddWindowStyles();
}

wxObject *wxColorIndicatorXmlHandler::DoCreateResource()
{
  XRC_MAKE_INSTANCE( ctrl, wxColorIndicator )

  ctrl->Create( m_parentAsWindow,
                GetID(),
                GetPosition(),
                GetSize(),
                GetStyle(),
                GetName() );

  SetupWindow( ctrl );

  return ctrl;
}

bool wxColorIndicatorXmlHandler::CanHandle( wxXmlNode *node )
{
  return IsOfClass( node, wxT( "wxColorIndicator" ) );
}
