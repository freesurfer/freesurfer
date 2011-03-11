/**
 * @file  wxColorIndicatorXmlHandler.h
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

#ifndef wxColorIndicatorXmlHandler_h
#define wxColorIndicatorXmlHandler_h


#include "wx/xrc/xmlres.h"



class WXDLLIMPEXP_XRC wxColorIndicatorXmlHandler : public wxXmlResourceHandler
{
  DECLARE_DYNAMIC_CLASS( wxColorIndicatorXmlHandler )


public:
  wxColorIndicatorXmlHandler();
  ~wxColorIndicatorXmlHandler()
  {}

  wxObject *DoCreateResource();

  bool CanHandle( wxXmlNode *node );

};

#endif


