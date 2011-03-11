/**
 * @file  Region2D.h
 * @brief Region2D data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:40 $
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

#ifndef Region2D_h
#define Region2D_h

#include <wx/wx.h>
#include "Broadcaster.h"
#include "Listener.h"

class RenderView2D;
class vtkRenderer;

class Region2D : public Broadcaster, public Listener
{
public:
  Region2D( RenderView2D* view );
  virtual ~Region2D();

  virtual void Offset( int nX, int nY ) = 0;
  
  virtual void UpdatePoint( int nIndex, int nX, int nY ) = 0;
  
  virtual bool Contains( int nX, int nY, int* nPointIndex = NULL ) = 0; 
  
  virtual void Highlight( bool bHighlight = true ) {}
  
  virtual void AppendProp( vtkRenderer* renderer ) = 0;
  
  virtual void Show( bool bshow );
  
  virtual void Update() {}
  
  virtual void UpdateStats();
  
  virtual void UpdateSlicePosition( int nPlane, double pos );
  
  virtual void GetWorldPoint( int nIndex, double* pt ) = 0;
  
  wxString GetShortStats()
  {
    return m_strShortStats;
  }
  
  wxArrayString GetLongStats()
  {
    return m_strsLongStats;
  }

protected:
  RenderView2D* m_view;
  wxString      m_strShortStats;
  wxArrayString m_strsLongStats;
};

#endif


