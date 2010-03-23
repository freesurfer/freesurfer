/**
 * @file  SurfaceLabel.h
 * @brief The common properties available to surface label
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/03/23 18:31:10 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2007-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef SurfaceLabel_h
#define SurfaceLabel_h

#include "Broadcaster.h"
#include "Listener.h"
#include <string>

extern "C"
{
#include "label.h"
}

class LayerSurface;

class SurfaceLabel  : public Broadcaster, public Listener
{
  friend class SurfaceLabelProperties;
public:
  SurfaceLabel ( LayerSurface* surf );
  ~SurfaceLabel ();

//  void SetSurface( LayerSurface* surf );
  
  const char* GetName();
  
  void SetName( const char* name );
  
  bool LoadLabel( const char* filename );
  
  void SetColor( double r, double g, double b );
  
  double* GetColor()
  {
    return m_rgbColor;
  }
  
  void MapLabel( unsigned char* colordata, int nVertexCount );
  
protected:
  virtual void DoListenToMessage ( std::string const iMessage, void* iData, void* sender );
  
private:
  LABEL*        m_label;  
  std::string   m_strName;
  LayerSurface* m_surface;
  double        m_rgbColor[3];
  bool          m_bTkReg;
};

#endif
