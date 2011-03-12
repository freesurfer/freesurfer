/**
 * @file  SurfaceRegionGroups.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:54 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2008-2009,
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

#ifndef SurfaceRegionGroups_h
#define SurfaceRegionGroups_h

#include <QObject>
#include <QColor>
#include <QList>

class LayerMRI;
class SurfaceRegion;

class SurfaceRegionGroups : public QObject
{
  Q_OBJECT
public:
  SurfaceRegionGroups( LayerMRI* owner );
  virtual ~SurfaceRegionGroups();

  QColor GetGroupColor( int nGroup );
  void SetGroupColor( int nGroup, const QColor& color );
  
  int GetGroupIdRange( SurfaceRegion* reg );

private:
  LayerMRI*       m_mri;
  QList<QColor>   m_colors;
};

#endif


