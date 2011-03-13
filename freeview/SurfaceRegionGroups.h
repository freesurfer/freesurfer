/**
 * @file  SurfaceRegionGroups.h
 * @brief Surface region from a surface selection in 3D view.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:18 $
 *    $Revision: 1.5 $
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


