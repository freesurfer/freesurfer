/**
 * @file  FSPointSet.h
 * @brief Base way points class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:47 $
 *    $Revision: 1.1 $
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

#ifndef FSPointSet_h
#define FSPointSet_h

#include <QObject>
#include "vtkMatrix4x4.h"
#include "CommonDataStruct.h"
#include <QList>

extern "C"
{
#include "label.h"
}


class FSVolume;
class wxWindow;
class wxCommandEvent;

struct WayPoint
{
  double pt[3];
  double value;
};

typedef QList<WayPoint> PointSet;

class FSPointSet : public QObject
{
public:
  FSPointSet( QObject* parent = NULL );
  virtual ~FSPointSet();

  bool ReadAsLabel( const QString& filename );
  bool ReadAsControlPoints( const QString& filename );
  bool WriteAsLabel( const QString& filename );
  bool WriteAsControlPoints( const QString& filename );
  
  static bool IsLabelFormat( const QString& filename );

  void UpdateLabel( PointSet& points_in, FSVolume* vol_ref );
  void LabelToPointSet( PointSet& points_out, FSVolume* vol_ref );

protected:
  // use label to save way points
  LABEL*   m_label;

};

#endif


