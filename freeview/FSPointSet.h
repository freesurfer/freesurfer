/**
 * @brief Base way points class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef FSPointSet_h
#define FSPointSet_h

#include <QObject>
#include "vtkMatrix4x4.h"
#include "CommonDataStruct.h"
#include <QList>
#include <QVariantMap>



#include "label.h"



class FSVolume;

struct ControlPoint
{
  double pt[3];
  double value;
  QVariantMap info; // enhanced fields
};

typedef QList<ControlPoint> PointSet;

class FSPointSet : public QObject
{
public:
  FSPointSet( QObject* parent = NULL );
  virtual ~FSPointSet();

  bool ReadAsLabel( const QString& filename );
  bool ReadAsControlPoints( const QString& filename );
  bool WriteAsLabel( const QString& filename );
  bool WriteAsControlPoints(const QString& filename);

  static bool IsLabelFormat( const QString& filename );

  void UpdateLabel( PointSet& points_in, FSVolume* vol_ref );
  void LabelToPointSet( PointSet& points_out, FSVolume* vol_ref );

  bool ReadFromStringAsControlPoints(const QString& content);
  QString WriteAsControlPointsToString();

  bool GetCentroidRASPosition(double* pos, FSVolume* ref_vol);

  double GetMinStat()
  {
    return m_dStatMin;
  }

  double GetMaxStat()
  {
    return m_dStatMax;
  }

protected:
  void UpdateStatRange();

  // use label to save way points
  LABEL*   m_label;
  double   m_dStatMin;
  double   m_dStatMax;
};

#endif


