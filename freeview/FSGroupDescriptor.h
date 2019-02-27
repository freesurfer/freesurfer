/**
 * @file  FSGroupDescriptor.h
 * @brief FSGD wrapper class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2016/02/03 21:38:19 $
 *    $Revision: 1.2 $
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

#ifndef FSGroupDescriptor_h
#define FSGroupDescriptor_h

#include <QObject>
#include <QStringList>
#include <QMap>
#include <QColor>



#include "fsgdf.h"
#undef isfinite  // conflicts with vnl


class FSVolume;

struct FSGDDataItem
{
  QString subject_id;
  int class_id;
  QList<double> variable_values;
  double measurement;
};

struct FSGDVariable
{
  QString label;
  double range[2];
};

struct FSGDClass
{
  QString label;
  QString marker;
  QColor  color;
};

class FSGroupDescriptor : public QObject
{
  Q_OBJECT
public:
  FSGroupDescriptor( QObject* parent );
  virtual ~FSGroupDescriptor();

  bool Read( const QString& filename );

  void UpdateData(int nVertex);

  FSGD*   m_fsgd;
  double  m_dXStart;
  double  m_dXDelta;
  QList<FSGDClass>    m_classes;
  QList<FSGDVariable> m_variables;
  QList<FSGDDataItem> m_data;
  int             m_nVertexNum;
  QString         m_title;
  QString         m_measureName;
  double m_dMeasurementRange[2];
};

#endif


