/**
 * @file  FSGroupDescriptor.h
 * @brief FSGD wrapper class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/03/29 20:35:50 $
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
 *
 */

#ifndef FSGroupDescriptor_h
#define FSGroupDescriptor_h

#include <QObject>
#include <QStringList>
#include <QMap>
#include <QColor>

extern "C"
{
#include "fsgdf.h"
}

class FSVolume;

struct FSGDDataItem
{
  QString subject_id;
  QString class_id;
  QList<double> variable_values;
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

  FSGD*   m_fsgd;
  double  m_dXStart;
  double  m_dXDelta;
  QList<FSGDDataItem> m_data;
  QList<FSGDVariable> m_variables;
  QList<FSGDClass>    m_classes;
  QString         m_title;

};

#endif


