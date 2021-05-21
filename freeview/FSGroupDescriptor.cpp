/**
 * @brief FSGD wrapper class that takes care of I/O and data conversion.
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


#include "FSGroupDescriptor.h"
#include <stdexcept>
#include "vtkImageData.h"
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QStringList>

FSGroupDescriptor::FSGroupDescriptor( QObject* parent ) : QObject( parent ),
  m_fsgd( NULL ),
  m_dXStart(0),
  m_dXDelta(1),
  m_nVertexNum(-1)
{}

FSGroupDescriptor::~FSGroupDescriptor()
{
  if ( m_fsgd )
  {
    ::gdfFree( &m_fsgd );
  }
}

bool FSGroupDescriptor::Read( const QString& filename )
{
  // first read XAxis start and delta
  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    return false;

  bool bOK = false;
  while (!file.atEnd())
  {
    QString line = file.readLine().trimmed();
    if (!line.isEmpty())
    {
      QStringList list = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
      if (list[0].toLower() == "xaxis" && list.size() == 3)
      {
        m_dXStart = list[1].toDouble(&bOK);
        if (bOK)
          m_dXDelta = list[2].toDouble(&bOK);
      }
    }
  }
  if (!bOK)
  {
    //    cout << "Could not find XAxis tag in file. Using default one." << endl;
    //    m_dXStart = 0;
    //    m_dXDelta = 1;
  }
  file.close();


  if ( m_fsgd )
  {
    ::gdfFree( &m_fsgd );
  }

  m_fsgd = ::gdfRead( filename.toLatin1().data(), 1 );
  if ( m_fsgd == NULL )
  {
    cerr << "gdfRead failed\n";
    return false;
  }
  float fMinValue, fMaxValue;
  MRIvalRange( m_fsgd->data, &fMinValue, &fMaxValue );
  m_dMeasurementRange[0] = fMinValue;
  m_dMeasurementRange[1] = fMaxValue;

  FSGDVariable gdv;
  for (int i = 0; i < m_fsgd->nvariables; i++)
  {
    gdv.label = m_fsgd->varlabel[i];
    m_variables << gdv;
  }

  QList<QColor> stockColors;
  stockColors << Qt::blue << Qt::red << Qt::green << Qt::yellow << Qt::cyan;
  QStringList stockMarkers;
  stockMarkers << "square" << "circle" << "plus" << "dot" << "asterisk" << "triangle" << "cross" << "diamond" ;
  for (int i = 0; i < m_fsgd->nclasses; i++)
  {
    FSGDClass gdc;
    gdc.label = m_fsgd->classlabel[i];
    gdc.marker = m_fsgd->classmarker[i];
    if (gdc.marker.isEmpty())
      gdc.marker = stockMarkers[i%stockMarkers.size()];
    gdc.color = QColor(m_fsgd->classcolor[i]);
    if (!gdc.color.isValid())
      gdc.color = stockColors[i%stockColors.size()];
    m_classes << gdc;
  }

  for (int i = 0; i < m_fsgd->ninputs; i++)
  {
    FSGDDataItem scd;
    scd.subject_id = m_fsgd->subjid[i];
    for (int n = 0; n < m_fsgd->nvariables; n++)
      scd.variable_values << m_fsgd->varvals[i][n];
    scd.class_id = m_fsgd->subjclassno[i];
    m_data << scd;
  }

  // find variable value range
  for (int i = 0; i < m_data.size(); i++)
  {
    if (i == 0)
    {
      for (int n = 0; n < m_variables.size(); n++)
      {
        m_variables[n].range[0] = m_variables[n].range[1]
            = m_data[i].variable_values[n];
      }
    }
    else
    {
      for (int n = 0; n < m_variables.size(); n++)
      {
        if (m_data[i].variable_values[n] < m_variables[n].range[0])
          m_variables[n].range[0] = m_data[i].variable_values[n];
        else if (m_data[i].variable_values[n] > m_variables[n].range[1])
          m_variables[n].range[1] = m_data[i].variable_values[n];
      }
    }
    //  qDebug() << m_data[i].class_id << m_data[i].subject_id << m_data[i].variable_values;
  }

  m_title = m_fsgd->title;
  m_measureName = m_fsgd->measname;
  //  UpdateData(0);

  return true;
}

void FSGroupDescriptor::UpdateData(int nVertex)
{
  if (nVertex >= 0 && nVertex < m_fsgd->data->width)
  {
    m_nVertexNum = nVertex;
    m_dMeasurementRange[0] = 1e10;
    m_dMeasurementRange[1] = -1e10;
    for (int i = 0; i < m_data.size(); i++)
    {
      float val;
      gdfGetNthSubjectMeasurement(m_fsgd, i, nVertex, 0, 0, &val);
      m_data[i].measurement = val;
      if (m_dMeasurementRange[0] > val)
        m_dMeasurementRange[0] = val;
      else if (m_dMeasurementRange[1] < val)
        m_dMeasurementRange[1] = val;
    }
  }
}
