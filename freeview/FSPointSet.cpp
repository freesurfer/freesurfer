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


#include "FSPointSet.h"
#include "FSVolume.h"
#include "vtkImageData.h"
#include "MyUtils.h"
#include <QFile>
#include <QTextStream>
#include <QDebug>

FSPointSet::FSPointSet( QObject* parent ) : QObject( parent ),
  m_label( NULL )
{}

FSPointSet::~FSPointSet()
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }
}

bool FSPointSet::IsLabelFormat( const QString& filename )
{
  QFile file(filename);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    return false;

  return file.readAll().contains("ascii label");
}

bool FSPointSet::ReadAsLabel( const QString& filename )
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  m_label = ::LabelRead( NULL, filename.toLatin1().data() );

  if ( m_label == NULL )
  {
    cerr << "LabelRead failed\n";
    return false;
  }

  UpdateStatRange();

  return true;
}

bool FSPointSet::ReadAsControlPoints( const QString& filename )
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  QFile file( filename );
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    cerr << qPrintable(file.errorString()) << "\n";;
    return false;
  }

  QTextStream in(&file);
  QString content = in.readAll();

  return ReadFromStringAsControlPoints(content);
}

bool FSPointSet::ReadFromStringAsControlPoints(const QString &content)
{
  QStringList ar = content.split("\n");
  int nCount = 0;
  QList<float> values;
  bool bRealRAS = true;
  for ( int i = 0; i < ar.size(); i++ )
  {
    QStringList subs = ar[i].split(" ", QString::SkipEmptyParts );
    if ( subs.size() > 2 )
    {
      for ( int j = 0; j < 3; j++ )
      {
        values.push_back( subs[j].toDouble() );
      }
      nCount ++;
    }
    else if (subs.size() > 1 &&
             subs[0].toLower() == "userealras" && subs[1] == "0")
      bRealRAS = false;
  }

  m_label = ::LabelAlloc( nCount, NULL, (char*)"" );
  m_label->n_points = nCount;
  if (bRealRAS)
  {
    m_label->coords = LABEL_COORDS_SCANNER_RAS;
    strncpy(m_label->space, "scanner", sizeof(m_label->space));
  }
  else
  {
    m_label->coords = LABEL_COORDS_TKREG_RAS;
    strncpy(m_label->space, "TkReg", sizeof(m_label->space));
  }
  for ( int i = 0; i < nCount; i++ )
  {
    m_label->lv[i].x = values[i*3];
    m_label->lv[i].y = values[i*3+1];
    m_label->lv[i].z = values[i*3+2];
    m_label->lv[i].vno = -1;
    m_label->lv[i].deleted = false;
    m_label->lv[i].stat = 0;
  }

  UpdateStatRange();

  return true;
}

bool FSPointSet::WriteAsLabel( const QString& filename )
{
  int err = ::LabelWrite( m_label, filename.toLatin1().data() );

  if ( err != 0 )
  {
    cerr << "Way Points Write failed\n";
  }

  return err == 0;
}

bool FSPointSet::WriteAsControlPoints( const QString& filename )
{
  QFile file( filename );
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  {
    QString strg = file.errorString();
    if (strg.isEmpty())
      cerr << "Can not open file for writing\n";
    else
      cerr << qPrintable(strg) << "\n";
    return false;
  }

  QTextStream out(&file);
  out << WriteAsControlPointsToString();

  return true;
}

QString FSPointSet::WriteAsControlPointsToString()
{
  QString strg;
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    strg += QString("%1 %2 %3\n").arg(m_label->lv[i].x).arg(m_label->lv[i].y).arg(m_label->lv[i].z);
  }
  strg += QString("info\nnumpoints %1\nuseRealRAS 1\n").arg( m_label->n_points );
  return strg;
}

void FSPointSet::UpdateLabel( PointSet& points_in, FSVolume* ref_vol )
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  int nCount = 0;
  double pos[3];
  QList<float> values;
  for ( int i = 0; i < points_in.size(); i++ )
  {
    ref_vol->TargetToRAS( points_in[i].pt, pos );
    ref_vol->RASToNativeRAS( pos, pos );
    values.push_back( pos[0] );
    values.push_back( pos[1] );
    values.push_back( pos[2] );
    values.push_back( points_in[i].value );
    nCount ++;
  }

  m_label = ::LabelAlloc( nCount, NULL, (char*)"" );
  ::LabelInit(m_label, ref_vol->GetMRI(), NULL, 0);
  ::LabelToScannerRAS(m_label, ref_vol->GetMRI(), m_label);

  m_label->n_points = nCount;
  for ( int i = 0; i < nCount; i++ )
  {
    m_label->lv[i].x = values[i*4];
    m_label->lv[i].y = values[i*4+1];
    m_label->lv[i].z = values[i*4+2];
    m_label->lv[i].vno = -1;
    m_label->lv[i].deleted = false;
    m_label->lv[i].stat = values[i*4+3];
  }
}

void FSPointSet::LabelToPointSet( PointSet& points_out, FSVolume* ref_vol )
{
  if ( !m_label )
  {
    //  cerr << "Label is empty\n";
    return;
  }

  ControlPoint wp;
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    wp.pt[0] = m_label->lv[i].x;
    wp.pt[1] = m_label->lv[i].y;
    wp.pt[2] = m_label->lv[i].z;
    wp.value = m_label->lv[i].stat;
    if ( m_label->coords == LABEL_COORDS_TKREG_RAS )
    {
      ref_vol->TkRegToNativeRAS( wp.pt, wp.pt );
    }
    else if (m_label->coords == LABEL_COORDS_VOXEL)
    {
      MRIvoxelToWorld(ref_vol->GetMRI(), wp.pt[0], wp.pt[1], wp.pt[2], wp.pt, wp.pt+1, wp.pt+2);
    }
    ref_vol->NativeRASToRAS( wp.pt, wp.pt );
    ref_vol->RASToTarget( wp.pt, wp.pt );
    points_out.push_back( wp );
  }
}

bool FSPointSet::GetCentroidRASPosition(double* pos, FSVolume* ref_vol)
{
  if (m_label && m_label->n_points > 0)
  {
    double x = 0, y = 0, z = 0;
    for ( int i = 0; i < m_label->n_points; i++ )
    {
      x += m_label->lv[i].x;
      y += m_label->lv[i].y;
      z += m_label->lv[i].z;
    }
    pos[0] = x / m_label->n_points;
    pos[1] = y / m_label->n_points;
    pos[2] = z / m_label->n_points;

    if ( m_label->coords == LABEL_COORDS_VOXEL )
    {
      MRIvoxelToWorld(ref_vol->GetMRI(), pos[0], pos[1], pos[2], pos, pos+1, pos+2);
    }
    else if (m_label->coords == LABEL_COORDS_TKREG_RAS)
    {
      ref_vol->TkRegToNativeRAS( pos, pos );
    }
    ref_vol->NativeRASToRAS( pos, pos );
    return true;
  }
  else
    return false;
}

void FSPointSet::UpdateStatRange()
{
  m_dStatMax = -1e10;
  m_dStatMin = 1e10;
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    if (m_label->lv[i].stat > m_dStatMax)
      m_dStatMax = m_label->lv[i].stat;
    if (m_label->lv[i].stat < m_dStatMin)
      m_dStatMin = m_label->lv[i].stat;
  }
  if (m_dStatMax <= m_dStatMin)
    m_dStatMax = m_dStatMin+1;
}
