/**
 * @file  FSLabel.h
 * @brief Base label class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/01/30 20:57:05 $
 *    $Revision: 1.20 $
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


#include "FSLabel.h"
#include <stdexcept>
#include "vtkImageData.h"
#include "FSVolume.h"
#include <QFile>
#include <QTextStream>
#include <vector>

using namespace std;

FSLabel::FSLabel( QObject* parent ) : QObject( parent ),
  m_label( NULL ),
  m_bTkReg( true )
{}

FSLabel::~FSLabel()
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }
}

bool FSLabel::LabelRead( const QString& filename )
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  m_label = ::LabelRead( NULL, filename.toAscii().data() );
  if ( m_label == NULL )
  {
    cerr << "LabelRead failed\n";
    return false;
  }

  QFile file( filename );
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    cerr << qPrintable(file.errorString()) << "\n";
    return false;
  }

  QTextStream in(&file);
  while (!in.atEnd())
  {
    QString line = in.readLine();
    if ( line.contains( "vox2ras=" ) &&
         !line.contains( "vox2ras=TkReg" ) )
    {
      m_bTkReg = false;
      break;
    }
  }

  return true;
}

void FSLabel::UpdateLabelFromImage( vtkImageData* rasImage,
                                    FSVolume* ref_vol )
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }

  int nCount = 0;
  int* dim = rasImage->GetDimensions();
  double* orig = rasImage->GetOrigin();
  double* vs = rasImage->GetSpacing();
  vector<float> values;
  float fvalue;
// int nProgressStep = ( 90 - event.GetInt() ) / 5;
  double pos[3];

  // ok. the following part looks tedious,
  // but it's 100 times faster than doing i, j, k 3d iterations!
  int nsize = dim[0] * dim[1] * dim[2];
  switch ( rasImage->GetScalarType() )
  {
  case VTK_UNSIGNED_CHAR:
  {
    unsigned char* p = (unsigned char*)rasImage->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue != 0 )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
        ref_vol->NativeRASToTkReg( pos, pos );
        values.push_back( pos[0] );
        values.push_back( pos[1] );
        values.push_back( pos[2] );
        values.push_back( fvalue );
        nCount ++;
      }
    }
  }
  break;
  case VTK_SHORT:
  {
    short* p = (short*)rasImage->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue != 0 )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
        ref_vol->NativeRASToTkReg( pos, pos );
        values.push_back( pos[0] );
        values.push_back( pos[1] );
        values.push_back( pos[2] );
        values.push_back( fvalue );
        nCount ++;
      }
    }
  }
  break;
  case VTK_FLOAT:
  {
    float* p = (float*)rasImage->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue != 0 )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
        ref_vol->NativeRASToTkReg( pos, pos );
        values.push_back( pos[0] );
        values.push_back( pos[1] );
        values.push_back( pos[2] );
        values.push_back( fvalue );
        nCount ++;
      }
    }
  }
  break;
  case VTK_LONG:
  {
    long* p = (long*)rasImage->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue != 0 )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
        ref_vol->NativeRASToTkReg( pos, pos );
        values.push_back( pos[0] );
        values.push_back( pos[1] );
        values.push_back( pos[2] );
        values.push_back( fvalue );
        nCount ++;
      }
    }
  }
  break;
  case VTK_INT:
  {
    int* p = (int*)rasImage->GetScalarPointer();
    for ( int i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue != 0 )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
        ref_vol->NativeRASToTkReg( pos, pos );
        values.push_back( pos[0] );
        values.push_back( pos[1] );
        values.push_back( pos[2] );
        values.push_back( fvalue );
        nCount ++;
      }
    }
  }
  break;
  }

  m_label = ::LabelAlloc( nCount, NULL, (char*)"" );
  m_label->n_points = nCount;
  for ( int i = 0; i < nCount; i++ )
  {
    m_label->lv[i].x = values[i*4];
    m_label->lv[i].y = values[i*4+1];
    m_label->lv[i].z = values[i*4+2];
    m_label->lv[i].vno = -1;
    m_label->lv[i].deleted = false;
    m_label->lv[i].stat = 1;
  }
}

void FSLabel::UpdateRASImage( vtkImageData* rasImage, FSVolume* ref_vol )
{
  if ( !m_label )
  {
    cerr << "Label is empty\n";
    return;
  }

  int n[3];
  double pos[3];
  int* dim = rasImage->GetDimensions();
  memset( rasImage->GetScalarPointer(),
          0,
          dim[0] * dim[1] * dim[2] * rasImage->GetScalarSize() );
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    pos[0] = m_label->lv[i].x;
    pos[1] = m_label->lv[i].y;
    pos[2] = m_label->lv[i].z;
    if ( m_bTkReg )
    {
      ref_vol->TkRegToNativeRAS( pos, pos );
    }
    ref_vol->NativeRASToRAS( pos, pos );
    ref_vol->RASToTargetIndex( pos, n );
    if ( n[0] >= 0 && n[0] < dim[0] && n[1] >= 0 && n[1] < dim[1] &&
         n[2] >= 0 && n[2] < dim[2] )
    {
      rasImage->SetScalarComponentFromFloat( n[0], n[1], n[2], 0, 1 );
    }
  }
}

bool FSLabel::LabelWrite( const QString& filename )
{
  int err = ::LabelWrite( m_label, filename.toAscii().data() );

  if ( err != 0 )
  {
    cerr << "LabelWrite failed\n";
  }

  return err == 0;
}
