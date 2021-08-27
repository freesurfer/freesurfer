/**
 * @brief Base label class that takes care of I/O and data conversion.
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


#include <stdexcept>
#include "vtkImageData.h"
#include "FSVolume.h"
#include "FSSurface.h"
#include "MyVTKUtils.h"
#include <QFile>
#include <QTextStream>
#include <vector>
#include <QDebug>
#include "FSLabel.h"
#include "vtkImageThreshold.h"

using namespace std;

FSLabel::FSLabel( QObject* parent, FSVolume* mri_template ) : QObject( parent ),
  m_label( NULL ), m_l2s( NULL ), m_mri_template(mri_template)
{
  m_dStatsRange[0] = 0;
  m_dStatsRange[1] = 1.0;
  m_label = ::LabelAlloc( 100, NULL, (char*)"" );
  if (mri_template)
  {
    ::LabelInit(m_label, mri_template->GetMRI(), NULL, CURRENT_VERTICES);
    LabelToScannerRAS(m_label, mri_template->GetMRI(), m_label);
  }
}

FSLabel::~FSLabel()
{
  if ( m_label )
  {
    ::LabelFree( &m_label );
  }
  foreach (LABEL* l, m_undoBuffer)
    ::LabelFree(&l);
  foreach (LABEL* l, m_redoBuffer)
    ::LabelFree(&l);
  if (m_l2s)
  {
    L2Sfree(&m_l2s);
  }
}

bool FSLabel::LabelRead( const QString& filename )
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
  if (m_label->coords != LABEL_COORDS_SCANNER_RAS && m_mri_template)
  {
    LabelToScannerRAS(m_label, m_mri_template->GetMRI(), m_label);
    cout << "Label coordinates are converted to scanner ras for " << qPrintable(filename) << endl;
  }

  if (m_label && m_label->n_points > 0)
  {
    double range[2];
    range[0] = range[1] = m_label->lv[0].stat;
    for (int i = 1; i < m_label->n_points; i++)
    {
      if (range[0] > m_label->lv[i].stat)
        range[0] = m_label->lv[i].stat;
      if (range[1] < m_label->lv[i].stat)
        range[1] = m_label->lv[i].stat;
    }

    m_dStatsRange[0] = range[0];
    m_dStatsRange[1] = range[1];
  }

  return true;
}

bool FSLabel::UpdateStatsRange(double val)
{
  if (m_dStatsRange[0] > val)
    m_dStatsRange[0] = val;
  else if (m_dStatsRange[1] < val)
    m_dStatsRange[1] = val;
  else
    return false;
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
  size_t nsize = ((size_t)dim[0]) * dim[1] * dim[2];
  switch ( rasImage->GetScalarType() )
  {
  case VTK_UNSIGNED_CHAR:
  {
    unsigned char* p = (unsigned char*)rasImage->GetScalarPointer();
    for ( size_t i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue >= m_dStatsRange[0] )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
//        ref_vol->NativeRASToTkReg( pos, pos );
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
    for ( size_t i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue >= m_dStatsRange[0] )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
//        ref_vol->NativeRASToTkReg( pos, pos );
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
    for ( size_t i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue >= m_dStatsRange[0] )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
//        ref_vol->NativeRASToTkReg( pos, pos );
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
    for ( size_t i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue >= m_dStatsRange[0] )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
//        ref_vol->NativeRASToTkReg( pos, pos );
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
    for ( size_t i = 0; i < nsize; i++ )
    {
      fvalue = p[i];
      if ( fvalue >= m_dStatsRange[0] )
      {
        pos[0] = (i%dim[0]) * vs[0] + orig[0];
        pos[1] = ( (i/dim[0])%dim[1] ) * vs[1] + orig[1];
        pos[2] = ( i/(dim[0]*dim[1]) ) * vs[2] + orig[2];
        ref_vol->TargetToRAS( pos, pos );
        ref_vol->RASToNativeRAS( pos, pos );
//        ref_vol->NativeRASToTkReg( pos, pos );
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
//  m_label->coords = LABEL_COORDS_TKREG_RAS;
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

void FSLabel::FillUnassignedVertices(FSSurface* surf, FSVolume* mri_template, int coords)
{
  LABEL* l = LabelSampleToSurface(surf->GetMRIS(), m_label, mri_template->GetMRI(), coords);
  ::LabelFree(&m_label);
  m_label = l;
  //    LabelFillUnassignedVertices(surf->GetMRIS(), m_label, coords);
}

void FSLabel::Initialize(FSVolume* ref_vol, FSSurface* surf, int coords)
{
  /*
  if (surf)
  {
    if (m_l2s)
      L2Sfree(&m_l2s);
    m_l2s = L2Salloc(1, "n/a");
    m_l2s->mri_template = ref_vol->GetMRI();
    m_l2s->surfs[0] = surf->GetMRIS();
    m_l2s->dmax = 3;
    m_l2s->hashres = 16;
    m_l2s->vol2surf = NULL;
    m_l2s->nhopsmax = 10;
    L2Sinit(m_l2s);
    L2SimportLabel(m_l2s, m_label, 1);
  }
  else
  */
    ::LabelInit(m_label, ref_vol->GetMRI(), surf?surf->GetMRIS():NULL, coords);
}

void FSLabel::UpdateRASImage( vtkImageData* rasImage, FSVolume* ref_vol, double threshold )
{
  if ( !m_label )
  {
    //   cerr << "Label is empty\n";
    return;
  }

  int n[3];
  double pos[3];
  int* dim = rasImage->GetDimensions();
  char* ptr = (char*)rasImage->GetScalarPointer();
  int scalar_type = rasImage->GetScalarType();
  memset( ptr,
          0,
          ((size_t)rasImage->GetScalarSize()) * dim[0] * dim[1] * dim[2]);
  if (m_dStatsRange[0] <= 0)
  {
    size_t nsize = ((size_t)dim[0])*dim[1]*dim[2];
    float* p = (float*)rasImage->GetScalarPointer();
    for (size_t i = 0; i < nsize; i++)
    {
      p[i] = m_dStatsRange[0]-1;
    }
  }

  QString error_strg;
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    if (!m_label->lv[i].deleted && (m_label->lv[i].stat >= threshold || m_dStatsRange[0] <= 0))
    {
      pos[0] = m_label->lv[i].x;
      pos[1] = m_label->lv[i].y;
      pos[2] = m_label->lv[i].z;
      if ( m_label->coords == LABEL_COORDS_VOXEL )
      {
        MRIvoxelToWorld(ref_vol->GetMRI(), pos[0], pos[1], pos[2], pos, pos+1, pos+2);
      }
      else if (m_label->coords == LABEL_COORDS_TKREG_RAS)
      {
        ref_vol->TkRegToNativeRAS( pos, pos );
      }
      ref_vol->NativeRASToRAS( pos, pos );
      ref_vol->RASToTargetIndex( pos, n );
      if ( n[0] >= 0 && n[0] < dim[0] && n[1] >= 0 && n[1] < dim[1] &&
           n[2] >= 0 && n[2] < dim[2] )
      {
        MyVTKUtils::SetImageDataComponent(ptr, dim, 1, n[0], n[1], n[2], 0, scalar_type,
            m_label->lv[i].stat >= threshold ? m_label->lv[i].stat : (m_dStatsRange[0]-1.0) );
      }
      else
      {
        error_strg = "Label coordinate out of bound";
      }
    }
  }
  if (!error_strg.isEmpty())
    cerr << qPrintable(error_strg) << endl;

  rasImage->Modified();
}

bool FSLabel::LabelWrite( const QString& filename )
{
  int err = ::LabelWrite( m_label, filename.toLatin1().data() );

  if ( err != 0 )
  {
    cerr << "LabelWrite failed\n";
  }

  return err == 0;
}

bool FSLabel::GetCentroidRASPosition(double* pos, FSVolume* ref_vol)
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

void FSLabel::GetStatsRange(double *range)
{
  range[0] = m_dStatsRange[0];
  range[1] = m_dStatsRange[1];
}

void FSLabel::EditVoxel(int nx, int ny, int nz, int coords, bool bAdd, int* vertices, int* pnum)
{
  if (bAdd)
    ::LabelAddVoxel(m_label, nx, ny, nz, coords, vertices, pnum);
  else
    ::LabelDeleteVoxel(m_label, nx, ny, nz, vertices, pnum);
//  Q_UNUSED(vertices);
//  Q_UNUSED(pnum);
//  qDebug() << L2SaddPoint(m_l2s, nx, ny, nz, bAdd?1:0);
}

bool FSLabel::HasUndo()
{
  return !m_undoBuffer.isEmpty();
}

bool FSLabel::HasRedo()
{
  return !m_redoBuffer.isEmpty();
}

void FSLabel::Undo()
{
  if (!m_undoBuffer.isEmpty())
  {
    LABEL* l = m_undoBuffer.last();
    LABEL* l2 = ::LabelCopy(m_label, NULL);
    ::LabelCopy(l, m_label);
    ::LabelFree(&l);
    m_undoBuffer.removeLast();
    m_redoBuffer << l2;
  }
}

void FSLabel::Redo()
{
  if (!m_redoBuffer.isEmpty())
  {
    LABEL* l = m_redoBuffer.last();
    LABEL* l2 = ::LabelCopy(m_label, NULL);
    ::LabelCopy(l, m_label);
    ::LabelFree(&l);
    m_redoBuffer.removeLast();
    m_undoBuffer << l2;
  }
}

void FSLabel::Clear()
{
  LabelFree(&m_label);
  m_label = ::LabelAlloc( 100, NULL, (char*)"" );
  if (m_mri_template)
  {
    ::LabelInit(m_label, m_mri_template->GetMRI(), NULL, CURRENT_VERTICES);
    ::LabelToScannerRAS(m_label, m_mri_template->GetMRI(), m_label);
  }
}

void FSLabel::SaveForUndo()
{
  LABEL* l = ::LabelCopy(m_label, NULL);
  m_undoBuffer << l;

  // clear redo buffer
  foreach (LABEL* l, m_redoBuffer)
    ::LabelFree(&l);
  m_redoBuffer.clear();
}
