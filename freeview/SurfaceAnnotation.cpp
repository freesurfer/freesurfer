/**
 * @file  SurfaceAnnotation.cxx
 * @brief Implementation for surface annotation.
 *
 * In 2D, the MRI is viewed as a single slice, and controls are
 * provided to change the color table and other viewing options. In
 * 3D, the MRI is viewed in three planes in 3D space, with controls to
 * move each plane axially.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/17 20:36:46 $
 *    $Revision: 1.25 $
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


#include "SurfaceAnnotation.h"
#include "vtkLookupTable.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkMath.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include <QFileInfo>
#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <QDebug>

SurfaceAnnotation::SurfaceAnnotation ( LayerSurface* surf ) :
  QObject( surf ),
  m_nIndices( NULL ),
  m_nOutlineIndices(NULL),
  m_nCenterVertices( NULL ),
  m_lut( NULL ),
  m_surface( surf ),
  m_bShowOutline(false),
  m_dOpacity(1.0)
{
}

SurfaceAnnotation::~SurfaceAnnotation ()
{
  Reset();
}

void SurfaceAnnotation::Reset()
{
  if ( m_nIndices )
  {
    delete[] m_nIndices;
  }

  if ( m_nCenterVertices )
  {
    delete[] m_nCenterVertices;
  }

  if ( m_nOutlineIndices )
  {
    delete[] m_nOutlineIndices;
  }

  m_nIndices = NULL;
  m_nOutlineIndices = NULL;
  m_lut = NULL;
  m_nCenterVertices = NULL;
}


bool SurfaceAnnotation::LoadAnnotation( const QString& fn )
{
  if ( m_surface )
  {
    m_strFilename = QFileInfo(fn).canonicalFilePath();
    MRIS* mris = m_surface->GetSourceSurface()->GetMRIS();
    m_nIndices = NULL;

    int ret;
    try {
      ret = MRISreadAnnotation( mris, fn.toAscii().data() );
    }
    catch (int return_code)
    {
      ret = return_code;
    }
    if ( ret != 0 )
    {
      cerr << "Could not load annotation from file " << qPrintable(fn) << ".\n";
      return false;
    }
    else
    {
      Reset();
      m_lut = mris->ct;
      m_nIndexSize = mris->nvertices;
      m_nIndices = new int[m_nIndexSize];

      // find valid annotations
      int n;
      QList<int> annotIndices;
      for ( int i = 0; i < m_nIndexSize; i++ )
      {
        if ( CTABfindAnnotation( m_lut, mris->vertices[i].annotation, &n ) == 0 &&
            n >= 0)
        {
          m_nIndices[i] = n;
          if (!annotIndices.contains(n))
            annotIndices << n;
        }
        else
          m_nIndices[i] = -1;
      }
      qSort(annotIndices);
      m_nAnnotations = annotIndices.size();
   //   CTABgetNumberOfValidEntries( m_lut, &m_nAnnotations );
      m_nCenterVertices = new int[m_nAnnotations];
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        m_nCenterVertices[i] = -1;
      }

      // convert annotations to lookup table indices
      // also find the "center vertex"
      double** pts = new double*[m_nAnnotations];
      int* vcount = new int[m_nAnnotations];
      memset( vcount, 0, sizeof( int )*m_nAnnotations );
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        pts[i] = new double[3];
        memset( pts[i], 0, sizeof( double ) * 3 );
      }

      for ( int i = 0; i < m_nIndexSize; i++ )
      {
        if ( m_nIndices[i] != -1 && mris->vertices[i].ripflag == 0 )
        {
          n = annotIndices.indexOf(m_nIndices[i]);
          vcount[n]++;
          pts[n][0] += mris->vertices[i].x;
          pts[n][1] += mris->vertices[i].y;
          pts[n][2] += mris->vertices[i].z;
          m_nCenterVertices[n] = i;
        }
      }

      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        if ( vcount[i] > 0 )
        {
          pts[i][0] /= vcount[i];
          pts[i][1] /= vcount[i];
          pts[i][2] /= vcount[i];

          int nVertex = m_surface->GetSourceSurface()->FindVertexAtSurfaceRAS( pts[i], NULL );
          if ( nVertex >=0 && m_nIndices[nVertex] == i )
          {
            m_nCenterVertices[i] = nVertex;
          }
        }
      }

      delete[] vcount;
      for ( int i = 0; i < m_nAnnotations; i++ )
      {
        delete[] pts[i];
      }
      delete[] pts;

      // build outline indices
      m_nOutlineIndices = new int[m_nIndexSize];
      memcpy(m_nOutlineIndices, m_nIndices, sizeof(int)*m_nIndexSize);
      for (int i = 0; i < m_nAnnotations; i++)
      {
          VERTEX *v;
          MRISclearMarks(mris);
          LABEL* label = MRISannotation_to_label(mris, annotIndices[i]);

          if (label)
          {
            LabelMarkSurface(label, mris);
            for (int n = 0 ; n < label->n_points ; n++)
            {
              if (label->lv[n].vno >= 0)
              {
                m_nOutlineIndices[label->lv[n].vno] = -1;
                v = &mris->vertices[label->lv[n].vno] ;
                if (v->ripflag)
                  continue;

                for (int m = 0 ; m < v->vnum ; m++)
                {
                  if (mris->vertices[v->v[m]].marked == 0)
                  {
                    m_nOutlineIndices[label->lv[n].vno] = m_nIndices[label->lv[n].vno];
                    break;
                  }
                }
              }
            }
            LabelFree(&label);
          }
      }
      return true;
    }
  }
  return false;
}

void SurfaceAnnotation::GetAnnotationPoint( int nIndex, double* pt_out )
{
  m_surface->GetTargetAtVertex( m_nCenterVertices[nIndex], pt_out );
}

QString SurfaceAnnotation::GetName()
{
  return m_strName;
}

void SurfaceAnnotation::SetName( const QString& name )
{
  m_strName = name;
}

int SurfaceAnnotation::GetIndexAtVertex( int nVertex )
{
  return m_nIndices[nVertex];
}

QString SurfaceAnnotation::GetAnnotationNameAtVertex( int nVertex )
{
  return GetAnnotationNameAtIndex( m_nIndices[nVertex] );
}

QString SurfaceAnnotation::GetAnnotationNameAtIndex( int nIndex )
{
  char name[128];
  int nValid = 0;
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex >= 0 && nIndex < nTotalCount )
  {
    CTABisEntryValid( m_lut, nIndex, &nValid );
  }
  if ( nValid && CTABcopyName( m_lut, nIndex, name, 128 ) == 0 )
  {
    return name;
  }

  return "";
}

void SurfaceAnnotation::GetAnnotationColorAtIndex( int nIndex, int* rgb )
{
  int nValid = 0;
  int nTotalCount = 0;
  CTABgetNumberOfTotalEntries( m_lut, &nTotalCount );
  if ( nIndex < nTotalCount )
  {
    CTABisEntryValid( m_lut, nIndex, &nValid );
  }
  if ( nValid )
  {
    CTABrgbAtIndexi( m_lut, nIndex, rgb, rgb+1, rgb+2 );
  }
}

void SurfaceAnnotation::SetShowOutline(bool bOutline)
{
  this->m_bShowOutline = bOutline;
}

void SurfaceAnnotation::MapAnnotationColor( unsigned char* colordata )
{
  int c[4];
  int* indices = (m_bShowOutline ? m_nOutlineIndices : m_nIndices);
  for ( int i = 0; i < m_nIndexSize; i++ )
  {
    if (indices[i] > 0)
    {
      if (CTABrgbAtIndexi( m_lut, indices[i], c, c+1, c+2 ) == 0) // no error
      {
        colordata[i*4] = ( int )( colordata[i*4] * ( 1 - m_dOpacity ) + c[0] * m_dOpacity );
        colordata[i*4+1] = ( int )( colordata[i*4+1] * ( 1 - m_dOpacity ) + c[1] * m_dOpacity );
        colordata[i*4+2] = ( int )( colordata[i*4+2] * ( 1 - m_dOpacity ) + c[2] * m_dOpacity );
      }
    }
  }
}
