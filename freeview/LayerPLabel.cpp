/**
 * @brief Layer class for p-label volumes.
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

#include "LayerPLabel.h"
#include "MainWindow.h"
#include "LayerPropertyMRI.h"
#include "MyUtils.h"
#include "MyVTKUtils.h"
#include "FSVolume.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkImageMapToColors.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "LUTDataHolder.h"
#include <vtkImageActor.h>
#include <vtkImageReslice.h>
#include <vtkAbstractTransform.h>
#include "vtkImageMapper3D.h"
#include <QFileInfo>
#include <QDebug>

LayerPLabel::LayerPLabel( LayerMRI* ref, QObject* parent ) : LayerMRI( ref, parent ),
  m_volumeTemp( NULL)
{
  m_strTypeNames.push_back( "PLabel" );

  SetEditable( false );
}

LayerPLabel::~LayerPLabel()
{
  if ( m_volumeTemp )
  {
    delete m_volumeTemp;
  }
}

bool LayerPLabel::LoadVolumeFiles()
{
  if ( m_sFilenames.isEmpty() )
  {
    return false;
  }

  m_sFilename = m_sFilenames[0];
  if ( !LayerMRI::LoadVolumeFromFile() )
  {
    return false;
  }

  if ( m_volumeTemp )
  {
    delete m_volumeTemp;
  }
  m_volumeTemp = NULL;

  LUTDataHolder* luts = MainWindow::GetMainWindow()->GetLUTData();
  COLOR_TABLE* ct = luts->GetColorTable( m_sLUT );
  if ( !ct )
  {
    cerr << "Can not find look up table.\n";
    return false;
  }

  m_imageData = vtkSmartPointer<vtkImageData>::New();
  m_imageIndex = vtkSmartPointer<vtkImageData>::New();
  for ( int i = 0; i < m_sFilenames.size(); i++ )
  {
    QString fn = QFileInfo( m_sFilenames[i] ).completeBaseName();
    if ( !m_sFilenamePrefix.isEmpty() )
    {
      fn = fn.mid( m_sFilenamePrefix.length() );
      SetName( m_sFilenamePrefix );
    }

    int r, g, b;
    if ( CTABrgbAtIndexi( ct, CTABentryNameToIndex( fn.toLatin1().data(), ct ), &r, &g, &b ) != 0 &&
         CTABrgbAtIndexi( ct, CTABentryNameToIndex( QString(fn).replace("-", "/").toLatin1().data(), ct ), &r, &g, &b ) != 0 &&
         CTABrgbAtIndexi( ct, CTABentryNameToIndex( QString(fn).replace("-", "_").toLatin1().data(), ct ), &r, &g, &b ) != 0)
    {
      cerr << "Can not find index for color name " << qPrintable(fn) << "\n";
      return false;
    }

    int color[4] = { r, g, b, 255 };
    m_volumeTemp = new FSVolume( m_volumeRef );
    if ( !m_volumeTemp->MRIRead( m_sFilenames[i],
                                 m_sRegFilename.size() > 0 ? m_sRegFilename : NULL ) )
    {
      cerr << "Can not load volume file " << qPrintable(m_sFilenames[i]) << "\n";
      return false;
    }

    vtkImageData* imageData = m_volumeTemp->GetImageOutput();
    if ( i == 0 ) // initialize m_imageData
    {
      m_imageData->SetDimensions( imageData->GetDimensions() );
      m_imageData->SetOrigin( imageData->GetOrigin() );
      m_imageData->SetSpacing( imageData->GetSpacing() );
      m_imageData->SetExtent( imageData->GetExtent() );
      m_imageIndex->SetDimensions( imageData->GetDimensions() );
      m_imageIndex->SetOrigin( imageData->GetOrigin() );
      m_imageIndex->SetSpacing( imageData->GetSpacing() );
      m_imageIndex->SetExtent( imageData->GetExtent() );
#if VTK_MAJOR_VERSION > 5
      m_imageData->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
      m_imageIndex->AllocateScalars(VTK_FLOAT, 2);
#else
      m_imageData->SetNumberOfScalarComponents(4);
      m_imageData->SetScalarTypeToUnsignedChar();
      m_imageData->AllocateScalars();
      m_imageIndex->SetNumberOfScalarComponents(2);
      m_imageData->SetScalarTypeToFloat();
      m_imageData->AllocateScalars();
#endif
    }

    int* dim = m_imageData->GetDimensions();
    char* ptr = (char*)m_imageData->GetScalarPointer();
    int scalar_type = m_imageData->GetScalarType();
    int n_frames = m_imageData->GetNumberOfScalarComponents();

    int* index_dim = m_imageIndex->GetDimensions();
    char* index_ptr = (char*)m_imageIndex->GetScalarPointer();
    int index_scalar_type = m_imageIndex->GetScalarType();
    int index_n_frames = m_imageIndex->GetNumberOfScalarComponents();

    int* img_dim = imageData->GetDimensions();
    char* img_ptr = (char*)imageData->GetScalarPointer();
    int img_scalar_type = imageData->GetScalarType();
    int img_n_frames = imageData->GetNumberOfScalarComponents();
    for ( int ni = 0; ni < dim[0]; ni++ )
    {
      for ( int nj = 0; nj < dim[1]; nj++ )
      {
        for ( int nk = 0; nk < dim[2]; nk++ )
        {
          float pvalue = MyVTKUtils::GetImageDataComponent(img_ptr, img_dim, img_n_frames, ni, nj, nk, 0, img_scalar_type );
          if ( i == 0)
          {
            MyVTKUtils::SetImageDataComponent(index_ptr, index_dim, index_n_frames, ni, nj, nk, 0, index_scalar_type, 0);
            MyVTKUtils::SetImageDataComponent(index_ptr, index_dim, index_n_frames, ni, nj, nk, 1, index_scalar_type, pvalue);
          }
          else
          {
            float old_pvalue = MyVTKUtils::GetImageDataComponent(index_ptr, index_dim, index_n_frames, ni, nj, nk, 0, index_scalar_type);
            if ( old_pvalue < pvalue)
            {
              MyVTKUtils::SetImageDataComponent(index_ptr, index_dim, index_n_frames, ni, nj, nk, 0, index_scalar_type, i);
              MyVTKUtils::SetImageDataComponent(index_ptr, index_dim, index_n_frames, ni, nj, nk, 1, index_scalar_type, pvalue);
            }
          }
          for ( int m = 0; m < 4; m++ )
          {
            float fvalue = 0;
            if ( i != 0 )
            {
              fvalue = MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, ni, nj, nk, m, scalar_type );
            }

            fvalue += color[m]*pvalue/255;
            if ( fvalue > 255 )
            {
              fvalue = 255;
            }
            MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, ni, nj, nk, m, scalar_type, fvalue );
          }
        }
      }
    }

    delete m_volumeTemp;
    m_volumeTemp = NULL;
  }

  InitializeActors();
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->GetMapper()->SetInputConnection( mReslice[i]->GetOutputPort() );
    m_sliceActor3D[i]->GetMapper()->SetInputConnection( mReslice[i]->GetOutputPort() );
  }
  return true;
}

void LayerPLabel::UpdateColorMap()
{
  // over-ride parent class, do nothing
}

double LayerPLabel::GetVoxelValue(double* pos)
{
  if ( m_imageIndex == NULL )
  {
    return 0;
  }

  vtkAbstractTransform* tr = mReslice[0]->GetResliceTransform();
  double pos_new[3];
  tr->TransformPoint( pos, pos_new );

  double* orig = m_imageData->GetOrigin();
  double* vsize = m_imageData->GetSpacing();
  int* ext = m_imageData->GetExtent();

  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos_new[i] - orig[i] ) / vsize[i] + 0.5 );
  }

  if ( n[0] < ext[0] || n[0] > ext[1] ||
       n[1] < ext[2] || n[1] > ext[3] ||
       n[2] < ext[4] || n[2] > ext[5] )
  {
    return 0;
  }
  else
  {
    return m_imageIndex->GetScalarComponentAsDouble( n[0], n[1], n[2], 1 );
  }
}

QString LayerPLabel::GetLabelName(double *pos)
{
  if ( m_imageIndex == NULL )
  {
    return "";
  }

  vtkAbstractTransform* tr = mReslice[0]->GetResliceTransform();
  double pos_new[3];
  tr->TransformPoint( pos, pos_new );

  double* orig = m_imageIndex->GetOrigin();
  double* vsize = m_imageIndex->GetSpacing();
  int* ext = m_imageIndex->GetExtent();

  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos_new[i] - orig[i] ) / vsize[i] + 0.5 );
  }

  int nIndex = 0;
  if ( n[0] < ext[0] || n[0] > ext[1] ||
       n[1] < ext[2] || n[1] > ext[3] ||
       n[2] < ext[4] || n[2] > ext[5] )
  {
    return "";
  }
  else
  {
    nIndex = (int)m_imageIndex->GetScalarComponentAsDouble( n[0], n[1], n[2], 0 );
    QString strg = m_sFilenames[nIndex].mid( m_sFilenamePrefix.length() );
    int n = strg.lastIndexOf('.');
    if (n >= 0 )
    {
      strg = strg.left(n);
    }
    return strg;
  }
}
