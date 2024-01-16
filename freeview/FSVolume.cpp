/**
 * @brief Base volume class that takes care of I/O and data conversion.
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

#include "FSVolume.h"
#include "MyUtils.h"
#include <stdexcept>
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkImageReslice.h"
#include "vtkMatrix4x4.h"
#include "vtkSmartPointer.h"
#include "vtkTransform.h"
#include "vtkImageChangeInformation.h"
#include "vtkMath.h"
#include "vtkTransform.h"
#include "MyVTKUtils.h"
#include "vtkMatrix4x4.h"
#include <QFileInfo>
#include <QFile>
#include <QTextStream>
#include <QDir>
#include <QDebug>
#include "ProgressCallback.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <QElapsedTimer>
#include "MigrationDefs.h"

#include "registerio.h"
#include "utils.h"
#include "macros.h"
#include "mrisegment.h"
#include "mri.h"
#include "mri2.h"


#define NUM_OF_HISTO_BINS 10000

using namespace std;

FSVolume::FSVolume( FSVolume* ref, QObject* parent ) : QObject( parent ),
  m_MRI( NULL ),
  m_MRITarget( NULL ),
  m_MRIRef( NULL ),
  m_MRIOrigTarget( NULL ),
  m_MRITemp( NULL ),
  m_matReg( NULL ),
  m_ctabEmbedded( NULL ),
  m_volumeRef( ref ),
  m_fMinValue( 0 ),
  m_fMaxValue( 1 ),
  m_bResampleToRAS( false ),
  m_bBoundsCacheDirty( true ),
  m_nInterpolationMethod( SAMPLE_NEAREST ),
  m_bConform( false ),
  m_bCrop( false ),
  m_bCropToOriginal( false ),
  m_histoCDF( NULL ),
  m_nHistoFrame(0),
  m_bValidHistogram(false),
  m_bSharedMRI(false),
  m_lta(NULL),
  m_bIgnoreHeader(false)
{
  m_imageData = NULL;
  if ( ref )
  {
    SetMRI( m_MRIRef, ref->m_MRI );
    m_bResampleToRAS = ref->m_bResampleToRAS;
  }
  strcpy( m_strOrientation, "RAS" );
  m_transform = vtkSmartPointer<vtkTransform>::New();
  m_transform->Identity();
  m_transform->PostMultiply();
}

FSVolume::~FSVolume()
{
  if ( m_MRI && !m_bSharedMRI )
  {
    ::MRIfree( &m_MRI );
  }

  if ( m_MRITarget )
  {
    ::MRIfree( &m_MRITarget );
  }

  if ( m_MRIRef )
  {
    ::MRIfree( &m_MRIRef );
  }

  if ( m_MRIOrigTarget )
  {
    ::MRIfree( &m_MRIOrigTarget );
  }

  if ( m_MRITemp )
  {
    ::MRIfree( &m_MRITemp );
  }

  if ( m_matReg )
  {
    ::MatrixFree( &m_matReg );
  }

  if ( m_ctabEmbedded )
  {
    ::CTABfree( &m_ctabEmbedded );
  }

  if (m_histoCDF)
  {
    ::HISTOfree(&m_histoCDF);
  }

  if (m_lta)
  {
    ::LTAfree(&m_lta);
  }
}

bool FSVolume::LoadMRI( const QString& filename, const QString& reg_filename )
{
//  if ( !reg_filename.isEmpty() && !m_MRIRef )
//  {
//    cerr << "Error: A target volume must be loaded first to apply registration matrix.\n";
//    return false;
//  }

  // save old header to release later so there is no gap where m_MRI becomes NULL during re-loading process
//  qDebug() << "Begin LoadMRI";
//  QElapsedTimer timer;
//  timer.start();
  MRI* tempMRI = m_MRI;
  try
  {
    m_MRI = ::MRIread( filename.toLatin1().data() );      // could be long process
  }
  catch (int ret)
  {
    return false;
  }

  if ( m_MRI == NULL )
  {
    cerr << "MRIread failed: Unable to read from " << qPrintable(filename) << "\n";
    if ( tempMRI )
    {
      m_MRI = tempMRI;
    }
    return false;
  }

  if (m_MRIRef && m_bIgnoreHeader)
    MRIcopyHeader( m_MRIRef, m_MRI );

  // if m_MRI successfully loaded, release old header.
  if ( tempMRI )
  {
    ::MRIfree( &tempMRI );
  }

  if ( m_MRI->ct != NULL )
  {
    m_ctabEmbedded = CTABdeepCopy( m_MRI->ct );
  }

  // update orientation index
  MRIdircosToOrientationString( m_MRI, m_strOrientation );

  /*
  if ( m_MRI->AutoAlign != NULL )
  {
    MATRIX* M = m_MRI->AutoAlign;
    cout << M->rptr[1][1] << " " << M->rptr[1][2] << " " << M->rptr[1][3] << " " << M->rptr[1][4] << endl
        << M->rptr[2][1] << " " << M->rptr[2][2] << " " << M->rptr[2][3] << " " << M->rptr[2][4] << endl
        << M->rptr[3][1] << " " << M->rptr[3][2] << " " << M->rptr[3][3] << " " << M->rptr[3][4] << endl
        << M->rptr[4][1] << " " << M->rptr[4][2] << " " << M->rptr[4][3] << " " << M->rptr[4][4] << endl;
  }
  */

  if (!m_MRIRef)
    SetMRI(m_MRIRef, m_MRI);

  // read registration matrix
  if ( !reg_filename.isEmpty() && !LoadRegistrationMatrix( reg_filename ) )
  {
    cerr << "Read registration failed\n";
    return false;
  }
//  qDebug() << timer.elapsed()/1000;
//  qDebug() << "begin MRIvalrange";
  MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );
  m_fMaxValueFull = m_fMaxValue;
//  qDebug() << timer.elapsed()/1000;
//  qDebug() << "begin UpdateHistoCDF";
  UpdateHistoCDF();
  if (m_bValidHistogram)
  {
    double val = GetHistoValueFromPercentile(0.95)*1.1; //+m_histoCDF->bin_size/2;
    if (m_fMaxValue > 10*val)
    {
      // abnormally high voxel value
      m_fMaxValue = val+m_fMaxValue/NUM_OF_HISTO_BINS/2;
      if (m_fMaxValue > 10*val)
      {
//        qDebug() << timer.elapsed()/1000;
        UpdateHistoCDF(0, m_fMaxValue, true);
        val = GetHistoValueFromPercentile(0.99)*1.1;
        m_fMaxValue = val+m_fMaxValue/NUM_OF_HISTO_BINS;
      }
    }
  }
//  qDebug() << timer.elapsed()/1000;
  return true;
}

bool FSVolume::MRIRead( const QString& filename, const QString& reg_filename )
{
#ifdef HAVE_OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1 ;
#endif
  int max_percent = 50;
  ::SetProgressCallback(ProgressCallback, 0, max_percent);
  if ( LoadMRI( filename, reg_filename ) )
  {
    this->CopyMatricesFromMRI();
    ::SetProgressCallback(ProgressCallback, max_percent, 100);

    if ( !this->MapMRIToImage() )
    {
      return false;
    }

    if ( m_volumeRef && m_volumeRef->m_MRIOrigTarget && !m_MRIOrigTarget )
    {
      m_MRIOrigTarget = CreateTargetMRI( m_MRI, m_volumeRef->m_MRIOrigTarget, false, m_bConform );
    }

    return true;
  }
  else
  {
    return false;
  }
}

bool FSVolume::CreateFromMRIData(MRI *mri)
{
  m_MRI = mri;
  m_bSharedMRI = true;
  MRIdircosToOrientationString( m_MRI, m_strOrientation );

  MRIvalRange( m_MRI, &m_fMinValue, &m_fMaxValue );
  if (m_fMaxValue == m_fMinValue)
    m_fMaxValue = m_fMinValue + 4;

  UpdateHistoCDF();

  this->CopyMatricesFromMRI();
  if ( !this->MapMRIToImage() )
  {
    return false;
  }

  if ( m_volumeRef && m_volumeRef->m_MRIOrigTarget && !m_MRIOrigTarget )
  {
    m_MRIOrigTarget = CreateTargetMRI( m_MRI, m_volumeRef->m_MRIOrigTarget, false, m_bConform );
  }
  return true;
}

bool FSVolume::Restore( const QString& filename, const QString& reg_filename )
{
  ::SetProgressCallback(ProgressCallback, 0, 100);
  if ( LoadMRI( filename, reg_filename ) )
  {
    // create m_MRITemp for save/rotate
    if ( m_MRITemp )
    {
      MRIfree( &m_MRITemp );
    }

    MRI* mri = MRIallocHeader( m_MRI->width,
                               m_MRI->height,
                               m_MRI->depth,
                               m_MRI->type,
                               m_MRI->nframes);
    MRIcopyHeader( m_MRI, mri );
    m_MRITemp = m_MRI;
    m_MRI = mri;

    return true;
  }
  else
  {
    cerr << "Restore failed.\n";
    return false;
  }
}

MATRIX* FSVolume::LoadRegistrationMatrix(const QString &filename, MRI *target, MRI *src, LTA** lta_out)
{
  QString ext = QFileInfo( filename ).suffix();
  MATRIX* matReg = NULL;
  if ( ext == "xfm" )  // MNI style
  {
    MATRIX* m = NULL;
    if ( regio_read_mincxfm( filename.toLatin1().data(), &m, NULL ) != 0 )
    {
      return NULL;
    }

    matReg = MRItkRegMtx( target, src, m );
    MatrixFree( &m );
  }
  else if ( ext == "mat" )  // fsl style
  {
    QFile file( filename );
    if ( !file.open( QIODevice::ReadOnly | QIODevice::Text ) )
    {
      cerr << qPrintable (file.errorString()) << "\n";;
      return NULL;
    }

    QTextStream in(&file);
    QString line = in.readLine();
    QStringList values;
    while ( !line.isNull() )
    {
      values += line.split( " ", MD_SkipEmptyParts );
      line = in.readLine();
    }
    if ( values.size() < 16 )
    {
      return NULL;
    }

    MATRIX* m = MatrixAlloc( 4, 4, MATRIX_REAL );
    for ( int i = 0; i < 16; i++ )
    {
      *MATRIX_RELT(m, (i/4)+1, (i%4)+1) = values[i].toDouble();
    }
    matReg = MRIfsl2TkReg( target, src, m );
    MatrixFree( &m );
  }
  else if ( ext == "dat" )  // tkregister style
  {
    char* subject = NULL;
    float inplaneres, betplaneres, intensity;
    int float2int;
    if ( regio_read_register( filename.toLatin1().data(),
                              &subject,
                              &inplaneres,
                              &betplaneres,
                              &intensity,
                              &matReg,
                              &float2int ) != 0 )
    {
      return NULL;
    }

    free( subject );
  }
  else  // LTA style & all possible other styles
  {
    TRANSFORM* FSXform = TransformRead( (char*)filename.toLatin1().data() );
    if ( FSXform == NULL )
    {
      return NULL;
    }
    LTA* lta = (LTA*) FSXform->xform;
    if ( lta->type != LINEAR_RAS_TO_RAS )
    {
      cout << "INFO: LTA input is not RAS to RAS...converting...\n";
      lta = LTAchangeType( lta, LINEAR_RAS_TO_RAS );
    }
    if ( lta->type != LINEAR_RAS_TO_RAS )
    {
      cerr << "ERROR: LTA input is not RAS to RAS\n";
      TransformFree( &FSXform );
      return NULL;
    }

    // Assume RAS2RAS and uses vox2ras from input volumes:
    // Note: This ignores the volume geometry in the LTA file.
    matReg = MRItkRegMtx( target, src, lta->xforms[0].m_L );
    if (lta_out)
    {
      *lta_out = LTAcopy(lta, NULL);
    }
    TransformFree( &FSXform );
  }

  return matReg;
}

// read in registration file and convert it to tkreg style
bool FSVolume::LoadRegistrationMatrix( const QString& filename )
{
  ClearRegistrationMatrix();
  m_matReg = LoadRegistrationMatrix(filename, m_MRIRef, m_MRI, &m_lta);
  return (m_matReg != NULL);
}

void FSVolume::ClearRegistrationMatrix()
{
  if ( m_matReg )
  {
    ::MatrixFree( &m_matReg );
  }
}

bool FSVolume::Create( FSVolume* src_vol, bool bCopyVoxelData, int data_type )
{
  if ( m_MRI )
  {
    ::MRIfree( &m_MRI );
  }
  if ( m_matReg )
  {
    ::MatrixFree( &m_matReg );
  }

  SetMRI( m_MRIRef, src_vol->GetMRI() );
  memcpy( m_RASToVoxelMatrix,
          src_vol->m_RASToVoxelMatrix,
          16 * sizeof( double ) );
  memcpy( m_VoxelToRASMatrix,
          src_vol->m_VoxelToRASMatrix,
          16 * sizeof( double ) );
  memcpy( m_MRIToImageMatrix,
          src_vol->m_MRIToImageMatrix,
          16 * sizeof( double ) );
  memcpy( m_VoxelToVoxelMatrix,
          src_vol->m_VoxelToVoxelMatrix,
          16 * sizeof( double ) );
  memcpy( m_RASToRASMatrix,
          src_vol->m_RASToRASMatrix,
          16 * sizeof( double ) );
  memcpy( m_RASToTkRegMatrix,
          src_vol->m_RASToTkRegMatrix,
          16 * sizeof( double ) );

  m_matReg = MatrixCopy( src_vol->m_matReg, NULL );

  m_fMinValue = src_vol->m_fMinValue;
  m_fMaxValue = src_vol->m_fMaxValue;
  m_bResampleToRAS = src_vol->m_bResampleToRAS;

  // SetOriginalOrigin( src_vol->m_fOriginalOrigin );

  if ( data_type == -1 )
  {
    data_type = src_vol->m_MRI->type;
  }

  try {
    m_MRI = MRIallocSequence( src_vol->m_MRI->width,
                              src_vol->m_MRI->height,
                              src_vol->m_MRI->depth,
                              data_type, 1 );
  }
  catch (int ret)
  {
    return false;
  }

  if ( NULL == m_MRI )
  {
    cerr << "Could not allocate new mri volume.\n";
    return false;
  }

  if ( bCopyVoxelData )
  {
    MRIcopy( src_vol->m_MRI, m_MRI );
  }

  SetMRITarget( src_vol->m_MRITarget );

  if ( src_vol->m_MRIOrigTarget )
  {
    MRI* mri = src_vol->m_MRIOrigTarget;
    m_MRIOrigTarget = MRIallocHeader( mri->width,
                                      mri->height,
                                      mri->depth,
                                      data_type,
                                      mri->nframes);
    MRIcopyHeader( mri, m_MRIOrigTarget );
  }

  // Copy the header from the template into the new mri.
  if ( !bCopyVoxelData )
  {
    MRIcopyHeader( src_vol->m_MRI, m_MRI );
  }

  if ( m_imageData == NULL )
  {
    m_imageData = vtkSmartPointer<vtkImageData>::New();
  }

  if ( bCopyVoxelData )
  {
    if ( data_type != src_vol->GetDataType() )
    {
      cerr << "Can not copy voxel data with different data type.\n";
    }
    else
    {
      m_imageData->DeepCopy( src_vol->m_imageData );
    }

    if ( !m_imageData->GetScalarPointer() )
    {
      cerr << "Unable to allocate voxel data.\n";
      return false;
    }
  }
  else
  {
    m_imageData->SetOrigin( src_vol->m_imageData->GetOrigin() );
    m_imageData->SetSpacing( src_vol->m_imageData->GetSpacing() );
    m_imageData->SetDimensions( src_vol->m_imageData->GetDimensions() );
#if VTK_MAJOR_VERSION > 5
    switch ( m_MRI->type )
    {
    case MRI_UCHAR:
      m_imageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
      break;
    case MRI_INT:
      m_imageData->AllocateScalars(VTK_INT, 1);
      break;
    case MRI_LONG:
      m_imageData->AllocateScalars(VTK_LONG, 1);
      break;
    case MRI_FLOAT:
      m_imageData->AllocateScalars(VTK_FLOAT, 1);
      break;
    case MRI_SHORT:
      m_imageData->AllocateScalars(VTK_SHORT, 1);
      break;
    case MRI_USHRT:
      m_imageData->AllocateScalars(VTK_UNSIGNED_SHORT, 1);
      break;
    default:
      break;
    }
#else
    m_imageData->SetNumberOfScalarComponents( 1 );
    switch ( m_MRI->type )
    {
    case MRI_UCHAR:
      m_imageData->SetScalarTypeToUnsignedChar();
      break;
    case MRI_INT:
      m_imageData->SetScalarTypeToInt();
      break;
    case MRI_LONG:
      m_imageData->SetScalarTypeToLong();
      break;
    case MRI_FLOAT:
      m_imageData->SetScalarTypeToFloat();
      break;
    case MRI_SHORT:
      m_imageData->SetScalarTypeToShort();
      break;
    case MRI_USHRT:
      m_imageData->SetScalarTypeToUnsignedShort();
      break;
    default:
      break;
    }
    m_imageData->AllocateScalars();
#endif

    char* ptr = ( char* )m_imageData->GetScalarPointer();
    int* nDim = m_imageData->GetDimensions();
    if ( !ptr )
    {
      cerr << "Unable to allocate voxel data.\n";
      return false;
    }
    if ( !bCopyVoxelData )
    {
      memset( ptr,
              0,
              ((size_t)m_imageData->GetScalarSize()) * nDim[0] * nDim[1] * nDim[2] );
    }
  }

  for ( int i = 0; i < 3; i++ )
  {
    m_strOrientation[i] = src_vol->m_strOrientation[i];
  }

  // Do not copy transform
  //  m_transform->DeepCopy( src_vol->m_transform );

  return true;
}

void FSVolume::SetMRI( MRI*& mri_out, MRI* mri_in )
{
  if ( mri_out )
  {
    ::MRIfree( &mri_out );
  }

  if ( mri_in )
  {
    mri_out = MRIallocHeader( mri_in->width,
                              mri_in->height,
                              mri_in->depth,
                              mri_in->type,
                              mri_in->nframes);
    MRIcopyHeader( mri_in, mri_out );
  }
  else
  {
    mri_out = NULL;
  }
}

void FSVolume::SetMRITarget( MRI* mri )
{
  if ( mri )
  {
    if ( m_MRITarget )
    {
      ::MRIfree( &m_MRITarget );
    }

    // must use source data type!
    m_MRITarget = MRIallocHeader( mri->width,
                                  mri->height,
                                  mri->depth,
                                  m_MRI->type,
                                  m_MRI->nframes);
    MRIcopyHeader( mri, m_MRITarget );
  }
}

vtkTransform* FSVolume::GetTransform()
{
  return m_transform;
}

MATRIX* FSVolume::GetTransformMatrixInRAS()
{
  vtkMatrix4x4* mat = m_transform->GetMatrix();
  MATRIX* m = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT(m, (i/4)+1, (i%4)+1) = mat->Element[i/4][i%4];
  }
  return m;
}

void FSVolume::ConvertTransformFromTargetToRAS(vtkMatrix4x4 *target_in, vtkMatrix4x4 *ras_out)
{
  MATRIX* m1 = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT(m1, (i/4)+1, (i%4)+1) = target_in->Element[i/4][i%4];
  }

  MATRIX* t2r = GetTargetToRASMatrix();
  MATRIX* r2t = MatrixInverse(t2r, NULL);
  MATRIX* m1_inverse = MatrixInverse(m1, NULL);
  MATRIX* m2 = MatrixMultiply(m1_inverse, r2t, NULL);
  MATRIX* m = MatrixMultiply(t2r, m2, NULL);
  for ( int i = 0; i < 16; i++ )
  {
    ras_out->Element[i/4][i%4] = *MATRIX_RELT(m, (i/4)+1, (i%4)+1);
  }
  MatrixFree(&r2t);
  MatrixFree(&t2r);
  MatrixFree(&m2);
  MatrixFree(&m1);
  MatrixFree(&m1_inverse);
  MatrixFree(&m);
}

bool FSVolume::SaveRegistration( const QString& filename )
{
  MATRIX* m = GetTransformMatrixInRAS();

  LINEAR_TRANSFORM *lt;
  VOL_GEOM srcG, dstG;
  LTA* lta;

  if (m_matReg)
  {
    MATRIX* reg = NULL;
    if (m_MRIRef)
    {
      LTA* temp = TransformRegDat2LTA(m_MRIRef, m_MRI, m_matReg);
      temp = LTAchangeType(temp, LINEAR_RAS_TO_RAS);
      reg = MatrixCopy(temp->xforms[0].m_L, NULL);
      LTAfree(&temp);
    }
    else
      reg = MatrixCopy(m_matReg, NULL);
    MATRIX* inv_reg = MatrixInverse(reg, NULL);
    MATRIX* md = MatrixMultiply(m, inv_reg, NULL);
    MatrixFree(&m);
    MatrixFree(&inv_reg);
    MatrixFree(&reg);
    m = md;
  }

  if (m_lta)
  {
    lta = LTAcopy(m_lta, NULL);
    lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
    lt = &lta->xforms[0];
    lt->m_L = m;
  }
  else
  {
    MATRIX* voxel_xform = MRIrasXformToVoxelXform( m_MRI, m_MRIRef, m, NULL );
    lta = LTAalloc(1, NULL);
    lt = &lta->xforms[0];
    lt->sigma = 1.0f ;
    lt->x0 = lt->y0 = lt->z0 = 0;
    lta->type = LINEAR_VOX_TO_VOX;
    lt->m_L = voxel_xform;
    MatrixFree(&m);
  }

  getVolGeom( m_MRI, &srcG );
  if (m_MRIRef)
    getVolGeom( m_MRIRef, &dstG );
  else
    getVolGeom( m_MRI, &dstG );
  lta->xforms[0].src = srcG;
  lta->xforms[0].dst = dstG;

  FILE* fp = fopen( filename.toLatin1().data(),"w" );
  bool ret = true;
  if( !fp )
  {
    cerr << "ERROR: cannot open for writing: " << qPrintable(filename) << "\n";
    ret = false;
  }
  else
  {
    LTAprint(fp, lta);
    fclose( fp );
  }

  LTAfree(&lta);

  return ret;
}

bool FSVolume::MRIWrite( const QString& filename, int nSampleMethod, bool resample )
{
  if ( !m_MRITemp )
  {
    cerr << "Volume not ready for save.\n";
    return false;
  }

  // check if transformation needed
  bool bTransformed = false;
  bool bRefTransformed = false;
  MATRIX* m = GetTransformMatrixInRAS();
  vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
  for ( int i = 0; i < 16; i++ )
  {
    mat->Element[i/4][i%4] = *MATRIX_RELT(m, (i/4)+1, (i%4)+1);
  }

  //  Do not call MatrixIsIdentity. It only checks the rotation part.
  //  if (MatrixIsIdentity(m))
  double scale = qMin(m_MRITemp->xsize, qMin(m_MRITemp->ysize, m_MRITemp->zsize));
  if (MyUtils::IsIdentity(mat->Element, scale))
  {
    if (m_volumeRef)
    {
      MatrixFree(&m);
      m = m_volumeRef->GetTransformMatrixInRAS();
      for ( int i = 0; i < 16; i++ )
      {
        mat->Element[i/4][i%4] = *MATRIX_RELT(m, (i/4)+1, (i%4)+1);
      }
      bRefTransformed = true;
    }
  }
  //  if (!MatrixIsIdentity(m))
  if ( !MyUtils::IsIdentity(mat->Element, scale) )
  {
    if ( resample ) // && MyUtils::IsOblique(mat->Element))
    {
      // find out the output voxel bounds
      float cornerFactor[3];
      double RAS[4] = {0, 0, 0, 1}, index[3];
      double indexBounds[6];
      indexBounds[0] = indexBounds[2] = indexBounds[4] = 1e10;
      indexBounds[1] = indexBounds[3] = indexBounds[5] = -1e10;
      for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
      {
        for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
        {
          for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
          {
            ::MRIvoxelToWorld( m_MRITemp,
                               cornerFactor[0]*(m_MRITemp->width-1),
                cornerFactor[1]*(m_MRITemp->height-1),
                cornerFactor[2]*(m_MRITemp->depth-1),
                &RAS[0], &RAS[1], &RAS[2] );
            mat->MultiplyPoint(RAS, RAS);
            ::MRIworldToVoxel( m_MRITemp,
                               RAS[0], RAS[1], RAS[2],
                &index[0], &index[1], &index[2] );

            if ( index[0] < indexBounds[0] )
            {
              indexBounds[0] = index[0];
            }
            if ( index[0] > indexBounds[1] )
            {
              indexBounds[1] = index[0];
            }
            if ( index[1] < indexBounds[2] )
            {
              indexBounds[2] = index[1];
            }
            if ( index[1] > indexBounds[3] )
            {
              indexBounds[3] = index[1];
            }
            if ( index[2] < indexBounds[4] )
            {
              indexBounds[4] = index[2];
            }
            if ( index[2] > indexBounds[5] )
            {
              indexBounds[5] = index[2];
            }
          }
        }
      }

      // calculate dimension of the converted target volume
      int dim[3];
      dim[0] = (int)( indexBounds[1] - indexBounds[0] + 1.0 );
      dim[1] = (int)( indexBounds[3] - indexBounds[2] + 1.0 );
      dim[2] = (int)( indexBounds[5] - indexBounds[4] + 1.0 );

      if (m_bCropToOriginal)
      {
        dim[0] = m_MRITemp->width;
        dim[1] = m_MRITemp->height;
        dim[2] = m_MRITemp->depth;
      }

      MRI* mri = NULL;
      try
      {
        mri = MRIallocSequence( dim[0], dim[1], dim[2], m_MRITemp->type, m_MRITemp->nframes );
      }
      catch (int ret)
      {
        return false;
      }

      if ( mri == NULL )
      {
        MatrixFree( &m );
        cerr << "Can not allocate mri.\n";
        return false;
      }

      MRIcopyHeader( m_MRITemp, mri );

      if (!m_bCropToOriginal)
      {
        double p0[3];
        ::MRIvoxelToWorld( m_MRITemp,
                           (int)indexBounds[0],
            (int)indexBounds[2],
            (int)indexBounds[4],
            &p0[0], &p0[1], &p0[2] );
        MRIp0ToCRAS( mri, p0[0], p0[1], p0[2] );
      }

      /*
      MATRIX* M = m;
      qDebug() << M->rptr[1][1] << " " << M->rptr[1][2] << " " << M->rptr[1][3] << " " << M->rptr[1][4]
          << M->rptr[2][1] << " " << M->rptr[2][2] << " " << M->rptr[2][3] << " " << M->rptr[2][4]
          << M->rptr[3][1] << " " << M->rptr[3][2] << " " << M->rptr[3][3] << " " << M->rptr[3][4]
          << M->rptr[4][1] << " " << M->rptr[4][2] << " " << M->rptr[4][3] << " " << M->rptr[4][4];
      */

      MATRIX* mi = MatrixIdentity(4, NULL);
      mri = MRIapplyRASlinearTransformInterp( m_MRITemp, mri, bRefTransformed?mi:m, nSampleMethod );
      MatrixFree(&mi);
      if ( !mri )
      {
        MatrixFree( &m );
        cerr << "MRIapplyRASlinearTransformInterp failed.\n";
        return false;
      }
      MRIfree( &m_MRITemp );
      m_MRITemp = mri;
    }
    else
    {
      // no resample, just modify the header
      MATRIX* old_v2r = extract_i_to_r(m_MRITemp);
      /*
      double p0[4] = {0, 0, 0, 1};
      ::MRIvoxelToWorld( m_MRITemp,
                         0, 0, 0,
                         &p0[0], &p0[1], &p0[2] );
      MATRIX* new_v2r = MatrixMultiply(m, old_v2r, 0);
      MRIsetVoxelToRasXform(m_MRITemp, new_v2r);
      mat->MultiplyPoint(p0, p0);
      MRIp0ToCRAS( m_MRITemp, p0[0], p0[1], p0[2] );
      */
      MATRIX* new_v2r = MatrixMultiply(m, old_v2r, 0);
      m_MRITemp->xsize = sqrt(new_v2r->rptr[1][1]*new_v2r->rptr[1][1] + new_v2r->rptr[2][1]*new_v2r->rptr[2][1] + new_v2r->rptr[3][1]*new_v2r->rptr[3][1]);
      m_MRITemp->ysize = sqrt(new_v2r->rptr[1][2]*new_v2r->rptr[1][2] + new_v2r->rptr[2][2]*new_v2r->rptr[2][2] + new_v2r->rptr[3][2]*new_v2r->rptr[3][2]);
      m_MRITemp->zsize = sqrt(new_v2r->rptr[1][3]*new_v2r->rptr[1][3] + new_v2r->rptr[2][3]*new_v2r->rptr[2][3] + new_v2r->rptr[3][3]*new_v2r->rptr[3][3]);

      MRIsetVoxelToRasXform(m_MRITemp, new_v2r);
      MatrixFree(&old_v2r);
      MatrixFree(&new_v2r);
    }
    MatrixFree( &m );
    bTransformed = true;

    // update ras2tkreg
    m = MRIgetRasToVoxelXform( m_MRITemp );
    for ( int i = 0; i < 16; i++ )
    {
      m_RASToVoxelMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
    }

    MATRIX* tkreg = MRIxfmCRS2XYZtkreg( m_MRITemp );
    MATRIX* m1 = MatrixMultiply( tkreg, m, NULL );
    for ( int i = 0; i < 16; i++ )
    {
      m_RASToTkRegMatrix[i] = (double) *MATRIX_RELT((m1),(i/4)+1,(i%4)+1);
    }
    MatrixFree( &m );
    MatrixFree( &tkreg );
    MatrixFree( &m1 );
  }

  // check if cropping is enabled
  // if so, calculate the extent to crop
  if ( m_bCrop )
  {
    int nmin[3] = { 10000000, 10000000, 1000000 };
    int nmax[3] = { -10000000, -10000000, -1000000 };
    int x = 0, y = 0, z = 0;
    double rx = 0, ry = 0, rz = 0;
    for ( int i = 0; i < 2; i++ )
    {
      for ( int j = 2; j < 4; j++ )
      {
        for ( int k = 4; k < 6; k++ )
        {
          this->TargetToRAS( m_dBounds[i], m_dBounds[j], m_dBounds[k], rx, ry, rz );
          this->RASToOriginalIndex( rx, ry, rz, x, y, z );
          if ( nmin[0] > x )
          {
            nmin[0] = x;
          }
          if ( nmin[1] > y )
          {
            nmin[1] = y;
          }
          if ( nmin[2] > z )
          {
            nmin[2] = z;
          }
          if ( nmax[0] < x )
          {
            nmax[0] = x;
          }
          if ( nmax[1] < y )
          {
            nmax[1] = y;
          }
          if ( nmax[2] < z )
          {
            nmax[2] = z;
          }
        }
      }
    }
    if ( nmin[0] < 0 )
    {
      nmin[0] = 0;
    }
    if ( nmin[1] < 0 )
    {
      nmin[1] = 0;
    }
    if ( nmin[2] < 0 )
    {
      nmin[2] = 0;
    }
    if ( nmax[0] >= m_MRITemp->width )
    {
      nmax[0] = m_MRITemp->width - 1;
    }
    if ( nmax[1] >= m_MRITemp->height )
    {
      nmax[1] = m_MRITemp->height - 1;
    }
    if ( nmax[2] >= m_MRITemp->depth )
    {
      nmax[2] = m_MRITemp->depth - 1;
    }

    int dx = nmax[0]-nmin[0]+1;
    int dy = nmax[1]-nmin[1]+1;
    int dz = nmax[2]-nmin[2]+1;
    if ( dx < 1 || dy < 1 || dz < 1 )
    {
      cerr << "Bad cropping range. \n";
      return false;
    }

    if ( nmin[0] != 0 || nmin[1] != 0 || nmin[2] != 0 ||
         dx != m_MRITemp->width || dy != m_MRITemp->height ||
         dz != m_MRITemp->depth )
    {
      MRI* mri = MRIextract( m_MRITemp, NULL, nmin[0], nmin[1], nmin[2],
          dx, dy, dz );
      if ( !mri )
      {
        cerr << "MRIextract failed.\n";
        return false;
      }
      MRIfree( &m_MRITemp );
      m_MRITemp = mri;
    }
  }

  // check if file is writable
  FILE* fp = fopen( filename.toLatin1().data(), "w" );
  if ( !fp )
  {
    cerr << "Failed to open file " << qPrintable(filename)
         << " for write. Please check if the directory exists and you have permission to write in that location.\n";
    MRIfree( &m_MRITemp );
    m_MRITemp = NULL;
    return false;
  }
  else
  {
    fclose( fp );
  }

  int err = 0;
  try
  {
    err = ::MRIwrite( m_MRITemp, filename.toLatin1().data() );
  }
  catch (int ret)
  {
    return false;
  }

  if ( err != 0 )
  {
    cerr << "MRIwrite failed\n";
  }
  else if ( bTransformed )    // save lta file only if volume was saved successfully
  {
    SaveRegistration( filename + ".lta" );
  }

  MRIfree( &m_MRITemp );
  m_MRITemp = NULL;

  return err == 0;
}

// if data_type < 0, use source data type
bool FSVolume::UpdateMRIFromImage( vtkImageData* rasImage, bool resampleToOriginal, int data_type )
{
  int nProgressStep = 5;

  MATRIX* vox2vox = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1) = m_VoxelToVoxelMatrix[i];
  }

  // create a target volume
  MRI* mri = NULL;
  try
  {
    mri = MRIallocSequence( m_MRITarget->width,
                            m_MRITarget->height,
                            m_MRITarget->depth,
                            data_type >= 0 ? data_type : m_MRITarget->type,
                            m_MRI->nframes );
  }
  catch (int ret)
  {
    return false;
  }

  if ( mri == NULL )
  {
    cout << "Can not allocate mri volume for buffering.\n";
    return false;
  }

  MRIcopyHeader( m_MRITarget, mri );

  int nProgress = 0;
  int nstart = global_progress_range[0];
  int nend = global_progress_range[1];
  char* ptr = (char*)rasImage->GetScalarPointer();
  int scalar_type = rasImage->GetScalarType();
  int* dim = rasImage->GetDimensions();
  int nNumberOfFrames = rasImage->GetNumberOfScalarComponents();
  int label_val = property("label_value").toInt();
  if ( true ) // mri->nframes > 1 )
  {
    global_progress_range[1] = nstart+(nend-nstart)*2/3;
    for ( int j = 0; j < mri->height; j++ )
    {
      for ( int k = 0; k < mri->depth; k++ )
      {
        for ( int i = 0; i < mri->width; i++ )
        {
          for ( int nFrame = 0; nFrame < mri->nframes; nFrame++ )
          {
            //            float val = rasImage->GetScalarComponentAsFloat(i, j, k, nFrame);
            float val = (float)MyVTKUtils::GetImageDataComponent(ptr, dim, nNumberOfFrames, i, j, k, nFrame, scalar_type);
            if (label_val > 0 && val != label_val)
              val = 0;
            switch ( mri->type )
            {
            case MRI_UCHAR:
              MRIseq_vox( mri, i, j, k, nFrame ) = (unsigned char)val;
              break;
            case MRI_INT:
              MRIIseq_vox( mri, i, j, k, nFrame ) = (int)val;
              break;
            case MRI_LONG:
              MRILseq_vox( mri, i, j, k, nFrame ) = (long)val;
              break;
            case MRI_FLOAT:
              MRIFseq_vox( mri, i, j, k, nFrame ) = val;
              break;
            case MRI_SHORT:
              MRISseq_vox( mri, i, j, k, nFrame ) = (short)val;
              break;
            case MRI_USHRT:
              MRIUSseq_vox( mri, i, j, k, nFrame ) = (unsigned short)val;
              break;
            default:
              break;
            }
          }
        }
      }
      if ( m_MRI->height >= 5 && j%(m_MRI->height/5) == 0 )
      {
        nProgress += nProgressStep;
        emit ProgressChanged( nProgress );
      }
      exec_progress_callback(j, mri->height, 0, 1);
    }
  }
  else
  {
    size_t bytes_per_slice = mri->bytes_per_vox * mri->height * mri->depth;
    global_progress_range[1] = nstart+(nend-nstart)/3;
    for ( int k = 0; k < mri->depth; k++ )
    {
      void* ptr = rasImage->GetScalarPointer( 0, 0, k );
      BUFTYPE* buf = &MRIseq_vox( mri, 0, 0, k, 0);
      memcpy( buf, ptr, bytes_per_slice );

      if ( mri->depth >= 5 && k%(mri->depth/5) == 0 )
      {
        nProgress += nProgressStep;
        emit ProgressChanged( nProgress );
      }
      exec_progress_callback(k, mri->depth, 0, 1);
    }
  }
  setProperty("label_value", 0);

  // create m_MRItemp for writing or rotation
  if ( m_MRITemp )
  {
    MRIfree( &m_MRITemp );
  }

  global_progress_range[0] = global_progress_range[1];
  global_progress_range[1] = nend;
  if ( resampleToOriginal )
  {
    try {
      m_MRITemp = MRIallocSequence( m_MRI->width,
                                    m_MRI->height,
                                    m_MRI->depth,
                                    data_type >= 0 ? data_type : m_MRI->type,
                                    m_MRI->nframes );
    } catch (int ret) {
      return false;
    }
    if ( m_MRITemp == NULL )
    {
      cout << "Can not allocate mri volume for buffering.\n";
      MRIfree( &mri );
      return false;
    }
    MRIcopyHeader( m_MRI, m_MRITemp );
    MRIvol2Vol( mri, m_MRITemp, vox2vox, m_nInterpolationMethod, 0 );
    MRIfree( &mri );
  }
  else
  {
    m_MRITemp = mri;
  }

  MatrixFree( &vox2vox );
  MRIvalRange( m_MRITemp, &m_fMinValue, &m_fMaxValue );

  return true;
}

double FSVolume::GetVoxelValue( int i, int j, int k, int frame )
{
  if ( i < 0 || i >= m_MRI->width ||
       j < 0 || j >= m_MRI->height ||
       k < 0 || k >= m_MRI->depth ||
       frame >= m_MRI->nframes )
  {
    return 0;
  }
  else
  {
    return MRIgetVoxVal( m_MRI, i, j, k, frame );
  }
}

int FSVolume::OriginalIndexToRAS( float iIdxX, float iIdxY, float iIdxZ,
                                  float& oRASX, float& oRASY, float& oRASZ )
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present.\n";
    return 1;
  }

  double ix, iy, iz, wx, wy, wz;
  int r;

  ix = iIdxX;
  iy = iIdxY;
  iz = iIdxZ;
  r = ::MRIvoxelToWorld( m_MRI, ix, iy, iz, &wx, &wy, &wz );
  oRASX = wx;
  oRASY = wy;
  oRASZ = wz;

  return r;
}

int FSVolume::RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                                   float& oIdxX, float& oIdxY, float& oIdxZ )
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present.\n";
    return 1;
  }

  double ix, iy, iz, wx, wy, wz;
  int r;

  if (m_matReg)
  {
    double v[4] = {iRASX, iRASY, iRASZ, 1};
    MATRIX* matNative = MRItkReg2Native(this->m_MRIRef, m_MRI, m_matReg);
    double m[16];
    for ( int i = 0; i < 16; i++ )
    {
      m[i] = (double) *MATRIX_RELT((matNative),(i/4)+1,(i%4)+1);
    }
    vtkMatrix4x4::Invert(m, m);
    vtkMatrix4x4::MultiplyPoint(m, v, v);
    MatrixFree(&matNative);
    wx = v[0];
    wy = v[1];
    wz = v[2];
  }
  else
  {
    wx = iRASX;
    wy = iRASY;
    wz = iRASZ;
  }
  r = ::MRIworldToVoxel( m_MRI, wx, wy, wz, &ix, &iy, &iz );
  oIdxX = ix;
  oIdxY = iy;
  oIdxZ = iz;

  return r;
}

int FSVolume::RASToOriginalIndex ( float iRASX, float iRASY, float iRASZ,
                                   int& oIdxX, int& oIdxY, int& oIdxZ )
{
  /*
  float ix, iy, iz;
  int r = RASToOriginalIndex(iRASX, iRASY, iRASZ, ix, iy, iz);
  if (r == 0) // no error
  {
    oIdxX = (int)( ix + 0.5 );
    oIdxY = (int)( iy + 0.5 );
    oIdxZ = (int)( iz + 0.5 );
  }
  return r;
  */
  double wx, wy, wz;
  int r;

  if (m_matReg)
  {
    double v[4] = {iRASX, iRASY, iRASZ, 1};
    MATRIX* matNative = MRItkReg2Native(this->m_MRIRef, m_MRI, m_matReg);
    double m[16];
    for ( int i = 0; i < 16; i++ )
    {
      m[i] = (double) *MATRIX_RELT((matNative),(i/4)+1,(i%4)+1);
    }
    vtkMatrix4x4::Invert(m, m);
    vtkMatrix4x4::MultiplyPoint(m, v, v);
    MatrixFree(&matNative);
    wx = v[0];
    wy = v[1];
    wz = v[2];
  }
  else
  {
    wx = iRASX;
    wy = iRASY;
    wz = iRASZ;
  }
  r = ::MRIworldToVoxelIndex(m_MRI, wx, wy, wz, &oIdxX, &oIdxY, &oIdxZ );

  return r;
}

bool FSVolume::RASToTalairach(const double *pos_in, double *pos_out)
{
  if (m_MRI->linear_transform)
  {
    double x, y, z;
    ::MRIworldToTalairach(m_MRI, pos_in[0], pos_in[1], pos_in[2], &x, &y, &z);
    pos_out[0] = x;
    pos_out[1] = y;
    pos_out[2] = z;
    return true;
  }
  else
    return false;
}

void FSVolume::TalairachToRAS(const double *pos_in, double *pos_out)
{
  double x, y, z, vx, vy, vz;
  ::MRItalairachToVoxel(m_MRI, pos_in[0], pos_in[1], pos_in[2], &vx, &vy, &vz);
  ::MRIvoxelToWorld(m_MRI, vx, vy, vz, &x, &y, &z);
  pos_out[0] = x;
  pos_out[1] = y;
  pos_out[2] = z;
}

MRI* FSVolume::CreateTargetMRI( MRI* src, MRI* refTarget, bool bAllocatePixel, bool bConform )
{
  MRI* mri = NULL;
  float cornerFactor[3];
  double RAS[3], index[3];
  double indexBounds[6];
  indexBounds[0] = indexBounds[2] = indexBounds[4] = 1e10;
  indexBounds[1] = indexBounds[3] = indexBounds[5] = -1e10;
  for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
  {
    for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
    {
      for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
      {
        ::MRIvoxelToWorld( src,
                           cornerFactor[0]*src->width-0.5,
            cornerFactor[1]*src->height-0.5,
            cornerFactor[2]*src->depth-0.5,
            &RAS[0], &RAS[1], &RAS[2] );
        ::MRIworldToVoxel( refTarget,
                           RAS[0], RAS[1], RAS[2],
            &index[0], &index[1], &index[2] );

        if ( index[0] < indexBounds[0] )
        {
          indexBounds[0] = index[0];
        }
        if ( index[0] > indexBounds[1] )
        {
          indexBounds[1] = index[0];
        }
        if ( index[1] < indexBounds[2] )
        {
          indexBounds[2] = index[1];
        }
        if ( index[1] > indexBounds[3] )
        {
          indexBounds[3] = index[1];
        }
        if ( index[2] < indexBounds[4] )
        {
          indexBounds[4] = index[2];
        }
        if ( index[2] > indexBounds[5] )
        {
          indexBounds[5] = index[2];
        }
      }
    }
  }

  if ( bConform )
  {
    int dim[3] = {refTarget->width, refTarget->height, refTarget->depth};
    //    for ( int i = 0; i < 3; i++ )
    //    {
    //      dim[i] = (int)(indexBounds[i*2+1] - indexBounds[i*2] + 0.5 );
    //    }
    if ( bAllocatePixel )
    {
      try {
        mri = MRIallocSequence( dim[0], dim[1], dim[2], src->type, src->nframes );
      } catch (int ret) {
        return NULL;
      }
    }
    else
    {
      mri = MRIallocHeader( dim[0], dim[1], dim[2], src->type, src->nframes );
    }

    if ( mri == NULL )
    {
      return NULL;
    }

    MRIcopyHeader( refTarget, mri );

    //    double p0[3];
    //    ::MRIvoxelToWorld( refTarget,
    //                       (int)(indexBounds[0]),
    //                       (int)(indexBounds[2]),
    //                       (int)(indexBounds[4]),
    //                       &p0[0], &p0[1], &p0[2] );
    //    MRIp0ToCRAS( mri, p0[0], p0[1], p0[2] );
    return mri;
  }
  else
  {
    // find out the voxelsize of the converted target volume
    MATRIX* mat = MRIgetVoxelToVoxelXform( src, refTarget );
    double m[16];
    for ( int i = 0; i < 16; i++ )
    {
      m[i] = (double) *MATRIX_RELT((mat),(i/4)+1,(i%4)+1);
    }
    double pixelSize[3] = { src->xsize, src->ysize, src->zsize };
    if ( fabs( m[0] ) > fabs( m[4] ) && fabs( m[0] ) > fabs( m[8] ) )
    {
      pixelSize[0] = src->xsize;
    }
    else if ( fabs( m[4] ) > fabs( m[8] ) )
    {
      pixelSize[1] = src->xsize;
    }
    else
    {
      pixelSize[2] = src->xsize;
    }

    if ( fabs( m[1] ) > fabs( m[5] ) && fabs( m[1] ) > fabs( m[9] ) )
    {
      pixelSize[0] = src->ysize;
    }
    else if ( fabs( m[5] ) > fabs( m[9] ) )
    {
      pixelSize[1] = src->ysize;
    }
    else
    {
      pixelSize[2] = src->ysize;
    }

    if ( fabs( m[2] ) > fabs( m[6] ) && fabs( m[2] ) > fabs( m[10] ) )
    {
      pixelSize[0] = src->zsize;
    }
    else if ( fabs( m[6] ) > fabs( m[10] ) )
    {
      pixelSize[1] = src->zsize;
    }
    else
    {
      pixelSize[2] = src->zsize;
    }

    // calculate dimension of the converted target volume
    int dim[3];
    dim[0] = (int)( ( indexBounds[1] - indexBounds[0] ) * refTarget->xsize / pixelSize[0] + 0.5 );
    dim[1] = (int)( ( indexBounds[3] - indexBounds[2] ) * refTarget->ysize / pixelSize[1] + 0.5 );
    dim[2] = (int)( ( indexBounds[5] - indexBounds[4] ) * refTarget->zsize / pixelSize[2] + 0.5 );

    if (qAbs(dim[0]-refTarget->width)%2 != 0)
      dim[0] = dim[0]+1;
    if (qAbs(dim[1]-refTarget->height)%2 != 0)
      dim[1] = dim[1]+1;
    if (qAbs(dim[2]-refTarget->depth)%2 != 0)
      dim[2] = dim[2]+1;

    if ( bAllocatePixel )
    {
      try {
        mri = MRIallocSequence( dim[0], dim[1], dim[2], src->type, src->nframes );
      } catch (int ret) {
        return NULL;
      }
    }
    else
    {
      mri = MRIallocHeader( dim[0], dim[1], dim[2], src->type, src->nframes );
    }

    if ( mri == NULL )
    {
      return NULL;
    }

    MRIcopyHeader( refTarget, mri );
    MRIsetResolution( mri, pixelSize[0], pixelSize[1], pixelSize[2] );
    double p0[3];
    ::MRIvoxelToWorld( refTarget,
                       indexBounds[0]+0.5*pixelSize[0]/refTarget->xsize,
        indexBounds[2]+0.5*pixelSize[1]/refTarget->ysize,
        indexBounds[4]+0.5*pixelSize[2]/refTarget->zsize,
        &p0[0], &p0[1], &p0[2] );
    MRIp0ToCRAS( mri, p0[0], p0[1], p0[2] );

    if (false) // disable this part
    {
      // make sure voxel boundaries are aligned
      double cpt[3], cpt_ext[3], dist[3], step_dist[3];
      dist[0] = refTarget->c_r - mri->c_r;
      dist[1] = refTarget->c_a - mri->c_a;
      dist[2] = refTarget->c_s - mri->c_s;
      step_dist[0] = qMin((float)pixelSize[0], refTarget->xsize)/2;
      step_dist[1] = qMin((float)pixelSize[1], refTarget->ysize)/2;
      step_dist[2] = qMin((float)pixelSize[2], refTarget->zsize)/2;

      for (int i = 0; i < 3; i++)
      {
        dist[i] = (qAbs(dist[i])/step_dist[i] - (int)(qAbs(dist[i])/step_dist[i]))*step_dist[i]*(dist[i]>=0?1:-1);
        if (dist[i] < -step_dist[i]/2)
            dist[i] += step_dist[i];
        else if (dist[i] > step_dist[i]/2)
            dist[i] -= step_dist[i];
      }

      mri->c_r += dist[0];
      mri->c_a += dist[1];
      mri->c_s += dist[2];

      MATRIX *tmp;
      tmp = extract_i_to_r( mri );
      AffineMatrixAlloc( &(mri->i_to_r__ ) );
      SetAffineMatrix( mri->i_to_r__, tmp );
      MatrixFree( &tmp );

      if( mri->r_to_i__ ) {
        MatrixFree(&mri->r_to_i__);
      }
      mri->r_to_i__ = extract_r_to_i(mri);
    }
  }

  return mri;
}

bool FSVolume::MapMRIToImage( bool do_not_create_image )
{
  // first create target MRI
  float bounds[6];
  double voxelSize[3];
  this->GetPixelSize( voxelSize );
  int dim[3];

  MRI* rasMRI = NULL;
  MATRIX* m = MatrixZero( 4, 4, NULL );
  if (m_matReg && m_MRIRef && m_volumeRef )
  {
    // if there is registration matrix, set target as the reference's target
    MRI* mri = m_volumeRef->m_MRITarget;
    try {
      rasMRI = MRIallocSequence( mri->width,
                                 mri->height,
                                 mri->depth,
                                 m_MRI->type,
                                 m_MRI->nframes );
    } catch (int ret) {
      return false;
    }

    if ( rasMRI == NULL )
    {
      cerr << "Can not allocate memory for volume transformation\n";
      MatrixFree( &m );
      return false;
    }
    MRIcopyHeader( mri, rasMRI );
  }
  else if ( m_bResampleToRAS && ( !m_volumeRef ) )
  {
    this->GetBounds( bounds );
    for ( int i = 0; i < 3; i++ )
    {
      dim[i] = (int) ( ( bounds[i*2+1] - bounds[i*2] ) / voxelSize[i] + 0.5 );
    }

    try {
      rasMRI = MRIallocSequence( dim[0], dim[1], dim[2],
          m_MRI->type, m_MRI->nframes );
    } catch (int ret) {
      return false;
    }

    if ( rasMRI == NULL )
    {
      cerr << "Can not allocate memory for volume transformation\n";
      MatrixFree( &m );
      return false;
    }
    MRIsetResolution( rasMRI, voxelSize[0], voxelSize[1], voxelSize[2] );

    *MATRIX_RELT( m, 1, 1 ) = voxelSize[0];
    *MATRIX_RELT( m, 2, 2 ) = voxelSize[1];
    *MATRIX_RELT( m, 3, 3 ) = voxelSize[2];
    *MATRIX_RELT( m, 4, 4 ) = 1;
    *MATRIX_RELT( m, 1, 4 ) = bounds[0];
    *MATRIX_RELT( m, 2, 4 ) = bounds[2];
    *MATRIX_RELT( m, 3, 4 ) = bounds[4];
    MRIsetVoxelToRasXform( rasMRI, m);
  }
  else
  {
    if ( !m_volumeRef )
    {
      double* rtv = this->GetVoxelToRASMatrix();
      int odim[3] = { m_MRI->width, m_MRI->height, m_MRI->depth };
      double n[4] = {0, 0, 0, 1}, ras_orig[4], ras0[4], ras1[4], ras2[4];
      vtkMatrix4x4::MultiplyPoint(rtv, n, ras_orig);
      n[0] = 1;
      vtkMatrix4x4::MultiplyPoint(rtv, n, ras0);
      n[0] = 0; n[1] = 1;
      vtkMatrix4x4::MultiplyPoint(rtv, n, ras1);
      n[1] = 0; n[2] = 1;
      vtkMatrix4x4::MultiplyPoint(rtv, n, ras2);
      double delta0[3], delta1[3], delta2[3];
      for (int i = 0; i < 3; i++)
      {
        delta0[i] = ras0[i]-ras_orig[i];
        delta1[i] = ras1[i]-ras_orig[i];
        delta2[i] = ras2[i]-ras_orig[i];
      }
      if ( fabs( delta0[0] ) >= fabs( delta0[1] ) &&
           fabs( delta0[0] ) >= fabs( delta0[2] ) )
      {
        *MATRIX_RELT( m, 1, 1 ) = ( delta0[0] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 1, 4 ) = ( delta0[0] > 0 ? 0 : odim[0] - 1 );
        dim[0] = odim[0];

        if (fabs(delta1[1]) >= fabs(delta1[2]))
        {
          *MATRIX_RELT( m, 2, 2 ) = ( delta1[1] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 2, 4 ) = ( delta1[1] > 0 ? 0 : odim[1] - 1 );
          dim[1] = odim[1];
          *MATRIX_RELT( m, 3, 3 ) = ( delta2[2] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 3, 4 ) = ( delta2[2] > 0 ? 0 : odim[2] - 1 );
          dim[2] = odim[2];
        }
        else
        {
          *MATRIX_RELT( m, 3, 2 ) = ( delta1[2] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 3, 4 ) = ( delta1[2] > 0 ? 0 : odim[1] - 1 );
          dim[2] = odim[1];
          *MATRIX_RELT( m, 2, 3 ) = ( delta2[1] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 2, 4 ) = ( delta2[1] > 0 ? 0 : odim[2] - 1 );
          dim[1] = odim[2];
        }
      }
      else if ( fabs( delta0[1] ) >= fabs( delta0[2] ) )
      {
        *MATRIX_RELT( m, 2, 1 ) = ( delta0[1] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 2, 4 ) = ( delta0[1] > 0 ? 0 : odim[0] - 1 );
        dim[1] = odim[0];
        if (fabs(delta1[0]) >= fabs(delta1[2]))
        {
          *MATRIX_RELT( m, 1, 2 ) = ( delta1[0] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 1, 4 ) = ( delta1[0] > 0 ? 0 : odim[1] - 1 );
          dim[0] = odim[1];
          *MATRIX_RELT( m, 3, 3 ) = ( delta2[2] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 3, 4 ) = ( delta2[2] > 0 ? 0 : odim[2] - 1 );
          dim[2] = odim[2];
        }
        else
        {
          *MATRIX_RELT( m, 1, 3 ) = ( delta2[0] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 1, 4 ) = ( delta2[0] > 0 ? 0 : odim[2] - 1 );
          dim[0] = odim[2];
          *MATRIX_RELT( m, 3, 2 ) = ( delta1[2] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 3, 4 ) = ( delta1[2] > 0 ? 0 : odim[1] - 1 );
          dim[2] = odim[1];
        }
      }
      else
      {
        *MATRIX_RELT( m, 3, 1 ) = ( delta0[2] > 0 ? 1 : -1 );
        *MATRIX_RELT( m, 3, 4 ) = ( delta0[2] > 0 ? 0 : odim[0] - 1 );
        dim[2] = odim[0];
        if (fabs(delta1[0]) >= fabs(delta1[1]))
        {
          *MATRIX_RELT( m, 1, 2 ) = ( delta1[0] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 1, 4 ) = ( delta1[0] > 0 ? 0 : odim[1] - 1 );
          dim[0] = odim[1];
          *MATRIX_RELT( m, 2, 3 ) = ( delta2[1] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 2, 4 ) = ( delta2[1] > 0 ? 0 : odim[2] - 1 );
          dim[1] = odim[2];
        }
        else
        {
          *MATRIX_RELT( m, 2, 2 ) = ( delta1[1] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 2, 4 ) = ( delta1[1] > 0 ? 0 : odim[1] - 1 );
          dim[1] = odim[1];
          *MATRIX_RELT( m, 1, 3 ) = ( delta2[0] > 0 ? 1 : -1 );
          *MATRIX_RELT( m, 1, 4 ) = ( delta2[0] > 0 ? 0 : odim[2] - 1 );
          dim[0] = odim[2];
        }
      }

      *MATRIX_RELT( m, 4, 4 ) = 1;

      try {
        rasMRI = MRIallocSequence( dim[0], dim[1], dim[2],
            m_MRI->type, m_MRI->nframes );
      } catch (int ret) {
        return false;
      }

      if ( rasMRI == NULL )
      {
        cerr << "Can not allocate memory for volume transformation\n";
        MatrixFree( &m );
        return false;
      }
      MRIsetResolution( rasMRI, voxelSize[0], voxelSize[1], voxelSize[2] );

      MATRIX* m1 = MRIgetVoxelToRasXform( m_MRI );
      if (!m1)
        m1 = MatrixIdentity(4, NULL);
      MATRIX* m_inv = MatrixInverse( m, NULL );
      MATRIX* m2 = MatrixMultiply( m1, m_inv, NULL );

      MRIsetVoxelToRasXform( rasMRI, m2 );

      MatrixFree( &m1 );
      MatrixFree( &m2 );
      MatrixFree( &m_inv );
    }
    else
    {
      rasMRI = CreateTargetMRI( m_MRI, m_volumeRef->m_MRITarget, true, m_bConform );
      if ( rasMRI == NULL )
      {
        cerr << "Can not allocate memory for volume transformation\n";
        MatrixFree( &m );
        return false;
      }
    }
  }
  MatrixFree( &m );

//  if ( m_matReg && !m_MRIRef )
//  {
//    cerr << "No target volume available! Cannot use registration matrix.\n";
//  }

  if ( m_matReg)
  {
    // registration matrix is always converted to tkReg style now
    // if ( m_nRegType == REG_TKREGISTER )
    {
      MATRIX *vox2vox, *Tmov, *invTmov, *Ttarg;
      Tmov = MRIxfmCRS2XYZtkreg( m_MRI );
      invTmov = MatrixInverse( Tmov, NULL );
      Ttarg = MRIxfmCRS2XYZtkreg( m_MRIRef );
      // vox2vox = invTmov*R*Ttarg
      vox2vox = MatrixMultiply( invTmov, m_matReg, NULL );
      MatrixMultiply( vox2vox, Ttarg, vox2vox );

      MATRIX* t2r = MRIgetVoxelToVoxelXform( rasMRI, m_MRIRef );
      MatrixMultiply( vox2vox, t2r, t2r );

      if (m_MRI->nframes == 3)
      {
        MATRIX* rot = MatrixAlloc(3, 3, MATRIX_REAL);
        double scale[4];
        for (int i = 1; i <= 3; i++)
        {
          scale[i] = qSqrt(m_matReg->rptr[1][i]*m_matReg->rptr[1][i] + m_matReg->rptr[2][i]*m_matReg->rptr[2][i] + m_matReg->rptr[3][i]*m_matReg->rptr[3][i]);
          if (scale[i] == 0)
            scale[i] = 1;
        }
        for (int i = 1; i <= 3; i++)
        {
          for (int j = 1; j <= 3; j++)
          {
            *MATRIX_RELT(rot, j, i) = m_matReg->rptr[i][j]/scale[j];
          }
        }
        MRIvol2VolR( m_MRI, rasMRI, t2r, m_nInterpolationMethod, 0, rot );
        MatrixFree(&rot);
      }
      else
        MRIvol2Vol( m_MRI, rasMRI, t2r, m_nInterpolationMethod, 0 );

      // copy vox2vox
      MatrixInverse( t2r, vox2vox );
      for ( int i = 0; i < 16; i++ )
      {
        m_VoxelToVoxelMatrix[i] =
            (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
      }

      MatrixFree( &vox2vox );
      MatrixFree( &Tmov );
      MatrixFree( &invTmov );
      MatrixFree( &Ttarg );
      MatrixFree( &t2r );
    }
  }
  else
  {
//    qDebug() << rasMRI->width << rasMRI->height << rasMRI->depth;

//    QElapsedTimer t; t.start();
//    qDebug() << "begin vol2vol";
    MRIvol2Vol( m_MRI, rasMRI, NULL, m_nInterpolationMethod, 0 );
//    qDebug() << "vol2vol time: " << t.elapsed()/1000;
    MATRIX* vox2vox = MRIgetVoxelToVoxelXform( m_MRI, rasMRI );
    for ( int i = 0; i < 16; i++ )
    {
      m_VoxelToVoxelMatrix[i] =
          (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
    }
    MatrixFree( &vox2vox );
  }

  SetMRITarget( rasMRI );
  UpdateRASToRASMatrix();

  if ( !do_not_create_image && !CreateImage( rasMRI ) )
  {
    return false;
  }

  // copy mri pixel data to vtkImage we will use for display
  CopyMRIDataToImage( rasMRI, m_imageData );

  // Need to recalc our bounds at some point.
  m_bBoundsCacheDirty = true;

  ::MRIfree( &rasMRI );

  return true;
}

void FSVolume::UpdateRASToRASMatrix()
{
  // compute ras2ras matrix
  MATRIX* r2r = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT( r2r, (i/4)+1, (i%4)+1 ) = m_VoxelToVoxelMatrix[i];
  }
  MATRIX* mMov = MRIgetVoxelToRasXform( m_MRI );
  MATRIX* mTarg = MRIgetVoxelToRasXform( m_MRITarget );
  MATRIX* mMov_inv = MatrixInverse( mMov, NULL );
  MatrixMultiply( r2r, mMov_inv, r2r );
  MatrixMultiply( mTarg, r2r, r2r );

  for ( int i = 0; i < 16; i++ )
  {
    m_RASToRASMatrix[i] = (double) *MATRIX_RELT((r2r),(i/4)+1,(i%4)+1);
  }

  MatrixFree( &r2r );
  MatrixFree( &mMov );
  MatrixFree( &mMov_inv );
  MatrixFree( &mTarg );
}

bool FSVolume::CreateImage( MRI* rasMRI )
{
  // first copy mri data to image
  vtkDataArray *scalars = NULL;
  vtkUnsignedCharArray  *ucharScalars = NULL;
  vtkIntArray           *intScalars = NULL;
  vtkShortArray         *shortScalars = NULL;
  vtkLongArray          *longScalars = NULL;
  vtkFloatArray         *floatScalars = NULL;
  vtkIdType cValues;
  //  int zElement=0;

  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present.\n";
    return false;
  }

  vtkIdType zX = rasMRI->width;
  vtkIdType zY = rasMRI->height;
  vtkIdType zZ = rasMRI->depth;
  vtkIdType zFrames = rasMRI->nframes;

  m_imageData = vtkSmartPointer<vtkImageData>::New();
  vtkImageData* imageData = m_imageData;

  // This object's output space is in voxel coordinates.
  imageData->SetDimensions( zX, zY, zZ );

  double origin[3] = { 0, 0, 0 };
  if ( m_bResampleToRAS && !m_volumeRef )
  {
    float bounds[6];
    GetBounds( bounds );
    origin[0] = bounds[0];
    origin[1] = bounds[2];
    origin[2] = bounds[4];
  }
  else if ( m_volumeRef )
  {
    double ras[3], cindex[3];

    ::MRIvoxelToWorld( rasMRI, 0., 0., 0., &ras[0], &ras[1], &ras[2] );
    ::MRIworldToVoxel( m_volumeRef->m_MRITarget,
                       ras[0], ras[1], ras[2],
        &cindex[0], &cindex[1], &cindex[2] );

    m_volumeRef->GetImageOutput()->GetOrigin( origin );
    for ( int i = 0; i < 3; i++ )
    {
      if ( fabs(cindex[i]-nint(cindex[i])) < 1e-4 )
      {
        cindex[i] = nint(cindex[i]);
      }
    }

    origin[0] += cindex[0] * m_volumeRef->m_MRITarget->xsize;
    origin[1] += cindex[1] * m_volumeRef->m_MRITarget->ysize;
    origin[2] += cindex[2] * m_volumeRef->m_MRITarget->zsize;

    //    qDebug() << m_volumeRef->m_MRITarget->x_r << m_volumeRef->m_MRITarget->y_r << m_volumeRef->m_MRITarget->z_r;
    //    qDebug() << m_volumeRef->m_MRITarget->x_a << m_volumeRef->m_MRITarget->y_a << m_volumeRef->m_MRITarget->z_a;
    //    qDebug() << m_volumeRef->m_MRITarget->x_s << m_volumeRef->m_MRITarget->y_s << m_volumeRef->m_MRITarget->z_s;

    //    qDebug() << rasMRI->c_r << rasMRI->c_a << rasMRI->c_s << m_volumeRef->m_MRITarget->c_r << m_volumeRef->m_MRITarget->c_a
    //             << m_volumeRef->m_MRITarget->c_s;
    //    qDebug() << cindex[0] << cindex[1] << cindex[2];

    /*
    double ras[3], ras2[3];
    ::MRIvoxelToWorld( rasMRI, 0., 0., 0., &ras[0], &ras[1], &ras[2] );
    ::MRIvoxelToWorld( m_volumeRef->m_MRITarget, 0., 0., 0., &ras2[0], &ras2[1], &ras2[2] );
    m_volumeRef->GetImageOutput()->GetOrigin( origin );
    for ( int i = 0; i < 3; i++ )
      origin[i] += ( ras[i] - ras2[i] );
    */
  }

  imageData->SetSpacing( rasMRI->xsize, rasMRI->ysize, rasMRI->zsize );
  imageData->SetOrigin( origin[0], origin[1], origin[2] );
  //  imageData->SetWholeExtent( 0, zX-1, 0, zY-1, 0, zZ-1 );
  imageData->SetDimensions(zX, zY, zZ);
  if (rasMRI->type == MRI_RGB)
    zFrames = 4;

  // create the scalars for all of the images. set the element size
  // for the data we will read.
#if VTK_MAJOR_VERSION > 5
  switch ( rasMRI->type )
  {
  case MRI_UCHAR:
  case MRI_RGB:
    imageData->AllocateScalars(VTK_UNSIGNED_CHAR, zFrames);
    break;
  case MRI_INT:
    imageData->AllocateScalars(VTK_INT, zFrames);
    break;
  case MRI_LONG:
    imageData->AllocateScalars(VTK_LONG, zFrames);
    break;
  case MRI_FLOAT:
    imageData->AllocateScalars(VTK_FLOAT, zFrames);
    break;
  case MRI_SHORT:
    imageData->AllocateScalars(VTK_SHORT, zFrames);
    break;
  case MRI_USHRT:
    imageData->AllocateScalars(VTK_UNSIGNED_SHORT, zFrames);
    break;
  default:
    return false;
  }
#else
  imageData->SetNumberOfScalarComponents(zFrames);
  switch ( rasMRI->type )
  {
  case MRI_UCHAR:
  case MRI_RGB:
    imageData->SetScalarTypeToUnsignedChar();
    break;
  case MRI_INT:
    imageData->SetScalarTypeToInt();
    break;
  case MRI_LONG:
    imageData->SetScalarTypeToLong();
    break;
  case MRI_FLOAT:
    imageData->SetScalarTypeToFloat();
    break;
  case MRI_SHORT:
    imageData->SetScalarTypeToShort();
    break;
  case MRI_USHRT:
    imageData->SetScalarTypeToUnsignedShort();
    break;
  default:
    return false;
  }
  imageData->AllocateScalars();
#endif

  return true;
}

bool FSVolume::ResizeRotatedImage( MRI* rasMRI, MRI* refTarget, vtkImageData* refImageData,
                                   double* rasPoint )
{
  // first copy mri data to image
  vtkDataArray *scalars = NULL;
  vtkUnsignedCharArray  *ucharScalars = NULL;
  vtkIntArray           *intScalars = NULL;
  vtkShortArray         *shortScalars = NULL;
  vtkLongArray          *longScalars = NULL;
  vtkFloatArray         *floatScalars = NULL;
  vtkIdType cValues;
  int zElement=0;

  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present.\n";
    return false;
  }

  vtkIdType zX = rasMRI->width;
  vtkIdType zY = rasMRI->height;
  vtkIdType zZ = rasMRI->depth;
  vtkIdType zFrames = rasMRI->nframes;

  m_imageData = vtkSmartPointer<vtkImageData>::New();
  vtkImageData* imageData = m_imageData;

  // This object's output space is in voxel coordinates.
  imageData->SetDimensions( zX, zY, zZ );

  double origin[3] = { 0, 0, 0 };
  refImageData->GetOrigin( origin );

  double vox[3], tvox[3];
  ::MRIworldToVoxel( rasMRI,
                     rasPoint[0], rasPoint[1], rasPoint[2],
      &vox[0], &vox[1], &vox[2] );
  ::MRIworldToVoxel( refTarget,
                     rasPoint[0], rasPoint[1], rasPoint[2],
      &tvox[0], &tvox[1], &tvox[2] );
  //  cout << "original: " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
  origin[0] -= ( vox[0] * rasMRI->xsize - tvox[0] * refTarget->xsize );
  origin[1] -= ( vox[1] * rasMRI->ysize - tvox[1] * refTarget->ysize );
  origin[2] -= ( vox[2] * rasMRI->zsize - tvox[2] * refTarget->zsize );

  //  cout << "new: " << origin[0] << " " << origin[1] << " " << origin[2] << endl; fflush(0);

  imageData->SetSpacing( rasMRI->xsize, rasMRI->ysize, rasMRI->zsize );
  imageData->SetOrigin( origin[0], origin[1], origin[2] );

  imageData->SetExtent( 0, zX-1, 0, zY-1, 0, zZ-1 );

  if (rasMRI->type == MRI_RGB)
    zFrames = 4;

  // create the scalars for all of the images. set the element size
  // for the data we will read.
#if VTK_MAJOR_VERSION > 5
  switch ( rasMRI->type )
  {
  case MRI_UCHAR:
  case MRI_RGB:
    imageData->AllocateScalars(VTK_UNSIGNED_CHAR, zFrames);
    break;
  case MRI_INT:
    imageData->AllocateScalars(VTK_INT, zFrames);
    break;
  case MRI_LONG:
    imageData->AllocateScalars(VTK_LONG, zFrames);
    break;
  case MRI_FLOAT:
    imageData->AllocateScalars(VTK_FLOAT, zFrames);
    break;
  case MRI_SHORT:
    imageData->AllocateScalars(VTK_SHORT, zFrames);
    break;
  case MRI_USHRT:
    imageData->AllocateScalars(VTK_UNSIGNED_SHORT, zFrames);
    break;
  default:
    return false;
  }
#else
  imageData->SetNumberOfScalarComponents(zFrames);
  switch ( rasMRI->type )
  {
  case MRI_UCHAR:
  case MRI_RGB:
    imageData->SetScalarTypeToUnsignedChar();
    break;
  case MRI_INT:
    imageData->SetScalarTypeToInt();
    break;
  case MRI_LONG:
    imageData->SetScalarTypeToLong();
    break;
  case MRI_FLOAT:
    imageData->SetScalarTypeToFloat();
    break;
  case MRI_SHORT:
    imageData->SetScalarTypeToShort();
    break;
  case MRI_USHRT:
    imageData->SetScalarTypeToUnsignedShort();
    break;
  default:
    break ;
  }
  imageData->AllocateScalars();
#endif

  return true;
}

bool FSVolume::Rotate( std::vector<RotationElement>& rotations,
                       int nSampleMethod )
{
  if ( rotations.size() == 0 )
  {
    return false;
  }

  if ( !m_MRITemp )
  {
    cerr << "Volume not ready for rotation\n";
    return false;
  }

  MRI* rasMRI = NULL;
  if ( rotations.size() == 0 || rotations[0].Plane == -1 )   // restore
  {
    if ( !m_MRIOrigTarget )  // try to restore but no where to restore
    {
      cout << "Can not restore because no original space was defined.\n";
      return false;
    }

    try {
      rasMRI = MRIallocSequence( m_MRIOrigTarget->width,
                                 m_MRIOrigTarget->height,
                                 m_MRIOrigTarget->depth,
                                 m_MRI->type,
                                 m_MRI->nframes );
    } catch (int ret) {
      return false;
    }

    if ( rasMRI == NULL )
    {
      cout << "Can not allocate memory for volume transformation.\n";
      return false;
    }
    MRIcopyHeader( m_MRIOrigTarget, rasMRI );
  }
  else    // rotate
  {
    MATRIX* m = MRIgetVoxelToRasXform( m_MRITarget );

    // calculate rotation matrix
    MATRIX* m_r = MatrixIdentity( 4, NULL ) ;
    for ( size_t i = 0; i < rotations.size(); i++ )
    {
      MATRIX* m_tmp = GetRotationMatrix( rotations[i].Plane,
                                         rotations[i].Angle,
                                         rotations[i].Point );
      MATRIX* m_tmp1 = MatrixMultiply( m_tmp, m_r, NULL );

      MatrixCopy( m_tmp1, m_r );
      MatrixFree( &m_tmp );
      MatrixFree( &m_tmp1 );
    }

    MATRIX* m_new = MatrixMultiply( m_r, m, NULL );

    rasMRI = MRIallocHeader( m_MRITarget->width,
                             m_MRITarget->height,
                             m_MRITarget->depth,
                             m_MRI->type,
                             m_MRI->nframes);
    MRIsetResolution( rasMRI,
                      m_MRITarget->xsize,
                      m_MRITarget->ysize,
                      m_MRITarget->zsize );
    MRIsetVoxelToRasXform( rasMRI, m_new );

    MatrixFree( &m );
    MatrixFree( &m_r );
    MatrixFree( &m_new );

    if ( !m_MRIOrigTarget )  // never rotated before, so save it as original
    {
      m_MRIOrigTarget = MRIallocHeader( m_MRITarget->width,
                                        m_MRITarget->height,
                                        m_MRITarget->depth,
                                        m_MRITarget->type,
                                        m_MRITarget->nframes);
      MRIcopyHeader( m_MRITarget, m_MRIOrigTarget );
    }
  }

  // change the field of view of the new target to full coverage of the image data
  if ( rotations[0].Plane != -1 )
  {
    MRI* mri = rasMRI;
    rasMRI = CreateTargetMRI( m_MRIOrigTarget, mri );
    ::MRIfree( &mri );
  }

  if ( nSampleMethod < 0 )
  {
    nSampleMethod = rotations[0].SampleMethod;
  }
  if ( rotations[0].Plane == -1 ) // restore
  {
    MRIvol2Vol( m_MRITemp, rasMRI, NULL, m_nInterpolationMethod, 0 );
  }
  else
  {
    MRIvol2Vol( m_MRITemp, rasMRI, NULL, nSampleMethod, 0 );
  }

  MRIfree( &m_MRITemp );
  m_MRITemp = NULL;

  // copy vox2vox
  MATRIX* vox2vox = MRIgetVoxelToVoxelXform( m_MRI, rasMRI );
  for ( int i = 0; i < 16; i++ )
  {
    m_VoxelToVoxelMatrix[i] = (double) *MATRIX_RELT((vox2vox),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &vox2vox );

  bool bUpdateImage;
  if ( rotations[0].Plane == -1 )
  {
    bUpdateImage = CreateImage( rasMRI );
  }
  else
    // copy mri pixel data to vtkImage we will use for display
  {
    bUpdateImage = ResizeRotatedImage( rasMRI, m_MRITarget, m_imageData, rotations[0].Point );
  }

  if ( !bUpdateImage )
  {
    return false;
  }

  CopyMRIDataToImage( rasMRI, m_imageData );

  SetMRITarget( rasMRI );
  UpdateRASToRASMatrix();

  // Need to recalc our bounds at some point.
  m_bBoundsCacheDirty = true;

  ::MRIfree( &rasMRI );

  m_bResampleToRAS = false;

  return true;
}

MATRIX* FSVolume::GetRotationMatrix( int nPlane, double angle, double* origin )
{
  // calculate rotation matrix
  double p0[3] = { 0, 0, 0 }, p1[3] = { 0, 0, 0 };
  p1[nPlane] = 1;
  TargetToRAS( p0, p0 );
  TargetToRAS( p1, p1 );

  vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
  tr->Translate( origin );
  tr->RotateWXYZ( angle, p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] );
  tr->Translate( -origin[0], -origin[1], -origin[2] );
  tr->Update();

  vtkSmartPointer<vtkMatrix4x4> m = vtkSmartPointer<vtkMatrix4x4>::New();
  tr->GetTranspose( m );

  MATRIX* m_r = MatrixZero( 4, 4, NULL );
  for ( int i = 0; i < 4; i++ )
    for ( int j = 0; j < 4; j++ )
    {
      *MATRIX_RELT( m_r, i+1, j+1 ) = m->GetElement( i, j );
    }

  return m_r;
}

void FSVolume::CopyMRIDataToImage( MRI* mri,
                                   vtkImageData* image )
{
  // Copy the slice data into the scalars.
  int zX = mri->width;
  int zY = mri->height;
  int zZ = mri->depth;
  int zFrames = mri->nframes;

  vtkIdType nTuple = 0;
  char* ptr = (char*)image->GetScalarPointer();
  int nProgressStep = 20;
  int nProgress = 0;
  for ( int nZ = 0; nZ < zZ; nZ++ )
  {
    for ( int nY = 0; nY < zY; nY++ )
    {
      for ( int nX = 0; nX < zX; nX++ )
      {
        if (mri->type == MRI_RGB)
        {
          int val = MRIIseq_vox(mri, nX, nY, nZ, 0);
          ptr[nTuple*4] = (val & 0x00ff);
          ptr[nTuple*4+1] = ((val >> 8) & 0x00ff);
          ptr[nTuple*4+2] = ((val >> 16) & 0x00ff);
          ptr[nTuple*4+3] = (char)255;
        }
        else
        {
          for ( int nFrame = 0; nFrame < zFrames; nFrame++ )
          {
            switch ( mri->type )
            {
            case MRI_UCHAR:
              ptr[nTuple*zFrames+nFrame] = MRIseq_vox( mri, nX, nY, nZ, nFrame );
              break;
            case MRI_INT:
              ((int*)ptr)[nTuple*zFrames+nFrame] = MRIIseq_vox( mri, nX, nY, nZ, nFrame );
              break;
            case MRI_LONG:
              ((long*)ptr)[nTuple*zFrames+nFrame] = MRILseq_vox( mri, nX, nY, nZ, nFrame );
              break;
            case MRI_FLOAT:
              ((float*)ptr)[nTuple*zFrames+nFrame] = MRIFseq_vox( mri, nX, nY, nZ, nFrame );
              break;
            case MRI_SHORT:
              ((short*)ptr)[nTuple*zFrames+nFrame] = MRISseq_vox( mri, nX, nY, nZ, nFrame );
              break;
            case MRI_USHRT:
              ((unsigned short*)ptr)[nTuple*zFrames+nFrame] = MRIUSseq_vox( mri, nX, nY, nZ, nFrame );
              break;
            default:
              break;
            }
          }
        }
        nTuple++;
      }
    }

    if ( nZ%(max(1, zZ/5)) == 0 )
    {
      nProgress += nProgressStep;
      emit ProgressChanged( nProgress );
    }
  }
}

vtkImageData* FSVolume::GetImageOutput()
{
  return m_imageData;
}

void FSVolume::CopyMatricesFromMRI ()
{
  if ( NULL == m_MRI )
  {
    return;
  }

  MATRIX* m = MRIgetVoxelToRasXform( m_MRI );
  if (m)
  {
    for ( int i = 0; i < 16; i++ )
    {
      m_VoxelToRASMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
    }
    MatrixFree( &m );
  }
  else
  {
    for ( int i = 0; i < 16; i++ )
    {
      m_VoxelToRASMatrix[i] = (i%5 ? 0:1);
    }
  }

  m = MRIgetRasToVoxelXform( m_MRI );
  if (m)
  {
    for ( int i = 0; i < 16; i++ )
    {
      m_RASToVoxelMatrix[i] = (double) *MATRIX_RELT((m),(i/4)+1,(i%4)+1);
    }

    MATRIX* tkreg = MRIxfmCRS2XYZtkreg( m_MRI );
    MATRIX* m1 = MatrixMultiply( tkreg, m, NULL );
    for ( int i = 0; i < 16; i++ )
    {
      m_RASToTkRegMatrix[i] = (double) *MATRIX_RELT((m1),(i/4)+1,(i%4)+1);
    }
    MatrixFree( &tkreg );
    MatrixFree( &m );
    MatrixFree( &m1 );
  }
  else
  {
    for ( int i = 0; i < 16; i++ )
    {
      m_RASToVoxelMatrix[i] = (i%5 ? 0:1);
      m_RASToTkRegMatrix[i] = (i%5 ? 0:1);
    }
  }
}

void FSVolume::GetBounds ( float oRASBounds[6] )
{
  if ( NULL == m_MRI )
  {

    oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
        oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;

    return;
  }

  if ( !m_bResampleToRAS )
  {
    //  if ( m_imageData.GetPointer() )
    if ( m_imageData != NULL )
    {
      double* origin = m_imageData->GetOrigin();
      int* dim = m_imageData->GetDimensions();
      double* vs = m_imageData->GetSpacing();

      for ( int i = 0; i < 3; i++ )
      {
        oRASBounds[i*2] = origin[i];
        oRASBounds[i*2+1] = origin[i] + dim[i] * vs[i];
      }
    }
    else if ( m_MRITarget )
    {
      MATRIX* m = MRIgetVoxelToRasXform( m_MRITarget );
      oRASBounds[0] = *MATRIX_RELT( m, 1, 4 );
      oRASBounds[2] = *MATRIX_RELT( m, 2, 4 );
      oRASBounds[4] = *MATRIX_RELT( m, 3, 4 );
      oRASBounds[1] = m_MRITarget->width * m_MRITarget->xsize + oRASBounds[0];
      oRASBounds[3] = m_MRITarget->height * m_MRITarget->ysize + oRASBounds[2];
      oRASBounds[5] = m_MRITarget->depth * m_MRITarget->zsize + oRASBounds[4];
      MatrixFree( &m );
    }
    else
    {
      cout << "Did not find bounds.\n";
    }
    return;
  }
  else if ( m_MRITarget )
  {
    MATRIX* m = MRIgetVoxelToRasXform( m_MRITarget );
    oRASBounds[0] = *MATRIX_RELT( m, 1, 4 );
    oRASBounds[2] = *MATRIX_RELT( m, 2, 4 );
    oRASBounds[4] = *MATRIX_RELT( m, 3, 4 );
    oRASBounds[1] = m_MRITarget->width * m_MRITarget->xsize + oRASBounds[0];
    oRASBounds[3] = m_MRITarget->height * m_MRITarget->ysize + oRASBounds[2];
    oRASBounds[5] = m_MRITarget->depth * m_MRITarget->zsize + oRASBounds[4];
    MatrixFree( &m );
    return;
  }
  else if ( m_bBoundsCacheDirty )
  {
    m_RASBounds[0] = m_RASBounds[2] = m_RASBounds[4] = 999999;
    m_RASBounds[1] = m_RASBounds[3] = m_RASBounds[5] = -999999;

    // For each corner, convert to RAS, and get the bounds.
    float cornerFactor[3];
    float RAS[3];
    for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
    {
      for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
      {
        for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
        {

          this->OriginalIndexToRAS( cornerFactor[0]*m_MRI->width,
              cornerFactor[1]*m_MRI->height,
              cornerFactor[2]*m_MRI->depth,
              RAS[0], RAS[1], RAS[2] );

          if ( RAS[0] < m_RASBounds[0] )
          {
            m_RASBounds[0] = RAS[0];
          }
          if ( RAS[0] > m_RASBounds[1] )
          {
            m_RASBounds[1] = RAS[0];
          }
          if ( RAS[1] < m_RASBounds[2] )
          {
            m_RASBounds[2] = RAS[1];
          }
          if ( RAS[1] > m_RASBounds[3] )
          {
            m_RASBounds[3] = RAS[1];
          }
          if ( RAS[2] < m_RASBounds[4] )
          {
            m_RASBounds[4] = RAS[2];
          }
          if ( RAS[2] > m_RASBounds[5] )
          {
            m_RASBounds[5] = RAS[2];
          }
        }
      }
    }

    m_bBoundsCacheDirty = false;
  }

  oRASBounds[0] = m_RASBounds[0];
  oRASBounds[1] = m_RASBounds[1];
  oRASBounds[2] = m_RASBounds[2];
  oRASBounds[3] = m_RASBounds[3];
  oRASBounds[4] = m_RASBounds[4];
  oRASBounds[5] = m_RASBounds[5];
}

void FSVolume::GetBounds( MRI* mri, float oRASBounds[6] )
{
  oRASBounds[0] = oRASBounds[2] = oRASBounds[4] = 1e10;
  oRASBounds[1] = oRASBounds[3] = oRASBounds[5] = -1e10;

  // For each corner, convert to RAS, and get the bounds.
  float cornerFactor[3];
  double RAS[3];
  for ( cornerFactor[2] = 0; cornerFactor[2] <= 1; cornerFactor[2]++ )
  {
    for ( cornerFactor[1] = 0; cornerFactor[1] <= 1; cornerFactor[1]++ )
    {
      for ( cornerFactor[0] = 0; cornerFactor[0] <= 1; cornerFactor[0]++ )
      {
        ::MRIvoxelToWorld( mri,
                           cornerFactor[0]*(mri->width-1),
            cornerFactor[1]*(mri->height-1),
            cornerFactor[2]*(mri->depth-1),
            &RAS[0], &RAS[1], &RAS[2] );

        if ( RAS[0] < oRASBounds[0] )
        {
          oRASBounds[0] = RAS[0];
        }
        if ( RAS[0] > oRASBounds[1] )
        {
          oRASBounds[1] = RAS[0];
        }
        if ( RAS[1] < oRASBounds[2] )
        {
          oRASBounds[2] = RAS[1];
        }
        if ( RAS[1] > oRASBounds[3] )
        {
          oRASBounds[3] = RAS[1];
        }
        if ( RAS[2] < oRASBounds[4] )
        {
          oRASBounds[4] = RAS[2];
        }
        if ( RAS[2] > oRASBounds[5] )
        {
          oRASBounds[5] = RAS[2];
        }
      }
    }
  }
}

int FSVolume::GetNumberOfFrames()
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present.\n";
    return 0;
  }

  return m_MRI->nframes;
}

void FSVolume::GetPixelSize( double* pixelSize )
{
  if ( m_MRI == NULL )
  {
    cerr << "No MRI is present.\n";
    return;
  }

  if ( !m_bResampleToRAS )
  {
    if ( m_MRITarget )
    {
      pixelSize[0] = m_MRITarget->xsize;
      pixelSize[1] = m_MRITarget->ysize;
      pixelSize[2] = m_MRITarget->zsize;

      return;
    }
  }
  else if ( m_matReg && m_MRIRef )
  {
    MRI* mri = m_volumeRef->m_MRITarget;
    if ( mri )
    {
      pixelSize[0] = mri->xsize;
      pixelSize[1] = mri->ysize;
      pixelSize[2] = mri->zsize;

      return;
    }
  }

  double* rtv = GetVoxelToRASMatrix();
  pixelSize[0] = m_MRI->xsize;
  pixelSize[1] = m_MRI->ysize;
  pixelSize[2] = m_MRI->zsize;
  double n[4] = {0, 0, 0, 1}, ras_orig[4], ras0[4], ras1[4], ras2[4];
  vtkMatrix4x4::MultiplyPoint(rtv, n, ras_orig);
  n[0] = 1;
  vtkMatrix4x4::MultiplyPoint(rtv, n, ras0);
  n[0] = 0; n[1] = 1;
  vtkMatrix4x4::MultiplyPoint(rtv, n, ras1);
  n[1] = 0; n[2] = 1;
  vtkMatrix4x4::MultiplyPoint(rtv, n, ras2);
  double delta0[3], delta1[3];
  for (int i = 0; i < 3; i++)
  {
    delta0[i] = ras0[i]-ras_orig[i];
    delta1[i] = ras1[i]-ras_orig[i];
    //    delta2[i] = ras2[i]-ras_orig[i];
  }
  double vs[3] = { m_MRI->xsize, m_MRI->ysize, m_MRI->zsize };
  if ( fabs( delta0[0] ) >= fabs( delta0[1] ) &&
       fabs( delta0[0] ) >= fabs( delta0[2] ) )
  {
    pixelSize[0] = vs[0];
    if (fabs(delta1[1]) >= fabs(delta1[2]))
    {
      pixelSize[1] = vs[1];
      pixelSize[2] = vs[2];
    }
    else
    {
      pixelSize[2] = vs[1];
      pixelSize[1] = vs[2];
    }
  }
  else if ( fabs( delta0[1] ) >= fabs( delta0[2] ) )
  {
    pixelSize[1] = vs[0];
    if (fabs(delta1[0]) >= fabs(delta1[2]))
    {
      pixelSize[0] = vs[1];
      pixelSize[2] = vs[2];
    }
    else
    {
      pixelSize[0] = vs[2];
      pixelSize[2] = vs[1];
    }
  }
  else
  {
    pixelSize[2] = vs[0];
    if (fabs(delta1[0]) >= fabs(delta1[1]))
    {
      pixelSize[0] = vs[1];
      pixelSize[1] = vs[2];
    }
    else
    {
      pixelSize[1] = vs[1];
      pixelSize[0] = vs[2];
    }
  }
}

void FSVolume::TargetToRAS( double x_in, double y_in, double z_in,
                            double& x_out, double& y_out, double& z_out )
{
  double pos[3] = { x_in, y_in, z_in };
  TargetToRAS( pos, pos );
  x_out = pos[0];
  y_out = pos[1];
  z_out = pos[2];
}

void FSVolume::TargetToRAS( const double* pos_in, double* pos_out )
{
  double pos[4] = { 0 };
  double vs[3] = {1,1,1};
  double origin[3] = {0};
  if (m_imageData != NULL)
  {
    m_imageData->GetSpacing(vs);
    m_imageData->GetOrigin(origin);
  }
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = ( pos_in[i] - origin[i] ) / vs[i];
  }

  double fpos[3];
  ::MRIvoxelToWorld( m_MRITarget,
                     (float)pos[0], (float)pos[1], (float)pos[2],
      &fpos[0], &fpos[1], &fpos[2] );
  //  cout << "out: " << fpos[0] << " " << fpos[1] << " " << fpos[2] << endl;
  for ( int i = 0; i < 3; i++ )
  {
    pos_out[i] = fpos[i];
  }
}

MATRIX* FSVolume::GetTargetToRASMatrix()
{
  double vs[3] = {1,1,1};
  double origin[3] = {0};
  if (m_imageData != NULL)
  {
    m_imageData->GetSpacing(vs);
    m_imageData->GetOrigin(origin);
  }
  MATRIX* m = MatrixZero( 4, 4, NULL );
  *MATRIX_RELT( m, 1, 1 ) = vs[0];
  *MATRIX_RELT( m, 2, 2 ) = vs[1];
  *MATRIX_RELT( m, 3, 3 ) = vs[2];
  *MATRIX_RELT( m, 4, 4 ) = 1;
  *MATRIX_RELT( m, 1, 4 ) = origin[0];
  *MATRIX_RELT( m, 2, 4 ) = origin[1];
  *MATRIX_RELT( m, 3, 4 ) = origin[2];
  MATRIX* m_invert = MatrixInverse( m, NULL );

  MATRIX* v2r = MRIgetVoxelToRasXform( m_MRITarget );
  MATRIX* t2r = MatrixMultiply( v2r, m_invert, NULL );
  MatrixFree( &m );
  MatrixFree( &v2r );
  MatrixFree( &m_invert );

  return t2r;
}

void FSVolume::TargetToRAS( const float* pos_in, float* pos_out )
{
  double pos[3] = { 0 };
  TargetToRAS( pos_in[0], pos_in[1], pos_in[2], pos[0], pos[1], pos[2] );
  pos_out[0] = pos[0];
  pos_out[1] = pos[1];
  pos_out[2] = pos[2];
}

void FSVolume::RASToTarget( const double* pos_in, double* pos_out )
{
  double pos[4] = { 0 };
  ::MRIworldToVoxel( m_MRITarget,
                     (float)pos_in[0],
      (float)pos_in[1],
      (float)pos_in[2],
      &pos[0],
      &pos[1],
      &pos[2] );
  double vs[3] = {1,1,1};
  double origin[3] = {0};
  if (m_imageData != NULL)
  {
    m_imageData->GetSpacing(vs);
    m_imageData->GetOrigin(origin);
  }
  for ( int i = 0; i < 3; i++ )
  {
    pos_out[i] = vs[i] * pos[i] + origin[i];
  }
}

void FSVolume::RASToTarget( const float* pos_in, float* pos_out )
{
  double p_in[3] = { pos_in[0], pos_in[1], pos_in[2] };
  double p_out[3];
  RASToTarget( p_in, p_out );
  pos_out[0] = p_out[0];
  pos_out[1] = p_out[1];
  pos_out[2] = p_out[2];
}

void FSVolume::RASToTargetIndex( const double* pos_in, int* index_out )
{
  double pos[3];
  RASToTarget( pos_in, pos );
  double vs[3] = {1,1,1};
  double origin[3] = {0};
  if (m_imageData != NULL)
  {
    m_imageData->GetSpacing(vs);
    m_imageData->GetOrigin(origin);
  }
  for ( int i = 0; i < 3; i++ )
  {
    index_out[i] = ( int )( ( pos[i] - origin[i] ) / vs[i] + 0.5 );
  }
}

void FSVolume::RASToNativeRAS( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  double m[16];
  vtkMatrix4x4::Invert( m_RASToRASMatrix, m );
  vtkMatrix4x4::MultiplyPoint( m, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}

void FSVolume::NativeRASToRAS( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  vtkMatrix4x4::MultiplyPoint( m_RASToRASMatrix, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}

void FSVolume::NativeRASToTkReg( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  vtkMatrix4x4::MultiplyPoint( m_RASToTkRegMatrix, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}

void FSVolume::TkRegToNativeRAS( const double* pos_in, double* pos_out )
{
  double p[4] = { pos_in[0], pos_in[1], pos_in[2], 1 };
  double m[16];
  vtkMatrix4x4::Invert( m_RASToTkRegMatrix, m );
  vtkMatrix4x4::MultiplyPoint( m, p, p );
  pos_out[0] = p[0];
  pos_out[1] = p[1];
  pos_out[2] = p[2];
}

void FSVolume::SetInterpolationMethod( int nMethod )
{
  m_nInterpolationMethod = nMethod;
}

void FSVolume::SetConform( bool bConform )
{
  m_bConform = bConform;
}

int FSVolume::GetDataType()
{
  if ( m_MRI )
  {
    return m_MRI->type;
  }
  else
  {
    return -1;
  }
}

void FSVolume::SetCroppingBounds( double* bounds )
{
  for ( int i = 0; i < 6; i++ )
  {
    m_dBounds[i] = bounds[i];
  }

  m_bCrop = true;
}

HISTOGRAM *MRIhistogramWithHighThreshold(
    MRI *mri, int nbins, HISTOGRAM *histo, float thresh, bool highThresh, int frame, float fmin_in, float fmax_in)
{
  int width, height, depth, z, x0, y0, z0, tid;
  float fmin, fmax;
#ifdef HAVE_OPENMP
  HISTOGRAM *histos[_MAX_FS_THREADS];
#else
  HISTOGRAM *histos[1];
#endif

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  x0 = y0 = z0 = 0;

  fmin = fmin_in;
  fmax = fmax_in;

  if (!nbins) nbins = nint(fmax - fmin + 1.0);

  if (!histo)
    histo = HISTOalloc(nbins);
  else {
    if (histo->nbins < nbins)
      HISTOrealloc(histo, nbins);
    else
      histo->nbins = nbins;
  }

#ifdef HAVE_OPENMP
  for (tid = 0; tid < _MAX_FS_THREADS; tid++) {
    histos[tid] = HISTOalloc(nbins);
    HISTOclear(histos[tid], histos[tid]);
    HISTOinit(histos[tid], nbins, fmin, fmax);
  }
#else
  histos[0] = histo;
#endif
  HISTOclear(histo, histo);
  HISTOinit(histo, nbins, fmin, fmax);

#ifdef HAVE_OPENMP
#pragma omp parallel for shared(histos, width, height, depth, fmin, fmax, frame, x0, y0, z0, histo)
#endif
  for (z = z0; z < depth; z++) {
    int y, x, tid;
    float val;
    for (y = y0; y < height; y++) {
      for (x = x0; x < width; x++) {
        val = MRIgetVoxVal(mri, x, y, z, frame);
        if (FZERO(val) || (highThresh && val > thresh)) continue;
//        val = MRIgetVoxVal(mri, x, y, z, frame);
#ifdef HAVE_OPENMP
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif
        HISTOaddSample(histos[tid], val, fmin, fmax);
      }
    }
  }

#ifdef HAVE_OPENMP
  for (tid = 0; tid < _MAX_FS_THREADS; tid++) {
    HISTOadd(histos[tid], histo, histo);
    HISTOfree(&histos[tid]);
  }
#endif

  return (histo);
}

void FSVolume::UpdateHistoCDF(int frame, float threshold, bool highThresh)
{
  float fMinValue, fMaxValue;
//  MRInonzeroValRange(m_MRI, &fMinValue, &fMaxValue);
  MRIvalRange(m_MRI, &fMinValue, &fMaxValue);
  if (fMinValue == fMaxValue)
  {
    m_bValidHistogram = false;
    return;
  }

  if (threshold < 0)
    threshold = fMinValue;

  /*
  HISTO* histo = HISTOinit(NULL, 1000, fMinValue, fMaxValue);

  for (int x = 0; x < m_MRI->width; x++)
    for (int y = 0; y < m_MRI->height; y++)
      for (int z = 0; z < m_MRI->depth; z++) {
        double val = MRIgetVoxVal(m_MRI, x, y, z, frame);
        if (FZERO(val) || (highThresh && val > threshold)) continue;
        HISTOaddSample(histo, val, 0, 0);
      }
      */
  HISTO* histo = MRIhistogramWithHighThreshold(m_MRI, 1000, NULL, threshold, highThresh, frame, fMinValue, fMaxValue);

  if (m_histoCDF)
    HISTOfree(&m_histoCDF);

  m_histoCDF = HISTOmakeCDF(histo, NULL);
  HISTOfree(&histo);
  if (!m_histoCDF)
  {
    cout << "Could not create HISTO" << endl;
    return;
  }

  m_bValidHistogram = true;
}

double FSVolume::GetHistoValueFromPercentile(double percentile, int frame)
{
  if (m_histoCDF)
  {
    if (m_nHistoFrame != frame)
    {
      UpdateHistoCDF(frame);
      m_nHistoFrame = frame;
    }
    int bin = HISTOfindBinWithCount(m_histoCDF, (float)percentile);
    return m_histoCDF->bins[bin];
//    return ::MRIfindPercentile(m_MRI, percentile, frame);
  }
  else
    return 0;
}

double FSVolume::GetHistoPercentileFromValue(double value, int frame)
{
  if (m_histoCDF)
  {
    if (m_nHistoFrame != frame)
    {
      UpdateHistoCDF(frame);
      m_nHistoFrame = frame;
    }
    return HISTOgetCount(m_histoCDF, (float)value);
  }
  else
    return 0;
}

void FSVolume::GetFrameValueRange(int frame, double *range)
{
  int      width, height, depth, x, y, z;
  float    fmin, fmax, *pf, val ;
  BUFTYPE  *pb ;

  MRI* mri = m_MRI;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  fmin = 1000000.0f ;
  fmax = -1000000.0f ;
  switch (mri->type)
  {
  case MRI_FLOAT:
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pf = &MRIFseq_vox(mri, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++)
        {
          val = *pf++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }
    break ;
  case MRI_INT:
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRIIseq_vox(mri, x, y, z, frame) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }
    break ;
  case MRI_SHORT:
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRISseq_vox(mri, x, y, z, frame) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }
    break ;
  case MRI_USHRT:
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRIUSseq_vox(mri, x, y, z, frame) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }
    break ;
  case MRI_UCHAR:
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pb = &MRIseq_vox(mri, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)*pb++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }
    break ;
  default:
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRIgetVoxVal(mri, x, y, z, frame) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }
    break ;
  }
  range[0] = fmin;
  range[1] = fmax;
}

bool FSVolume::Segment(int min_label_index, int max_label_index, int min_number_of_voxels)
{
  if (!UpdateMRIFromImage(m_imageData))
    return false;

  MRI_SEGMENTATION *mseg ;
  mseg = MRIsegment(m_MRITemp, min_label_index, max_label_index) ;
  MRIeraseSmallSegments(mseg, m_MRITemp, min_number_of_voxels) ;
  MRIsegmentFree(&mseg);
  MRI* temp = m_MRI;
  m_MRI = m_MRITemp;
  MapMRIToImage(true);
  m_MRI = temp;
  MRIfree(&m_MRITemp);

  return true;
}
