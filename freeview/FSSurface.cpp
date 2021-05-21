/**
 * @brief Base surface class that takes care of I/O and data conversion.
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

#include "FSSurface.h"
#include <stdexcept>
#include "vtkShortArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkIntArray.h"
#include "vtkSmartPointer.h"
#include "vtkImageReslice.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "vtkImageChangeInformation.h"
#include "vtkPolyData.h"
#include "vtkTubeFilter.h"
#include "vtkMath.h"
#include "vtkLine.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkSelectEnclosedPoints.h"
#include "vtkDelaunay3D.h"
#include "vtkUnstructuredGrid.h"
#include "FSVolume.h"
#include "MyUtils.h"
#include <QFileInfo>
#include <QDebug>
#include <QDir>

#include "mri_identify.h"

using namespace std;

FSSurface::FSSurface( FSVolume* ref, QObject* parent ) : QObject( parent ),
  m_MRIS( NULL ),
  m_MRISTarget( NULL ),
  m_bBoundsCacheDirty( true ),
  m_bCurvatureLoaded( false ),
  m_nActiveSurface( SurfaceMain ),
  m_volumeRef( ref ),
  m_nActiveVector( -1 ),
  m_bSharedMRIS(false),
  m_dMaxSegmentLength(10.0),
  m_bIgnoreVG(false)
{
  m_polydata = vtkSmartPointer<vtkPolyData>::New();
  m_polydataVector = vtkSmartPointer<vtkPolyData>::New();
  m_polydataVertices = vtkSmartPointer<vtkPolyData>::New();
  m_polydataWireframes = vtkSmartPointer<vtkPolyData>::New();
  m_polydataTarget = vtkSmartPointer<vtkPolyData>::New();

  for ( int i = 0; i < 3; i++ )
  {
    m_polydataVector2D[i] = vtkSmartPointer<vtkPolyData>::New();
    m_polydataVertex2D[i] = vtkSmartPointer<vtkPolyData>::New();
  }

  for ( int i = 0; i < NUM_OF_VSETS; i++ )
  {
    m_fVertexSets[i] = NULL;
    m_fNormalSets[i] = NULL;
    m_bSurfaceLoaded[i] = false;
    m_HashTable[i] = NULL;
  }
  m_fSmoothedNormal = NULL;

  m_targetToRasMatrix[0] = 1;
  m_targetToRasMatrix[1] = 0;
  m_targetToRasMatrix[2] = 0;
  m_targetToRasMatrix[3] = 0;
  m_targetToRasMatrix[4] = 0;
  m_targetToRasMatrix[5] = 1;
  m_targetToRasMatrix[6] = 0;
  m_targetToRasMatrix[7] = 0;
  m_targetToRasMatrix[8] = 0;
  m_targetToRasMatrix[9] = 0;
  m_targetToRasMatrix[10] = 1;
  m_targetToRasMatrix[11] = 0;
  m_targetToRasMatrix[12] = 0;
  m_targetToRasMatrix[13] = 0;
  m_targetToRasMatrix[14] = 0;
  m_targetToRasMatrix[15] = 1;
  if (ref)
  {
    MATRIX* mat = ref->GetTargetToRASMatrix();
    for ( int i = 0; i < 16; i++ )
    {
      m_targetToRasMatrix[i] = (double) *MATRIX_RELT((mat),(i/4)+1,(i%4)+1);
    }
  }
  m_targetToRasTransform = vtkSmartPointer<vtkTransform>::New();
  m_targetToRasTransform->SetMatrix( m_targetToRasMatrix );
}

FSSurface::~FSSurface()
{
  if ( m_MRIS && !m_bSharedMRIS )
  {
    ::MRISfree( &m_MRIS );
  }

  if ( m_MRISTarget )
  {
    ::MRISfree( &m_MRISTarget );
  }

  for ( int i = 0; i < NUM_OF_VSETS; i++ )
  {
    if ( m_fNormalSets[i] )
    {
      delete[] m_fNormalSets[i];
    }

    if ( m_fVertexSets[i] )
    {
      delete[] m_fVertexSets[i];
    }

    if ( m_HashTable[i] )
    {
      MHTfree( &m_HashTable[i] );
    }
  }

  if (m_fSmoothedNormal)
    delete[] m_fSmoothedNormal;

  for ( size_t i = 0; i <  m_vertexVectors.size(); i++ )
  {
    delete[] m_vertexVectors[i].data;
  }
  m_vertexVectors.clear();
}

bool FSSurface::MRISRead( const QString& filename,
                          const QString& vector_filename,
                          const QString& patch_filename,
                          const QString& target_filename,
                          const QString& sphere_filename,
                          const QStringList& sup_files )
{
  if ( m_MRIS )
  {
    ::MRISfree( &m_MRIS );
  }

//  try {
//    m_MRIS = ::MRISread( filename.toLatin1().data() );
//  }
//  catch (int ret)
//  {
//    return false;
//  }

  m_MRIS = ::MRISread( filename.toLatin1().data() );
  if ( m_MRIS == NULL )
  {
    cerr << "MRISread failed\n";
    return false;
  }
  else
  {
    if (!sphere_filename.isEmpty())
      MRISreadCanonicalCoordinates(m_MRIS, qPrintable(sphere_filename));
    return InitializeData(vector_filename, patch_filename, target_filename, sup_files);
  }
}

bool FSSurface::CreateFromMRIS(MRIS *mris)
{
  m_MRIS = mris;
  m_bSharedMRIS = true;
  return InitializeData();
}

bool FSSurface::LoadPatch(const QString &filename)
{
  if (::MRISreadPatchNoRemove(m_MRIS, filename.toLatin1().data() ) == 0 )
  {
    RipFaces();
    UpdateHashTable();
    UpdatePolyData();
    SaveVertices( m_MRIS, SurfaceMain );
    SaveNormals ( m_MRIS, SurfaceMain );
    RestoreVertices( m_MRIS, SurfaceMain );
    RestoreNormals( m_MRIS, SurfaceMain );
    return true;
  }
  else
    return false;
}

bool FSSurface::InitializeData(const QString &vector_filename,
                               const QString &patch_filename,
                               const QString &target_filename,
                               const QStringList &sup_files)
{
  // backup ripflags
  m_originalRipflags.clear();
  for (int i = 0; i < m_MRIS->nvertices; i++)
    m_originalRipflags << m_MRIS->vertices[i].ripflag;

  if ( !patch_filename.isEmpty() )
  {
    if ( ::MRISreadPatchNoRemove(m_MRIS, patch_filename.toLatin1().data() ) != 0 )
    {
      cerr << "Can not load patch file " << qPrintable(patch_filename) << "\n";
    }
  }

  // Get some info from the MRIS. This can either come from the volume
  // geometry data embedded in the surface; this is done for newer
  // surfaces. Or it can come from the source information in the
  // transform. We use it to get the RAS center offset for the
  // surface->RAS transform.
  m_SurfaceToRASMatrix[0] = 1;
  m_SurfaceToRASMatrix[1] = 0;
  m_SurfaceToRASMatrix[2] = 0;
  m_SurfaceToRASMatrix[3] = 0;
  m_SurfaceToRASMatrix[4] = 0;
  m_SurfaceToRASMatrix[5] = 1;
  m_SurfaceToRASMatrix[6] = 0;
  m_SurfaceToRASMatrix[7] = 0;
  m_SurfaceToRASMatrix[8] = 0;
  m_SurfaceToRASMatrix[9] = 0;
  m_SurfaceToRASMatrix[10] = 1;
  m_SurfaceToRASMatrix[11] = 0;
  m_SurfaceToRASMatrix[12] = 0;
  m_SurfaceToRASMatrix[13] = 0;
  m_SurfaceToRASMatrix[14] = 0;
  m_SurfaceToRASMatrix[15] = 1;

  m_bValidVolumeGeometry = false;
  if ( m_MRIS->vg.valid && !m_bIgnoreVG )
  {
    MRI* tmp = MRIallocHeader(m_MRIS->vg.width, m_MRIS->vg.height, m_MRIS->vg.depth, MRI_UCHAR, 1);
    useVolGeomToMRI(&m_MRIS->vg, tmp);
    MATRIX* vox2rasScanner = MRIxfmCRS2XYZ(tmp, 0);
    MATRIX* vo2rasTkReg = MRIxfmCRS2XYZtkreg(tmp);
    if (vo2rasTkReg)
    {
      MATRIX* vox2rasTkReg_inv = MatrixInverse( vo2rasTkReg, NULL );
      if (vox2rasTkReg_inv)
      {
        MATRIX* M = MatrixMultiply( vox2rasScanner, vox2rasTkReg_inv, NULL );
        for ( int i = 0; i < 16; i++ )
        {
          m_SurfaceToRASMatrix[i] =
              (double) *MATRIX_RELT( M, (i/4)+1, (i%4)+1 );
        }
        MatrixFree( &vox2rasTkReg_inv );
        MatrixFree( &M );
        m_bValidVolumeGeometry = true;
      }
      else
        cout << "Warning: MatrixInverse failed." << endl;
    }

    MRIfree( &tmp );
    MatrixFree( &vox2rasScanner );
    MatrixFree( &vo2rasTkReg );
  }
  else if ( !m_MRIS->useRealRAS && m_volumeRef )
  {
    MATRIX* m = RASFromSurfaceRAS_( m_volumeRef->GetMRI());
    for ( int i = 0; i < 16; i++ )
    {
      m_SurfaceToRASMatrix[i] =
          (double) *MATRIX_RELT( m, (i/4)+1, (i%4)+1 );
    }
    MatrixFree(&m);
  }

  // Make our transform object and set the matrix.
  m_SurfaceToRASTransform = vtkSmartPointer<vtkTransform>::New();
  m_SurfaceToRASTransform->SetMatrix( m_SurfaceToRASMatrix );

  // Make the hash table. This makes it with v->x,y,z.
  UpdateHashTable(0, CURRENT_VERTICES);
  UpdatePolyData();

  // Save main vertices and normals
  SaveVertices( m_MRIS, SurfaceMain );
  SaveNormals ( m_MRIS, SurfaceMain );
  m_bSurfaceLoaded[SurfaceMain] = true;

  //  if ( patch_filename.isEmpty() )
  {
    if (sup_files.contains("white"))
      LoadSurface ( "white",    SurfaceWhite );
    if (sup_files.contains("pial"))
      LoadSurface ( "pial",     SurfacePial );
    if (sup_files.contains("orig"))
      LoadSurface ( "orig",     SurfaceOriginal );
    if (sup_files.contains("inflated"))
      LoadSurface ( "inflated", SurfaceInflated );
    if (sup_files.contains("orig.nofix"))
      LoadSurface ( "orig.nofix",    SurfaceWhite );
  }

  RestoreVertices( m_MRIS, SurfaceMain );
  RestoreNormals( m_MRIS, SurfaceMain );

  if ( !target_filename.isEmpty() )
  {
    LoadTargetSurface( target_filename );
  }

  if ( !vector_filename.isEmpty() )
  {
    LoadVectors ( vector_filename );
  }

  UpdateSmoothedNormals();

  QFileInfo fi(m_MRIS->fname.data());
  if (QFileInfo(fi.absoluteDir(), fi.completeBaseName() + ".curv").exists())
    LoadCurvature();

  return true;
}

vtkTransform* FSSurface::GetSurfaceToRasTransform()
{
    return m_SurfaceToRASTransform;
}

void FSSurface::UpdateHashTable(int nSet, int coord)
{
  if (m_HashTable[nSet])
    MHTfree( &m_HashTable[nSet] );
  double max_spacing;
  int max_vno;
  MRIScomputeVertexSpacingStats(m_MRIS, NULL, NULL, &max_spacing, NULL, &max_vno, coord);
  m_HashTable[nSet] = MHTcreateVertexTable_Resolution( m_MRIS, coord, max_spacing/2 );
}

void FSSurface::LoadTargetSurface( const QString& filename )
{
  if ( m_MRISTarget )
  {
    ::MRISfree( &m_MRISTarget );
  }

  try
  {
    m_MRISTarget = ::MRISread( filename.toLatin1().data() );
  }
  catch (int ret)
  {
    m_MRISTarget = NULL;
  }

  if ( m_MRISTarget == NULL )
  {
    cerr << "MRISread failed. Can not load target surface.\n";
    return;
  }

  UpdatePolyData( m_MRISTarget, m_polydataTarget );
}

bool FSSurface::MRISWrite( const QString& filename )
{
  if ( m_MRIS == NULL )
  {
    cerr << "No MRIS to write.\n";
    return false;
  }

  int ret = 0;
  try
  {
    ret = ::MRISwrite( m_MRIS, filename.toLatin1().data() );
  }
  catch (int ret)
  {
    return false;
  }

  return ( ret == 0 );
}

bool FSSurface::MRISReadVectors( const QString& filename )
{
  return LoadVectors( filename );
}

bool FSSurface::LoadSurface( const QString& filename, int nSet )
{
  if ( ::MRISreadVertexPositions( m_MRIS, filename.toLatin1().data() ) != 0 )
  {
    cerr << "could not load surface from " << qPrintable(filename) << "\n";
    m_bSurfaceLoaded[nSet] = false;
    return false;
  }
  else
  {
    if (filename == "white")
      ::MRISsaveVertexPositions(m_MRIS, WHITE_VERTICES);

    UpdateHashTable(nSet, (filename == "white" ? WHITE_VERTICES:CURRENT_VERTICES));
    ComputeNormals();
    SaveVertices( m_MRIS, nSet );
    SaveNormals ( m_MRIS, nSet );
    m_bSurfaceLoaded[nSet] = true;
    return true;
  }
}

bool FSSurface::LoadCurvature( const QString& filename )
{
  if ( ::MRISreadCurvatureFile( m_MRIS, (char*)(filename.isEmpty() ? "curv" : filename.toLatin1().data() ) ) != 0 )
  {
    cerr << "could not read curvature from " << qPrintable(filename.isEmpty() ? "curv" : filename) << "\n";
    m_bCurvatureLoaded = false;
    return false;
  }
  else
  {
    int cVertices = m_MRIS->nvertices;

    vtkSmartPointer<vtkFloatArray> curvs = vtkSmartPointer<vtkFloatArray>::New();
    curvs->Allocate( cVertices );
    curvs->SetNumberOfComponents( 1 );
    curvs->SetName( "Curvature" );

    for ( int vno = 0; vno < cVertices; vno++ )
    {
      curvs->InsertNextValue( m_MRIS->vertices[vno].curv );
    }
    m_polydata->GetPointData()->SetScalars( curvs );
    m_polydataWireframes->GetPointData()->SetScalars( curvs );
    m_bCurvatureLoaded = true;

    return true;
  }
}

bool FSSurface::LoadOverlay( const QString& filename, const QString& fn_reg,
                             float** data_out, int* nvertices_out, int* nframes_out,
                             bool bUseSecondHalfData )
{
  //    int mritype = mri_identify((char*)( filename.toLatin1().data() ));
  //    qDebug() << "mritype " << mritype;
  MRI* mriheader = MRIreadHeader(filename.toLatin1().data(), MRI_VOLUME_TYPE_UNKNOWN);

  if (mriheader && mriheader->width*mriheader->height*mriheader->depth != m_MRIS->nvertices &&
      mriheader->width*mriheader->height*mriheader->depth*mriheader->nframes != m_MRIS->nvertices )
  {
    if (mriheader->height == 1 && mriheader->depth == 1)
    {
      // likely wrong overlay data
      printf("Number of vertices in overlay data (%d or %d) does not match with surface (%d).\n",
             mriheader->width*mriheader->height*mriheader->depth,
             mriheader->width*mriheader->height*mriheader->depth*mriheader->nframes,
             m_MRIS->nvertices);
      return false;
    }

    // try load as volume
    MRI* mri = MRIread(filename.toLatin1().data());
    if (!mri)
    {
      cerr << "could not read overlay data from " << qPrintable(filename) << "\n";
      return false;
    }

    // if there is registration file, read it
    MATRIX* tkregMat = MRIxfmCRS2XYZtkreg(mri);
    MATRIX* ras2vox_tkreg = MatrixInverse(tkregMat, NULL);
    MatrixFree(&tkregMat);
    if (!fn_reg.isEmpty())
    {
      MRI* tmp = MRIallocHeader(m_MRIS->vg.width, m_MRIS->vg.height, m_MRIS->vg.depth, MRI_UCHAR, 1);
      useVolGeomToMRI(&m_MRIS->vg, tmp);
      MATRIX* matReg = FSVolume::LoadRegistrationMatrix(fn_reg, tmp, mri);
      if (matReg == NULL)
      {
        cerr << "could not read registration data from " << qPrintable(fn_reg) << "\n";
      }
      else
      {
        ras2vox_tkreg = MatrixMultiply(ras2vox_tkreg, matReg, NULL);
      }
      MatrixFree(&matReg);
      MRIfree(&tmp);
    }
    double m[16];
    for ( int i = 0; i < 16; i++ )
    {
      m[i] = (double) *MATRIX_RELT((ras2vox_tkreg),(i/4)+1,(i%4)+1);
    }

    float* data = new float[mri->nframes*m_MRIS->nvertices];
    if (!data)
      return false;
    for ( int i = 0; i < m_MRIS->nvertices; i++ )
    {
      double v[4] = { m_MRIS->vertices[i].x, m_MRIS->vertices[i].y, m_MRIS->vertices[i].z, 1 };
      vtkMatrix4x4::MultiplyPoint(m, v, v);
      int nx = (int)(v[0]+0.5);
      int ny = (int)(v[1]+0.5);
      int nz = (int)(v[2]+0.5);
      if (nx >= 0 && nx < mri->width && ny >= 0 && ny < mri->height && nz >= 0 && nz < mri->depth)
      {
        for (int j = 0; j < mri->nframes; j++)
          data[j*m_MRIS->nvertices + i] = ::MRIgetVoxVal(mri, nx, ny, nz, j);
      }
    }

    *data_out = data;
    *nframes_out = mri->nframes;
    *nvertices_out = m_MRIS->nvertices;
    MRIfree(&mri);
    MatrixFree(&ras2vox_tkreg);
  }
  else
  {
    MRI* mri = MRIread(filename.toLatin1().data());
    if (!mri)
    {
      cerr << "could not read overlay data from " << qPrintable(filename) << "\n";
      return false;
    }
    int nPerFrame = mri->width*mri->height*mri->depth;
    int nframes = nPerFrame*mri->nframes / m_MRIS->nvertices;
    float* data = new float[nPerFrame*mri->nframes];
    for (int nx = 0; nx < mri->width; nx++)
    {
      for (int ny = 0; ny < mri->height; ny++)
      {
        for (int nz = 0; nz < mri->depth; nz++)
        {
          for (int nk = 0; nk < mri->nframes; nk++)
          {
            data[nk*nPerFrame + nz*mri->height*mri->width + ny*mri->width + nx]
                = ::MRIgetVoxVal(mri, nx, ny, nz, nk);
          }
        }
      }
    }
    if (bUseSecondHalfData)
    {
      float* data2 = new float[m_MRIS->nvertices];
      memcpy(data2, data+m_MRIS->nvertices, sizeof(float)*m_MRIS->nvertices);
      delete[] data;
      *data_out = data2;
      *nframes_out = 1;
    }
    else
    {
      *data_out = data;
      *nframes_out = nframes;
    }
    *nvertices_out = m_MRIS->nvertices;
  }

  return true;
}

/*
bool FSSurface::LoadOverlay( const QString& filename )
{
  // user read curvature routine because read values does not handle filename properly
  int nCount = m_MRIS->nvertices;
  for ( int i = 0; i < nCount; i++ )
    m_MRIS->vertices[i].val = m_MRIS->vertices[i].curv;
  float fMin = m_MRIS->min_curv;
  float fMax = m_MRIS->max_curv;

  if ( ::MRISreadCurvatureFile( m_MRIS, (char*)( filename.toLatin1().data()) ) != 0 )
  {
    cerr << "could not read overlay data from " << filename.toLatin1().data() << "\n";
    return false;
  }
  else
  {
    // move curv to val and restore real curv
    float ftemp;
    for ( int i = 0; i < nCount; i++ )
    {
      ftemp = m_MRIS->vertices[i].curv;
      m_MRIS->vertices[i].curv = m_MRIS->vertices[i].val;
      m_MRIS->vertices[i].val = ftemp;
    }
    m_MRIS->min_curv = fMin;
    m_MRIS->max_curv = fMax;

    return true;
  }
}
*/

bool FSSurface::LoadVectors( const QString& filename )
{
  MRI* mri = ::MRIread( filename.toLatin1().data() );

  VertexVectorItem vector;
  vector.name = QFileInfo( filename ).fileName();
  if ( mri )
  {
    if ( mri->type != MRI_FLOAT )
    {
      cerr << "Vector volume must be in type of MRI_FLOAT.\n";
    }
    else if ( SaveVertices( mri, vector.data ) )
    {
      m_vertexVectors.push_back( vector );
      m_nActiveVector = m_vertexVectors.size() - 1;
      UpdateVectors();
      ::MRIfree( &mri );
      cout << "vector data loaded for surface from " << qPrintable(filename) << "\n";
      return true;
    }
  }
  else if ( ::MRISreadVertexPositions( m_MRIS, (char*)filename.toLatin1().data() ) == 0 )
  {
    // compute vectors
    if ( ComputeVectors( m_MRIS, vector.data ) )
    {
      m_vertexVectors.push_back( vector );
      m_nActiveVector = m_vertexVectors.size() - 1;

      // restore original vertices in m_MRIS because MRISreadVertexPositions changed it!
      RestoreVertices( m_MRIS, m_nActiveSurface );

      UpdateVectors();

      cout << "vector data loaded for surface from " << qPrintable(filename) << "\n";
      return true;
    }
  }

  if ( mri )
  {
    ::MRIfree( &mri );
  }
  return false;
}

bool FSSurface::ComputeVectors( MRIS* mris, VertexItem*& buffer )
{
  int nvertices = mris->nvertices;
  VERTEX *v;
  for ( int vno = 0; vno < nvertices; vno++ )
  {
    v = &mris->vertices[vno];
    v->x -= m_fVertexSets[m_nActiveSurface][vno].x;
    v->y -= m_fVertexSets[m_nActiveSurface][vno].y;
    v->z -= m_fVertexSets[m_nActiveSurface][vno].z;
  }
  return SaveVertices( mris, buffer );
}

void FSSurface::SaveVertices( MRIS* mris, int nSet )
{
  if ( !mris || nSet >= NUM_OF_VSETS )
  {
    return;
  }

  SaveVertices( mris, m_fVertexSets[nSet] );
}

bool FSSurface::SaveVertices( MRIS* mris, VertexItem*& buffer )
{
  int nvertices = mris->nvertices;
  VERTEX *v;

  if ( buffer == NULL )
  {
    buffer = new VertexItem[nvertices];
    if ( !buffer )
    {
      cerr << "Can not allocate memory for vertex sets.\n";
      return false;
    }
  }
  for ( int vno = 0; vno < nvertices; vno++ )
  {
    v = &mris->vertices[vno];
    buffer[vno].x = v->x;
    buffer[vno].y = v->y;
    buffer[vno].z = v->z;
  }
  return true;
}


bool FSSurface::SaveVertices( MRI* mri, VertexItem*& buffer )
{
  int nvertices = mri->width;

  if ( buffer == NULL )
  {
    buffer = new VertexItem[nvertices];
    if ( !buffer )
    {
      cerr << "Can not allocate memory for vertex sets.\n";
      return false;
    }
  }
  for ( int vno = 0; vno < nvertices; vno++ )
  {
    buffer[vno].x = MRIFseq_vox( mri, vno, 0, 0, 0 );
    buffer[vno].y = MRIFseq_vox( mri, vno, 0, 0, 1 );
    buffer[vno].z = MRIFseq_vox( mri, vno, 0, 0, 2 );
  }
  return true;
}

void FSSurface::RestoreVertices( MRIS* mris, int nSet )
{
  if ( !mris || nSet >= NUM_OF_VSETS || m_fVertexSets[nSet] == NULL)
  {
    return;
  }

  int nvertices = mris->nvertices;
  VERTEX *v;

  for ( int vno = 0; vno < nvertices; vno++ )
  {
    v = &mris->vertices[vno];
    v->x = m_fVertexSets[nSet][vno].x;
    v->y = m_fVertexSets[nSet][vno].y;
    v->z = m_fVertexSets[nSet][vno].z;
  }
}

void FSSurface::SaveNormals( MRIS* mris, int nSet )
{
  if ( !mris || nSet >= NUM_OF_VSETS )
  {
    return;
  }

  int nvertices = mris->nvertices;
  VERTEX *v;

  if ( m_fNormalSets[nSet] == NULL )
  {
    m_fNormalSets[nSet] = new VertexItem[nvertices];
    if ( !m_fNormalSets[nSet] )
    {
      cerr << "Can not allocate memory for normal sets.\n";
      return;
    }
  }
  for ( int vno = 0; vno < nvertices; vno++ )
  {
    v = &mris->vertices[vno];
    m_fNormalSets[nSet][vno].x = v->nx;
    m_fNormalSets[nSet][vno].y = v->ny;
    m_fNormalSets[nSet][vno].z = v->nz;
  }
}

void FSSurface::RestoreNormals( MRIS* mris, int nSet )
{
  if ( !mris || nSet >= NUM_OF_VSETS || m_fNormalSets[nSet] == NULL)
  {
    return;
  }

  int nvertices = mris->nvertices;
  VERTEX *v;

  for ( int vno = 0; vno < nvertices; vno++ )
  {
    v = &mris->vertices[vno];
    v->nx = m_fNormalSets[nSet][vno].x;
    v->ny = m_fNormalSets[nSet][vno].y;
    v->nz = m_fNormalSets[nSet][vno].z;
  }
}

void FSSurface::GetNormalAtVertex( int nVertex, double* vec_out )
{
  vec_out[0] = m_fNormalSets[m_nActiveSurface][nVertex].x;
  vec_out[1] = m_fNormalSets[m_nActiveSurface][nVertex].y;
  vec_out[2] = m_fNormalSets[m_nActiveSurface][nVertex].z;
}

void FSSurface::UpdatePolyData()
{
  UpdatePolyData( m_MRIS, m_polydata, m_polydataVertices, m_polydataWireframes, true );
}

void FSSurface::UpdatePolyData( MRIS* mris,
                                vtkPolyData* polydata,
                                vtkPolyData* polydata_verts,
                                vtkPolyData* polydata_wireframe,
                                bool create_segs)
{
  // Allocate all our arrays.
  int cVertices = mris->nvertices;
  int cFaces = mris->nfaces;

  vtkSmartPointer<vtkPoints> newPoints =
      vtkSmartPointer<vtkPoints>::New();
  newPoints->Allocate( cVertices );

  vtkSmartPointer<vtkCellArray> newPolys =
      vtkSmartPointer<vtkCellArray>::New();
  newPolys->Allocate( newPolys->EstimateSize(cFaces,VERTICES_PER_FACE) );

  vtkSmartPointer<vtkFloatArray> newNormals =
      vtkSmartPointer<vtkFloatArray>::New();
  newNormals->Allocate( cVertices );
  newNormals->SetNumberOfComponents( 3 );
  newNormals->SetName( "Normals" );;

  vtkSmartPointer<vtkCellArray> verts;
  if ( polydata_verts )
  {
    verts = vtkSmartPointer<vtkCellArray>::New();
    verts->Allocate( cVertices );
  }

  // Go through the surface and copy the vertex and normal for each
  // vertex. We need to transform them from surface RAS into standard
  // RAS.
  bool bInvertNormal = (vtkMatrix4x4::Determinant(m_SurfaceToRASMatrix) < 0);
  float point[3], normal[3], surfaceRAS[3];
  for ( int vno = 0; vno < cVertices; vno++ )
  {
    surfaceRAS[0] = mris->vertices[vno].x;
    surfaceRAS[1] = mris->vertices[vno].y;
    surfaceRAS[2] = mris->vertices[vno].z;
    this->ConvertSurfaceToRAS( surfaceRAS, point );
    m_targetToRasTransform->GetInverse()->TransformPoint(point, point);
    newPoints->InsertNextPoint( point );

    normal[0] = mris->vertices[vno].nx;
    normal[1] = mris->vertices[vno].ny;
    normal[2] = mris->vertices[vno].nz;
    float orig[3] = { 0, 0, 0 };
    this->ConvertSurfaceToRAS( orig, orig );
    this->ConvertSurfaceToRAS( normal, normal );
    m_targetToRasTransform->GetInverse()->TransformPoint(orig, orig);
    m_targetToRasTransform->GetInverse()->TransformPoint(normal, normal);

    for (int i = 0; i < 3; i++)
      normal[i] = bInvertNormal?(orig[i] - normal[i]):(normal[i] - orig[i]);
    vtkMath::Normalize(normal);
    newNormals->InsertNextTuple( normal );

    if ( polydata_verts )
    {
      vtkIdType n = vno;
      verts->InsertNextCell( 1, &n );
    }
  }

  // Go through and add the face indices.
  vtkSmartPointer<vtkCellArray> lines;
  if ( polydata_wireframe )
  {
    lines = vtkSmartPointer<vtkCellArray>::New();
    lines->Allocate( cFaces * 6 );
  }
  vtkIdType face[VERTICES_PER_FACE];
  float dx = 0, dy = 0, dz = 0;
  for ( int fno = 0; fno < cFaces; fno++ )
  {
    if ( mris->faces[fno].ripflag == 0 )
    {
      face[0] = mris->faces[fno].v[0];
      face[1] = mris->faces[fno].v[1];
      face[2] = mris->faces[fno].v[2];
      newPolys->InsertNextCell( 3, face );
      if ( polydata_wireframe )
      {
        lines->InsertNextCell( 2, face );
        lines->InsertNextCell( 2, face+1 );
        vtkIdType t[2] = { face[0], face[2] };
        lines->InsertNextCell( 2, t );
      }
      if (create_segs)
      {
        dx = qMax(dx, qAbs(mris->vertices[face[0]].x - mris->vertices[face[1]].x));
        dx = qMax(dx, qAbs(mris->vertices[face[1]].x - mris->vertices[face[2]].x));
        dx = qMax(dx, qAbs(mris->vertices[face[0]].x - mris->vertices[face[2]].x));
        dy = qMax(dy, qAbs(mris->vertices[face[0]].y - mris->vertices[face[1]].y));
        dy = qMax(dy, qAbs(mris->vertices[face[1]].y - mris->vertices[face[2]].y));
        dy = qMax(dy, qAbs(mris->vertices[face[0]].y - mris->vertices[face[2]].y));
        dz = qMax(dz, qAbs(mris->vertices[face[0]].z - mris->vertices[face[1]].z));
        dz = qMax(dz, qAbs(mris->vertices[face[1]].z - mris->vertices[face[2]].z));
        dz = qMax(dz, qAbs(mris->vertices[face[0]].z - mris->vertices[face[2]].z));
      }
    }
  }
  if (create_segs)
  {
    m_dMaxSegmentLength = sqrt(dx*dx+dy*dy+dz*dz);
    //  qDebug() << m_dMaxSegmentLength;
  }

  polydata->SetPoints( newPoints );
  polydata->GetPointData()->SetNormals( newNormals );
  newPolys->Squeeze(); // since we've estimated size; reclaim some space
  polydata->SetPolys( newPolys );
  if ( polydata_verts )
  {
    polydata_verts->SetPoints( newPoints );
    polydata_verts->SetVerts( verts );
  }
  if ( polydata_wireframe )
  {
    polydata_wireframe->SetPoints( newPoints );
    polydata_wireframe->SetLines( lines );
  }
}

void FSSurface::UpdateSmoothedNormals()
{
  if ( m_fSmoothedNormal == NULL )
  {
    m_fSmoothedNormal = new VertexItem[m_MRIS->nvertices];
    if ( !m_fSmoothedNormal )
    {
      cerr << "Can not allocate memory for normal sets.\n";
      return;
    }
  }
  MRISsmoothSurfaceNormals(m_MRIS, 50);
  int nvertices = m_MRIS->nvertices;
  float normal[3];
  for ( int vno = 0; vno < nvertices; vno++ )
  {
    normal[0] = m_MRIS->vertices[vno].nx;
    normal[1] = m_MRIS->vertices[vno].ny;
    normal[2] = m_MRIS->vertices[vno].nz;
    float orig[3] = { 0, 0, 0 };
    this->ConvertSurfaceToRAS( orig, orig );
    this->ConvertSurfaceToRAS( normal, normal );
    m_targetToRasTransform->GetInverse()->TransformPoint(orig, orig);
    m_targetToRasTransform->GetInverse()->TransformPoint(normal, normal);

    for (int i = 0; i < 3; i++)
      normal[i] = normal[i] - orig[i];
    vtkMath::Normalize(normal);
    m_fSmoothedNormal[vno].x = normal[0];
    m_fSmoothedNormal[vno].y = normal[1];
    m_fSmoothedNormal[vno].z = normal[2];
  }
  RestoreNormals(m_MRIS, SurfaceMain);
}

void FSSurface::UpdateVerticesAndNormals()
{
  // Allocate all our arrays.
  int cVertices = m_MRIS->nvertices;

  vtkSmartPointer<vtkPoints> newPoints =
      vtkSmartPointer<vtkPoints>::New();
  newPoints->Allocate( cVertices );

  vtkSmartPointer<vtkFloatArray> newNormals =
      vtkSmartPointer<vtkFloatArray>::New();
  newNormals->Allocate( cVertices );
  newNormals->SetNumberOfComponents( 3 );
  newNormals->SetName( "Normals" );

  // Go through the surface and copy the vertex and normal for each
  // vertex. We have to transform them from surface RAS into normal
  // RAS.
  float point[3], normal[3], surfaceRAS[3];
  for ( int vno = 0; vno < cVertices; vno++ )
  {
    surfaceRAS[0] = m_MRIS->vertices[vno].x;
    surfaceRAS[1] = m_MRIS->vertices[vno].y;
    surfaceRAS[2] = m_MRIS->vertices[vno].z;
    this->ConvertSurfaceToRAS( surfaceRAS, point );
    /*
    if ( m_volumeRef )
    {
      m_volumeRef->RASToTarget( point, point );
    }
    */
    m_targetToRasTransform->GetInverse()->TransformPoint(point, point);
    newPoints->InsertNextPoint( point );

    normal[0] = m_MRIS->vertices[vno].nx;
    normal[1] = m_MRIS->vertices[vno].ny;
    normal[2] = m_MRIS->vertices[vno].nz;
    newNormals->InsertNextTuple( normal );
  }

  m_polydata->SetPoints( newPoints );
  m_polydata->GetPointData()->SetNormals( newNormals );
  m_polydataVertices->SetPoints( newPoints );
  m_polydataWireframes->SetPoints( newPoints );
  //  m_polydata->Update();

  // if vector data exist
  UpdateVectors();
}

void FSSurface::UpdateVectors()
{
  if ( HasVectorSet() && m_nActiveVector >= 0 )
  {
    VertexItem* vectors = m_vertexVectors[m_nActiveVector].data;
    int cVertices = m_MRIS->nvertices;
    vtkPoints* oldPoints = m_polydata->GetPoints();
    vtkFloatArray* normals = vtkFloatArray::SafeDownCast(m_polydata->GetPointData()->GetNormals());
    float point[3], surfaceRAS[3];
    vtkIdType n = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> testPoints = vtkSmartPointer<vtkPoints>::New();
    double v_normal[3], v[3];
    for ( int vno = 0; vno < cVertices; vno++ )
    {
      if ( vectors[vno].x != 0 || vectors[vno].y != 0 || vectors[vno].z != 0 )
      {
        surfaceRAS[0] = m_fVertexSets[m_nActiveSurface][vno].x + vectors[vno].x;
        surfaceRAS[1] = m_fVertexSets[m_nActiveSurface][vno].y + vectors[vno].y;
        surfaceRAS[2] = m_fVertexSets[m_nActiveSurface][vno].z + vectors[vno].z;
        this->ConvertSurfaceToRAS( surfaceRAS, point );
        /*
        if ( m_volumeRef )
        {
          m_volumeRef->RASToTarget( point, point );
        }
        */
        m_targetToRasTransform->GetInverse()->TransformPoint(point, point);
        double* p0 = oldPoints->GetPoint( vno );
        if (normals)
        {
          // flip the vector if it is pointing inward
          normals->GetTuple(vno, v_normal);
          v[0] = point[0] - p0[0];
          v[1] = point[1] - p0[1];
          v[2] = point[2] - p0[2];
          vtkMath::Normalize(v);
          if (vtkMath::Dot(v_normal, v) < 0)
          {
            point[0] = p0[0]-(point[0]-p0[0]);
            point[1] = p0[1]-(point[1]-p0[1]);
            point[2] = p0[2]-(point[2]-p0[2]);
          }
        }
        points->InsertNextPoint( p0 );
        points->InsertNextPoint( point );

        /*
        testPoints->InsertNextPoint(p0[0]+(point[0]-p0[0])/5,
                                    p0[1]+(point[1]-p0[1])/5,
                                    p0[2]+(point[2]-p0[2])/5);
        */

        verts->InsertNextCell( 1, &n );

        lines->InsertNextCell( 2 );
        lines->InsertCellPoint( n++ );
        lines->InsertCellPoint( n++ );
      }
    }
    /*
    vtkSmartPointer<vtkPolyData> testPolydata = vtkSmartPointer<vtkPolyData>::New();
    testPolydata->SetPoints(testPoints);
    vtkSmartPointer<vtkDelaunay3D> selectEnclosedPoints = vtkSmartPointer<vtkDelaunay3D>::New();
    selectEnclosedPoints->SetInput(testPolydata);
    selectEnclosedPoints->Update();

    double p0[3], p1[3];
    double pcoords[3], weights[3];
    vtkIdType cellId;
    int subId;
    for(vtkIdType i = 0; i < testPoints->GetNumberOfPoints(); i++)
    {
        cellId = selectEnclosedPoints->GetOutput()->FindCell(testPoints->GetPoint(i), NULL, 0, .0001,
                                            subId, pcoords, weights);
        if (cellId >= 0)
        {
            points->GetPoint(i*2, p0);
            points->GetPoint(i*2+1, p1);
            points->SetPoint(i*2+1, p0[0]-(p1[0]-p0[0]),
                             p0[1]-(p1[1]-p0[1]),
                             p0[2]-(p1[2]-p0[2]));
        }
    }
    */

    m_polydataVector->SetPoints( points );
    m_polydataVector->SetLines( lines );
    m_polydataVector->SetVerts( verts );
  }
}

void FSSurface::UpdateVector2D( int nPlane, double slice_pos, vtkPolyData* contour_polydata )
{
  if ( HasVectorSet() && m_nActiveVector >= 0 )
  {
    VertexItem* vectors = m_vertexVectors[m_nActiveVector].data;
    int cVertices = m_MRIS->nvertices;
    int cFaces = m_MRIS->nfaces;

    // first figure out what vertices crossing the plane
    unsigned char* mask = new unsigned char[cVertices];
    memset( mask, 0, cVertices );
    vtkPoints* oldPoints = m_polydata->GetPoints();
    double pt_a[3], pt_b[3];
    for ( int fno = 0; fno < cFaces; fno++ )
    {
      if ( m_MRIS->faces[fno].ripflag == 0 )
      {
        const int * np = m_MRIS->faces[fno].v;
        vtkIdType lines[3][2] = { {np[0], np[1]}, {np[1], np[2]}, {np[2], np[0]} };
        for ( int i = 0; i < 3; i++ )
        {
          oldPoints->GetPoint( lines[i][0], pt_a );
          oldPoints->GetPoint( lines[i][1], pt_b );
          if ( (pt_a[nPlane] >= slice_pos && pt_b[nPlane] <= slice_pos) ||
               (pt_a[nPlane] <= slice_pos && pt_b[nPlane] >= slice_pos) )
          {
            mask[lines[i][0]] = 1;
            mask[lines[i][1]] = 1;
          }
        }
      }
    }

    // build vector actor
    vtkIdType n = 0;
    double point[3], surfaceRAS[3];
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
    double old_pt[3];
    vtkPoints* contour_pts = NULL;
    vtkCellArray* contour_lines = NULL;
    if ( contour_polydata )
    {
      contour_pts = contour_polydata->GetPoints();
      contour_lines = contour_polydata->GetLines();
    }

    vtkPolyData* target_polydata = NULL;
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    if ( m_polydataTarget->GetPoints() && m_polydataTarget->GetPoints()->GetNumberOfPoints() > 0 )
    {
      vtkSmartPointer<vtkPlane> slicer = vtkSmartPointer<vtkPlane>::New();
      double pos[3] = { 0, 0, 0 };
      pos[nPlane] = slice_pos;
      slicer->SetOrigin( pos );
      slicer->SetNormal( (nPlane==0), (nPlane==1), (nPlane==2) );
#if VTK_MAJOR_VERSION > 5
      cutter->SetInputData( m_polydataTarget );
#else
      cutter->SetInput( m_polydataTarget );
#endif
      cutter->SetCutFunction( slicer );
      cutter->Update();
      target_polydata = cutter->GetOutput();
    }

    for ( int vno = 0; vno < cVertices; vno++ )
    {
      if ( mask[vno] )
      {
        oldPoints->GetPoint( vno, old_pt );
        surfaceRAS[0] = m_fVertexSets[m_nActiveSurface][vno].x + vectors[vno].x;
        surfaceRAS[1] = m_fVertexSets[m_nActiveSurface][vno].y + vectors[vno].y;
        surfaceRAS[2] = m_fVertexSets[m_nActiveSurface][vno].z + vectors[vno].z;
        this->ConvertSurfaceToRAS( surfaceRAS, point );
        /*
        if ( m_volumeRef )
        {
          m_volumeRef->RASToTarget( point, point );
        }
        */
        m_targetToRasTransform->GetInverse()->TransformPoint(point, point);

        if ( contour_pts )
        {
          double new_pt[3];
          ProjectVectorPoint2D( old_pt, contour_pts, contour_lines, new_pt );

          for ( int i = 0; i < 3; i++ )
          {
            point[i] += (new_pt[i] - old_pt[i] );
          }
          points->InsertNextPoint( new_pt );

          if ( target_polydata && target_polydata->GetPoints() )
          {
            ProjectVectorPoint2D( point, target_polydata->GetPoints(), target_polydata->GetLines(), point );
          }
          points->InsertNextPoint( point );
        }
        else
        {
          points->InsertNextPoint( old_pt );
          points->InsertNextPoint( point );
        }

        verts->InsertNextCell( 1, &n );
        lines->InsertNextCell( 2 );
        lines->InsertCellPoint( n++ );
        lines->InsertCellPoint( n++ );
      }
    }
    m_polydataVector2D[nPlane]->SetPoints( points );
    m_polydataVector2D[nPlane]->SetLines( lines );
    m_polydataVector2D[nPlane]->SetVerts( verts );
    delete[] mask;
  }
}

bool FSSurface::ProjectVectorPoint2D( double* pt_in,
                                      vtkPoints* contour_pts,
                                      vtkCellArray* contour_lines,
                                      double* pt_out )
{
  // first find the closest point on the contour
  double pt0[3];
  double dist2 = 1e10;
  int n0 = 0;
  double* old_pt = pt_in;

  //  cout << contour_pts << " " << contour_lines << "\n"; fflush(0);
  //  cout << contour_pts->GetNumberOfPoints() << " " <<
  //      contour_lines->GetNumberOfCells() << "\n";
  for ( int i = 0; i < contour_pts->GetNumberOfPoints(); i++ )
  {
    contour_pts->GetPoint( i, pt0 );
    double d = vtkMath::Distance2BetweenPoints( old_pt, pt0 );
    if ( d < dist2 )
    {
      dist2 = d;
      n0 = i;
    }
  }

  // find the closest projection on the contour
  double cpt1[3], cpt2[3], t1, t2, p0[3], p1[3], p2[3];
  int n1 = -1, n2 = -1;
  contour_lines->InitTraversal();
  vtkIdType npts = 0, *cellpts;
  while ( contour_lines->GetNextCell( npts, cellpts ) )
  {
    if ( cellpts[0] == n0 )
    {
      if ( n1 < 0 )
      {
        n1 = cellpts[1];
      }
      else
      {
        n2 = cellpts[1];
      }
    }
    else if ( cellpts[1] == n0 )
    {
      if ( n1 < 0 )
      {
        n1 = cellpts[0];
      }
      else
      {
        n2 = cellpts[0];
      }
    }
    if ( n1 >= 0 && n2 >= 0 )
    {
      break;
    }
  }

  if ( n1 >= 0 && n2 >= 0 )
  {
    contour_pts->GetPoint( n0, p0 );
    contour_pts->GetPoint( n1, p1 );
    contour_pts->GetPoint( n2, p2 );
    double d1 = vtkLine::DistanceToLine( old_pt, p0, p1, t1, cpt1 );
    double d2 = vtkLine::DistanceToLine( old_pt, p0, p2, t2, cpt2 );
    if ( d1 < d2 )
    {
      pt_out[0] = cpt1[0];
      pt_out[1] = cpt1[1];
      pt_out[2] = cpt1[2];
    }
    else
    {
      pt_out[0] = cpt2[0];
      pt_out[1] = cpt2[1];
      pt_out[2] = cpt2[2];
    }
    return true;
  }
  else
  {
    return false;
  }
}

void FSSurface::GetVectorAtVertex( int nVertex, double* vec_out, int nVector )
{
  int nv = nVector;
  if ( nv < 0 )
  {
    nv = m_nActiveVector;
  }

  VertexItem* vectors = m_vertexVectors[nv].data;
  vec_out[0] = vectors[nVertex].x;
  vec_out[1] = vectors[nVertex].y;
  vec_out[2] = vectors[nVertex].z;
}

bool FSSurface::SetActiveSurface( int nIndex )
{
  if ( nIndex == m_nActiveSurface )
  {
    return false;
  }

  m_nActiveSurface = nIndex;

  RestoreVertices( m_MRIS, nIndex );

  if ( m_fNormalSets[nIndex] == NULL )
  {
    ComputeNormals();
    SaveNormals( m_MRIS, nIndex );
  }
  else
  {
    RestoreNormals( m_MRIS, nIndex );
  }

  UpdateVerticesAndNormals();

  return true;
}

bool FSSurface::SetActiveVector( int nIndex )
{
  if ( nIndex == m_nActiveVector )
  {
    return false;
  }

  m_nActiveVector = nIndex;

  UpdateVectors();
  return true;
}

void FSSurface::NormalFace(int fac, int n, float *norm)
{
  int n0,n1;
  FACE *f;
  MRIS* mris = m_MRIS;
  float v0[3],v1[3];

  n0 = (n==0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n==VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
  Normalize(v0);
  Normalize(v1);
  norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
  norm[1] = v1[0]*v0[2] - v0[0]*v1[2];
  norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];
}

float FSSurface::TriangleArea(int fac, int n)
{
  int n0,n1;
  FACE *f;
  MRIS* mris = m_MRIS;
  float v0[3],v1[3],d1,d2,d3;

  n0 = (n==0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n==VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
  d1 = -v1[1]*v0[2] + v0[1]*v1[2];
  d2 = v1[0]*v0[2] - v0[0]*v1[2];
  d3 = -v1[0]*v0[1] + v0[0]*v1[1];
  return sqrt(d1*d1+d2*d2+d3*d3)/2;
}

void FSSurface::Normalize( float v[3] )
{
  float d;

  d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d>0)
  {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }
}

void FSSurface::ComputeNormals()
{
  if ( !m_MRIS )
  {
    return;
  }

  MRIS* mris = m_MRIS;
  int k,n;
  VERTEX *v;
  VERTEX_TOPOLOGY* vt;
  FACE *f;
  float norm[3],snorm[3];

  for (k=0; k<mris->nfaces; k++)
  {
    if (mris->faces[k].ripflag)
    {
      f = &mris->faces[k];
      for (n=0; n<VERTICES_PER_FACE; n++)
      {
        mris->vertices[f->v[n]].border = TRUE;
      }
    }
  }
  for (k=0; k<mris->nvertices; k++)
  {
    v = &mris->vertices[k];
    vt = &mris->vertices_topology[k];
    if (!mris->vertices[k].ripflag)
    {
      snorm[0]=snorm[1]=snorm[2]=0;
      v->area = 0;
      for (n=0; n<vt->num; n++)
        if (!mris->faces[vt->f[n]].ripflag)
        {
          NormalFace(vt->f[n],vt->n[n],norm);
          snorm[0] += norm[0];
          snorm[1] += norm[1];
          snorm[2] += norm[2];
          v->area += TriangleArea(vt->f[n],vt->n[n]);
          /* Note: overest. area by 2! */
        }
      Normalize( snorm );

      if (v->origarea<0)
      {
        v->origarea = v->area;
      }

      v->nx = snorm[0];
      v->ny = snorm[1];
      v->nz = snorm[2];
    }
  }
}

void FSSurface::ConvertSurfaceToRAS ( float iX, float iY, float iZ,
                                      float& oX, float& oY, float& oZ ) const
{
  float surface[3];
  float ras[3];

  surface[0] = iX;
  surface[1] = iY;
  surface[2] = iZ;

  this->ConvertSurfaceToRAS( surface, ras );

  oX = ras[0];
  oY = ras[1];
  oZ = ras[2];
}

void FSSurface::ConvertSurfaceToRAS ( double iX, double iY, double iZ,
                                      double& oX, double& oY, double& oZ ) const
{
  double surface[3];
  double ras[3];

  surface[0] = iX;
  surface[1] = iY;
  surface[2] = iZ;

  this->ConvertSurfaceToRAS( surface, ras );

  oX = ras[0];
  oY = ras[1];
  oZ = ras[2];
}

void FSSurface::ConvertRASToSurface ( float iX, float iY, float iZ,
                                      float& oX, float& oY, float& oZ ) const
{
  float ras[3];
  float surface[3];

  ras[0] = iX;
  ras[1] = iY;
  ras[2] = iZ;

  this->ConvertRASToSurface( ras, surface );

  oX = surface[0];
  oY = surface[1];
  oZ = surface[2];
}

void FSSurface::ConvertRASToSurface ( double iX, double iY, double iZ,
                                      double& oX, double& oY, double& oZ ) const
{
  double ras[3];
  double surface[3];

  ras[0] = iX;
  ras[1] = iY;
  ras[2] = iZ;

  this->ConvertRASToSurface( ras, surface );

  oX = surface[0];
  oY = surface[1];
  oZ = surface[2];
}

void FSSurface::ConvertSurfaceToRAS ( float const iSurf[3], float oRAS[3] ) const
{
  m_SurfaceToRASTransform->TransformPoint( iSurf, oRAS );
}

void FSSurface::ConvertSurfaceToRAS ( double const iSurf[3], double oRAS[3] ) const
{
  m_SurfaceToRASTransform->TransformPoint( iSurf, oRAS );
}

void FSSurface::ConvertRASToSurface ( float const iRAS[3], float oSurf[3] ) const
{
  m_SurfaceToRASTransform->GetInverse()->TransformPoint( iRAS, oSurf );
}

void FSSurface::ConvertRASToSurface ( double const iRAS[3], double oSurf[3] ) const
{
  m_SurfaceToRASTransform->GetInverse()->TransformPoint( iRAS, oSurf );
}

void FSSurface::ConvertTargetToRAS(const double iTarget[], double oRAS[]) const
{
  m_targetToRasTransform->TransformPoint(iTarget, oRAS);
}

void FSSurface::ConvertRASToTarget(const double iRAS[], double oTarget[]) const
{
  m_targetToRasTransform->GetInverse()->TransformPoint(iRAS, oTarget);
}

void FSSurface::GetBounds ( float oRASBounds[6] )
{
  if ( NULL == m_MRIS )
  {
    oRASBounds[0] = oRASBounds[1] = oRASBounds[2] =
        oRASBounds[3] = oRASBounds[4] = oRASBounds[5] = 0;
    return;
  }

  if ( m_bBoundsCacheDirty )
  {
    m_RASBounds[0] = m_RASBounds[2] = m_RASBounds[4] = 999999;
    m_RASBounds[1] = m_RASBounds[3] = m_RASBounds[5] = -999999;

    // Find the bounds.
    for ( int vno = 0; vno < m_MRIS->nvertices; vno++ )
    {

      // Translate to actual RAS coords.
      float rasX, rasY, rasZ;
      this->ConvertSurfaceToRAS( m_fVertexSets[SurfaceMain][vno].x,
                                 m_fVertexSets[SurfaceMain][vno].y,
                                 m_fVertexSets[SurfaceMain][vno].z,
                                 rasX, rasY, rasZ );

      if ( rasX < m_RASBounds[0] )
      {
        m_RASBounds[0] = rasX;
      }
      if ( rasX > m_RASBounds[1] )
      {
        m_RASBounds[1] = rasX;
      }
      if ( rasY < m_RASBounds[2] )
      {
        m_RASBounds[2] = rasY;
      }
      if ( rasY > m_RASBounds[3] )
      {
        m_RASBounds[3] = rasY;
      }
      if ( rasZ < m_RASBounds[4] )
      {
        m_RASBounds[4] = rasZ;
      }
      if ( rasZ > m_RASBounds[5] )
      {
        m_RASBounds[5] = rasZ;
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


int FSSurface::GetNumberOfVertices () const
{
  if ( m_MRIS )
  {
    return m_MRIS->nvertices;
  }
  else
  {
    return 0;
  }
}


int FSSurface::FindVertexAtRAS ( float const iRAS[3], float* oDistance, int surface_type )
{
  float surf[3];
  this->ConvertRASToSurface( iRAS, surf );

  return this->FindVertexAtSurfaceRAS( surf, oDistance, surface_type );
}

int FSSurface::FindVertexAtRAS ( double const iRAS[3], double* oDistance, int surface_type )
{
  double surf[3];
  this->ConvertRASToSurface( iRAS, surf );

  return this->FindVertexAtSurfaceRAS( surf, oDistance, surface_type );
}

int FSSurface::FindVertexAtSurfaceRAS ( float const iSurfaceRAS[3], float* oDistance, int surface_type )
{
  VERTEX v;
  v.x = iSurfaceRAS[0];
  v.y = iSurfaceRAS[1];
  v.z = iSurfaceRAS[2];
  float distance;
  int nSurface = (surface_type < 0 ? m_nActiveSurface : surface_type);
  int nClosestVertex =
      MHTfindClosestVertexNoXYZ( m_HashTable[nSurface], m_MRIS, v.x,v.y,v.z, &distance );

  if ( -1 == nClosestVertex )
  {
    // cerr << "No vertices found.\n";
    return -1;
  }

  if ( NULL != oDistance )
  {
    *oDistance = distance;
  }

  return nClosestVertex;
}

int FSSurface::FindVertexAtSurfaceRAS ( double const iSurfaceRAS[3], double* oDistance, int surface_type )
{
  VERTEX v;
  v.x = static_cast<float>(iSurfaceRAS[0]);
  v.y = static_cast<float>(iSurfaceRAS[1]);
  v.z = static_cast<float>(iSurfaceRAS[2]);
  float distance;
  int nSurface = (surface_type < 0 ? m_nActiveSurface : surface_type);
  int nClosestVertex = MHTfindClosestVertexNoXYZ( m_HashTable[nSurface], m_MRIS, v.x,v.y,v.z, &distance );
  if ( -1 == nClosestVertex )
  {
    // cerr << "No vertices found.";
  }

  if ( NULL != oDistance )
  {
    *oDistance = static_cast<double>(distance);
  }

  return nClosestVertex;
}

bool FSSurface::GetRASAtVertex ( int inVertex, float ioRAS[3], int surface_type )
{
  float surfaceRAS[3];
  if ( this->GetSurfaceRASAtVertex( inVertex, surfaceRAS, surface_type ) )
  {
    this->ConvertSurfaceToRAS( surfaceRAS, ioRAS );
    return true;
  }
  else
  {
    return false;
  }
}

bool FSSurface::GetRASAtVertex ( int inVertex, double ioRAS[3], int surface_type )
{
  Q_UNUSED(surface_type);
  double surfaceRAS[3];
  if ( this->GetSurfaceRASAtVertex( inVertex, surfaceRAS ) )
  {
    this->ConvertSurfaceToRAS( surfaceRAS, ioRAS );
    return true;
  }
  else
  {
    return false;
  }
}

bool FSSurface::GetSurfaceRASAtVertex ( int inVertex, float ioRAS[3], int surface_type )
{
  if ( m_MRIS == NULL )
    //  throw runtime_error( "GetRASAtVertex: m_MRIS was NULL" );
  {
    return false;
  }

  if ( inVertex < 0 || inVertex >= m_MRIS->nvertices )
    //  throw runtime_error( "GetRASAtVertex: inVertex was invalid" );
  {
    return false;
  }

  int nSurface = (surface_type < 0 ? m_nActiveSurface : surface_type);
  if ( nSurface >= 0 && m_fVertexSets[nSurface] != NULL )
  {
    ioRAS[0] = m_fVertexSets[nSurface][inVertex].x;
    ioRAS[1] = m_fVertexSets[nSurface][inVertex].y;
    ioRAS[2] = m_fVertexSets[nSurface][inVertex].z;
  }
  else
  {
    ioRAS[0] = m_MRIS->vertices[inVertex].x;
    ioRAS[1] = m_MRIS->vertices[inVertex].y;
    ioRAS[2] = m_MRIS->vertices[inVertex].z;
  }

  return true;
}

bool FSSurface::GetSurfaceRASAtVertex ( int inVertex, double ioRAS[3], int surface_type )
{
  if ( m_MRIS == NULL )
    //  throw runtime_error( "GetRASAtVertex: m_MRIS was NULL" );
  {
    return false;
  }

  if ( inVertex < 0 || inVertex >= m_MRIS->nvertices )
    //  throw runtime_error( "GetRASAtVertex: inVertex was invalid" );
  {
    return false;
  }

  int nSurface = (surface_type < 0 ? m_nActiveSurface : surface_type);
  if ( nSurface >= 0 && m_fVertexSets[nSurface] != NULL )
  {
    ioRAS[0] = m_fVertexSets[nSurface][inVertex].x;
    ioRAS[1] = m_fVertexSets[nSurface][inVertex].y;
    ioRAS[2] = m_fVertexSets[nSurface][inVertex].z;
  }
  else
  {
    ioRAS[0] = static_cast<double>(m_MRIS->vertices[inVertex].x);
    ioRAS[1] = static_cast<double>(m_MRIS->vertices[inVertex].y);
    ioRAS[2] = static_cast<double>(m_MRIS->vertices[inVertex].z);
  }

  return true;
}

QString FSSurface::GetVectorSetName( int nSet )
{
  if ( nSet >= 0 )
  {
    return m_vertexVectors[nSet].name;
  }
  else
  {
    return "";
  }
}

double FSSurface::GetCurvatureValue( int nVertex )
{
  return m_MRIS->vertices[nVertex].curv;
}

void FSSurface::Reposition( FSVolume *volume, int target_vno, double target_val, int nsize, double sigma, int flags )
{
  MRISsaveVertexPositions( m_MRIS, INFLATED_VERTICES );
  float fval = (float)target_val;
  MRISrepositionSurface( m_MRIS, volume->GetMRI(), &target_vno, &fval, 1, nsize, sigma, flags );
  PostEditProcess();
}

void FSSurface::Reposition( FSVolume *volume, int target_vno, double* coord, int nsize, double sigma, int flags )
{
  MRISsaveVertexPositions( m_MRIS, INFLATED_VERTICES );
  MRISrepositionSurfaceToCoordinate( m_MRIS, volume->GetMRI(), target_vno, coord[0], coord[1], coord[2], nsize, sigma, flags );
  PostEditProcess();
}

void FSSurface::RepositionSmooth(int vert_n, int nbhd_size, int nsmoothing_steps)
{
  MRISsaveVertexPositions( m_MRIS, INFLATED_VERTICES );
  for (int vno = 0 ; vno < m_MRIS->nvertices ; vno++)
  {
    m_MRIS->vertices[vno].ripflag = 1;
  }
  for (int n = 0 ; n < m_MRIS->vertices_topology[vert_n].vnum ; n++)
  {
    m_MRIS->vertices[m_MRIS->vertices_topology[vert_n].v[n]].ripflag = 0;
  }
  m_MRIS->vertices[vert_n].ripflag = 0;
  MRISerodeRipped(m_MRIS, nbhd_size);
  MRISaverageVertexPositions(m_MRIS, nsmoothing_steps);
  MRISunrip(m_MRIS);
  PostEditProcess();
}

bool FSSurface::Smooth(int nMethod, int niters, double lambda, double K_bp)
{
  MRISsaveVertexPositions( m_MRIS, INFLATED_VERTICES );
  if (nMethod == 0 ) // Taubin
  {
    double mu = (1.0)/((K_bp)-1.0/lambda);
    if (MRIStaubinSmooth(m_MRIS, niters, lambda, mu, TAUBIN_UNIFORM_WEIGHTS) != 0)
      return false;
  }
  else if (nMethod == 1) // standard
  {
    if (MRISaverageVertexPositions(m_MRIS, niters) != 0)
      return false;
  }
  PostEditProcess();
  return true;
}

void FSSurface::RemoveIntersections()
{
  MRISremoveIntersections(m_MRIS);
  PostEditProcess();
}

void FSSurface::UpdateCoords()
{
  SaveVertices( m_MRIS, m_nActiveSurface );
  ComputeNormals();
  SaveNormals( m_MRIS, m_nActiveSurface );
  UpdateVerticesAndNormals();
}

void FSSurface::PostEditProcess()
{
  UpdateHashTable(m_nActiveSurface);
  UpdateCoords();
}

void FSSurface::RepositionVertex(int vno, double *coord)
{
  MRISsaveVertexPositions( m_MRIS, INFLATED_VERTICES );
  VERTEX *v = &m_MRIS->vertices[vno];
  v->x = coord[0];
  v->y = coord[1];
  v->z = coord[2];
  PostEditProcess();
}

void FSSurface::UndoReposition()
{
  MRISrestoreVertexPositions( m_MRIS, INFLATED_VERTICES );
  SaveVertices( m_MRIS, m_nActiveSurface );
  ComputeNormals();
  SaveNormals( m_MRIS, m_nActiveSurface );
  UpdateVerticesAndNormals();
}

bool FSSurface::FindPath(int* vert_vno, int num_vno,
                         int* path, int* path_length )
{
  MRIS* mris = m_MRIS;
  bool flag2d = false;
  int max_path_length = 100;
  int cur_vert_vno;
  int src_vno;
  int dest_vno;
  int vno;
  char* check;
  float* dist;
  int* pred;
  char done;
  VERTEX* v;
  VERTEX_TOPOLOGY* vt;
  VERTEX* u;
  float closest_dist;
  int closest_vno;
  int neighbor;
  int neighbor_vno;
  float dist_uv;
  int path_vno;
  //  int num_path = 0;
  int num_checked;
  //  float vu_x, vu_y, vu_z;

  dist = (float*) calloc (mris->nvertices, sizeof(float));
  pred = (int*) calloc (mris->nvertices, sizeof(int));
  check = (char*) calloc (mris->nvertices, sizeof(char));
  num_checked = 0;
  (*path_length) = 0;

  for (cur_vert_vno = 0; cur_vert_vno < num_vno-1; cur_vert_vno++)
  {

    /* clear everything */
    for (vno = 0; vno < mris->nvertices; vno++)
    {
      dist[vno] = 999999;
      pred[vno] = -1;
      check[vno] = FALSE;
    }

    /* Set src and dest */
    src_vno = vert_vno[cur_vert_vno+1];
    dest_vno = vert_vno[cur_vert_vno];

    /* make sure both are in range. */
    if (src_vno < 0 || src_vno >= mris->nvertices ||
        dest_vno < 0 || dest_vno >= mris->nvertices)
      continue;

    if (src_vno == dest_vno)
      continue;

    /* pull the src vertex in. */
    dist[src_vno] = 0;
    pred[src_vno] = vno;
    check[src_vno] = TRUE;

    done = FALSE;
    while (!done)
    {
      /* find the vertex with the shortest edge. */
      closest_dist = 999999;
      closest_vno = -1;
      for (vno = 0; vno < mris->nvertices; vno++)
        if (check[vno])
          if (dist[vno] < closest_dist)
          {
            closest_dist = dist[vno];
            closest_vno = vno;
          }
      v = &(mris->vertices[closest_vno]);
      vt = &(mris->vertices_topology[closest_vno]);
      check[closest_vno] = FALSE;

      /* if this is the dest node, we're done. */
      if (closest_vno == dest_vno)
      {
        done = TRUE;
      }
      else
      {
        /* relax its neighbors. */
        for (neighbor = 0; neighbor < vt->vnum; neighbor++)
        {
          neighbor_vno = vt->v[neighbor];
          u = &(mris->vertices[neighbor_vno]);

          /* calc the vector from u to v. */
          //          vu_x = u->x - v->x;
          //          vu_y = u->y - v->y;
          //          if (flag2d)
          //            vu_z = 0;
          //          else
          //            vu_z = u->z - v->z;

          /* recalc the weight. */
          if (flag2d)
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
                           ((v->y - u->y) * (v->y - u->y)));
          else
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
                           ((v->y - u->y) * (v->y - u->y)) +
                           ((v->z - u->z) * (v->z - u->z)));

          /* if this is a new shortest path, update the predecessor,
                             weight, and add it to the list of ones to check next. */
          if (dist_uv + dist[closest_vno] < dist[neighbor_vno])
          {
            pred[neighbor_vno] = closest_vno;
            dist[neighbor_vno] = dist_uv + dist[closest_vno];
            check[neighbor_vno] = TRUE;
          }
        }
      }
      num_checked++;
      if ((num_checked % 100) == 0)
      {

      }
    }

    /* add the predecessors from the dest to the src to the path. */
    path_vno = dest_vno;
    path[(*path_length)++] = dest_vno;
    while (pred[path_vno] != src_vno &&
           (*path_length) < max_path_length )
    {
      path[(*path_length)++] = pred[path_vno];
      path_vno = pred[path_vno];
    }
  }

  free (dist);
  free (pred);
  free (check);

  return true;
}

void FSSurface::RipFaces()
{
  MRIS* mris = m_MRIS;
  int n,k;
  FACE *f;

  for (k=0;k<mris->nfaces;k++)
    mris->faces[k].ripflag = FALSE;
  for (k=0;k<mris->nfaces;k++)
  {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
      if (mris->vertices[f->v[n]].ripflag)
        f->ripflag = TRUE;
  }
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].border = FALSE;
  for (k=0;k<mris->nfaces;k++)
    if (mris->faces[k].ripflag)
    {
      f = &mris->faces[k];
      for (n=0;n<VERTICES_PER_FACE;n++)
        mris->vertices[f->v[n]].border = TRUE;
    }
}

QVector<int> FSSurface::FloodFillFromSeed(int seed_vno)
{
  MRIS* mris = m_MRIS;
  char* filled;
  int num_filled_this_iter;
  int num_filled;
  int iter;
  int min_vno, max_vno, step_vno;
  int vno;
  //  int this_label = 0;
  int neighbor_index;
  int neighbor_vno;
  VERTEX* v;
  VERTEX_TOPOLOGY* vt;
  VERTEX* neighbor_v;
  //  float fvalue = 0;
  //  float seed_curv = 0;
  //  float seed_fvalue = 0;
  //  int new_index;
  //  int num_labels_found, found_label_index;
  //  int skip;
  int count;

  QVector<int> filled_verts;

  if (seed_vno < 0 || seed_vno >= mris->nvertices)
    return filled_verts;

  /* init filled array. */
  filled = (char*) calloc (mris->nvertices, sizeof(char));
  memset(filled, 0, sizeof(char)*mris->nvertices);

  /* start with the seed filled.*/
  filled[seed_vno] = TRUE;
  filled_verts << seed_vno;

  /* find seed values for some conditions. */
  //  if (params->dont_cross_label)
  //    this_label = labl_selected_label;
  //  if (params->dont_cross_cmid)
  //    seed_curv = mris->vertices[seed_vno].curv;
  //  if (params->dont_cross_fthresh)
  //    sclv_get_value (&mris->vertices[seed_vno],
  //                    sclv_current_field, &seed_fvalue);

  /* while we're still filling stuff in a pass... */
  num_filled_this_iter = 1;
  num_filled = 0;
  iter = 0;
  while (num_filled_this_iter > 0)
  {
    num_filled_this_iter = 0;

    /* switch between iterating forward and backwards. */
    if ((iter%2)==0)
    {
      min_vno = 0;
      max_vno = mris->nvertices-1;
      step_vno = 1;
    }
    else
    {
      min_vno = mris->nvertices-1;
      max_vno = 0;
      step_vno = -1;
    }

    /* for each vertex, if it's filled, check its neighbors. for the
       rules that are up-to-and-including, make the check on this
       vertex. for the rules that are up-to-and-not-including, check
       on the neighbor. */
    for (vno = min_vno; vno != max_vno; vno += step_vno)
    {
      if (filled[vno])
      {

        /* check the neighbors... */
        v = &mris->vertices[vno];
        vt = &mris->vertices_topology[vno];

        /* if this vert is ripped, move on. */
        if (v->ripflag)
        {
          continue;
        }

        /* if we're not crossing paths, check if this is a
           path. if so, move on. */
        //        if (params->dont_cross_path &&
        //            path_is_vertex_on_path (vno))
        //        {
        //          continue;
        //        }

        /* if we're not crossing the cmid, see if the cmid at this
           vertex is on the other side of the cmid as the seed
           point. if so, move on. */
        //        if (params->dont_cross_cmid &&
        //            ((seed_curv <= cmid && v->curv > cmid) ||
        //             (seed_curv >= cmid && v->curv < cmid)))
        //        {
        //          continue;
        //        }

        for (neighbor_index = 0;
             neighbor_index < vt->vnum;
             neighbor_index++)
        {
          neighbor_vno = vt->v[neighbor_index];
          neighbor_v = &mris->vertices[neighbor_vno] ;

          /* if the neighbor is filled, move on. */
          if (filled[neighbor_vno])
            continue;

          if (neighbor_v->ripflag)
            continue;

          /* if we're not crossing labels, check if the label at
             this vertex is the same as the one at the seed. if not,
             move on. */
          //          if (params->dont_cross_label ||
          //              params->dont_fill_unlabeled)
          //          {

          //            labl_find_label_by_vno (neighbor_vno, 0,
          //                                    label_index_array,
          //                                    LABL_MAX_LABELS,
          //                                    &num_labels_found);
          //            if (num_labels_found > 0 &&
          //                params->dont_cross_label)
          //            {
          //              skip = 0;
          //              for (found_label_index = 0;
          //                   found_label_index < num_labels_found;
          //                   found_label_index++)
          //              {
          //                if (label_index_array[found_label_index] !=
          //                    this_label)
          //                {
          //                  skip = 1;
          //                  break;
          //                }
          //              }
          //              if (skip) continue;
          //            }
          //            if (num_labels_found == 0 &&
          //                params->dont_fill_unlabeled)
          //            {
          //              continue;
          //            }
          //          }

          /* if we're not crossing the fthresh, make sure this
             point is above it, or, if our initial functional
             value was negative, make sure it's not above
             -fthresh. if not, move on. */
          //          if (params->dont_cross_fthresh)
          //          {
          //            sclv_get_value (neighbor_v, sclv_current_field, &fvalue);
          //            if ((fthresh != 0 &&
          //                 seed_fvalue > 0 &&
          //                 fvalue < fthresh) ||
          //                (fthresh != 0 &&
          //                 seed_fvalue < 0 &&
          //                 fvalue > -fthresh) ||
          //                (fthresh == 0 && (fvalue * seed_fvalue < 0)))
          //            {
          //              continue;
          //            }
          //          }

          /* mark this vertex as filled. */
          filled[neighbor_vno] = TRUE;
          filled_verts << neighbor_vno;
          num_filled_this_iter++;
          num_filled++;
        }
      }
    }

    iter++;
  }

  /* mark all filled vertices. */
  for (vno = 0; vno < mris->nvertices; vno++ )
  {
    mris->vertices[vno].ripflag = (!filled[vno]);
  }

  free (filled);

  return filled_verts;
}

QVector<int> FSSurface::MakeCutLine(const QVector<int>& verts)
{
  QVector<int> old;
  MRIS* mris = m_MRIS;
  for (int i = 0; i < verts.size(); i++)
  {
    VERTEX *v;
    int vno = verts[i];
    v = &mris->vertices[vno];
    if (v->ripflag != 1)
      old << vno;
    v->ripflag = 1;
  }
  return old;
}

void FSSurface::ClearCuts(const QVector<int> &verts)
{
  MRIS* mris = m_MRIS;
  if (verts.isEmpty())
  {
    for (int i = 0; i < mris->nvertices; i++)
    {
      VERTEX *v;
      v = &mris->vertices[i];
      v->ripflag = m_originalRipflags[i];
    }
  }
  else
  {
    for (int i = 0; i < verts.size(); i++)
    {
      VERTEX *v;
      int vno = verts[i];
      v = &mris->vertices[vno];
      v->ripflag = 0;
    }
  }
}

bool FSSurface::SaveTransform(vtkTransform *t_in, const QString &filename)
{
  if (!m_bValidVolumeGeometry)
  {
    qDebug() << "Surface does not contain valid volume information. Cannot save transformation.";
    return false;
  }

  vtkSmartPointer<vtkTransform> surf2targ = vtkSmartPointer<vtkTransform>::New();
  surf2targ->DeepCopy(m_SurfaceToRASTransform);
  surf2targ->PostMultiply();
  vtkSmartPointer<vtkMatrix4x4> mat = m_targetToRasTransform->GetMatrix();
  mat->Invert();
  surf2targ->Concatenate(mat);

  vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
  t->DeepCopy(surf2targ);
  t->PostMultiply();
  t->Concatenate(t_in->GetMatrix());
  mat = surf2targ->GetMatrix();
  mat->Invert();
  t->Concatenate(mat);

  mat = t->GetMatrix();
  MATRIX* m = MatrixAlloc( 4, 4, MATRIX_REAL );
  for ( int i = 0; i < 16; i++ )
  {
    *MATRIX_RELT(m, (i/4)+1, (i%4)+1) = mat->Element[i/4][i%4];
  }

  MRI* mri = MRIallocHeader(m_MRIS->vg.width, m_MRIS->vg.height, m_MRIS->vg.depth, MRI_UCHAR, 1);
  useVolGeomToMRI(&m_MRIS->vg, mri);

  MATRIX* voxel_xform = MRIrasXformToVoxelXform( mri, NULL, m, NULL );
  LTA* lta = LTAalloc(1, NULL);
  LINEAR_TRANSFORM *lt = &lta->xforms[0];
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0;
  lta->type = LINEAR_VOX_TO_VOX;
  lt->m_L = voxel_xform;

  VOL_GEOM srcG, dstG;
  getVolGeom( mri, &srcG );
  getVolGeom( mri, &dstG );
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
