/**
 * @brief Layer class for MRI volume.
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

#include "LayerMRI.h"
#include "MainWindow.h"
#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkImageReslice.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkCubeSource.h"
#include "vtkVolume.h"
#include "vtkImageThreshold.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkContourFilter.h"
#include "vtkTubeFilter.h"
#include "vtkImageCast.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkAppendPolyData.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkMarchingSquares.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkSimpleLabelEdgeFilter.h"
#include "vtkImageMask.h"
#include "vtkImageResample.h"
#include "vtkPolyDataWriter.h"
#include "vtkMath.h"
#include "vtkImageThreshold.h"
#include "vtkImageShiftScale.h"
#include "vtkImageMapper3D.h"
#include "MyUtils.h"
#include "MyVTKUtils.h"
#include "FSVolume.h"
#include "LayerPropertyMRI.h"
//#include "BuildContourThread.h"
#include "Contour2D.h"
#include "SurfaceRegion.h"
#include "SurfaceRegionGroups.h"
#include "ThreadBuildContour.h"
#include <QtGlobal>
#include <QFile>
#include <QDebug>
#include "ProgressCallback.h"
#include "LayerMRIWorkerThread.h"
#include "vtkImageFlip.h"
#include "LayerSurface.h"
#include "vtkImageResample.h"
#include "vtkImageExtractComponents.h"
#include "vtkMaskPoints.h"
#include <QVariantMap>
#include "LayerROI.h"
#include <QFileInfo>
#include "GeoSWorker.h"
#include "BrushProperty.h"
#include "vtkImageResliceMapper.h"
#include "vtkSTLWriter.h"
#include "vtkImageMathematics.h"
#include "vtkImageExtractComponents.h"
#include "vtkImageAppendComponents.h"
#include "Region3D.h"
#include "LayerPointSet.h"
#include <QTimer>

#include "utils.h"
#include "geos.h"


#define IMAGE_RESAMPLE_FACTOR     4.0     // must be multiples of 2

LayerMRI::LayerMRI( LayerMRI* ref, QObject* parent ) : LayerVolumeBase( parent ),
  m_volumeSource( NULL),
  m_volumeRef( ref ? ref->GetSourceVolume() : NULL ),
  m_bResampleToRAS( false ),
  m_bReorient( false ),
  m_nSampleMethod( SAMPLE_NEAREST ),
  m_bConform( false ),
  m_bWriteResampled( true ),
  m_currentSurfaceRegion( NULL ),
  m_current3DRegion(NULL),
  m_nGotoLabelSlice(-1),
  m_nGotoLabelOrientation(-1),
  m_layerMask(NULL),
  m_correlationSurface(NULL),
  m_bIgnoreHeader(false)
{
  m_strTypeNames.push_back( "Supplement" );
  m_strTypeNames.push_back( "MRI" );
  m_sPrimaryType = "MRI";
  
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  for ( int i = 0; i < 3; i++ )
  {
    // m_nSliceNumber[i] = 0;
    m_sliceActor2D[i] = vtkImageActor::New();
    m_sliceActor3D[i] = vtkImageActor::New();
    /*
    m_sliceActor2D[i]->GetProperty()->SetAmbient( 1 );
    m_sliceActor2D[i]->GetProperty()->SetDiffuse( 0 );
    m_sliceActor3D[i]->GetProperty()->SetAmbient( 1 );
    m_sliceActor3D[i]->GetProperty()->SetDiffuse( 0 );
    */
    m_sliceActor2D[i]->InterpolateOff();
    m_sliceActor3D[i]->InterpolateOff();

    m_glyphActor2D[i] = vtkActor::New();
    m_glyphActor3D[i] = vtkActor::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    m_glyphActor2D[i]->SetMapper( mapper );
    m_glyphActor3D[i]->SetMapper( mapper2 );
    m_vectorDotActor2D[i] = vtkActor::New();
    m_vectorDotActor2D[i]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    m_vectorDotActor2D[i]->SetProperty( m_vectorDotActor2D[i]->MakeProperty() );
    m_vectorDotActor2D[i]->GetProperty()->SetPointSize(3*ratio);
    m_vectorDotActor2D[i]->GetProperty()->SetInterpolationToFlat();
    m_projectionMapActor[i] = vtkImageActor::New();
    m_projectionMapActor[i]->VisibilityOff();
#if VTK_MAJOR_VERSION > 5
    m_sliceActor2D[i]->ForceOpaqueOn();
    m_sliceActor3D[i]->ForceOpaqueOn();
    m_projectionMapActor[i]->ForceOpaqueOn();
    m_glyphActor2D[i]->ForceOpaqueOn();
    m_glyphActor3D[i]->ForceOpaqueOn();
    m_vectorDotActor2D[i]->ForceOpaqueOn();
#endif
  }
  
  m_actorContour = vtkSmartPointer<vtkActor>::New();
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION > 5
  mapper->SetInputData( vtkSmartPointer<vtkPolyData>::New() );
  m_actorContour->ForceOpaqueOn();
#else
  mapper->SetInput( vtkSmartPointer<vtkPolyData>::New() );
#endif
  m_actorContour->SetMapper( mapper );

  m_actorCurrentContour = NULL;
  
  m_propVolume = vtkSmartPointer<vtkVolume>::New();
  
  m_nThreadID = 0;
  m_surfaceRegionGroups = new SurfaceRegionGroups( this );
  
  private_buf1_3x3 = new double*[3];
  private_buf2_3x3 = new double*[3];
  for ( int i = 0; i < 3; i++ )
  {
    private_buf1_3x3[i] = new double[3];
    private_buf2_3x3[i] = new double[3];
  }
  
  mProperty = new LayerPropertyMRI( this );
  ConnectProperty();
  
  qRegisterMetaType< IntList >( "IntList" );
  m_worker = new LayerMRIWorkerThread(this);
  connect(m_worker, SIGNAL(LabelInformationReady()), this, SLOT(OnLabelInformationReady()));

  connect(this, SIGNAL(Modified()), this, SLOT(UpdateLabelInformation()));
  
  QVariantMap map = MainWindow::GetMainWindow()->GetDefaultSettings();
  if (map["Smoothed"].toBool())
  {
    GetProperty()->SetTextureSmoothing(1);
    UpdateTextureSmoothing();
  }

  m_geos = NULL;
}

LayerMRI::~LayerMRI()
{
  if (m_worker->isRunning())
    m_worker->Abort();
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->Delete();
    m_sliceActor3D[i]->Delete();
    m_glyphActor2D[i]->Delete();
    m_glyphActor3D[i]->Delete();
    m_vectorDotActor2D[i]->Delete();
    m_projectionMapActor[i]->Delete();
  }

  if ( m_sFilename.size() > 0 )
  {
    GetProperty()->SaveSettings( m_sFilename );
  }
  
  if ( m_volumeSource )
  {
    delete m_volumeSource;
  }
  
  for ( int i = 0; i < m_surfaceRegions.size(); i++ )
  {
    delete m_surfaceRegions[i];
  }
  
  for ( int i = 0; i < 3; i++ )
  {
    delete[] private_buf1_3x3[i];
    delete[] private_buf2_3x3[i];
  }
  delete[] private_buf1_3x3;
  delete[] private_buf2_3x3;
  
  QList<int> keys = m_labelActors.keys();
  foreach (int i, keys)
    m_labelActors[i]->Delete();
  m_labelActors.clear();

  if (m_geos)
    delete m_geos;
}

void LayerMRI::ConnectProperty()
{
  LayerPropertyMRI* p = GetProperty();
  connect( p, SIGNAL(ColorMapChanged()), this, SLOT(UpdateColorMap()) );
  connect( p, SIGNAL(ContourChanged()), this, SLOT(UpdateContour()) );
  connect( p, SIGNAL(ContourColorChanged()), this, SLOT(UpdateContourColor()) );
  connect( p, SIGNAL(ContourShown(bool)), this, SLOT(ShowContour()) );
  connect( p, SIGNAL(ContourSmoothIterationChanged(int)), this, SLOT(RebuildContour()));
  connect( p, SIGNAL(ContourNeedsRebuild()), this, SLOT(RebuildContour()));
  connect( p, SIGNAL(DisplayModeChanged()), this, SLOT(UpdateDisplayMode()) );
  connect( p, SIGNAL(LabelOutlineChanged(bool)), this, SLOT(UpdateLabelOutline()) );
  connect( p, SIGNAL(OpacityChanged(double)), this, SLOT(UpdateOpacity()) );
  //  connect( p, SIGNAL(ResliceInterpolationChanged()), this, SLOT(UpdateResliceInterpolation()) );
  connect( p, SIGNAL(TextureSmoothingChanged()), this, SLOT(UpdateTextureSmoothing()) );
  connect( p, SIGNAL(UpSampleMethodChanged(int)), this, SLOT(UpdateUpSampleMethod()) );
  connect( this, SIGNAL(SurfaceRegionAdded()), this, SIGNAL(ActorChanged()));
  connect( this, SIGNAL(SurfaceRegionRemoved()), this, SIGNAL(ActorChanged()));
  connect( this, SIGNAL(Region3DAdded()), this, SIGNAL(ActorChanged()));
  connect( this, SIGNAL(Region3DRemoved()), this, SIGNAL(ActorChanged()));
  connect( p, SIGNAL(ProjectionMapChanged()), this, SLOT(UpdateProjectionMap()));
  connect( this, SIGNAL(ActiveFrameChanged(int)), this, SLOT(UpdateContour()));
  connect( p, SIGNAL(LabelContourChanged(int)), this, SLOT(OnLabelContourChanged(int)));
  connect( p, SIGNAL(VectorLineWidthChanged(double)), this, SLOT(UpdateVectorLineWidth(double)));
  connect( p, SIGNAL(VectorSkipChanged(int)), SLOT(UpdateVectorActor()));
}

void LayerMRI::SetResampleToRAS( bool bResample )
{
  m_bResampleToRAS = bResample;
}

QString LayerMRI::GetOrientationString()
{
  if ( m_volumeSource )
  {
    return m_volumeSource->GetOrientationString();
  }
  else
  {
    return "RAS";
  }
}

void LayerMRI::SetConform( bool bConform )
{
  m_bConform = bConform;
}

void LayerMRI::ResetRef()
{
  m_volumeSource->ResetRef();
}

void swap_nifti_header(struct nifti_1_header *hdr)
{
  int i;

  hdr->sizeof_hdr = swapInt(hdr->sizeof_hdr);

  for (i = 0; i < 8; i++) hdr->dim[i] = swapShort(hdr->dim[i]);

  hdr->intent_p1 = swapFloat(hdr->intent_p1);
  hdr->intent_p2 = swapFloat(hdr->intent_p2);
  hdr->intent_p3 = swapFloat(hdr->intent_p3);
  hdr->intent_code = swapShort(hdr->intent_code);
  hdr->datatype = swapShort(hdr->datatype);
  hdr->bitpix = swapShort(hdr->bitpix);
  hdr->slice_start = swapShort(hdr->slice_start);

  for (i = 0; i < 8; i++) hdr->pixdim[i] = swapFloat(hdr->pixdim[i]);

  hdr->vox_offset = swapFloat(hdr->vox_offset);
  hdr->scl_slope = swapFloat(hdr->scl_slope);
  hdr->scl_inter = swapFloat(hdr->scl_inter);
  hdr->slice_end = swapShort(hdr->slice_end);
  hdr->cal_max = swapFloat(hdr->cal_max);
  hdr->cal_min = swapFloat(hdr->cal_min);
  hdr->slice_duration = swapFloat(hdr->slice_duration);
  hdr->toffset = swapFloat(hdr->toffset);
  hdr->qform_code = swapShort(hdr->qform_code);
  hdr->sform_code = swapShort(hdr->sform_code);
  hdr->quatern_b = swapFloat(hdr->quatern_b);
  hdr->quatern_c = swapFloat(hdr->quatern_c);
  hdr->quatern_d = swapFloat(hdr->quatern_d);
  hdr->qoffset_x = swapFloat(hdr->qoffset_x);
  hdr->qoffset_y = swapFloat(hdr->qoffset_y);
  hdr->qoffset_z = swapFloat(hdr->qoffset_z);

  for (i = 0; i < 4; i++) hdr->srow_x[i] = swapFloat(hdr->srow_x[i]);

  for (i = 0; i < 4; i++) hdr->srow_y[i] = swapFloat(hdr->srow_y[i]);

  for (i = 0; i < 4; i++) hdr->srow_z[i] = swapFloat(hdr->srow_z[i]);

  return;
}

bool LayerMRI::LoadVolumeFromFile()
{
  if ( m_volumeSource )
  {
    delete m_volumeSource;
  }
  
  m_volumeSource = new FSVolume( m_volumeRef );
  if (m_volumeRef)
    connect(m_volumeRef, SIGNAL(destroyed(QObject*)), this, SLOT(ResetRef()));
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  m_volumeSource->SetConform( m_bConform );
  m_volumeSource->SetInterpolationMethod( m_nSampleMethod );
  m_volumeSource->SetIgnoreHeader(m_bIgnoreHeader);
  
  if ( !m_volumeSource->MRIRead( m_sFilename.toLatin1().data(),
                                 m_sRegFilename.size() > 0 ? m_sRegFilename.toLatin1().data() : NULL ) )
  {
    return false;
  }
  
  ParseSubjectName(m_sFilename);
  InitializeVolume();
  InitializeActors();
  
  GetProperty()->SetVolumeSource( m_volumeSource );
  GetProperty()->RestoreSettings( m_sFilename );

  //  int* dim = m_imageData->GetDimensions();
  //  qDebug() << dim[0] << dim[1] << dim[2];
  
  if (m_nGotoLabelOrientation >= 0)
    m_nGotoLabelSlice = this->GoToLabel(m_nGotoLabelOrientation, m_strGotoLabelName);

  UpdateNiftiHeader();

  if (GetDataType() == MRI_RGB)
    GetProperty()->SetDisplayRGB(true);

  return true;
}

void LayerMRI::UpdateNiftiHeader()
{
  QFileInfo fi(m_sFilename);
  if (fi.suffix() == "nii" || fi.completeSuffix().contains("nii.gz"))
  {
    znzFile fp = znzopen(qPrintable(m_sFilename), "r", fi.suffix() == "gz");
    if (fp)
    {
      znzread(&m_niftiHeader, sizeof(nifti_1_header), 1, fp);
      if (m_niftiHeader.sizeof_hdr != 348)
        swap_nifti_header(&m_niftiHeader);
    }
  }
}

bool LayerMRI::LoadVolumeTransform()
{
  if (!m_volumeSource->LoadRegistrationMatrix(m_sRegFilename))
  {
    cerr << "Could not load transformation from " << qPrintable(m_sRegFilename) << endl;
    return false;
  }
  m_volumeSource->MapMRIToImage();
  InitializeVolume();
  InitializeActors();
  return true;
}

void LayerMRI::UnloadVolumeTransform()
{
  m_sRegFilename.clear();
  m_volumeSource->ClearRegistrationMatrix();
  m_volumeSource->MapMRIToImage();
  InitializeVolume();
  InitializeActors();
}

bool LayerMRI::CreateFromMRIData(void *mri_ptr)
{
  MRI* mri = (MRI*)mri_ptr;
  if ( m_volumeSource )
  {
    delete m_volumeSource;
  }
  m_volumeSource = new FSVolume( m_volumeRef );
  if (m_volumeRef)
    connect(m_volumeRef, SIGNAL(destroyed(QObject*)), this, SLOT(ResetRef()));
  m_volumeSource->SetResampleToRAS( m_bResampleToRAS );
  m_volumeSource->SetConform( m_bConform );
  m_volumeSource->SetInterpolationMethod( m_nSampleMethod );
  m_volumeSource->SetIgnoreHeader(m_bIgnoreHeader);
  if (!m_volumeSource->CreateFromMRIData(mri))
    return false;
  
  ParseSubjectName(m_sFilename);
  InitializeVolume();
  InitializeActors();
  
  GetProperty()->SetVolumeSource( m_volumeSource );
  GetProperty()->RestoreSettings( m_sFilename );
  
  if (m_nGotoLabelOrientation >= 0)
    m_nGotoLabelSlice = this->GoToLabel(m_nGotoLabelOrientation, m_strGotoLabelName);
  
  return true;
}

bool LayerMRI::Create( LayerMRI* mri, bool bCopyVoxelData, int data_type, int voxel_option )
{
  if ( m_volumeSource )
  {
    delete m_volumeSource;
  }
  
  m_volumeSource = new FSVolume( mri->m_volumeSource );
  if ( !m_volumeSource->Create( mri->m_volumeSource, bCopyVoxelData, data_type ) )
  {
    return false;
  }
  
  m_bResampleToRAS = mri->m_bResampleToRAS;
  m_bReorient = mri->m_bReorient;
  m_imageDataRef = mri->GetImageData();
  if ( m_imageDataRef != NULL )
  {
    SetWorldOrigin( mri->GetWorldOrigin() );
    SetWorldVoxelSize( mri->GetWorldVoxelSize() );
    SetWorldSize( mri->GetWorldSize() );
    
    m_imageData = m_volumeSource->GetImageOutput();
    
    int* dim = m_imageData->GetDimensions();
    char* ptr = (char*)m_imageData->GetScalarPointer();
    int scalar_type = m_imageData->GetScalarType();
    int len = qMin(dim[0], qMin(dim[1], dim[2]))/3;
    int len2 = len*len;
    if (voxel_option == 0)    // sphere
    {
      int c[3] = {dim[0]/2, dim[1]/2, dim[2]/2};
      double val = 100;
      for (int i = c[0]-len; i <= c[0]+len; i++)
      {
        for (int j = c[1]-len; j <= c[1]+len; j++)
        {
          for (int k = c[2]-len; k <= c[2]+len; k++)
          {
            if ((i-c[0])*(i-c[0])+(j-c[1])*(j-c[1])+(k-c[2])*(k-c[2]) < len2)
              MyVTKUtils::SetImageDataComponent(ptr, dim, 1, i, j, k, 0, scalar_type, val);
          }
        }
      }
    }
    else if (voxel_option == 1)  // cube
    {
      int c[3] = {dim[0]/2, dim[1]/2, dim[2]/2};
      double val = 100;
      for (int i = c[0]-len; i <= c[0]+len; i++)
      {
        for (int j = c[1]-len; j <= c[1]+len; j++)
        {
          for (int k = c[2]-len; k <= c[2]+len; k++)
          {
            MyVTKUtils::SetImageDataComponent(ptr, dim, 1, i, j, k, 0, scalar_type, val);
          }
        }
      }
    }
    else if (voxel_option == 2) // mask
    {
      double w = mri->GetProperty()->GetWindow(),
          l = mri->GetProperty()->GetLevel();
      /*
      double dMin = l - w/2, dMax = l + w/2;
      for (int i = 0; i < dim[0]; i++)
      {
        for (int j = 0; j < dim[1]; j++)
        {
          for (int k = 0; k < dim[2]; k++)
          {
            double val = m_imageDataRef->GetScalarComponentAsDouble(i, j, k, mri->GetActiveFrame());
            if (val >= dMin && val <= dMax)
              m_imageData->SetScalarComponentFromDouble(i, j, k, 0, 1);
          }
        }
      }
      */
      vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
      threshold->ThresholdByUpper(l-w/2);
      threshold->SetInValue(1);
      threshold->SetOutValue(0);
      threshold->ReplaceInOn();
      threshold->ReplaceOutOn();
      threshold->SetOutputScalarType(m_imageData->GetScalarType());
#if VTK_MAJOR_VERSION > 5
      threshold->SetInputData(m_imageDataRef);
#else
      threshold->SetInput(m_imageDataRef);
#endif
      threshold->Update();
      m_imageData->DeepCopy(threshold->GetOutput());
    }
    
    InitializeActors();
    
    GetProperty()->SetVolumeSource( m_volumeSource );
    
    m_sFilename = "";
    
    if ( bCopyVoxelData )
    {
      GetProperty()->CopySettings( mri->GetProperty() );
      SetModified();
    }
    else if (voxel_option == 2) // mask
    {
      SetModified();
      GetProperty()->SetWindowLevel(0, 0.5);
    }
  }
  
  return true;
}

void LayerMRI::SetReorient( bool bReorient )
{
  m_bReorient = bReorient;
}

void LayerMRI::SetSampleMethod( int nSampleMethod )
{
  m_nSampleMethod = nSampleMethod;
  if (m_volumeSource)
    m_volumeSource->SetInterpolationMethod(nSampleMethod);
}

bool LayerMRI::SaveVolume()
{
  if ( m_sFilename.isEmpty() || m_imageData == NULL )
  {
    return false;
  }
  
  ::SetProgressCallback(ProgressCallback, 0, 60);
  m_volumeSource->setProperty("label_value", property("label_value"));
  if ( !m_volumeSource->UpdateMRIFromImage( m_imageData, !m_bReorient ) )
  {
    m_volumeSource->setProperty("label_value", 0);
    return false;
  }
  
  // now first save a copy of the old file if requested
  if ( MainWindow::GetMainWindow()->GetSaveCopy() && property("label_fn").toString().isEmpty())
  {
    QString new_fn = m_sFilename + "~";
    QFile::remove( new_fn );
    QFile::rename( m_sFilename, new_fn );
  }
  
  QString fn = m_sFilename;
  if (!property("label_fn").toString().isEmpty())
    fn = property("label_fn").toString();

  ::SetProgressCallback(ProgressCallback, 60, 100);
  int nSampleMethod = GetProperty()->GetResliceInterpolation();
  bool bSaved = m_volumeSource->MRIWrite( fn.toLatin1().data(),
                                          nSampleMethod,
                                          m_bWriteResampled);
  m_bModified = !bSaved;
  setProperty("saved_name", fn);
  setProperty("label_fn", "");
  setProperty("label_value", 0);
  
  return bSaved;
}

bool LayerMRI::IsTransformed()
{
  vtkMatrix4x4* mat = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() )->GetMatrix();
  return !MyUtils::IsIdentity( mat->Element, qMin(m_dWorldVoxelSize[0], qMin(m_dWorldVoxelSize[1], m_dWorldVoxelSize[2])) );
}

void LayerMRI::DoRestore()
{
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  slice_tr->Identity();
  ras_tr->Identity();
  for ( int i = 0; i < 3; i++ )
  {
    mReslice[i]->SetInterpolationModeToNearestNeighbor();
    mReslice[i]->Modified();
  }
}

void LayerMRI::DoTransform(double *m, int sample_method)
{
  if ( GetProperty()->GetColorMap() == LayerPropertyMRI::LUT || sample_method == SAMPLE_NEAREST )
    GetProperty()->SetResliceInterpolation(SAMPLE_NEAREST);
  else
    GetProperty()->SetResliceInterpolation(sample_method);
  
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
  mat->DeepCopy(m);
  vtkSmartPointer<vtkMatrix4x4> mat_inv = vtkSmartPointer<vtkMatrix4x4>::New();
  mat_inv->DeepCopy(mat);
  mat_inv->Invert();
  slice_tr->Concatenate(mat_inv); // reslice uses inverse transformation matrix
  
  // also record transformation in RAS space
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  MATRIX* t2ras = GetSourceVolume()->GetTargetToRASMatrix();
  double m_t2r[16], m_a[16];
  for ( int i = 0; i < 16; i++ )
  {
    m_t2r[i] = (double) *MATRIX_RELT((t2ras),(i/4)+1,(i%4)+1);
  }
  vtkMatrix4x4::Multiply4x4(m_t2r, m, m_a);
  vtkMatrix4x4::Invert(m_t2r, m_t2r);
  vtkMatrix4x4::Multiply4x4(m_a, m_t2r, m_a);
  mat->DeepCopy(m_a);
  ras_tr->Concatenate(mat);
  
  MatrixFree(&t2ras);
  for ( int i = 0; i < 3; i++ )
    mReslice[i]->Modified();
}

void LayerMRI::DoTransform(int sample_method)
{
  Q_UNUSED(sample_method);
  bool bTransformed = IsTransformed();
  UpdateResliceInterpolation();
  vtkTransform* slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  // also record transformation in RAS space
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  slice_tr->Identity();
  ras_tr->Identity();
  
  double pos[3];
  if (bTransformed || m_bUseRotationCenter)
  {
    pos[0] = m_dRotationCenter[0];
    pos[1] = m_dRotationCenter[1];
    pos[2] = m_dRotationCenter[2];
  }
  else
  {
    if (m_bRotateAroundCenter)
    {
      pos[0] = m_dWorldOrigin[0] + m_dWorldSize[0]/2;
      pos[1] = m_dWorldOrigin[1] + m_dWorldSize[1]/2;
      pos[2] = m_dWorldOrigin[2] + m_dWorldSize[2]/2;
    }
    else
    {
      pos[0] = m_dSlicePosition[0];
      pos[1] = m_dSlicePosition[1];
      pos[2] = m_dSlicePosition[2];
    }
    m_dRotationCenter[0] = pos[0];
    m_dRotationCenter[1] = pos[1];
    m_dRotationCenter[2] = pos[2];
  }
  
  // flip first
  for (int i = 0; i < 3; i++)
  {
    if (m_bFlip[i])
    {
      double dTargetPoint[3] = { pos[0], pos[1], pos[2] };
      double m[16] = { 1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1 };
      m[i*4+i] = -1;
      double ras[3];
      this->TargetToRAS(pos, ras);
      slice_tr->Translate( dTargetPoint[0], dTargetPoint[1], dTargetPoint[2] );
      slice_tr->Concatenate(m);
      slice_tr->Translate( -dTargetPoint[0], -dTargetPoint[1], -dTargetPoint[2] );
      
      ras_tr->Translate( -ras[0], -ras[1], -ras[2] );
      for (int j = 0; j < 16; j++)
        m[j] = -m[j];
      ras_tr->Concatenate(m);
      ras_tr->Translate( ras );
    }
  }
  
  // then scale
  double cpt[3], target_cpt[3];
  GetRASCenter( cpt );
  RASToTarget( cpt, target_cpt );
  slice_tr->Translate( target_cpt[0], target_cpt[1], target_cpt[2] );
  slice_tr->Scale( 1.0/m_dScale[0], 1.0/m_dScale[1], 1.0/m_dScale[2] );
  slice_tr->Translate( -target_cpt[0], -target_cpt[1], -target_cpt[2] );
  
  ras_tr->Translate( -cpt[0], -cpt[1], -cpt[2] );
  ras_tr->Scale( m_dScale );
  ras_tr->Translate( cpt[0], cpt[1], cpt[2] );
  
  // rotate
  for ( int i = 0; i < 3; i++ )
  {
    double v[3] = { 0, 0, 0 };
    v[i] = 1;
    double dTargetPoint[3] = { pos[0], pos[1], pos[2] };
    double ras[3];
    this->TargetToRAS(pos, ras);
    
    slice_tr->Translate( dTargetPoint[0], dTargetPoint[1], dTargetPoint[2] );
    slice_tr->RotateWXYZ( -m_dRotate[i], v );
    slice_tr->Translate( -dTargetPoint[0], -dTargetPoint[1], -dTargetPoint[2] );
    
    // record transformation in RAS space
    ras_tr->Translate( -ras[0], -ras[1], -ras[2] );
    double vp[3];
    vp[0] = dTargetPoint[0] + v[0];
    vp[1] = dTargetPoint[1] + v[1];
    vp[2] = dTargetPoint[2] + v[2];
    TargetToRAS( vp, vp );
    v[0] = vp[0] - ras[0];
    v[1] = vp[1] - ras[1];
    v[2] = vp[2] - ras[2];
    ras_tr->RotateWXYZ( m_dRotate[i], v );
    ras_tr->Translate( ras );
  }
  
  // translate
  slice_tr->Translate( -m_dTranslate[0], -m_dTranslate[1], -m_dTranslate[2] );
  ras_tr->Translate( m_dTranslate );
  
  // new approach, directly convert to ras_tr
  ras_tr->Identity();
  m_volumeSource->ConvertTransformFromTargetToRAS(slice_tr->GetMatrix(), ras_tr->GetMatrix());
  
  for ( int i = 0; i < 3; i++ )
  {
    mReslice[i]->Modified();
  }
}

bool LayerMRI::DoRotate( std::vector<RotationElement>& rotations )
{
  if ( GetProperty()->GetColorMap() == LayerPropertyMRI::LUT || rotations[0].SampleMethod == SAMPLE_NEAREST)
    GetProperty()->SetResliceInterpolation(SAMPLE_NEAREST);
  else
    GetProperty()->SetResliceInterpolation(rotations[0].SampleMethod);
  
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  // also record transformation in RAS space
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  
  for ( size_t i = 0; i < rotations.size(); i++ )
  {
    double v[3] = { 0, 0, 0 };
    v[rotations[i].Plane] = 1;
    double dTargetPoint[3];
    RASToTarget( rotations[i].Point, dTargetPoint );
    
    slice_tr->Translate( dTargetPoint[0], dTargetPoint[1], dTargetPoint[2] );
    slice_tr->RotateWXYZ( -rotations[i].Angle, v );
    slice_tr->Translate( -dTargetPoint[0], -dTargetPoint[1], -dTargetPoint[2] );
    
    // record transformation in RAS space
    ras_tr->Translate( -rotations[i].Point[0], -rotations[i].Point[1], -rotations[i].Point[2] );
    double vp[3];
    vp[0] = dTargetPoint[0] + v[0];
    vp[1] = dTargetPoint[1] + v[1];
    vp[2] = dTargetPoint[2] + v[2];
    TargetToRAS( vp, vp );
    v[0] = vp[0] - rotations[i].Point[0];
    v[1] = vp[1] - rotations[i].Point[1];
    v[2] = vp[2] - rotations[i].Point[2];
    ras_tr->RotateWXYZ( rotations[i].Angle, v );
    ras_tr->Translate( rotations[i].Point[0], rotations[i].Point[1], rotations[i].Point[2] );
  }
  
  for ( int i = 0; i < 3; i++ )
  {
    mReslice[i]->Modified();
  }
  
  return true;
}

void LayerMRI::DoTranslate( double* offset )
{
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  slice_tr->Translate( -offset[0], -offset[1], -offset[2] );
  ras_tr->Translate( offset );
  for ( int i = 0; i < 3; i++ )
  {
    mReslice[i]->Modified();
  }
}

void LayerMRI::DoScale( double* scale, int nSampleMethod )
{
  if ( GetProperty()->GetColorMap() == LayerPropertyMRI::LUT || nSampleMethod == SAMPLE_NEAREST )
    GetProperty()->SetResliceInterpolation(SAMPLE_NEAREST);
  else
    GetProperty()->SetResliceInterpolation(nSampleMethod);
  
  vtkSmartPointer<vtkTransform> slice_tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
  vtkTransform* ras_tr = m_volumeSource->GetTransform();
  double cpt[3], target_cpt[3];
  GetRASCenter( cpt );
  RASToTarget( cpt, target_cpt );
  slice_tr->Translate( target_cpt[0], target_cpt[1], target_cpt[2] );
  slice_tr->Scale( 1.0/scale[0], 1.0/scale[1], 1.0/scale[2] );
  slice_tr->Translate( -target_cpt[0], -target_cpt[1], -target_cpt[2] );
  
  ras_tr->Translate( -cpt[0], -cpt[1], -cpt[2] );
  ras_tr->Scale( scale );
  ras_tr->Translate( cpt[0], cpt[1], cpt[2] );
  
  for ( int i = 0; i < 3; i++ )
  {
    mReslice[i]->Modified();
  }
}

void LayerMRI::InitializeVolume()
{
  if ( m_volumeSource == NULL )
  {
    return;
  }
  
  FSVolume* source = m_volumeSource;
  
  float RASBounds[6];
  source->GetBounds( RASBounds );
  m_dWorldOrigin[0] = RASBounds[0];
  m_dWorldOrigin[1] = RASBounds[2];
  m_dWorldOrigin[2] = RASBounds[4];
  source->GetPixelSize( m_dWorldVoxelSize );
  
  m_dWorldSize[0] = ( ( int )( (RASBounds[1] - RASBounds[0]) / m_dWorldVoxelSize[0] ) ) * m_dWorldVoxelSize[0];
  m_dWorldSize[1] = ( ( int )( (RASBounds[3] - RASBounds[2]) / m_dWorldVoxelSize[1] ) ) * m_dWorldVoxelSize[1];
  m_dWorldSize[2] = ( ( int )( (RASBounds[5] - RASBounds[4]) / m_dWorldVoxelSize[2] ) ) * m_dWorldVoxelSize[2];
  //  qDebug() << RASBounds[0] << RASBounds[1] << RASBounds[2] << RASBounds[3] << RASBounds[4] << RASBounds[5];
  //  qDebug() << m_dWorldSize[0] << m_dWorldSize[1] << m_dWorldSize[2];
  
  m_imageData = source->GetImageOutput();
}


void LayerMRI::InitializeActors()
{
  vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
  tr->Identity();
  m_actorContour->SetUserTransform( tr->GetLinearInverse() );
  for ( int i = 0; i < 3; i++ )
  {
    // The reslice object just takes a slice out of the volume.
    //
    mReslice[i] = vtkSmartPointer<vtkImageReslice>::New();
#if VTK_MAJOR_VERSION > 5
    mReslice[i]->SetInputData( m_imageData );
#else
    mReslice[i]->SetInput( m_imageData );
#endif
    mReslice[i]->BorderOn();
    mReslice[i]->SetResliceTransform( tr );
    mReslice[i]->AutoCropOutputOn();
    
    // This sets us to extract slices.
    mReslice[i]->SetOutputDimensionality( 2 );
    
    // This will change depending what orienation we're in.
    mReslice[i]->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 1, 0,
                                                 0, 0, 1 );
    
    // This will change to select a different slice.
    mReslice[i]->SetResliceAxesOrigin( 0, 0, 0 );
    mReslice[i]->SetInterpolationModeToNearestNeighbor();
    
    //
    // Image to colors using color table.
    //
    mColorMap[i] = vtkSmartPointer<vtkImageMapToColors>::New();
    mColorMap[i]->SetLookupTable( GetProperty()->GetGrayScaleTable() );
    mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
    mColorMap[i]->SetOutputFormatToRGBA();
    mColorMap[i]->PassAlphaToOutputOff();
    
    //
    // Prop in scene with plane mesh and texture.
    //
    m_sliceActor2D[i]->GetMapper()->SetInputConnection(mColorMap[i]->GetOutputPort());
    m_sliceActor3D[i]->GetMapper()->SetInputConnection( mColorMap[i]->GetOutputPort());
    
    mEdgeFilter[i] = vtkSmartPointer<vtkSimpleLabelEdgeFilter>::New();
    mResample[i] = vtkSmartPointer<vtkImageReslice>::New();
    //    mResample[i]->SetAxisMagnificationFactor( 0, IMAGE_RESAMPLE_FACTOR );
    //    mResample[i]->SetAxisMagnificationFactor( 1, IMAGE_RESAMPLE_FACTOR );
    //    mResample[i]->SetAxisMagnificationFactor( 2, IMAGE_RESAMPLE_FACTOR );
    mResample[i]->SetInterpolationModeToNearestNeighbor();
    
    // Set ourselves up.
    this->OnSlicePositionChanged( i );
  }
  
  this->blockSignals( true );
  //  this->UpdateResliceInterpolation();
  this->UpdateTextureSmoothing();
  this->UpdateOpacity();
  this->UpdateColorMap();
  this->blockSignals( false );
}

void LayerMRI::UpdateOpacity()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetOpacity( GetProperty()->GetOpacity() );
    m_sliceActor3D[i]->SetOpacity( GetProperty()->GetOpacity() );
    m_projectionMapActor[i]->SetOpacity(GetProperty()->GetOpacity());
  }
  m_actorContour->GetProperty()->SetOpacity( GetProperty()->GetOpacity() );
  QList<int> keys = m_labelActors.keys();
  foreach (int i, keys)
    m_labelActors[i]->GetProperty()->SetOpacity(GetProperty()->GetOpacity());
  emit ActorUpdated();
}

void LayerMRI::UpdateColorMap()
{
  assert( GetProperty() );
  
  if ( !mColorMap[0].GetPointer() )
  {
    return;
  }
  
  for ( int i = 0; i < 3; i++ )
  {
    mColorMap[i]->SetActiveComponent( m_nActiveFrame );
  }
  
  for ( int i = 0; i < 3; i++ )
  {
    mColorMap[i]->SetLookupTable( GetProperty()->GetActiveLookupTable() );
    if (mColorMapMaxProjection[i].GetPointer() != NULL)
      mColorMapMaxProjection[i]->SetLookupTable(GetProperty()->GetActiveLookupTable());
  }
  
  m_actorContour->GetMapper()->SetLookupTable( GetProperty()->GetActiveLookupTable() );
  emit ActorUpdated();
  
  if (this->m_nAvailableLabels.isEmpty() || m_listLabelCenters.isEmpty())
    UpdateLabelInformation();
}

void LayerMRI::UpdateLabelInformation()
{
  if (GetProperty()->GetColorMap() == LayerPropertyMRI::LUT)
  {
    if (m_worker->isRunning())
    {
      m_worker->Abort();
      QTimer::singleShot(100, m_worker, SLOT(start()));
    }
    else
      m_worker->start();
  }
}

void LayerMRI::UpdateResliceInterpolation ()
{
  assert( GetProperty() );
  
  for ( int i = 0; i < 3; i++ )
  {
    if ( mReslice[i].GetPointer() )
    {
      switch (GetProperty()->GetResliceInterpolation())
      {
      case SAMPLE_TRILINEAR:
        mReslice[i]->SetInterpolationModeToLinear();
        break;
      case SAMPLE_CUBIC_BSPLINE:
        mReslice[i]->SetInterpolationModeToCubic();
        break;
      default:
        mReslice[i]->SetInterpolationModeToNearestNeighbor();
        break;
      }
    }
  }
  emit ActorUpdated();
}

void LayerMRI::UpdateTextureSmoothing ()
{
  assert( GetProperty() );
  
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetInterpolate( GetProperty()->GetTextureSmoothing() );
    m_sliceActor3D[i]->SetInterpolate( GetProperty()->GetTextureSmoothing() );
  }
  emit ActorUpdated();
}

void LayerMRI::UpdateContour( int nSegValue )
{
  if ( this->GetProperty()->GetShowAsContour() )
  {
    UpdateContourActor( nSegValue );
  }
}

void LayerMRI::UpdateContourActor( int nSegValue )
{
  // Generate a new thread id before creating the thread. so that mainwindow will be able to determine
  // if a build contour result is already expired, by comparing the returned id and current id. If they
  // are different, it means a new thread is rebuilding the contour
  m_nThreadID++;
  ThreadBuildContour* thread = new ThreadBuildContour(this);
  connect(thread, SIGNAL(Finished(int)), this, SLOT(OnContourThreadFinished(int)));
  thread->BuildContour( this, nSegValue, m_nThreadID );
  emit IsoSurfaceUpdating();
}

// Contour mapper is ready, attach it to the actor
void LayerMRI::OnContourThreadFinished(int thread_id)
{
  if (m_nThreadID == thread_id)
  {
    if (GetProperty()->GetShowAsLabelContour())
    {
      QList<int> labels = m_labelActorsTemp.keys();
      if (!labels.isEmpty())
      {
        foreach (int n, labels)
        {
          m_labelActors[n] = m_labelActorsTemp[n];
#if VTK_MAJOR_VERSION > 5
          m_labelActors[n]->ForceTranslucentOn();
#endif
          m_labelActors[n]->GetMapper()->SetLookupTable( GetProperty()->GetLUTTable() );
        }
        OnLabelContourChanged();
        emit ActorChanged();
      }
    }
    else if( m_actorContourTemp.GetPointer() && m_actorContourTemp->GetMapper() )
    {
      m_actorContour->SetMapper( m_actorContourTemp->GetMapper() );
      UpdateContourColor();
      emit ActorChanged();
    }
    emit IsoSurfaceUpdated();
  }
}

void LayerMRI::ShowContour()
{
  if ( !m_actorContourTemp.GetPointer() && this->GetProperty()->GetShowAsContour() )
  {
    UpdateContour();
  }
  emit ActorChanged();
}

void LayerMRI::UpdateContourColor()
{
  if ( GetProperty()->GetContourUseImageColorMap() || GetProperty()->GetShowAsLabelContour())
  {
    m_actorContour->GetMapper()->ScalarVisibilityOn();
  }
  else
  {
    m_actorContour->GetMapper()->ScalarVisibilityOff();
    m_actorContour->GetProperty()->SetColor( GetProperty()->GetContourColor() );
  }
  UpdateColorMap();
}

void LayerMRI::UpdateVolumeRendering()
{
  /*
  if ( GetProperty()->GetShowAsContour() )
  {
    MyUtils::BuildVolume( GetImageData(),
                          GetProperty()->GetContourMinThreshold(),
                          GetProperty()->GetContourMaxThreshold(),
                          m_propVolume );
  }
  */
}

void LayerMRI::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  if ( GetProperty()->GetDisplayVector() || GetProperty()->GetDisplayTensor() )
  {
    renderer->AddViewProp( m_glyphActor2D[nPlane] );
    if (GetProperty()->GetVectorRepresentation() == LayerPropertyMRI::VR_Direction_Line)
      renderer->AddViewProp(m_vectorDotActor2D[nPlane]);
  }
  else
  {
    renderer->AddViewProp(m_sliceActor2D[nPlane]);
    renderer->AddViewProp(m_projectionMapActor[nPlane]);
  }
}

void LayerMRI::Remove2DProps( vtkRenderer* renderer, int nPlane )
{
  if ( GetProperty()->GetDisplayVector() || GetProperty()->GetDisplayTensor() )
  {
    renderer->RemoveViewProp( m_glyphActor2D[nPlane] );
    renderer->RemoveViewProp(m_vectorDotActor2D[nPlane]);
  }
  else
  {
    renderer->RemoveViewProp( m_sliceActor2D[nPlane] );
  }
}

void LayerMRI::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  bool bContour = GetProperty()->GetShowAsContour();
  for ( int i = 0; i < 3; i++ )
  {
    if ( !bContour && ( bSliceVisibility == NULL || bSliceVisibility[i] ) )
    {
      if ( GetProperty()->GetDisplayVector() || GetProperty()->GetDisplayTensor() )
      {
        renderer->AddViewProp( m_glyphActor3D[i] );
      }
      else
      {
        renderer->AddViewProp( m_sliceActor3D[i] );
      }
    }
  }
  
  if ( bContour )
  {
    if (GetProperty()->GetShowAsLabelContour())
    {
      QList<int> keys = m_labelActors.keys();
      foreach (int i, keys)
        renderer->AddViewProp(m_labelActors[i]);
    }
    else
    {
      renderer->AddViewProp( m_actorContour );
    }

    for ( int i = 0; i < m_surfaceRegions.size(); i++ )
    {
      m_surfaceRegions[i]->AppendProps( renderer );
    }

    for (int i = 0; i < m_3DRegions.size(); i++)
    {
      m_3DRegions[i]->AppendProps( renderer );
    }
  }
}

void LayerMRI::SetSlicePositionToWorldCenter()
{
  if ( m_volumeSource == NULL )
  {
    return;
  }
  
  // Get some values from the MRI.
  double pos[3];
  for ( int i = 0; i < 3; i++ )
  {
    pos[i] = ((int)( m_dWorldSize[i]/2/m_dWorldVoxelSize[i] ) + 0.0 ) * m_dWorldVoxelSize[i] + m_dWorldOrigin[i];
  }
  
  SetSlicePosition( pos );
}

void LayerMRI::OnSlicePositionChanged( int nPlane )
{
  if ( m_volumeSource == NULL || nPlane < 0 || nPlane > 2)
  {
    return;
  }
  
  assert( GetProperty() );
  
  // display slice image
  vtkSmartPointer<vtkMatrix4x4> matrix =
      vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  double origin[3] = {0};
  m_imageData->GetOrigin(origin);
  switch ( nPlane )
  {
  case 0:
    m_sliceActor2D[0]->PokeMatrix( matrix );
    m_sliceActor2D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
    m_sliceActor2D[0]->RotateX( 90 );
    m_sliceActor2D[0]->RotateY( 90 );
    m_sliceActor3D[0]->PokeMatrix( matrix );
    m_sliceActor3D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
    m_sliceActor3D[0]->RotateX( 90 );
    m_sliceActor3D[0]->RotateY( 90 );
    m_projectionMapActor[0]->PokeMatrix( matrix );
    m_projectionMapActor[0]->SetPosition( m_dSlicePosition[0], origin[1], origin[2] );
    m_projectionMapActor[0]->RotateX( 90 );
    m_projectionMapActor[0]->RotateY( 90 );
    
    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
    mReslice[0]->SetResliceAxesDirectionCosines( 0, 1, 0,
                                                 0, 0, 1,
                                                 1, 0, 0 );
    mReslice[0]->SetResliceAxesOrigin( m_dSlicePosition[0], 0, 0  );
    mReslice[0]->Modified();
    break;
  case 1:
    m_sliceActor2D[1]->PokeMatrix( matrix );
    m_sliceActor2D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
    m_sliceActor2D[1]->RotateX( 90 );
    m_sliceActor3D[1]->PokeMatrix( matrix );
    m_sliceActor3D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
    m_sliceActor3D[1]->RotateX( 90 );
    m_projectionMapActor[1]->PokeMatrix( matrix );
    m_projectionMapActor[1]->SetPosition( origin[0], m_dSlicePosition[1], origin[2] );
    m_projectionMapActor[1]->RotateX( 90 );
    
    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
    mReslice[1]->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 0, 1,
                                                 0, 1, 0 );
    mReslice[1]->SetResliceAxesOrigin( 0, m_dSlicePosition[1], 0 );
    mReslice[1]->Modified();
    break;
  case 2:
    m_sliceActor2D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
    // m_sliceActor2D[2]->RotateY( 180 );
    m_sliceActor3D[2]->SetPosition( 0, 0, m_dSlicePosition[2] );
    // m_sliceActor3D[2]->RotateY( 180 );
    m_projectionMapActor[2]->SetPosition( origin[0], origin[1], m_dSlicePosition[2] );
    
    mReslice[2]->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 1, 0,
                                                 0, 0, 1 );
    mReslice[2]->SetResliceAxesOrigin( 0, 0, m_dSlicePosition[2] );
    mReslice[2]->Modified();
    break;
  }
  // display 4D data as vector
  if ( GetProperty()->GetDisplayVector() )
  {
    UpdateVectorActor( nPlane );
  }
  else if ( GetProperty()->GetDisplayTensor() )
  {
    UpdateTensorActor( nPlane );
  }
  else if ( /*GetProperty()->GetColorMap() == LayerPropertyMRI::LUT &&*/ GetProperty()->GetShowLabelOutline() )
  {
    UpdateLabelOutline();
  }
}

void LayerMRI::UpdateDisplayMode()
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( GetProperty()->GetDisplayTensor() &&
         GetProperty()->GetTensorRepresentation() == LayerPropertyMRI::TR_Ellipsoid )
    {
      m_glyphActor2D[i]->GetProperty()->SetInterpolationToGouraud();
      m_glyphActor3D[i]->GetProperty()->SetInterpolationToGouraud();
    }
    else
    {
      m_glyphActor2D[i]->GetProperty()->SetInterpolationToFlat();
      m_glyphActor3D[i]->GetProperty()->SetInterpolationToFlat();
    }
    if (GetProperty()->GetDisplayRGB())
    {
      vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
      vtkSmartPointer<vtkImageThreshold> upper = vtkSmartPointer<vtkImageThreshold>::New();
      upper->ThresholdByUpper(0);
      upper->SetOutValue(0);
      upper->SetReplaceOut(1);
      upper->SetInputConnection(mReslice[i]->GetOutputPort());
      vtkSmartPointer<vtkImageThreshold> lower = vtkSmartPointer<vtkImageThreshold>::New();
      lower->ThresholdByLower(255);
      lower->SetOutValue(255);
      lower->SetReplaceOut(1);
      lower->SetInputConnection(upper->GetOutputPort());
      cast->SetInputConnection(lower->GetOutputPort());
      cast->SetOutputScalarTypeToUnsignedChar();
      m_sliceActor2D[i]->GetMapper()->SetInputConnection( cast->GetOutputPort() );
      m_sliceActor3D[i]->GetMapper()->SetInputConnection( cast->GetOutputPort() );
    }
    else
    {
      m_sliceActor2D[i]->GetMapper()->SetInputConnection( mColorMap[i]->GetOutputPort() );
      m_sliceActor3D[i]->GetMapper()->SetInputConnection( mColorMap[i]->GetOutputPort() );
    }
  }
  if ( GetProperty()->GetDisplayVector() )
  {
    this->UpdateVectorActor();
  }
  else if ( GetProperty()->GetDisplayTensor() )
  {
    this->UpdateTensorActor();
  }
  else
  {
    emit ActorChanged();
  }
}

void LayerMRI::SetVisible( bool bVisible )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetVisibility( bVisible && !GetProperty()->GetShowProjectionMap() );
    m_sliceActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_glyphActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_glyphActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_vectorDotActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
    
    m_projectionMapActor[i]->SetVisibility( bVisible && GetProperty()->GetShowProjectionMap() );
  }
  m_actorContour->SetVisibility( bVisible ? 1 : 0 );
  if (GetProperty()->GetShowAsContour() && GetProperty()->GetShowAsLabelContour())
    OnLabelContourChanged();
  LayerVolumeBase::SetVisible(bVisible);
}

bool LayerMRI::IsVisible()
{
  return (m_sliceActor2D[0]->GetVisibility() || m_projectionMapActor[0]->GetVisibility());
}

double LayerMRI::GetVoxelValue( double* pos )
{
  if ( m_imageData == NULL )
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
    if (GetCorrelationSurface())
      return m_imageRawDisplay->GetScalarComponentAsDouble( n[0], n[1], n[2], 0 );
    else if (m_layerMask)
      return m_imageDataBackup->GetScalarComponentAsDouble( n[0], n[1], n[2], m_nActiveFrame );
    else
      return m_imageData->GetScalarComponentAsDouble( n[0], n[1], n[2], m_nActiveFrame );
  }
}

double LayerMRI::GetVoxelValueByOriginalIndex( int i, int j, int k, int frame )
{
  if (frame < 0)
    frame = m_nActiveFrame;
  return m_volumeSource->GetVoxelValue( i, j, k, frame );
}

// trilinear interpolated value
double LayerMRI::GetSampledVoxelValueByRAS(double* ras, int frame)
{
  if (frame < 0)
    frame = m_nActiveFrame;
  MRI* mri = m_volumeSource->GetMRI();
  float fx, fy, fz;
  m_volumeSource->RASToOriginalIndex(ras[0], ras[1], ras[2], fx, fy, fz);
  double val = -1;
  if (GetNumberOfFrames() > 1)
  {
    float fval;
    MRIsampleSeqVolume(mri, fx, fy, fz, &fval, frame, frame);
    val = fval;
  }
  else
    MRIsampleVolume(mri, fx, fy, fz, &val);
  return val;
}

std::vector<double> LayerMRI::GetSampledVoxelValues(std::vector<std::vector<double> > &line3d, int frame)
{
  std::vector<double> vals;
  for (size_t i = 0; i < line3d.size(); i++)
  {
    double pt[3] = { line3d[i][0], line3d[i][1], line3d[i][2] };
    this->TargetToRAS(pt, pt);
    vals.push_back(GetSampledVoxelValueByRAS(pt, frame));
  }
  return vals;
}

std::vector<double> LayerMRI::GetMeanSegmentValues(std::vector<std::vector<double> > &line3d, int frame)
{
  double* spacing = this->m_imageData->GetSpacing();
  double voxel_length = qMin(spacing[0], qMin(spacing[1], spacing[2]));
  std::vector<double> vals;
  for (size_t i = 0; i < line3d.size()-1; i++)
  {
    double pt1[3] = { line3d[i][0], line3d[i][1], line3d[i][2] };
    double pt2[3] = { line3d[i+1][0], line3d[i+1][1], line3d[i+1][2] };
    this->TargetToRAS(pt1, pt1);
    this->TargetToRAS(pt2, pt2);
    double dist = sqrt(vtkMath::Distance2BetweenPoints(pt1, pt2));
    if (dist == 0)
      dist = 1.0;   // not going to happen
    double v[3];
    for (int j = 0; j < 3; j++)
      v[j] = (pt2[j] - pt1[j]) / dist;
    int n = (int)(dist/voxel_length/2.0+0.5);
    if (n < 1)
      n = 1;
    dist = dist/n;
    double sum = 0;
    for (int j = 0; j <= n; j++)
    {
      double pt[3];
      pt[0] = pt1[0] + v[0]*j*dist;
      pt[1] = pt1[1] + v[1]*j*dist;
      pt[2] = pt1[2] + v[2]*j*dist;
      double val = GetSampledVoxelValueByRAS(pt, frame);
      sum += val;
      if (j != 0 && j != n)
        sum += val;
    }
    vals.push_back(sum/n/2);
  }
  return vals;
}

QList<double> LayerMRI::GetVoxelValueByOriginalIndexAllFrames(int i, int j, int k)
{
  QList<double> list;
  for (int frame = 0; frame < GetNumberOfFrames(); frame++)
    list << m_volumeSource->GetVoxelValue( i, j, k, frame );
  return list;
}

void LayerMRI::GetVoxelValueByOriginalIndexAllFrames(int i, int j, int k, float* buffer)
{
  for (int frame = 0; frame < GetNumberOfFrames(); frame++)
    buffer[frame] = m_volumeSource->GetVoxelValue( i, j, k, frame );
}

void LayerMRI::SetModified()
{
  mReslice[0]->Modified();
  mReslice[1]->Modified();
  mReslice[2]->Modified();
  
  LayerVolumeBase::SetModified();
}

QString LayerMRI::GetLabelName( double value )
{
  int nIndex = (int)value;
  if ( GetProperty()->GetColorMap() == LayerPropertyMRI::LUT )
  {
    COLOR_TABLE* ct = GetProperty()->GetLUTCTAB();
    if ( !ct )
    {
      return "";
    }
    char name[128];
    int nValid = 0;
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( ct, &nTotalCount );
    if ( nIndex > 0 && nIndex < nTotalCount )
    {
      CTABisEntryValid( ct, nIndex, &nValid );
      if ( nValid && CTABcopyName( ct, nIndex, name, 128 ) == 0 )
      {
        return name;
      }
    }
  }
  
  return "";
}

void LayerMRI::RemapPositionToRealRAS( const double* pos_in, double* pos_out )
{
  m_volumeSource->TargetToRAS( pos_in, pos_out );
}

void LayerMRI::RemapPositionToRealRAS( double x_in, double y_in, double z_in,
                                       double& x_out, double& y_out, double& z_out )
{
  m_volumeSource->TargetToRAS( x_in, y_in, z_in, x_out, y_out, z_out );
}

void LayerMRI::RASToTarget( const double* pos_in, double* pos_out )
{
  m_volumeSource->RASToTarget( pos_in, pos_out );
}

void LayerMRI::NativeRASToTkReg( const double* pos_in, double* pos_out )
{
  m_volumeSource->NativeRASToTkReg( pos_in, pos_out );
}

void LayerMRI::TkRegToNativeRAS( const double* pos_in, double* pos_out )
{
  m_volumeSource->TkRegToNativeRAS( pos_in, pos_out );
}

int LayerMRI::GetNumberOfFrames()
{
  if ( m_imageData )
  {
    return m_imageData->GetNumberOfScalarComponents();
  }
  else
  {
    return 1;
  }
}

void LayerMRI::SetActiveFrame( int nFrame )
{
  if ( nFrame != m_nActiveFrame && nFrame >= 0 && nFrame < this->GetNumberOfFrames() )
  {
    m_nActiveFrame = nFrame;
    m_listLabelCenters.clear();
    GetProperty()->UpdateActiveFrame(nFrame);
    UpdateColorMap();
    emit ActiveFrameChanged( nFrame );
    emit ActorUpdated();
  }
}

bool LayerMRI::HasProp( vtkProp* prop )
{
  if ( GetProperty()->GetShowAsContour() )
  {
    if ( m_actorContour.GetPointer() == prop )
    {
      return true;
    }
    else if (GetProperty()->GetShowAsLabelContour())
    {
      QList<int> ids = m_labelActors.keys();
      foreach (int id, ids)
      {
        if (m_labelActors[id] == prop)
          return true;
      }
    }

    for ( int i = 0; i < m_surfaceRegions.size(); i++ )
    {
      if ( m_surfaceRegions[i]->GetMeshActor() == prop )
      {
        return true;
      }
    }
    return false;
  }
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
      if ( m_sliceActor3D[i] == prop )
      {
        return true;
      }
    }
    return false;
  }
}

void LayerMRI::RASToOriginalIndex( const double* pos, int* n )
{
  m_volumeSource->RASToOriginalIndex( (float)(pos[0]), (float)(pos[1]), (float)(pos[2]),
      n[0], n[1], n[2] );
}

void LayerMRI::RASToOriginalIndex(const double *pos, double *n_out)
{
  float fPos[3] = {(float)pos[0], (float)pos[1], (float)pos[2]};
  float fout0, fout1, fout2;
  m_volumeSource->RASToOriginalIndex(fPos[0], fPos[1], fPos[2], fout0, fout1, fout2);
  n_out[0] = fout0;
  n_out[1] = fout1;
  n_out[2] = fout2;
}

void LayerMRI::OriginalIndexToRAS( const int* n, double* pos )
{
  float x, y, z;
  m_volumeSource->OriginalIndexToRAS( n[0], n[1], n[2], x, y, z );
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
}

void LayerMRI::OriginalVoxelToRAS(const double *vcoord, double *pos)
{
  float x, y, z;
  m_volumeSource->OriginalIndexToRAS( vcoord[0], vcoord[1], vcoord[2], x, y, z );
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
}

void LayerMRI::TargetIndexToOriginalIndex(const int *n_in, int *n_out)
{
  double pos[3];
  double* orig = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for (int i = 0; i < 3; i++)
  {
    pos[i] = orig[i] + voxel_size[i]*n_in[i];
  }
  TargetToRAS(pos, pos);
  RASToOriginalIndex(pos, n_out);
}

void LayerMRI::UpdateVectorLineWidth(double val)
{
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  for ( int i = 0; i < 3; i++ )
  {
    m_glyphActor2D[i]->GetProperty()->SetLineWidth(val*ratio);
    m_vectorDotActor2D[i]->GetProperty()->SetPointSize(val*(val>1?2:3)*ratio);
  }
  if (GetProperty()->GetVectorRepresentation() != LayerPropertyMRI::VR_Bar)
  {
    emit ActorUpdated();
  }
  else
  {
    UpdateVectorActor();
  }
}

void LayerMRI::ReorderColorComponent(unsigned char *c)
{
}

void LayerMRI::UpdateVectorActor()
{
  this->blockSignals( true );
  double val = GetProperty()->GetVectorLineWidth();
  double ratio = 1;
#if VTK_MAJOR_VERSION > 7
  ratio = MainWindow::GetMainWindow()->devicePixelRatio();
#endif
  for ( int i = 0; i < 3; i++ )
  {
    m_glyphActor2D[i]->GetProperty()->SetLineWidth(val*ratio);
    m_vectorDotActor2D[i]->GetProperty()->SetPointSize(val*(val>1?2:3)*ratio);
  }
  for ( int i = 0; i < 3; i++ )
  {
    UpdateVectorActor( i );
  }
  this->blockSignals( false );
  
  emit ActorChanged();
}

void LayerMRI::UpdateVectorActor( int nPlane )
{
  UpdateVectorActor( nPlane, m_imageDataBackup.GetPointer()?m_imageDataBackup:m_imageData );
}

void LayerMRI::UpdateVectorActor( int nPlane, vtkImageData* imagedata, vtkImageData* scaledata, vtkImageData* brightnessData )
{
  double* pos = GetSlicePosition();
  double* orig = imagedata->GetOrigin();
  double* voxel_size = imagedata->GetSpacing();
  int* dim = imagedata->GetDimensions();
  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos[i] - orig[i] ) / voxel_size[i] + 0.5 );
  }
  
  //  vtkPolyData* polydata = vtkPolyDataMapper::SafeDownCast( m_vectorActor2D[nPlane]->GetMapper() )->GetInput();
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  scalars->SetNumberOfComponents( 4 );
  polydata->SetPoints( points );
  polydata->SetLines( lines );
  polydata->GetPointData()->SetScalars( scalars );
  if ( n[0] < 0 || n[0] >= dim[0] ||
       n[1] < 0 || n[1] >= dim[1] ||
       n[2] < 0 || n[2] >= dim[2] )
  {
#if VTK_MAJOR_VERSION > 5
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInputData( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInputData( polydata );
#else
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
#endif
    return;
  }
  
  int nCnt = 0;
  bool bNormalizeVector = GetProperty()->GetNormalizeVector();
  int nFrames = imagedata->GetNumberOfScalarComponents();
  if (nFrames == 6)
    bNormalizeVector = false;
  double scale_overall = GetProperty()->GetVectorScale();
  int nVectorRep = GetProperty()->GetVectorRepresentation();
  if ( nVectorRep == LayerPropertyMRI::VR_Bar )
  {
    vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION > 5
    tube->SetInputData( polydata );
#else
    tube->SetInput( polydata );
#endif
    tube->SetNumberOfSides( 4 );
    tube->SetRadius( qMin( qMin( voxel_size[0], voxel_size[1] ), voxel_size[2] ) / 8 * GetProperty()->GetVectorLineWidth());
    tube->CappingOn();
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInputConnection( tube->GetOutputPort() );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInputConnection( tube->GetOutputPort() );
  }
  else
  {
#if VTK_MAJOR_VERSION > 5
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInputData( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInputData( polydata );
#else
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
#endif
    
    if (nVectorRep == LayerPropertyMRI::VR_Direction_Line)
    {
      vtkSmartPointer<vtkMaskPoints> pts = vtkSmartPointer<vtkMaskPoints>::New();
      pts->GenerateVerticesOn();
      pts->SetOnRatio(2);
#if VTK_MAJOR_VERSION > 5
      pts->SetInputData(polydata);
#else
      pts->SetInput(polydata);
#endif
      vtkPolyDataMapper::SafeDownCast( m_vectorDotActor2D[nPlane]->GetMapper() )->SetInputConnection( pts->GetOutputPort() );
    }
  }
  
  QString orient = GetOrientationString();
  int flag_assign[3] = { 0, 1, 2 };
  int flag_sign[3] = {1, 1, 1};
  QString default_orient = "RASLPI";
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      if (orient[i] == default_orient[j] || orient[i] == default_orient[j+3])
      {
        flag_assign[i] = j;
        flag_sign[i] = (orient[i] == default_orient[j]?1:-1);
      }
    }
  }
  
  unsigned char c[4] = { 0, 0, 0, 255 };
  double scale = scale_overall;
  scale *= GetProperty()->GetVectorDisplayScale();
  scale_overall = scale;
  
  char* ptr = (char*)imagedata->GetScalarPointer();
  int scalar_type = imagedata->GetScalarType();
  char* scale_ptr = NULL;
  int scale_scalar_type = 0;
  int* scale_dim = NULL;
  if (scaledata)
  {
    scale_ptr = (char*)scaledata->GetScalarPointer();
    scale_dim = scaledata->GetDimensions();
    scale_scalar_type = scaledata->GetScalarType();
  }
  char* brightness_ptr = NULL;
  int brightness_scalar_type = 0;
  int* brightness_dim = NULL;
  int brightness_nframes = 1;
  if (brightnessData)
  {
    brightness_ptr = (char*)brightnessData->GetScalarPointer();
    brightness_dim = brightnessData->GetDimensions();
    brightness_scalar_type = brightnessData->GetScalarType();
    brightness_nframes = brightnessData->GetNumberOfScalarComponents();
  }

  double actor_pos[3] = {0,0,0};
  actor_pos[nPlane] = voxel_size[nPlane]*(nPlane==2?-dim[nPlane]:dim[nPlane])/2;
  m_glyphActor2D[nPlane]->SetPosition(actor_pos);
  m_vectorDotActor2D[nPlane]->SetPosition(actor_pos);
  int nSkip = GetProperty()->GetVectorSkip()+1;
  char* mask_ptr = NULL;
  int mask_scalar_type;
  int mask_frames;
  if (m_layerMask)
  {
    mask_ptr = (char*)m_layerMask->GetImageData()->GetScalarPointer();
    mask_scalar_type = m_layerMask->GetImageData()->GetScalarType();
    mask_frames = m_layerMask->GetNumberOfFrames();
  }
  double dNormTh = GetProperty()->GetVectorNormThreshold();
  if (nFrames == 6)
    scale *= 2;

  double brightness = 1;
  switch ( nPlane )
  {
  case 0:
    for ( int i = 0; i < dim[1]; i+=nSkip )
    {
      for ( int j = 0; j < dim[2]; j+=nSkip )
      {
        if (mask_ptr && MyVTKUtils::GetImageDataComponent(mask_ptr, dim, mask_frames, n[0], i, j, 0, mask_scalar_type) < m_dMaskThreshold)
          continue;
        double v[3], v2[3] = {0}, vn[3], v_temp[3], v_temp2[3];
        double* vp = v;
        double pt[3];
        v[0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n[0], i, j, 0, scalar_type );
        v[1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n[0], i, j, 1, scalar_type );
        v[2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n[0], i, j, 2, scalar_type );
        if (scaledata)
          scale = MyVTKUtils::GetImageDataComponent(scale_ptr, scale_dim, 1, n[0], i, j, 0, scale_scalar_type ) * scale_overall;
        if (brightnessData)
          brightness = qMin(1.0, MyVTKUtils::GetImageDataComponent(brightness_ptr, brightness_dim, brightness_nframes, n[0], i, j, 0, brightness_scalar_type));
        if (nFrames == 6)
        {
          v2[0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n[0], i, j, 3, scalar_type );
          v2[1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n[0], i, j, 4, scalar_type );
          v2[2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n[0], i, j, 5, scalar_type );

          vn[0] = v2[0] - v[0];
          vn[1] = v2[1] - v[1];
          vn[2] = v2[2] - v[2];
          vp = vn;
        }
        if (vtkMath::Norm(vp) < dNormTh)
          continue;
        
        if (v[0] != 0 || v[1] != 0 || v[2] != 0)
        {
          if (bNormalizeVector)
            vtkMath::Normalize(v);

          if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_X )
          {
            v[0] = -v[0];
            v2[0] = -v2[0];
          }
          else if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_Y )
          {
            v[1] = -v[1];
            v2[1] = -v2[1];
          }
          else if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_Z )
          {
            v[2] = -v[2];
            v2[2] = -v2[2];
          }

          for (int k = 0; k < 3; k++)
          {
            v_temp[k] = v[k];
            v_temp2[k] = v2[k];
          }
          for (int k = 0; k < 3; k++)
          {
            v[flag_assign[k]] = v_temp[k]*flag_sign[k];
            v2[flag_assign[k]] = v_temp2[k]*flag_sign[k];
          }
          
          pt[0] = orig[0] + voxel_size[0] * n[0];
          pt[1] = orig[1] + voxel_size[1] * i;
          pt[2] = orig[2] + voxel_size[2] * j;
          lines->InsertNextCell( 2 );
          points->InsertNextPoint( pt[0] + scale * v[0],
              pt[1] + scale * v[1],
              pt[2] + scale * v[2] );
          if (nFrames == 6)
          {
            points->InsertNextPoint( pt[0] + scale * v2[0],
                pt[1] + scale * v2[1],
                pt[2] + scale * v2[2] );
          }
          else
          {
            points->InsertNextPoint( pt[0] - scale * v[0],
                pt[1] - scale * v[1],
                pt[2] - scale * v[2] );
          }
          lines->InsertCellPoint( nCnt++ );
          lines->InsertCellPoint( nCnt++ );
          
          if (nVectorRep == LayerPropertyMRI::VR_Direction_Line)
          {
            GetColorWheelColor(v, nPlane, c);
          }
          else
          {
            if (!bNormalizeVector)
              vtkMath::Normalize( v );  // normalize v for color computing
            c[0] = (int)(fabs( v[0] *255*brightness ) );
            c[1] = (int)(fabs( v[1] *255*brightness ) );
            c[2] = (int)(fabs( v[2] *255*brightness ) );
            ReorderColorComponent(c);
          }
#if VTK_MAJOR_VERSION > 5
          scalars->InsertNextTypedTuple( c );
          scalars->InsertNextTypedTuple( c );
#else
          scalars->InsertNextTupleValue( c );
          scalars->InsertNextTupleValue( c );
#endif
        }
      }
    }
    break;
  case 1:
    for ( int i = 0; i < dim[0]; i+=nSkip )
    {
      for ( int j = 0; j < dim[2]; j+=nSkip )
      {
        if (mask_ptr && MyVTKUtils::GetImageDataComponent(mask_ptr, dim, mask_frames, i, n[1], j, 0, mask_scalar_type) < m_dMaskThreshold)
          continue;
        double v[3], v2[3] = {0}, vn[3], v_temp[3], v_temp2[3];
        double* vp = v;
        double pt[3];
        v[0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, n[1], j, 0, scalar_type );
        v[1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, n[1], j, 1, scalar_type );
        v[2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, n[1], j, 2, scalar_type );
        if (scaledata)
          scale = MyVTKUtils::GetImageDataComponent(scale_ptr, scale_dim, 1, i, n[1], j, 0, scale_scalar_type ) * scale_overall;
        if (brightnessData)
          brightness = qMin(1.0, MyVTKUtils::GetImageDataComponent(brightness_ptr, brightness_dim, brightness_nframes, i, n[1], j, 0, brightness_scalar_type));
        if (nFrames == 6)
        {
          v2[0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, n[1], j, 3, scalar_type );
          v2[1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, n[1], j, 4, scalar_type );
          v2[2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, n[1], j, 5, scalar_type );

          vn[0] = v2[0] - v[0];
          vn[1] = v2[1] - v[1];
          vn[2] = v2[2] - v[2];
          vp = vn;
        }
        if (vtkMath::Norm(vp) < dNormTh)
          continue;

        if (v[0] != 0 || v[1] != 0 || v[2] != 0)
        {
          if (bNormalizeVector)
            vtkMath::Normalize(v);

          if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_X )
          {
            v[0] = -v[0];
            v2[0] = -v2[0];
          }
          else if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_Y )
          {
            v[1] = -v[1];
            v2[1] = -v2[1];
          }
          else if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_Z )
          {
            v[2] = -v[2];
            v2[2] = -v2[2];
          }

          for (int k = 0; k < 3; k++)
          {
            v_temp[k] = v[k];
            v_temp2[k] = v2[k];
          }
          for (int k = 0; k < 3; k++)
          {
            v[flag_assign[k]] = v_temp[k]*flag_sign[k];
            v2[flag_assign[k]] = v_temp2[k]*flag_sign[k];
          }
          
          pt[0] = orig[0] + voxel_size[0] * i;
          pt[1] = orig[1] + voxel_size[1] * n[1];
          pt[2] = orig[2] + voxel_size[2] * j;
          lines->InsertNextCell( 2 );
          points->InsertNextPoint( pt[0] + scale * v[0],
              pt[1] + scale * v[1],
              pt[2] + scale * v[2] );
          if (nFrames == 6)
          {
            points->InsertNextPoint( pt[0] + scale * v2[0],
                pt[1] + scale * v2[1],
                pt[2] + scale * v2[2] );
          }
          else
          {
            points->InsertNextPoint( pt[0] - scale * v[0],
                pt[1] - scale * v[1],
                pt[2] - scale * v[2] );
          }
          lines->InsertCellPoint( nCnt++ );
          lines->InsertCellPoint( nCnt++ );
          
          if (nVectorRep == LayerPropertyMRI::VR_Direction_Line)
          {
            GetColorWheelColor(v, nPlane, c);
          }
          else
          {
            if (!bNormalizeVector)
              vtkMath::Normalize( v );  // normalize v for color computing
            c[0] = (int)(fabs( v[0] *255*brightness ) );
            c[1] = (int)(fabs( v[1] *255*brightness ) );
            c[2] = (int)(fabs( v[2] *255*brightness ) );
            ReorderColorComponent(c);
          }
#if VTK_MAJOR_VERSION > 5
          scalars->InsertNextTypedTuple( c );
          scalars->InsertNextTypedTuple( c );
#else
          scalars->InsertNextTupleValue( c );
          scalars->InsertNextTupleValue( c );
#endif
        }
      }
    }
    break;
  case 2:
    for ( int i = 0; i < dim[0]; i+=nSkip )
    {
      for ( int j = 0; j < dim[1]; j+=nSkip )
      {
        if (mask_ptr && MyVTKUtils::GetImageDataComponent(mask_ptr, dim, mask_frames, i, j, n[2], 0, mask_scalar_type) < m_dMaskThreshold)
          continue;
        double v[3], v2[3] = {0}, vn[3], v_temp[3], v_temp2[3];
        double* vp = v;
        double pt[3];
        v[0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, n[2], 0, scalar_type );
        v[1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, n[2], 1, scalar_type );
        v[2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, n[2], 2, scalar_type );
        if (scaledata)
          scale = MyVTKUtils::GetImageDataComponent(scale_ptr, scale_dim, 1, i, j, n[2], 0, scale_scalar_type) * scale_overall;
        if (brightnessData)
          brightness = qMin(1.0, MyVTKUtils::GetImageDataComponent(brightness_ptr, brightness_dim, brightness_nframes, i, j, n[2], 0, brightness_scalar_type));
        if (nFrames == 6)
        {
          v2[0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, n[2], 3, scalar_type );
          v2[1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, n[2], 4, scalar_type );
          v2[2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, n[2], 5, scalar_type );

          vn[0] = v2[0] - v[0];
          vn[1] = v2[1] - v[1];
          vn[2] = v2[2] - v[2];
          vp = vn;
        }
        if (vtkMath::Norm(vp) < dNormTh)
          continue;

        if (v[0] != 0 || v[1] != 0 || v[2] != 0)
        {
          if (bNormalizeVector)
            vtkMath::Normalize(v);

          if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_X )
          {
            v[0] = -v[0];
            v2[0] = -v2[0];
          }
          else if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_Y )
          {
            v[1] = -v[1];
            v2[1] = -v2[1];
          }
          else if ( GetProperty()->GetVectorInversion() == LayerPropertyMRI::VI_Z )
          {
            v[2] = -v[2];
            v2[2] = -v2[2];
          }

          for (int k = 0; k < 3; k++)
          {
            v_temp[k] = v[k];
            v_temp2[k] = v2[k];
          }
          for (int k = 0; k < 3; k++)
          {
            v[flag_assign[k]] = v_temp[k]*flag_sign[k];
            v2[flag_assign[k]] = v_temp2[k]*flag_sign[k];
          }
          
          pt[0] = orig[0] + voxel_size[0] * i;
          pt[1] = orig[1] + voxel_size[1] * j;
          pt[2] = orig[2] + voxel_size[2] * n[2];
          lines->InsertNextCell( 2 );
          points->InsertNextPoint( pt[0] + scale * v[0],
              pt[1] + scale * v[1],
              pt[2] + scale * v[2] );
          if (nFrames == 6)
          {
            points->InsertNextPoint( pt[0] + scale * v2[0],
                pt[1] + scale * v2[1],
                pt[2] + scale * v2[2] );
          }
          else
          {
            points->InsertNextPoint( pt[0] - scale * v[0],
                pt[1] - scale * v[1],
                pt[2] - scale * v[2] );
          }
          lines->InsertCellPoint( nCnt++ );
          lines->InsertCellPoint( nCnt++ );
          
          if (nVectorRep == LayerPropertyMRI::VR_Direction_Line)
          {
            GetColorWheelColor(v, nPlane, c);
          }
          else
          {
            if (!bNormalizeVector)
              vtkMath::Normalize( v );  // normalize v for color computing
            c[0] = (int)(fabs( v[0] *255*brightness ) );
            c[1] = (int)(fabs( v[1] *255*brightness ) );
            c[2] = (int)(fabs( v[2] *255*brightness ) );
            ReorderColorComponent(c);
          }
#if VTK_MAJOR_VERSION > 5
          scalars->InsertNextTypedTuple( c );
          scalars->InsertNextTypedTuple( c );
#else
          scalars->InsertNextTupleValue( c );
          scalars->InsertNextTupleValue( c );
#endif
        }
      }
    }
    break;
  default:
    break;
  }
  
  emit ActorUpdated();
}

void LayerMRI::GetColorWheelColor(double *v, int nPlane, unsigned char *c_out)
{
  double x = v[0], y = v[1];
  if (nPlane == 0)
    x = v[2];
  else if (nPlane == 1)
    y = v[2];
  double v1[2] = {x, y}, v2[2] = {1, 0};
  vtkMath::Normalize2D(v1);
  double angle = acos(vtkMath::Dot2D(v1, v2))/(2*vtkMath::Pi());
  if (y < 0)
    angle = 1 - angle;
  double hsv[3] = {angle, 1, 1}, rgb[3];
  vtkMath::HSVToRGB(hsv, rgb);
  c_out[0] = (int)(rgb[0]*255);
  c_out[1] = (int)(rgb[1]*255);
  c_out[2] = (int)(rgb[2]*255);
}

void LayerMRI::UpdateTensorActor()
{
  this->blockSignals( true );
  for ( int i = 0; i < 3; i++ )
  {
    UpdateTensorActor( i );
  }
  this->blockSignals( false );
  
  emit ActorChanged();
}

void LayerMRI::UpdateTensorActor( int nPlane, vtkImageData* imagedata_in )
{
  vtkImageData* imagedata = imagedata_in;
  if ( !imagedata )
  {
    imagedata = m_imageData;
  }
  
  double* pos = GetSlicePosition();
  double* orig = imagedata->GetOrigin();
  double* voxel_size = imagedata->GetSpacing();
  int* dim = imagedata->GetDimensions();
  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos[i] - orig[i] ) / voxel_size[i] + 0.5 );
  }
  
  if ( n[0] < 0 || n[0] >= dim[0] ||
       n[1] < 0 || n[1] >= dim[1] ||
       n[2] < 0 || n[2] >= dim[2] )
  {
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
#if VTK_MAJOR_VERSION > 5
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInputData( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInputData( polydata );
#else
    vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
    vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
#endif
    return;
  }
  
  double scale = qMin( qMin( voxel_size[0], voxel_size[1] ), voxel_size[2] ) * 1.0;
  
  vtkSmartPointer<vtkPolyDataAlgorithm> objsource;
  if ( GetProperty()->GetTensorRepresentation() == LayerPropertyMRI::TR_Ellipsoid )
  {
    objsource = vtkSmartPointer<vtkSphereSource>::New();
  }
  else
  {
    objsource = vtkSmartPointer<vtkCubeSource>::New();
  }
  objsource->Update();
  vtkPolyData* srcpolydata = objsource->GetOutput();
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  scalars->SetNumberOfComponents( 4 );
  //  srcpolydata->GetPointData()->SetNormals( NULL );    // remove normals
  vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
  double pt[3];
  int nSkip = 1;
  switch ( nPlane )
  {
  case 0:
    for ( int i = 0; i < dim[1]; i+=nSkip )
    {
      for ( int j = 0; j < dim[2]; j+=nSkip )
      {
        pt[0] = orig[0] + voxel_size[0] * n[0];
        pt[1] = orig[1] + voxel_size[1] * i;
        pt[2] = orig[2] + voxel_size[2] * j;
        BuildTensorGlyph( imagedata, n[0], i, j, pt, scale, srcpolydata, scalars, append ) ;
      }
    }
    break;
  case 1:
    for ( int i = 0; i < dim[0]; i+=nSkip )
    {
      for ( int j = 0; j < dim[2]; j+=nSkip )
      {
        pt[0] = orig[0] + voxel_size[0] * i;
        pt[1] = orig[1] + voxel_size[1] * n[1];
        pt[2] = orig[2] + voxel_size[2] * j;
        BuildTensorGlyph( imagedata, i, n[1], j, pt, scale, srcpolydata, scalars, append );
      }
    }
    break;
  case 2:
    for ( int i = 0; i < dim[0]; i+=nSkip )
    {
      for ( int j = 0; j < dim[1]; j+=nSkip )
      {
        pt[0] = orig[0] + voxel_size[0] * i;
        pt[1] = orig[1] + voxel_size[1] * j;
        pt[2] = orig[2] + voxel_size[2] * n[2];
        BuildTensorGlyph( imagedata, i, j, n[2], pt, scale, srcpolydata, scalars, append );
      }
    }
    break;
  default:
    break;
  }
  append->Update();
  vtkPolyData* polydata = append->GetOutput();
  polydata->GetPointData()->SetScalars( scalars );
#if VTK_MAJOR_VERSION > 5
  vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInputData( polydata );
  vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInputData( polydata );
#else
  vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() )->SetInput( polydata );
  vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() )->SetInput( polydata );
#endif
  emit ActorUpdated();
}

void LayerMRI::BuildTensorGlyph( vtkImageData* imagedata,
                                 int i, int j, int k,
                                 double* pt, double scale,
                                 vtkPolyData* sourcepolydata,
                                 vtkUnsignedCharArray* scalars,
                                 vtkPolyDataAlgorithm* a)
{
  double** D = private_buf1_3x3;
  double** v = private_buf2_3x3;
  double w[3];
  char* ptr = (char*)imagedata->GetScalarPointer();
  int* dim = imagedata->GetDimensions();
  int scalar_type = imagedata->GetScalarType();
  int nFrames = imagedata->GetNumberOfScalarComponents();
  D[0][0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 0, scalar_type );
  D[0][1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 1, scalar_type );
  D[0][2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 2, scalar_type );
  D[1][0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 3, scalar_type );
  D[1][1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 4, scalar_type );
  D[1][2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 5, scalar_type );
  D[2][0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 6, scalar_type );
  D[2][1] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 7, scalar_type );
  D[2][2] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, 8, scalar_type );
  if ( vtkMath::Jacobi( D, w, v ) )
  {
    v[1][0] = -v[1][0];         // by default invert Y !!
    v[1][1] = -v[1][1];
    v[1][2] = -v[1][2];
    if ( GetProperty()->GetTensorInversion() == LayerPropertyMRI::VI_X )
    {
      v[0][0] = -v[0][0];
      v[0][1] = -v[0][1];
      v[0][2] = -v[0][2];
    }
    else if ( GetProperty()->GetTensorInversion() == LayerPropertyMRI::VI_Y )
    {
      v[1][0] = -v[1][0];
      v[1][1] = -v[1][1];
      v[1][2] = -v[1][2];
    }
    else if ( GetProperty()->GetTensorInversion() == LayerPropertyMRI::VI_Z )
    {
      v[2][0] = -v[2][0];
      v[2][1] = -v[2][1];
      v[2][2] = -v[2][2];
    }
    
    // make the vectors in right hand coordinate
    double v0[3] = { v[0][0], v[1][0], v[2][0] };
    double v1[3] = { v[0][1], v[1][1], v[2][1] };
    double v2[3];
    vtkMath::Cross( v0, v1, v2 );
    v[0][2] = v2[0];
    v[1][2] = v2[1];
    v[2][2] = v2[2];
    
    vtkMath::Normalize( w );
    //   double w_sum = 1;//fabs(w[0]) + fabs(w[1]) + fabs(w[2]);
    w[0] = fabs(w[0]*scale);
    w[1] = fabs(w[1]*scale);
    w[2] = fabs(w[2]*scale);
    
    vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
    tr->Identity();
    tr->Translate( pt );
    double m[16];
    memset( m, 0, sizeof(double)*16 );
    m[15] = 1;
    
    m[0] = v[0][0];
    m[1] = v[0][1];
    m[2] = v[0][2];
    
    m[4] = v[1][0];
    m[5] = v[1][1];
    m[6] = v[1][2];
    
    
    m[8] = v[2][0];
    m[9] = v[2][1];
    m[10]= v[2][2];
    
    tr->Concatenate(m);
    tr->Scale( w );
    //    tr->RotateZ( -90 );
    unsigned char c[4];
    c[0] = (int)(fabs( v[0][0] *255 ) );
    c[1] = (int)(fabs( v[1][0] *255 ) );
    c[2] = (int)(fabs( v[2][0] *255 ) );
    c[3] = 255;
    int nPts = sourcepolydata->GetPoints()->GetNumberOfPoints();
    for ( int i = 0; i < nPts; i++ )
    {
#if VTK_MAJOR_VERSION > 5
      scalars->InsertNextTypedTuple( c );
#else
      scalars->InsertNextTupleValue( c );
#endif
    }
    
    vtkSmartPointer<vtkTransformPolyDataFilter> filter =
        vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    filter->SetTransform( tr );
#if VTK_MAJOR_VERSION > 5
    filter->SetInputData( sourcepolydata );
    a->AddInputData( filter->GetOutput() );
#else
    filter->SetInput( sourcepolydata );
    a->AddInput( filter->GetOutput() );
#endif
  }
}

void LayerMRI::GetRASCenter( double* rasPt )
{
  MRI* mri = m_volumeSource->GetMRITarget();
  ::MRIvoxelToWorld( mri,
                     mri->width / 2.0 - 0.5,
                     mri->height / 2.0 - 0.5,
                     mri->depth / 2.0 - 0.5,
                     &rasPt[0], &rasPt[1], &rasPt[2] );
}


void LayerMRI::UpdateVoxelValueRange( double dValue )
{
  GetProperty()->UpdateValueRange( dValue );
}

// Get voxel value range of a selected rectangle region defined by pt0, pt1
bool LayerMRI::GetVoxelValueRange( const double* pt0, const double* pt1, int nPlane, double* range_out )
{
  double* orig = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int* dim = m_imageData->GetDimensions();
  int scalar_type = m_imageData->GetScalarType();
  
  if ( nPlane < 0 ) // || nPlane >= dim[nPlane] )
  {
    return false;
  }
  
  // find the index range of the selection
  int n0[3] = { 10000000, 10000000, 10000000 },
      n1[3] = { -10000000, -10000000, -10000000 };
  n0[nPlane] = n1[nPlane] = (int)( ( pt0[nPlane] - orig[nPlane] ) / voxel_size[nPlane] + 0.5 );
  for ( int i = 0; i < 3; i++ )
  {
    if ( i != nPlane )
    {
      int p0 = (int)( ( pt0[i] - orig[i] ) / voxel_size[i] + 0.5 );
      int p1 = (int)( ( pt1[i] - orig[i] ) / voxel_size[i] + 0.5 );
      p0 = qMax( 0, qMin( dim[i]-1, p0 ) );
      p1 = qMax( 0, qMin( dim[i]-1, p1 ) );
      n0[i] = qMin( p0, qMin( p1, n0[i] ) );
      n1[i] = qMax( p0, qMax( p1, n0[i] ) );
      
    }
  }
  
  int nActiveComp = GetActiveFrame();
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int nFrames = m_imageData->GetNumberOfScalarComponents();
  range_out[0] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n0[0], n0[1], n0[2], nActiveComp, scalar_type);
  range_out[1] = range_out[0];
  
  for ( int i = n0[0]; i <= n1[0]; i++ )
  {
    for ( int j = n0[1]; j <= n1[1]; j++ )
    {
      for ( int k = n0[2]; k <= n1[2]; k++ )
      {
        double value = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, nActiveComp, scalar_type);
        if ( range_out[0] > value )
        {
          range_out[0] = value;
        }
        if ( range_out[1] < value )
        {
          range_out[1] = value;
        }
      }
    }
  }
  
  return true;
}

// Get rectangle region stats
bool LayerMRI::GetVoxelStatsRectangle( const double* pt0, const double* pt1, int nPlane, double* mean_out, double* sd_out, int* cnt_out )
{
  double* orig = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int* dim = m_imageData->GetDimensions();
  int scalar_type = m_imageData->GetScalarType();
  
  if ( nPlane < 0 ) //|| nPlane >= dim[nPlane] )
  {
    return false;
  }
  
  // find the index range of the selection
  int n0[3] = { 10000000, 10000000, 10000000 },
      n1[3] = { -10000000, -10000000, -10000000 };
  n0[nPlane] = n1[nPlane] = (int)( ( pt0[nPlane] - orig[nPlane] ) / voxel_size[nPlane] + 0.5 );
  for ( int i = 0; i < 3; i++ )
  {
    if ( i != nPlane )
    {
      int p0 = (int)( ( pt0[i] - orig[i] ) / voxel_size[i] + 0.5 );
      int p1 = (int)( ( pt1[i] - orig[i] ) / voxel_size[i] + 0.5 );
      p0 = qMax( 0, qMin( dim[i]-1, p0 ) );
      p1 = qMax( 0, qMin( dim[i]-1, p1 ) );
      n0[i] = qMin( p0, qMin( p1, n0[i] ) );
      n1[i] = qMax( p0, qMax( p1, n1[i] ) );
    }
  }
  
  int nActiveComp = GetActiveFrame();
  double dMean = 0;
  int nCount = 0;
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int nFrames = m_imageData->GetNumberOfScalarComponents();
  for ( int i = n0[0]; i <= n1[0]; i++ )
  {
    for ( int j = n0[1]; j <= n1[1]; j++ )
    {
      for ( int k = n0[2]; k <= n1[2]; k++ )
      {
        dMean += MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, nActiveComp, scalar_type );
        nCount++;
      }
    }
  }
  if ( nCount > 0 )
  {
    *mean_out = dMean / nCount;
  }
  
  if ( sd_out )
  {
    double sd = 0;
    for ( int i = n0[0]; i <= n1[0]; i++ )
    {
      for ( int j = n0[1]; j <= n1[1]; j++ )
      {
        for ( int k = n0[2]; k <= n1[2]; k++ )
        {
          double value = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, nActiveComp, scalar_type );
          sd += ( value-(*mean_out) ) * ( value-(*mean_out) );
        }
      }
    }
    if (nCount > 1)
    {
      *sd_out = sqrt( sd / (nCount-1) );
    }
    else
    {
      *sd_out = 0;
    }
  }
  
  if ( cnt_out )
  {
    *cnt_out = nCount;
  }
  
  return true;
}

bool LayerMRI::GetVoxelStats(QVector<int> &indices, double *mean_out, double *sd_out)
{
  int nActiveComp = GetActiveFrame();
  double dMean = 0;
  int nCount = 0;
  int* dim = m_imageData->GetDimensions();
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int scalar_type = m_imageData->GetScalarType();
  int nFrames = m_imageData->GetNumberOfScalarComponents();
  for ( int n = 0; n < indices.size(); n+=3 )
  {
    int i = indices[n];
    int j = indices[n+1];
    int k = indices[n+2];
    if (i >= 0 && i < dim[0] && j >= 0 && j < dim[1] && k >= 0 && k < dim[2])
    {
      dMean += MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, nActiveComp, scalar_type );
      nCount++;
    }
  }
  if ( nCount > 0 )
  {
    *mean_out = dMean / nCount;
  }
  
  if ( sd_out )
  {
    double sd = 0;
    for ( int n = 0; n < indices.size(); n+=3 )
    {
      int i = indices[n];
      int j = indices[n+1];
      int k = indices[n+2];
      if (i >= 0 && i < dim[0] && j >= 0 && j < dim[1] && k >= 0 && k < dim[2])
      {
        double value = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, nActiveComp, scalar_type );
        sd += ( value-(*mean_out) ) * ( value-(*mean_out) );
      }
    }
    if (nCount > 1)
    {
      *sd_out = sqrt( sd / (nCount-1) );
    }
    else
    {
      *sd_out = 0;
    }
  }
  return true;
}

bool LayerMRI::GetVoxelStatsByTargetRAS(QVector<float> &coords, double* mean_out, double *sd_out)
{
  double* orig = m_imageData->GetOrigin();
  double* vsize = m_imageData->GetSpacing();
  
  QVector<int> indices;
  for (int i = 0; i < coords.size(); i+=3)
  {
    indices << (int)( ( coords[i] - orig[0] ) / vsize[0] + 0.5 )
        << (int)( ( coords[i+1] - orig[1] ) / vsize[1] + 0.5 )
        << (int)( ( coords[i+2] - orig[2] ) / vsize[2] + 0.5 );
  }
  return GetVoxelStats(indices, mean_out, sd_out);
}

// memory allocated for indice_out and value_out need to be freed outside of this function!
bool LayerMRI::GetVoxelsOnLine( const double* pt0, const double* pt1, int nPlane, int*& indice_out, double*& value_out, int* cnt_out )
{
  double* orig = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  int* dim = m_imageData->GetDimensions();
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int scalar_type = m_imageData->GetScalarType();
  int nFrames = m_imageData->GetNumberOfScalarComponents();
  
  if ( nPlane < 0 )// || nPlane >= dim[nPlane] )
  {
    return false;
  }
  
  // find the index range of the selection
  int n0[3], n1[3];
  n0[nPlane] = n1[nPlane] = (int)( ( pt0[nPlane] - orig[nPlane] ) / voxel_size[nPlane] + 0.5 );
  for ( int i = 0; i < 3; i++ )
  {
    if ( i != nPlane )
    {
      int p0 = (int)( ( pt0[i] - orig[i] ) / voxel_size[i] + 0.5 );
      int p1 = (int)( ( pt1[i] - orig[i] ) / voxel_size[i] + 0.5 );
      n0[i] = qMax( 0, qMin( dim[i]-1, p0 ) );
      n1[i] = qMax( 0, qMin( dim[i]-1, p1 ) );
    }
  }
  
  std::vector<int> indices = GetVoxelIndicesBetweenPoints( n0, n1 );
  std::vector<double> values;
  
  int nActiveComp = GetActiveFrame();
  double dMean = 0;
  int nCount = 0;
  for ( size_t i = 0; i < indices.size(); i += 3 )
  {
    double value = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, indices[i], indices[i+1], indices[i+2], nActiveComp, scalar_type );
    dMean += value;
    values.push_back( value );
    nCount++;
  }
  
  indice_out = new int[nCount*3];
  value_out = new double[nCount];
  double pt[3], ras[3];
  int nIndex[3];
  for ( int i = 0; i < nCount; i++ )
  {
    pt[0] = indices[i*3]*voxel_size[0] + orig[0];
    pt[1] = indices[i*3+1]*voxel_size[1] + orig[1];
    pt[2] = indices[i*3+2]*voxel_size[2] + orig[2];
    RemapPositionToRealRAS( pt, ras );
    RASToOriginalIndex( ras, nIndex );
    indice_out[i*3]   = nIndex[0];
    indice_out[i*3+1] = nIndex[1];
    indice_out[i*3+2] = nIndex[2];
    value_out[i] = values[i];
  }
  
  *cnt_out = nCount;
  return true;
}

std::vector<int> LayerMRI::GetVoxelIndicesBetweenPoints( int* n0, int* n1 )
{
  std::vector<int> indices;
  if ( n1[0] == n0[0] && n1[1] == n0[1] && n1[2] == n0[2] )
  {
    indices.push_back( n0[0] );
    indices.push_back( n0[1] );
    indices.push_back( n0[2] );
  }
  else if ( fabs(n0[0]-n1[0]) <= 1 && fabs(n0[1]-n1[1]) <= 1 && fabs(n0[2]-n1[2]) <= 1 )
  {
    indices.push_back( n0[0] );
    indices.push_back( n0[1] );
    indices.push_back( n0[2] );
    indices.push_back( n1[0] );
    indices.push_back( n1[1] );
    indices.push_back( n1[2] );
  }
  else
  {
    int n[3];
    for ( int i = 0; i < 3; i++ )
    {
      n[i] = (int)( (n0[i]+n1[i]) / 2.0 + 0.5 );
    }
    
    indices = GetVoxelIndicesBetweenPoints( n0, n );
    std::vector<int> indices1 = GetVoxelIndicesBetweenPoints( n, n1 );
    for ( size_t i = 3; i < indices1.size(); i++ )
    {
      indices.push_back( indices1[i] );
    }
  }
  
  return indices;
}

void LayerMRI::ResetWindowLevel()
{
  double range[2];
  m_imageData->GetScalarRange( range );
  GetProperty()->SetMinMaxGrayscaleWindow( range[0], range[1] );
  GetProperty()->SetMinMaxGenericThreshold( range[0], range[1] );
  GetProperty()->SetHeatScale( range[0], (range[0]+range[1])/2, range[1] );
}

int LayerMRI::GetDataType()
{
  return ( m_volumeSource ? m_volumeSource->GetDataType() : -1 );
}

COLOR_TABLE* LayerMRI::GetEmbeddedColorTable()
{
  return ( m_volumeSource ? m_volumeSource->GetEmbeddedColorTable(): NULL );
}

void LayerMRI::SnapToVoxelCenter( const double* pt_in, double* pt_out )
{
  if ( m_imageData == NULL )
  {
    pt_out[0] = pt_in[0];
    pt_out[1] = pt_in[1];
    pt_out[2] = pt_in[2];
    return;
  }
  
  double* orig = m_imageData->GetOrigin();
  double* vsize = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    pt_out[i] = ( (int)( (pt_in[i] - orig[i])/vsize[i] + 0.5 ) ) * vsize[i] + orig[i];
  }
}

void LayerMRI::UpdateLabelOutline()
{
  if ( /*GetProperty()->GetColorMap() == LayerPropertyMRI::LUT &&*/ GetProperty()->GetShowLabelOutline() )
  {
    double* vsize = m_imageData->GetSpacing();
    for ( int i = 0; i < 3; i++ )
    {
      //      mResample[i]->SetAxisMagnificationFactor( 0, IMAGE_RESAMPLE_FACTOR );
      //      mResample[i]->SetAxisMagnificationFactor( 1, IMAGE_RESAMPLE_FACTOR );
      //      mResample[i]->SetAxisMagnificationFactor( 2, IMAGE_RESAMPLE_FACTOR );
      double pos[3] = { vsize[0]/IMAGE_RESAMPLE_FACTOR/2, vsize[1]/IMAGE_RESAMPLE_FACTOR/2, vsize[2]/IMAGE_RESAMPLE_FACTOR/2 };
      mResample[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
      mResample[i]->SetOutputSpacing(vsize[0]/IMAGE_RESAMPLE_FACTOR, vsize[1]/IMAGE_RESAMPLE_FACTOR, vsize[2]/IMAGE_RESAMPLE_FACTOR);
      mResample[i]->SetInterpolationModeToNearestNeighbor();
      mEdgeFilter[i]->SetInputConnection( mResample[i]->GetOutputPort() );
      mColorMap[i]->SetInputConnection( mEdgeFilter[i]->GetOutputPort() );
      pos[i] = m_dSlicePosition[i];
      m_sliceActor2D[i]->SetPosition( pos );
      m_sliceActor3D[i]->SetPosition( pos );
    }
  }
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
      double pos[3] = { 0, 0, 0};
      mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
      pos[i] = m_dSlicePosition[i];
      m_sliceActor2D[i]->SetPosition( pos );
      m_sliceActor3D[i]->SetPosition( pos );
    }
  }
  emit ResampleFactorChanged();
  emit ActorUpdated();
}

void LayerMRI::UpdateUpSampleMethod()
{
  switch ( GetProperty()->GetUpSampleMethod() )
  {
  case LayerPropertyMRI::UM_NearestNeighbor:
    for ( int i = 0; i < 3; i++ )
    {
      mResample[i]->SetInterpolationModeToNearestNeighbor();
    }
    break;
  case LayerPropertyMRI::UM_Linear:
    for ( int i = 0; i < 3; i++ )
    {
      mResample[i]->SetInterpolationModeToLinear();
    }
    break;
  case LayerPropertyMRI::UM_Cubic:
    for ( int i = 0; i < 3; i++ )
    {
      mResample[i]->SetInterpolationModeToCubic();
    }
    break;
  default:
    break;
  }
  if ( GetProperty()->GetUpSampleMethod() == LayerPropertyMRI::UM_None )
  {
    for ( int i = 0; i < 3; i++ )
    {
      mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
    }
  }
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
      mResample[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
      mColorMap[i]->SetInputConnection( mResample[i]->GetOutputPort() );
      if ( !GetProperty()->GetShowLabelOutline() )
      {
        //        mResample[i]->SetAxisMagnificationFactor( 0, IMAGE_RESAMPLE_FACTOR/2 );
        //        mResample[i]->SetAxisMagnificationFactor( 1, IMAGE_RESAMPLE_FACTOR/2 );
        //        mResample[i]->SetAxisMagnificationFactor( 2, IMAGE_RESAMPLE_FACTOR/2 );
      }
    }
  }
  
  emit ResampleFactorChanged();
  emit ActorUpdated();
}


void LayerMRI::GetCurrentLabelStats(int nPlane, float *label_out, int *count_out, float *area_out,
                                    LayerMRI *underlying_mri, double *mean_out, double *sd_out)
{
  if ( !m_imageData || nPlane < 0 || nPlane > 2 )
  {
    return;
  }
  
  double* origin = m_imageData->GetOrigin();
  int* dim = m_imageData->GetDimensions();
  double* pos = GetSlicePosition();
  double vs[3];
  m_imageData->GetSpacing( vs );
  
  int n[3];
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = (int)( ( pos[i] - origin[i] ) / vs[i]+0.5 );
  }
  
  float fLabel = 0;
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int scalar_type = m_imageData->GetScalarType();
  int nFrames = m_imageData->GetNumberOfScalarComponents();
  if ( n[0] >= 0 && n[0] < dim[0] && n[1] >= 0 && n[1] < dim[1] && n[2] >= 0 && n[2] < dim[2] )
  {
    fLabel = (float)MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, n[0], n[1], n[2], m_nActiveFrame, scalar_type );
  }
  
  int cnt = 0;
  int ext[3][2] = { { 0, dim[0]-1 }, {0, dim[1]-1}, {0, dim[2]-1} };
  ext[nPlane][0] = ext[nPlane][1] = n[nPlane];
  //  QList<int> indices;
  QVector<float> coords;
  for ( int i = ext[0][0]; i <= ext[0][1]; i++ )
  {
    for ( int j = ext[1][0]; j <= ext[1][1]; j++ )
    {
      for ( int k = ext[2][0]; k <= ext[2][1]; k++ )
      {
        if ( MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, i, j, k, m_nActiveFrame, scalar_type ) == fLabel )
        {
          cnt++;
          //        indices << i << j << k;
          coords << i*vs[0] + origin[0] << j*vs[1] + origin[1] << k*vs[2] + origin[2];
        }
      }
    }
  }
  vs[nPlane] = 1.0;
  
  *label_out = fLabel;
  *count_out = cnt;
  *area_out = cnt*vs[0]*vs[1]*vs[2];
  
  if (underlying_mri)
    underlying_mri->GetVoxelStatsByTargetRAS(coords, mean_out, sd_out);
}

vtkImageData* LayerMRI::GetSliceImageData( int nPlane )
{
  return mReslice[nPlane]->GetOutput();
}

bool LayerMRI::FloodFillByContour2D( double* ras, Contour2D* c2d )
{
  int nPlane = c2d->GetPlane();
  vtkImageData* image = c2d->GetThresholdedImage();
  if (!image)
    return false;

  vtkImageData* original_image = GetSliceImageData( nPlane );
  int* nDim = image->GetDimensions();      // 2D image
  int n[3];
  double* origin = m_imageData->GetOrigin();
  double* voxel_size = m_imageData->GetSpacing();
  for ( int i = 0; i < 3; i++ )
  {
    n[i] = ( int )( ( ras[i] - origin[i] ) / voxel_size[i] + 0.5 );
  }
  int nx = nDim[0], ny = nDim[1], x = 0, y = 0;
  switch ( nPlane )
  {
  case 0:
    x = n[1];
    y = n[2];
    break;
  case 1:
    x = n[0];
    y = n[2];
    break;
  case 2:
    x = n[0];
    y = n[1];
    break;
  }
  
  char* mask = new char[nDim[0]*nDim[1]];
  memset( mask, 0, nDim[0]*nDim[1] );
  char* ptr = (char*)image->GetScalarPointer();
  int scalar_type = image->GetScalarType();
  int n_frames = image->GetNumberOfScalarComponents();
  int* orig_dim = original_image->GetDimensions();
  char* orig_ptr = (char*)original_image->GetScalarPointer();
  int orig_scalar_type = original_image->GetScalarType();
  int orig_n_frames = original_image->GetNumberOfScalarComponents();
  double dVoxelValue = MyVTKUtils::GetImageDataComponent(orig_ptr, orig_dim, orig_n_frames, x, y, 0, 0, orig_scalar_type );
  double dMaskValue = MyVTKUtils::GetImageDataComponent(ptr, nDim, n_frames, x, y, 0, 0, scalar_type );
  for ( int i = 0; i < nDim[0]; i++ )
  {
    for ( int j = 0; j < nDim[1]; j++ )
    {
      if ( MyVTKUtils::GetImageDataComponent(orig_ptr, orig_dim, orig_n_frames, i, j, 0, 0, orig_scalar_type ) == dVoxelValue &&
           MyVTKUtils::GetImageDataComponent(ptr, nDim, n_frames, i, j, 0, 0, scalar_type ) == dMaskValue )
      {
        mask[j*nDim[0]+i] = 1;
      }
    }
  }
  
  int nFillValue = 2;
  MyUtils::FloodFill( mask, x, y, nx, ny, nFillValue, 0 );
  
  nDim = m_imageData->GetDimensions();
  ptr = (char*)m_imageData->GetScalarPointer();
  scalar_type = m_imageData->GetScalarType();
  int nActiveComp = this->GetActiveFrame();
  int cnt = 0;
  switch ( nPlane )
  {
  case 0:
    for ( int i = 0; i < nx; i++ )
    {
      for ( int j = 0; j < ny; j++ )
      {
        if ( mask[j*nx+i] == nFillValue )
        {
          MyVTKUtils::SetImageDataComponent(ptr, nDim, n_frames, n[nPlane], i, j, nActiveComp, scalar_type, m_fFillValue );
          cnt++;
        }
      }
    }
    break;
  case 1:
    for ( int i = 0; i < nx; i++ )
    {
      for ( int j = 0; j < ny; j++ )
      {
        if ( mask[j*nx+i] == nFillValue )
        {
          MyVTKUtils::SetImageDataComponent(ptr, nDim, n_frames, i, n[nPlane], j, nActiveComp, scalar_type, m_fFillValue );
          cnt++;
        }
      }
    }
    break;
  case 2:
    for ( int i = 0; i < nx; i++ )
    {
      for ( int j = 0; j < ny; j++ )
      {
        if ( mask[j*nx+i] == nFillValue )
        {
          MyVTKUtils::SetImageDataComponent(ptr, nDim, n_frames, i, j, n[nPlane], nActiveComp, scalar_type, m_fFillValue);
          cnt++;
        }
      }
    }
    break;
  }
  SetModified();
  emit ActorUpdated();
  return true;
}

SurfaceRegion* LayerMRI::CreateNewSurfaceRegion( double* pt, vtkProp* prop )
{
  m_actorCurrentContour = vtkActor::SafeDownCast(prop);
  if (!m_actorCurrentContour)
    return NULL;

  SurfaceRegion* r = new SurfaceRegion( this );
  connect( r, SIGNAL(ColorChanged(QColor)), this, SIGNAL(ActorUpdated()), Qt::UniqueConnection);
  r->SetInput( vtkPolyData::SafeDownCast( m_actorCurrentContour->GetMapper()->GetInput() ) );
  r->AddPoint( pt );
  r->SetId( m_surfaceRegions.size() + 1 );
  m_surfaceRegions.push_back( r );
  int nGroup = 1;
  if ( m_currentSurfaceRegion )
  {
    m_currentSurfaceRegion->Highlight( false );
    nGroup = m_currentSurfaceRegion->GetGroup();
  }
  m_currentSurfaceRegion = r;
  r->SetGroup( nGroup );
  emit SurfaceRegionAdded();
  return r;
}

void LayerMRI::AddSurfaceRegionLoopPoint( double* pt )
{
  if ( m_currentSurfaceRegion )
  {
    m_currentSurfaceRegion->AddPoint( pt );
    emit ActorUpdated();
  }
}

void LayerMRI::CloseSurfaceRegion()
{
  if ( m_currentSurfaceRegion )
  {
    if ( !m_currentSurfaceRegion->Close() )
    {
      qDebug() << "failed to close region";
      DeleteCurrentSurfaceRegion();
    }
    emit ActorUpdated();
  }
}

SurfaceRegion* LayerMRI::SelectSurfaceRegion( double* pos )
{
  for ( int i = 0; i < m_surfaceRegions.size(); i++ )
  {
    if ( m_surfaceRegions[i]->HasPoint( pos ) )
    {
      if ( m_currentSurfaceRegion )
      {
        m_currentSurfaceRegion->Highlight( false );
      }
      m_currentSurfaceRegion = m_surfaceRegions[i];
      m_currentSurfaceRegion->Highlight( true );
      return m_currentSurfaceRegion;
    }
  }
  
  return NULL;
}

SurfaceRegion* LayerMRI::SelectSurfaceRegion( int nId )
{
  for ( int i = 0; i < m_surfaceRegions.size(); i++ )
  {
    if ( m_surfaceRegions[i]->GetId() == nId )
    {
      if ( m_currentSurfaceRegion )
      {
        m_currentSurfaceRegion->Highlight( false );
      }
      m_currentSurfaceRegion = m_surfaceRegions[i];
      m_currentSurfaceRegion->Highlight( true );
      return m_currentSurfaceRegion;
    }
  }
  return NULL;
}

bool LayerMRI::DeleteCurrentSurfaceRegion()
{
  if ( m_currentSurfaceRegion )
  {
    for ( int i = 0; i < m_surfaceRegions.size(); i++ )
    {
      if ( m_surfaceRegions[i] == m_currentSurfaceRegion )
      {
        m_surfaceRegions.erase( m_surfaceRegions.begin() + i );
        delete m_currentSurfaceRegion;
        m_currentSurfaceRegion = NULL;
        ResetSurfaceRegionIds();
        emit SurfaceRegionRemoved();
        return true;
      }
    }
  }
  return false;
}

void LayerMRI::ResetSurfaceRegionIds()
{
  for ( int i = 0; i < (int)m_surfaceRegions.size(); i++ )
  {
    m_surfaceRegions[i]->SetId( i+1 );
  }
}

bool LayerMRI::SaveAllSurfaceRegions( const QString& fn )
{
  FILE* fp = fopen( fn.toLatin1().data(), "w" );
  if ( !fp )
  {
    return false;
  }
  
  if (m_surfaceRegions.size() == 0)
  {
    cerr << "No surface regions to save.\n";
    return false;
  }
  
  bool ret = SurfaceRegion::WriteHeader( fp, this, m_surfaceRegions.size() );
  for ( int i = 0; i < m_surfaceRegions.size(); i++ )
  {
    if ( !m_surfaceRegions[i]->WriteBody( fp ) )
    {
      ret = false;
    }
  }
  fclose( fp );
  return ret;
}

bool LayerMRI::LoadSurfaceRegions( const QString& fn )
{
  FILE* fp = fopen( fn.toLatin1().data(), "r" );
  if ( !fp )
  {
    cerr << "Can not open file " << qPrintable(fn) << endl;
    return false;
  }
  int nNum = 0;
  char ch[1000];
  float dTh_low, dTh_high;
  fscanf( fp, "VOLUME_PATH %s\nVOLUME_THRESHOLD %f %f\nSURFACE_REGIONS %d", ch, &dTh_low, &dTh_high, &nNum );
  if ( nNum == 0 )
  {
    return false;
  }
  
  bool bSuccess = true;
  vtkActor* actor = m_actorCurrentContour;
  if (!actor)
    actor = m_actorContour;
  for ( int i = 0; i < nNum; i++ )
  {
    SurfaceRegion* r = new SurfaceRegion( this );
    connect( r, SIGNAL(ColorChanged(QColor)), this, SIGNAL(ActorUpdated()), Qt::UniqueConnection);
    if ( r->Load( fp ) )
    {
      r->SetInput( vtkPolyData::SafeDownCast( actor->GetMapper()->GetInput() ) );
      m_surfaceRegions.push_back( r );
      r->Highlight( false );
    }
    else
    {
      fclose( fp );
      bSuccess = false;
      break;
    }
  }
  
  ResetSurfaceRegionIds();
  if ( bSuccess )
  {
    GetProperty()->SetContourThreshold( dTh_low, dTh_high );
    GetProperty()->SetShowAsContour( true );
  }
  
  emit SurfaceRegionAdded();
  return bSuccess;
}

void LayerMRI::SetCroppingBounds( double* bounds )
{
  m_volumeSource->SetCroppingBounds( bounds );
}

void LayerMRI::SetCropToOriginal(bool bCropToOriginal)
{
  m_volumeSource->SetCropToOriginal(bCropToOriginal);
}

void LayerMRI::GetDisplayBounds( double* bounds )
{
  m_imageData->GetBounds( bounds );
  if ( mReslice[0].GetPointer() && mReslice[0]->GetAutoCropOutput() )
  {
    double d[6];
    m_imageData->GetBounds( d );
    vtkTransform* tr = vtkTransform::SafeDownCast( mReslice[0]->GetResliceTransform() );
    double pt[3];
    for ( int i = 0; i < 2; i++ )
    {
      for ( int j = 0; j < 2; j++ )
      {
        for ( int k = 0; k < 2; k++ )
        {
          pt[0] = d[i];
          pt[1] = d[2+j];
          pt[2] = d[4+k];
          tr->GetLinearInverse()->TransformPoint( pt, pt );
          if ( pt[0] < bounds[0] )
          {
            bounds[0] = pt[0];
          }
          if ( pt[0] > bounds[1] )
          {
            bounds[1] = pt[0];
          }
          if ( pt[1] < bounds[2] )
          {
            bounds[2] = pt[1];
          }
          if ( pt[1] > bounds[3] )
          {
            bounds[3] = pt[1];
          }
          if ( pt[2] < bounds[4] )
          {
            bounds[4] = pt[2];
          }
          if ( pt[2] > bounds[5] )
          {
            bounds[5] = pt[2];
          }
        }
      }
    }
  }
}

bool LayerMRI::SaveRegistration( const QString& filename )
{
  return m_volumeSource->SaveRegistration( filename );
}

struct LabelStatsPrivate
{
  int id;
  int count;
  std::vector<double> values;
};

void LayerMRI::GetLabelStats( LayerMRI* label, int nPlane,
                              std::vector<int>& ids,
                              std::vector<int>& numbers,
                              std::vector<double>& means,
                              std::vector<double>& stds )
{
  vtkImageData* mri_image = mReslice[nPlane]->GetOutput();
  vtkImageData* label_image = label->mReslice[nPlane]->GetOutput();
  int*    mri_dim = mri_image->GetDimensions();
  double* mri_vs = mri_image->GetSpacing();
  double* mri_orig = mri_image->GetOrigin();
  int mri_frames = mri_image->GetNumberOfScalarComponents();
  int*    label_dim = label_image->GetDimensions();
  double* label_vs = label_image->GetSpacing();
  double* label_orig = label_image->GetOrigin();
  char* mri_ptr = (char*)mri_image->GetScalarPointer();
  int   mri_scalar_type = mri_image->GetScalarType();
  char* label_ptr = (char*)label_image->GetScalarPointer();
  int   label_scalar_type = label_image->GetScalarType();
  
  // first find all label ids
  std::vector<LabelStatsPrivate> labels;
  for ( int i = 0; i < label_dim[0]; i++ )
  {
    for ( int j = 0; j < label_dim[1]; j++ )
    {
      int nId = (int)MyVTKUtils::GetImageDataComponent(label_ptr, label_dim, 1, i, j, 0, 0, label_scalar_type );
      if ( nId > 0 )
      {
        int mi = (int)( ( i*label_vs[0] + label_orig[0] - mri_orig[0] ) / mri_vs[0] );
        int mj = (int)( ( j*label_vs[1] + label_orig[1] - mri_orig[1] ) / mri_vs[1] );
        if ( mi >= 0 && mi < mri_dim[0] && mj >= 0 && mj < mri_dim[1] )
        {
          bool bFound = false;
          for ( size_t n = 0; n < labels.size(); n++ )
          {
            if ( nId == labels[n].id )
            {
              labels[n].values.push_back( MyVTKUtils::GetImageDataComponent(mri_ptr, mri_dim, mri_frames, mi, mj, 0, 0, mri_scalar_type ) );
              labels[n].count++;
              bFound = true;
              break;
            }
            else if ( nId < labels[n].id )
            {
              LabelStatsPrivate l;
              l.id = nId;
              l.count = 1;
              l.values.push_back( MyVTKUtils::GetImageDataComponent(mri_ptr, mri_dim, mri_frames, mi, mj, 0, 0, mri_scalar_type ) );
              labels.insert( labels.begin() + n, l );
              bFound = true;
              break;
            }
          }
          if ( !bFound )
          {
            LabelStatsPrivate l;
            l.id = nId;
            l.count = 1;
            l.values.push_back( MyVTKUtils::GetImageDataComponent(mri_ptr, mri_dim, mri_frames, mi, mj, 0, 0, mri_scalar_type ) );
            labels.push_back( l );
          }
        }
      }
    }
  }
  for ( size_t n = 0; n < labels.size(); n++ )
  {
    if ( labels[n].id > 0 )
    {
      ids.push_back( labels[n].id );
      numbers.push_back( labels[n].count );
      double dvalue = 0;
      for ( int i = 0; i < labels[n].count; i++ )
      {
        dvalue += labels[n].values[i];
      }
      double dmean = dvalue / labels[n].count;
      means.push_back( dmean );
      double sd = 0;
      for ( int i = 0; i < labels[n].count; i++ )
      {
        sd += (labels[n].values[i]-dmean)*(labels[n].values[i]-dmean);
      }
      if ( labels[n].count > 1 )
      {
        sd = sqrt( sd / ( labels[n].count-1 ) );
      }
      else
      {
        sd = 0;
      }
      stds.push_back( sd );
    }
  }
}

bool LayerMRI::SaveContourToFile(const QString &fn)
{
  MATRIX* mat = m_volumeSource->GetTargetToRASMatrix();
  double m[16];
  for ( int i = 0; i < 16; i++ )
  {
    m[i] = (double) *MATRIX_RELT((mat),(i/4)+1,(i%4)+1);
  }
  MatrixFree( &mat );
  vtkSmartPointer<vtkMatrix4x4> vmat = vtkSmartPointer<vtkMatrix4x4>::New();
  vmat->DeepCopy( m );
  vtkSmartPointer<vtkTransform> tr = vtkSmartPointer<vtkTransform>::New();
  tr->SetMatrix( vmat );
  vtkSmartPointer<vtkTransformPolyDataFilter> filter =
      vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  filter->SetTransform( tr );
#if VTK_MAJOR_VERSION > 5
  filter->SetInputData( vtkPolyDataMapper::SafeDownCast( m_actorContour->GetMapper())->GetInput() );
#else
  filter->SetInput( vtkPolyDataMapper::SafeDownCast( m_actorContour->GetMapper())->GetInput() );
#endif
  filter->Update();
  vtkWriter* writer;
  QFileInfo fi(fn);
  if (fi.suffix().toLower() == "stl")
  {
    writer = vtkSTLWriter::New();
    vtkSTLWriter::SafeDownCast(writer)->SetFileName(fn.toLatin1().constData());
  }
  else
  {
    writer = vtkPolyDataWriter::New();
    vtkPolyDataWriter::SafeDownCast(writer)->SetFileName(fn.toLatin1().constData());
  }
#if VTK_MAJOR_VERSION > 5
  writer->SetInputData( filter->GetOutput() );
#else
  writer->SetInput( filter->GetOutput() );
#endif
  bool ret = writer->Write();
  writer->Delete();
  return ret;
}

int LayerMRI::GoToLabel(int orientation, const QString& label_name)
{
  bool bOK;
  int nLabel = label_name.toInt(&bOK);
  if (!bOK)
  {
    COLOR_TABLE* ctab = GetProperty()->GetLUTCTAB ();
    if (ctab)
    {
      CTABfindName(ctab, qPrintable(label_name), &nLabel);
    }
    else
    {
      cerr << "Did not find the label name in the color table.";
      return -1;
    }
  }
  if (nLabel >= 0)
  {
    return ::MRIfindSliceWithMostStructure(m_volumeSource->GetMRI(), orientation, nLabel);
  }
  else
  {
    return -1;
  }
}

void LayerMRI::ReplaceVoxelValue(double orig_value, double new_value, int nPlane)
{
  this->SaveForUndo(-1);
  int* dim = m_imageData->GetDimensions();
  size_t range[3][2];
  range[0][0] = range[1][0] = range[2][0] = 0;
  range[0][1] = dim[0]-1;
  range[1][1] = dim[1]-1;
  range[2][1] = dim[2]-1;
  if (nPlane >= 0)
  {
    double* pos = GetSlicePosition();
    double* orig = m_imageData->GetOrigin();
    double* voxel_size = m_imageData->GetSpacing();
    int n[3];
    for ( int i = 0; i < 3; i++ )
    {
      n[i] = (int)( ( pos[i] - orig[i] ) / voxel_size[i] + 0.5 );
    }
    if (n[nPlane] >= 0 && n[nPlane] < dim[nPlane])
    {
      range[nPlane][0] = range[nPlane][1] = n[nPlane];
    }
    else
      range[nPlane][1] = range[nPlane][0]-1;
  }
  char* ptr = (char*)m_imageData->GetScalarPointer();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  for (size_t i = range[0][0]; i <= range[0][1]; i++)
  {
    for (size_t j = range[1][0]; j <= range[1][1]; j++)
    {
      for (size_t k = range[2][0]; k <= range[2][1]; k++)
      {
        double val = MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, i, j, k, m_nActiveFrame, scalar_type);
        if (val == orig_value)
        {
          MyVTKUtils::SetImageDataComponent(ptr, dim, n_frames, i, j, k, m_nActiveFrame, scalar_type, new_value);
        }
      }
    }
  }
  SetModified();
  emit ActorUpdated();
}

void LayerMRI::UpdateProjectionMap()
{
  if (GetProperty()->GetShowProjectionMap())
  {
    int m_dim[3];
    m_imageData->GetDimensions(m_dim);
    vtkSmartPointer<vtkImageData> images[3];
    float* ptrs[3];
    for (int i = 0; i < 3; i++)
    {
      vtkSmartPointer<vtkImageData> image = images[i] = vtkSmartPointer<vtkImageData>::New();
      int dim[3];
      m_imageData->GetDimensions(dim);
      image->SetSpacing(m_imageData->GetSpacing());
      dim[i] = 1;
      image->SetDimensions(dim);
#if VTK_MAJOR_VERSION > 5
      image->AllocateScalars(VTK_FLOAT, 1);
#else
      image->SetNumberOfScalarComponents(1);
      image->SetScalarTypeToFloat();
      image->AllocateScalars();
#endif
      float* ptr = ( float* )image->GetScalarPointer();
      memset(ptr, 0, ((size_t)sizeof(float))*dim[0]*dim[1]*dim[2]);
      ptrs[i] = ptr;
    }
    vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
#if VTK_MAJOR_VERSION > 5
    cast->SetInputData(m_imageData);
#else
    cast->SetInput(m_imageData);
#endif
    cast->SetOutputScalarTypeToFloat();
    cast->Update();
    vtkImageData* new_image = cast->GetOutput();
    float* ptr = (float*)new_image->GetScalarPointer();
    int nType = GetProperty()->GetProjectionMapType();
    int nRange[6];
    GetProperty()->GetProjectionMapRange(nRange);
    for (int i = 0; i < 3; i++)
    {
      if (nRange[i*2+1] < 0)
        nRange[i*2+1] = m_dim[i]-1;
    }
    for (qlonglong x = 0; x < m_dim[0]; x++)
    {
      for (qlonglong y = 0; y < m_dim[1]; y++)
      {
        for (qlonglong z = 0; z < m_dim[2]; z++)
        {
          float val = ptr[z*m_dim[0]*m_dim[1]+y*m_dim[0]+x];
          if (nType == LayerPropertyMRI::PM_Maximum)
          {
            if (x >= nRange[0] && x <= nRange[1] && ptrs[0][z*m_dim[1]+y] < val)
              ptrs[0][z*m_dim[1]+y] = val;
            if (y >= nRange[2] && y <= nRange[3] && ptrs[1][z*m_dim[0]+x] < val)
              ptrs[1][z*m_dim[0]+x] = val;
            if (z >= nRange[4] && z <= nRange[5] && ptrs[2][y*m_dim[0]+x] < val)
              ptrs[2][y*m_dim[0]+x] = val;
          }
          else if (nType == LayerPropertyMRI::PM_Mean)
          {
            ptrs[0][z*m_dim[1]+y] += val/(nRange[1]-nRange[0]);
            ptrs[1][z*m_dim[0]+x] += val/(nRange[3]-nRange[2]);
            ptrs[2][y*m_dim[0]+x] += val/(nRange[5]-nRange[4]);
          }
        }
      }
    }
    for (int i = 0; i < 3; i++)
    {
      vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
#if VTK_MAJOR_VERSION > 5
      reslice->SetInputData(images[i]);
#else
      reslice->SetInput(images[i]);
#endif
      reslice->BorderOff();
      //  reslice->SetResliceTransform( tr );
      reslice->SetOutputDimensionality( 2 );
      switch (i)
      {
      case 0:
        reslice->SetResliceAxesDirectionCosines( 0, 1, 0,
                                                 0, 0, 1,
                                                 1, 0, 0 );
        break;
      case 1:
        reslice->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 0, 1,
                                                 0, 1, 0 );
        break;
      case 2:
        reslice->SetResliceAxesDirectionCosines( 1, 0, 0,
                                                 0, 1, 0,
                                                 0, 0, 1 );
        break;
      }
      reslice->SetResliceAxesOrigin( 0, 0, 0 );
      if (true) // this->m_projectionMapActor[i]->GetInput() == NULL)
      {
        mColorMapMaxProjection[i] = vtkSmartPointer<vtkImageMapToColors>::New();
        mColorMapMaxProjection[i]->SetInputConnection(reslice->GetOutputPort());
        mColorMapMaxProjection[i]->SetLookupTable(GetProperty()->GetActiveLookupTable());
        m_projectionMapActor[i]->GetMapper()->SetInputConnection(mColorMapMaxProjection[i]->GetOutputPort());
        m_projectionMapActor[i]->InterpolateOff();
      }
    }
  }
  SetVisible(IsVisible());
  emit ActorUpdated();
}

double LayerMRI::GetTR()
{
  return m_volumeSource->GetMRI()->tr;
}

bool LayerMRI::SaveIsoSurface(const QString &fn)
{
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  vtkPolyData* polydata = vtkPolyData::SafeDownCast(m_actorContour->GetMapper()->GetInput());
  if (polydata)
  {
#if VTK_MAJOR_VERSION > 5
    writer->SetInputData( polydata );
#else
    writer->SetInput( polydata );
#endif
    writer->SetFileName( qPrintable(fn) );
    return writer->Write();
  }
  else
    return false;
}

bool LayerMRI::HasReg()
{
  return GetSourceVolume()->GetRegMatrix();
}

void LayerMRI::SetMaskLayer(LayerMRI *layer_mask)
{
  m_layerMask = layer_mask;
  if (GetProperty()->GetDisplayVector() || GetProperty()->GetDisplayTensor())
  {
    UpdateDisplayMode();
    GetProperty()->EmitChangeSignal();
    return;
  }

  vtkImageData* source = this->GetImageData();
  if (layer_mask == NULL)
  {
    if (m_imageDataBackup.GetPointer())
      source->DeepCopy(m_imageDataBackup);
  }
  else
  {
    vtkImageData* mask = layer_mask->GetImageData();
    if (!m_imageDataBackup.GetPointer())
    {
      m_imageDataBackup = vtkSmartPointer<vtkImageData>::New();
      m_imageDataBackup->DeepCopy(source);
    }

    double range[2];
    mask->GetScalarRange(range);
    if (range[1] <= 0)
      range[1] = 1;
    
    if (m_mapMaskThresholds.contains(layer_mask))
      m_dMaskThreshold = m_mapMaskThresholds[layer_mask];
    else
      m_dMaskThreshold = (range[0]+range[1])/5.0;
    
    double s1[3], s2[3];
    source->GetSpacing(s1);
    mask->GetSpacing(s2);
    if (mask->GetNumberOfScalarComponents() > 1)
    {
      vtkSmartPointer<vtkImageAppendComponents> append = vtkSmartPointer<vtkImageAppendComponents>::New();
      for (int i = 0; i < GetNumberOfFrames(); i++)
      {
        int nFrame = i;
        if (nFrame >= mask->GetNumberOfScalarComponents())
          nFrame = mask->GetNumberOfScalarComponents()-1;
        vtkSmartPointer<vtkImageExtractComponents> extractor = vtkSmartPointer<vtkImageExtractComponents>::New();
        extractor->SetInputData(mask);
        extractor->SetComponents(nFrame);
        vtkSmartPointer<vtkImageReslice> resampler = vtkSmartPointer<vtkImageReslice>::New();
        vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
        vtkSmartPointer<vtkImageMask> mask_filter = vtkSmartPointer<vtkImageMask>::New();
        resampler->SetInputConnection(extractor->GetOutputPort());
        resampler->SetOutputSpacing(s1);
        resampler->SetInterpolationModeToNearestNeighbor();
        threshold->ThresholdByUpper(m_dMaskThreshold+1e-12);
        threshold->SetInputConnection(resampler->GetOutputPort());
        threshold->ReplaceInOn();
        threshold->ReplaceOutOn();
        threshold->SetInValue(1);
        threshold->SetOutValue(0);
        threshold->SetOutputScalarType(VTK_UNSIGNED_CHAR);
        threshold->Update();
        extractor = vtkSmartPointer<vtkImageExtractComponents>::New();
        extractor->SetInputData(m_imageDataBackup);
        extractor->SetComponents(i);
        if (m_imageDataBackup->GetScalarType() == VTK_UNSIGNED_CHAR)
        {
          vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
          cast->SetOutputScalarType(VTK_SHORT);
          cast->SetInputConnection(extractor->GetOutputPort());
          cast->Update();
          mask_filter->SetImageInputData(cast->GetOutput());
        }
        else
        {
          extractor->Update();
          mask_filter->SetImageInputData(extractor->GetOutput());
        }
        mask_filter->SetMaskInputData(threshold->GetOutput());
        mask_filter->SetMaskedOutputValue(VTK_SHORT_MIN);
        append->AddInputConnection(0, mask_filter->GetOutputPort());
      }
      append->Update();
      source->DeepCopy(append->GetOutput());
    }
    else
    {
      vtkSmartPointer<vtkImageReslice> resampler = vtkSmartPointer<vtkImageReslice>::New();
      vtkSmartPointer<vtkImageMask> mask_filter = vtkSmartPointer<vtkImageMask>::New();
      vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
#if VTK_MAJOR_VERSION > 5
      resampler->SetInputData(mask);
#else
      resampler->SetInput(mask);
#endif
      resampler->SetOutputSpacing(s1);
      resampler->SetInterpolationModeToNearestNeighbor();
      threshold->ThresholdByUpper(m_dMaskThreshold);
      threshold->SetInputConnection(resampler->GetOutputPort());
      threshold->ReplaceInOn();
      threshold->ReplaceOutOn();
      threshold->SetInValue(1);
      threshold->SetOutValue(0);
      threshold->SetOutputScalarType(VTK_UNSIGNED_CHAR);
      threshold->Update();
#if VTK_MAJOR_VERSION > 5
      if (m_imageDataBackup->GetScalarType() == VTK_UNSIGNED_CHAR)
      {
        vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
        cast->SetOutputScalarType(VTK_SHORT);
        cast->SetInputData(m_imageDataBackup);
        cast->Update();
        mask_filter->SetImageInputData(cast->GetOutput());
      }
      else
        mask_filter->SetImageInputData(m_imageDataBackup);
      mask_filter->SetMaskInputData(threshold->GetOutput());
#else
      mask_filter->SetImageInputData(m_imageDataBackup);
      mask_filter->SetMaskInputData(threshold->GetOutput());
#endif
      mask_filter->SetMaskedOutputValue(VTK_SHORT_MIN);
      mask_filter->Update();
      source->DeepCopy(mask_filter->GetOutput());
    }
  }
  GetProperty()->OnColorMapChanged();
  emit ActorUpdated();
  GetProperty()->EmitChangeSignal();
}

void LayerMRI::Threshold(int frame, LayerMRI* src, int src_frame, double th_low, double th_high,
                         bool replace_in, double in_value, bool replace_out, double out_value)
{
  if (!m_imageDataBackup.GetPointer())
  {
    m_imageDataBackup = vtkSmartPointer<vtkImageData>::New();
    m_imageDataBackup->DeepCopy(GetImageData());
  }
  if (frame == -1)
  {
    for (int i = 0; i < GetNumberOfFrames(); i++)
    {
      Threshold(i, src, src_frame, th_low, th_high, replace_in, in_value, replace_out, out_value);
    }
  }
  else
  {
    vtkSmartPointer<vtkImageThreshold> threshold = vtkSmartPointer<vtkImageThreshold>::New();
    threshold->ThresholdBetween(th_low, th_high);
    vtkImageData* image = src->GetImageData();
    if (m_imageDataBackup.GetPointer() && src == this)
    {
      image = m_imageDataBackup;
    }
    if (image->GetNumberOfScalarComponents() > 1)
    {
      vtkSmartPointer<vtkImageExtractComponents> extract = vtkSmartPointer<vtkImageExtractComponents>::New();
#if VTK_MAJOR_VERSION > 5
      extract->SetInputData(image);
#else
      extract->SetInput(image);
#endif
      extract->SetComponents(src_frame);
      threshold->SetInputConnection(extract->GetOutputPort());
    }
    else
    {
#if VTK_MAJOR_VERSION > 5
      threshold->SetInputData(image);
#else
      threshold->SetInput(image);
#endif
    }
    threshold->SetOutputScalarTypeToChar();
    threshold->SetInValue(1.0);
    threshold->SetOutValue(0);
    threshold->ReplaceInOn();
    threshold->ReplaceOutOn();
    threshold->Update();
    vtkSmartPointer<vtkImageData> src_image = threshold->GetOutput();
    vtkImageData* target_image = GetImageData();
    char* src_ptr = (char*)src_image->GetScalarPointer();
    int *dim = src_image->GetDimensions();
    char* target_ptr = (char*)target_image->GetScalarPointer();
    int nFrames = target_image->GetNumberOfScalarComponents();
    int nBytes = target_image->GetScalarSize();
    char* backup_ptr = (char*)m_imageDataBackup->GetScalarPointer();
    char in_value_ptr[16], out_value_ptr[16];
    switch (target_image->GetScalarType())
    {
    case VTK_CHAR:
    case VTK_SIGNED_CHAR:
    case VTK_UNSIGNED_CHAR:
      in_value_ptr[0] = (char)in_value;
      out_value_ptr[0] = (char)out_value;
      break;
    case VTK_SHORT:
    {
      short* ptr = (short*)in_value_ptr;
      ptr[0] = (short)in_value;
      ptr = (short*)out_value_ptr;
      ptr[0] = (short)out_value;
    }
      break;
    case VTK_UNSIGNED_SHORT:
    {
      unsigned short* ptr = (unsigned short*)in_value_ptr;
      ptr[0] = (unsigned short)in_value;
      ptr = (unsigned short*)out_value_ptr;
      ptr[0] = (unsigned short)out_value;
    }
      break;
    case VTK_INT:
    {
      int* ptr = (int*)in_value_ptr;
      ptr[0] = (int)in_value;
      ptr = (int*)out_value_ptr;
      ptr[0] = (int)out_value;
    }
      break;
    case VTK_UNSIGNED_INT:
    {
      unsigned int* ptr = (unsigned int*)in_value_ptr;
      ptr[0] = (unsigned int)in_value;
      ptr = (unsigned int*)out_value_ptr;
      ptr[0] = (unsigned int)out_value;
    }
      break;
    case VTK_FLOAT:
    {
      float* ptr = (float*)in_value_ptr;
      ptr[0] = (float)in_value;
      ptr = (float*)out_value_ptr;
      ptr[0] = (float)out_value;
    }
      break;
    case VTK_DOUBLE:
    {
      double* ptr = (double*)in_value_ptr;
      ptr[0] = (double)in_value;
      ptr = (double*)out_value_ptr;
      ptr[0] = (double)out_value;
    }
      break;
    }
    
    for (qlonglong i = 0; i < dim[0]; i++)
    {
      for (qlonglong j = 0; j < dim[1]; j++)
      {
        for (qlonglong k = 0; k < dim[2]; k++)
        {
          qlonglong n = k*dim[0]*dim[1] + j*dim[0] + i;
          qlonglong offset = (n*nFrames+frame)*nBytes;
          if (src_ptr[n] < 1)
          {
            if (replace_out)
              memcpy(target_ptr + offset, out_value_ptr, nBytes);
            else
              memcpy(target_ptr + offset, backup_ptr + offset, nBytes);
          }
          else
          {
            if (replace_in)
              memcpy(target_ptr + offset, in_value_ptr, nBytes);
            else
              memcpy(target_ptr + offset, backup_ptr + offset, nBytes);
          }
        }
      }
    }
    SetModified();
    emit ActorUpdated();
    GetProperty()->EmitChangeSignal();
  }
}

bool LayerMRI::Segment(int min_label_index, int max_label_index, int min_num_of_voxels)
{
  if (!m_imageDataBackup.GetPointer())
  {
    m_imageDataBackup = vtkSmartPointer<vtkImageData>::New();
    m_imageDataBackup->DeepCopy(GetImageData());
  }
  
  if (!m_volumeSource->Segment(min_label_index, max_label_index, min_num_of_voxels))
  {
    return false;
  }
  SetModified();
  emit ActorUpdated();
  GetProperty()->EmitChangeSignal();
  return true;
}

void LayerMRI::RestoreFromBackup()
{
  if (m_imageDataBackup.GetPointer())
  {
    m_imageData->DeepCopy(m_imageDataBackup);
    emit ActorUpdated();
  }
}

void LayerMRI::SetCorrelationSurface(LayerSurface *surf)
{
  if (m_correlationSurface)
    disconnect(m_correlationSurface, 0, this, 0);
  
  m_correlationSurface = surf;
  if (m_correlationSurface)
  {
    this->SetActiveFrame(0);
    if (!m_imageRawDisplay.GetPointer())
    {
      m_imageRawDisplay = vtkSmartPointer<vtkImageData>::New();
      m_imageRawDisplay->SetDimensions(m_imageData->GetDimensions());
      m_imageRawDisplay->SetExtent(m_imageData->GetExtent());
      m_imageRawDisplay->SetSpacing(m_imageData->GetSpacing());
      m_imageRawDisplay->SetOrigin(m_imageData->GetOrigin());
#if VTK_MAJOR_VERSION > 5
      m_imageRawDisplay->AllocateScalars(VTK_FLOAT, 1);
#else
      m_imageRawDisplay->SetNumberOfScalarComponents(1);
      m_imageRawDisplay->SetScalarType(VTK_FLOAT);
      m_imageRawDisplay->AllocateScalars();
#endif
      GetProperty()->SetWindowLevel(1, 0);
      GetProperty()->SetHeatScale(0, 0.5, 1);
    }
    for ( int i = 0; i < 3; i++ )
    {
#if VTK_MAJOR_VERSION > 5
      mReslice[i]->SetInputData( m_imageRawDisplay );
#else
      mReslice[i]->SetInput( m_imageRawDisplay );
#endif
    }
    connect(m_correlationSurface, SIGNAL(CurrentVertexChanged(int)), this, SLOT(UpdateSurfaceCorrelationData()));
  }
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
#if VTK_MAJOR_VERSION > 5
      mReslice[i]->SetInputData( m_imageData );
#else
      mReslice[i]->SetInput( m_imageData );
#endif
    }
  }
  emit ActorUpdated();
  emit CorrelationSurfaceChanged(surf);
}

void LayerMRI::UpdateSurfaceCorrelationData()
{
  if (m_correlationSurface && m_correlationSurface->GetCurrentVertex() >= 0)
  {
    int nFrames = GetNumberOfFrames();
    float* buffer = new float[nFrames];
    
    m_correlationSurface->GetCorrelationOverlayDataAtVertex(m_correlationSurface->GetCurrentVertex(), buffer, nFrames);
    float* x;
    int dim[3];
    m_imageRawDisplay->GetDimensions(dim);
    float* inPixel = static_cast<float*>(m_imageData->GetScalarPointer());
    float* outPixel = static_cast<float*>(m_imageRawDisplay->GetScalarPointer());
    for (qlonglong i = 0; i < dim[0]; i++)
    {
      for (qlonglong j = 0; j < dim[1]; j++)
      {
        for (qlonglong k = 0; k < dim[2]; k++)
        {
          qlonglong nOffset = k*dim[0]*dim[1] + j*dim[0] + i;
          x = inPixel + nOffset*nFrames;
          bool masked = true;
          for (int n = 0; n < nFrames; n++)
          {
            if (x[n] != 0)
            {
              masked = false;
              break;
            }
          }
          if (!masked)
            outPixel[nOffset] = MyUtils::CalculateCorrelationCoefficient(buffer, x, nFrames);
          else
            outPixel[nOffset] = 0;
        }
      }
    }
    delete[] buffer;
    m_imageRawDisplay->Modified();
    emit ActorUpdated();
  }
}

double LayerMRI::GetHistoValueFromPercentile(double percentile)
{
  return m_volumeSource ? m_volumeSource->GetHistoValueFromPercentile(percentile, GetActiveFrame()) : 0;
}

double LayerMRI::GetHistoPercentileFromValue(double value)
{
  return m_volumeSource ? m_volumeSource->GetHistoPercentileFromValue(value, GetActiveFrame()) : 0;
}

bool LayerMRI::HasValidHistogram()
{
  return m_volumeSource ? m_volumeSource->HasValidHistogram() : false;
}

void LayerMRI::UpdateMRIToImage()
{
  m_volumeSource->MapMRIToImage(true);
  GetProperty()->SetVolumeSource(m_volumeSource);
  for (int i = 0; i < 3; i++)
    mReslice[i]->Modified();
  emit ActorUpdated();
}

bool LayerMRI::GetLayerLabelCenter(double val, double *pos_out)
{
  if (m_listLabelCenters.contains((int)val))
  {
    QList<double> center = m_listLabelCenters[(int)val];
    pos_out[0] = center[0];
    pos_out[1] = center[1];
    pos_out[2] = center[2];
    return true;
  }
  else
    return false;
}

bool LayerMRI::IsWindowAdjustable()
{
  return IsVisible() && GetProperty()->GetOpacity() > 0 && GetProperty()->GetColorMap() != LayerPropertyMRI::LUT &&
      GetProperty()->GetColorMap() != LayerPropertyMRI::DirectionCoded && !GetProperty()->GetDisplayVector();
}

bool LayerMRI::IsObscuring()
{
  return IsVisible() && GetProperty()->GetOpacity() == 1 && GetProperty()->GetColorMap() == LayerPropertyMRI::Grayscale;
}

void LayerMRI::OnLabelContourChanged(int n)
{
  QList<int> labels = GetProperty()->GetSelectedLabels();
  QList<int> keys = m_labelActors.keys();
  if (n >= 0 && keys.contains(n))
  {
    keys.clear();
    keys << n;
  }
  bool bVisible = IsVisible();
  foreach (int i, keys)
  {
    m_labelActors[i]->SetVisibility((bVisible && labels.contains(i))?1:0);
  }
  emit ActorUpdated();
}

void LayerMRI::RebuildContour()
{
  QList<int> keys = m_labelActors.keys();
  foreach (int i, keys)
  {
    m_labelActors[i]->Delete();
  }
  m_labelActors.clear();
  UpdateContour();
}

void LayerMRI::OnLabelInformationReady()
{
  if (GetProperty()->GetColorMap() == LayerPropertyMRI::LUT &&
      GetProperty()->GetShowAsContour() && m_labelActors.isEmpty())
  {
    UpdateContour();
  }
  
  emit LabelStatsReady();
}

void LayerMRI::SetMaskThreshold(double val)
{
  m_dMaskThreshold = val;
  if (m_layerMask)
    m_mapMaskThresholds[m_layerMask] = val;
  
  SetMaskLayer(m_layerMask);
}

VOXEL_LIST* LabelToVoxelList(MRI* mri, LABEL *area)
{
  double xd, yd, zd;

  VOXEL_LIST* vlist = VLSTalloc(area->n_points);
  vlist->nvox = 0;
  for (int n = 0; n < area->n_points; n++)
  {
    MRIworldToVoxel(mri, area->lv[n].x, area->lv[n].y, area->lv[n].z, &xd, &yd, &zd);
    VLSTadd(vlist, nint(xd), nint(yd), nint(zd), xd, yd, zd);
  }
  return vlist;
}

bool LayerMRI::GeodesicSegmentation(LayerMRI* seeds, double lambda, int wsize, double max_dist, double smoothing, LayerMRI *mask, double max_foreground_dist)
{
  if (!m_geos)
  {
    m_geos = new GeoSWorker;
    connect(m_geos, SIGNAL(ComputeFinished(double)), this, SIGNAL(GeodesicSegmentationFinished(double)));
    connect(m_geos, SIGNAL(Progress(double)), this, SIGNAL(GeodesicSegmentationProgress(double)));
  }

  m_geos->Compute((LayerMRI*)m_propertyBrush->GetReferenceLayer(), this, seeds, (int)max_dist, smoothing, mask, mask?mask->GetFillValue():-1, (int)max_foreground_dist);
  return true;
}

void LayerMRI::GeodesicSegmentationAbort()
{
  if (m_geos)
    m_geos->Abort();
}

void LayerMRI::GeodesicSegmentationApply(LayerMRI *filled)
{
  if (!m_geos)
    m_geos = new GeoSWorker;

  connect(m_geos, SIGNAL(ApplyFinished()), this, SIGNAL(GeodesicSegmentationApplied()), Qt::UniqueConnection);
  m_geos->Apply(this, filled);
}

void LayerMRI::GetVolumeInfo(int *dim, double *voxel_size)
{
  vtkImageData* image = m_volumeSource->GetImageOutput();
  image->GetDimensions(dim);
  image->GetSpacing(voxel_size);
}

QVector<double> LayerMRI::GetVoxelList(int nVal, bool bForce)
{
  if (!bForce && m_voxelLists.contains(nVal))
    return m_voxelLists[nVal];

  QVector<double> vlist;
  int* dim = m_imageData->GetDimensions();
  int scalar_type = m_imageData->GetScalarType();
  int n_frames = m_imageData->GetNumberOfScalarComponents();
  char* ptr = (char*)m_imageData->GetScalarPointer();
  double* origin = m_imageData->GetOrigin();
  double* vs = m_imageData->GetSpacing();
  for ( int k = 0; k < dim[2]; k++ )
  {
    for ( int j = 0; j < dim[1]; j++ )
    {
      for ( int i = 0; i < dim[0]; i++ )
      {
        int val = (int)MyVTKUtils::GetImageDataComponent(ptr, dim, n_frames, i, j, k, m_nActiveFrame, scalar_type);
        if (val == nVal)
        {
          vlist << i*vs[0] + origin[0] << j*vs[1] + origin[1] << k*vs[2] + origin[2];
        }
      }
    }
  }
  return vlist;
}

QVariantMap LayerMRI::GetTimeSeriesInfo()
{
  QVariantMap info;
  if (m_niftiHeader.sizeof_hdr == 348)
  {
    int t = XYZT_TO_TIME(m_niftiHeader.xyzt_units);
    if (t == NIFTI_UNITS_SEC)
      info["unit"] = "s";
    else if (t == NIFTI_UNITS_MSEC)
      info["unit"] = "ms";
    else if (t == NIFTI_UNITS_USEC)
      info["unit"] = "µs";
    else if (t == NIFTI_UNITS_HZ)
      info["unit"] = "Hz";
    else if (t == NIFTI_UNITS_PPM)
      info["unit"] = "ppm";
    else if (t == NIFTI_UNITS_RADS)
      info["unit"] = "rad/s";

    info["offset"] = m_niftiHeader.toffset;
    info["tr"] = m_niftiHeader.pixdim[4];
  }
  else
  {
    info["unit"] = "msec";
    info["offset"] = 0;
    info["tr"] = m_volumeSource->GetMRI()->tr;
  }
  return info;
}

QString LayerMRI::GetGeoSegErrorMessage()
{
  return m_geos?m_geos->GetErrorMessage():"";
}

bool LayerMRI::ExportLabelStats(const QString &fn)
{
  QFile file(fn);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    return false;

  QTextStream out(&file);
  out << "Label,Name,Count,Volume (mm3)\n";
  double* vs = m_imageData->GetSpacing();
  QList<int> labels = m_nAvailableLabels;
  qSort(labels);
  for (int i = 0; i < labels.size(); i++)
  {
    QVector<double> list = GetVoxelList(labels[i]);
    if (!list.isEmpty())
      out << QString("%1,%4,%2,%3\n").arg(labels[i])
             .arg(list.size()/3).arg(list.size()/3*vs[0]*vs[1]*vs[2]).arg(GetLabelName(labels[i]));
  }
  return true;
}

QList<vtkActor*> LayerMRI::GetContourActors(bool bVisibleOnly)
{
  QList<vtkActor*> actors;
  if (GetProperty()->GetShowAsLabelContour())
  {
    QList<int> keys = m_labelActors.keys();
    foreach (int n, keys)
    {
      if (!bVisibleOnly || m_labelActors[n]->GetVisibility())
        actors << m_labelActors[n];
    }
  }
  else
    actors << m_actorContour;
  return actors;
}

int LayerMRI::GetLabelCount(int nVal)
{
  if (m_labelVoxelCounts.contains(nVal))
    return m_labelVoxelCounts[nVal];
  else
    return 0;
}

Region3D* LayerMRI::CreateNew3DRegion( double* pt, vtkProp* prop )
{
  m_actorCurrentContour = vtkActor::SafeDownCast(prop);
  if (!m_actorCurrentContour)
    return NULL;

  Region3D* r = new Region3D( this );
  connect( r, SIGNAL(ColorChanged(QColor)), this, SIGNAL(ActorUpdated()), Qt::UniqueConnection);
  r->AddPoint( pt );
  m_3DRegions << r;
  int nGroup = 1;
  if ( m_current3DRegion )
    m_current3DRegion->Highlight( false );

  m_current3DRegion = r;
  emit Region3DAdded();
  return r;
}

bool LayerMRI::DeleteCurrent3DRegion()
{
  if ( m_current3DRegion )
  {
    for ( int i = 0; i < m_3DRegions.size(); i++ )
    {
      if ( m_3DRegions[i] == m_current3DRegion )
      {
        m_3DRegions.erase( m_3DRegions.begin() + i );
        m_current3DRegion->deleteLater();
        m_current3DRegion = NULL;
        emit Region3DRemoved();
        return true;
      }
    }
  }
  return false;
}

void LayerMRI::DeleteAll3DRegions()
{
  for ( int i = 0; i < m_3DRegions.size(); i++ )
  {
    m_3DRegions[i]->deleteLater();
  }
  m_3DRegions.clear();
  m_current3DRegion = NULL;
  emit Region3DRemoved();
}

void LayerMRI::Add3DRegionPoint( double* pt )
{
  if ( m_current3DRegion )
  {
    m_current3DRegion->AddPoint( pt );
    emit ActorUpdated();
  }
}

Region3D* LayerMRI::Select3DRegion( double* pos, double dist )
{
  for ( int i = 0; i < m_3DRegions.size(); i++ )
  {
    if ( m_3DRegions[i]->HasPoint( pos, dist ) )
    {
      if ( m_current3DRegion )
      {
        m_current3DRegion->Highlight( false );
      }
      m_current3DRegion = m_3DRegions[i];
      m_current3DRegion->Highlight( true );
      return m_current3DRegion;
    }
  }

  return NULL;
}

bool LayerMRI::SaveAll3DRegions( const QString& fn )
{
  FILE* fp = fopen( fn.toLatin1().data(), "w" );
  if ( !fp )
  {
    return false;
  }

  if (m_3DRegions.size() == 0)
  {
    cerr << "No surface regions to save.\n";
    return false;
  }

  bool ret = Region3D::WriteHeader( fp, this, m_3DRegions.size() );
  for ( int i = 0; i < m_3DRegions.size(); i++ )
  {
    if ( !m_3DRegions[i]->WriteBody( fp ) )
    {
      ret = false;
    }
  }
  fclose( fp );
  return ret;
}

bool LayerMRI::Load3DRegions( const QString& fn )
{
  FILE* fp = fopen( fn.toLatin1().data(), "r" );
  if ( !fp )
  {
    cerr << "Can not open file " << qPrintable(fn) << endl;
    return false;
  }
  int nNum = 0;
  char ch[1000];
  float dTh_low, dTh_high;
  fscanf( fp, "VOLUME_PATH %s\nVOLUME_THRESHOLD %f %f\nNUM_OF_REGIONS %d", ch, &dTh_low, &dTh_high, &nNum );
  if ( nNum == 0 )
  {
    return false;
  }

  bool bSuccess = true;
  vtkActor* actor = m_actorCurrentContour;
  if (!actor)
    actor = m_actorContour;
  for ( int i = 0; i < nNum; i++ )
  {
    Region3D* r = new Region3D( this );
    connect( r, SIGNAL(ColorChanged(QColor)), this, SIGNAL(ActorUpdated()), Qt::UniqueConnection);
    if ( r->Load( fp ) )
    {
      m_3DRegions << r;
      r->Highlight( false );
    }
    else
    {
      fclose( fp );
      bSuccess = false;
      break;
    }
  }

  emit Region3DAdded();
  return bSuccess;
}


void LayerMRI::Close3DRegion()
{
  if ( m_current3DRegion )
  {
    m_current3DRegion->Close();
    emit ActorUpdated();
  }
}

void LayerMRI::UpdateVoxelsByPointSet(LayerPointSet *ps, int nPlane)
{
  vtkPoints* pts = ps->GetSplinedPoints();
  if (pts)
  {
    double pt1[3], pt2[3];
    double* pos = GetSlicePosition();
    for (int i = 0; i < pts->GetNumberOfPoints()-1; i++)
    {
      pts->GetPoint(i, pt1);
      pts->GetPoint(i+1, pt2);
      pt1[nPlane] = pt2[nPlane] = pos[nPlane];
      this->SetVoxelByRAS(pt1, pt2, nPlane, true, true);
    }
  }
}

void LayerMRI::LocateLocalMaximumAtRAS(double* ras_in, double dx, double dy, double dz, double* ras_out, double sigma, double dist_in_vox)
{
  double ras[3];
  TargetToRAS(ras_in, ras);
  RASToOriginalIndex(ras, ras);
  double pt0[3] = {0, 0, 0}, pt1[3] = {dx, dy, dz};
  TargetToRAS(pt0, pt0);
  RASToOriginalIndex(pt0, pt0);
  TargetToRAS(pt1, pt1);
  RASToOriginalIndex(pt1, pt1);
  double r = sqrt(vtkMath::Distance2BetweenPoints(pt0, pt1));
  dx = (pt1[0] - pt0[0])/r;
  dy = (pt1[1] - pt0[1])/r;
  dz = (pt1[2] - pt0[2])/r;

  double dMax = 0;
  double x = ras[0], y = ras[1], z = ras[2];
  double x_out = x, y_out = y, z_out = z;
  ::MRIsampleVolumeDerivativeScale(m_volumeSource->GetMRI(), x, y, z, dx, dy, dz, &dMax, sigma);
  dMax = qAbs(dMax);
  for (double d = 0; d <= dist_in_vox; d += 0.25)
  {
    double mag;
    x = ras[0] + d*dx;
    y = ras[1] + d*dy;
    z = ras[2] + d*dz;
    ::MRIsampleVolumeDerivativeScale(m_volumeSource->GetMRI(), x, y, z, dx, dy, dz, &mag, sigma);
    if (qAbs(mag) > dMax)
    {
      dMax = qAbs(mag);
      x_out = x;
      y_out = y;
      z_out = z;
    }
    x = ras[0] - d*dx;
    y = ras[1] - d*dy;
    z = ras[2] - d*dz;
    ::MRIsampleVolumeDerivativeScale(m_volumeSource->GetMRI(), x, y, z, dx, dy, dz, &mag, sigma);
    if (qAbs(mag) > dMax)
    {
      dMax = qAbs(mag);
      x_out = x;
      y_out = y;
      z_out = z;
    }
  }

  ras[0] = x_out;
  ras[1] = y_out;
  ras[2] = z_out;
  OriginalVoxelToRAS(ras, ras);
  RASToTarget(ras, ras_out);
}
