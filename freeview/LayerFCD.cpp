/**
 * @file  LayerFCD.cpp
 * @brief Support for FCD functionality
 *
 */
/*
 * Original Author: Ruopeng Wang 
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/17 20:36:46 $
 *    $Revision: 1.15 $
 *
 * Copyright © 2014 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include "LayerFCD.h"
#include "LayerPropertyFCD.h"
#include "vtkRenderer.h"
#include "vtkImageReslice.h"
#include "vtkImageActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapToColors.h"
#include "vtkImageFlip.h"
#include "vtkTransform.h"
#include "vtkPlaneSource.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTexture.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkRGBAColorTransferFunction.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include "vtkFreesurferLookupTable.h"
#include "vtkImageChangeInformation.h"
#include "LayerPropertyROI.h"
#include "LayerMRI.h"
#include "LayerPropertyMRI.h"
#include "LayerFCDWorkerThread.h"
#include "FSVolume.h"
#include "ProgressCallback.h"
#include "LayerSurface.h"
#include "LayerPropertySurface.h"
#include <QDebug>
#include <QTimer>

LayerFCD::LayerFCD(LayerMRI* layerMRI, 
                   QObject *parent) : LayerVolumeBase(parent),
  m_fcd(NULL),
  m_mri_difference(NULL)
{
  m_strTypeNames.push_back( "FCD" );
  m_sPrimaryType = "FCD";
  mProperty = new LayerPropertyFCD(this);
  for ( int i = 0; i < 3; i++ )
  {
    m_dSlicePosition[i] = 0;
    m_sliceActor2D[i] = vtkSmartPointer<vtkImageActor>::New();
    m_sliceActor3D[i] = vtkSmartPointer<vtkImageActor>::New();
    m_sliceActor2D[i]->InterpolateOff();
    m_sliceActor3D[i]->InterpolateOff();
  }

  mProperty = new LayerPropertyFCD( this );

  m_layerSource = layerMRI;
  if (m_layerSource)
  {
    InitializeData();
  }

  LayerPropertyFCD* p = GetProperty();
  connect( p, SIGNAL(ColorMapChanged()), this, SLOT(UpdateColorMap()) );
  connect( p, SIGNAL(OpacityChanged(double)), this, SLOT(UpdateOpacity()) );
  connect( p, SIGNAL(ThicknessThresholdChanged(double)), 
           this, SLOT(Recompute()));
  connect( p, SIGNAL(SigmaChanged(double)), this, SLOT(Recompute()));
  connect( p, SIGNAL(MinAreaChanged(int)), this, SLOT(Recompute()));

  // pre allocate MRIs & surfaces
  m_mri_norm = new LayerMRI(NULL);
  m_mri_flair = new LayerMRI(NULL);
  m_mri_t2 = new LayerMRI(NULL);
  m_mri_aseg = new LayerMRI(NULL);
  m_mri_difference = new LayerMRI(NULL);
  connect(m_mri_norm, SIGNAL(destroyed()), 
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  connect(m_mri_flair, SIGNAL(destroyed()), 
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  connect(m_mri_t2, SIGNAL(destroyed()), 
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  connect(m_mri_aseg, SIGNAL(destroyed()), 
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  connect(m_mri_difference, SIGNAL(destroyed()),
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);

  m_surf_lh = new LayerSurface(NULL);
  m_surf_rh = new LayerSurface(NULL);
  connect(m_surf_lh, SIGNAL(destroyed()),
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  connect(m_surf_rh, SIGNAL(destroyed()),
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  m_surf_lh_pial = new LayerSurface(NULL);
  m_surf_rh_pial = new LayerSurface(NULL);
  connect(m_surf_lh_pial, SIGNAL(destroyed()),
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  connect(m_surf_rh_pial, SIGNAL(destroyed()),
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);

  m_surf_lh_sphere_d1 = new LayerSurface(NULL);
  m_surf_rh_sphere_d1 = new LayerSurface(NULL);
  connect(m_surf_lh_sphere_d1, SIGNAL(destroyed()),
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);
  connect(m_surf_rh_sphere_d1, SIGNAL(destroyed()),
          this, SLOT(OnLayerDestroyed()), Qt::UniqueConnection);

  m_worker = new LayerFCDWorkerThread(this);
  connect(m_worker, SIGNAL(started()), this, SIGNAL(StatusChanged()));
  connect(m_worker, SIGNAL(finished()), this, SIGNAL(StatusChanged()));
  connect(m_worker, SIGNAL(terminated()), this, SIGNAL(StatusChanged()));
}

LayerFCD::~LayerFCD()
{
  if (m_fcd)
  {
    // if layers still exist, do not free mri and mris 
    // because they are still being shared.
    if (m_mri_norm)
    {
      m_fcd->mri_norm = NULL;
    }
    if (m_mri_flair)
    {
      m_fcd->mri_flair = NULL;
    }
    if (m_mri_t2)
    {
      m_fcd->mri_t2 = NULL;
    }
    if (m_mri_aseg)
    {
      m_fcd->mri_aseg = NULL;
    }
    if (m_mri_difference)
    {
      m_fcd->mri_thickness_difference = NULL;
    }
    if (m_surf_lh)
    {
      m_fcd->mris_lh = NULL;
    }
    if (m_surf_rh)
    {
      m_fcd->mris_rh = NULL;
    }
    if (m_surf_lh_pial)
    {
      m_fcd->mris_lh_pial = NULL;
    }
    if (m_surf_rh_pial)
    {
      m_fcd->mris_rh_pial = NULL;
    }
    if (m_surf_lh_sphere_d1)
    {
        m_fcd->mris_lh_sphere_d1 = NULL;
    }
    if (m_surf_rh_sphere_d1)
    {
        m_fcd->mris_rh_sphere_d1 = NULL;
    }

    FCDfree(&m_fcd);
  }
}

bool LayerFCD::IsBusy()
{
  return m_worker->isRunning();
}

void LayerFCD::InitializeData()
{
  if ( m_layerSource )
  {
    GetProperty()->SetLookupTable(m_layerSource->GetProperty()->GetLUTTable());
    m_imageDataRef = m_layerSource->GetImageData();
    SetWorldOrigin( m_layerSource->GetWorldOrigin() );
    SetWorldVoxelSize( m_layerSource->GetWorldVoxelSize() );
    SetWorldSize( m_layerSource->GetWorldSize() );

    m_imageData = vtkSmartPointer<vtkImageData>::New();
    // m_imageData->DeepCopy( m_layerSource->GetRASVolume() );

    m_imageData->SetNumberOfScalarComponents( 1 );
    m_imageData->SetScalarTypeToUnsignedChar();
    m_imageData->SetOrigin( GetWorldOrigin() );
    m_imageData->SetSpacing( GetWorldVoxelSize() );
    m_imageData->SetDimensions( 
      ( int )( m_dWorldSize[0] / m_dWorldVoxelSize[0] + 0.5 ),
      ( int )( m_dWorldSize[1] / m_dWorldVoxelSize[1] + 0.5 ),
      ( int )( m_dWorldSize[2] / m_dWorldVoxelSize[2] + 0.5 ) );
    m_imageData->AllocateScalars();
    void* ptr = m_imageData->GetScalarPointer();
    int* nDim = m_imageData->GetDimensions();
    // cout << nDim[0] << ", " << nDim[1] << ", " << nDim[2] << endl;
    memset( ptr, 0, m_imageData->GetScalarSize() * 
            nDim[0] * nDim[1] * nDim[2] );
    InitializeActors();
  }
}

bool LayerFCD::Load(const QString &subdir, const QString &subject)
{
  m_sSubjectDir = subdir;
  if (m_sSubjectDir[m_sSubjectDir.length()-1] == '/')
  {
    m_sSubjectDir.resize(m_sSubjectDir.size()-1);
  }
  m_sSubject = subject;
  SetFileName(m_sSubjectDir + "/" + m_sSubject);

  return LoadFromFile();
}

bool LayerFCD::LoadFromFile()
{
  ::SetProgressCallback(ProgressCallback, 0, 50);
  try
  {
    m_fcd = ::FCDloadData(m_sSubjectDir.toAscii().data(), 
                          m_sSubject.toAscii().data());
  }
  catch (int ret)
  {
    return false;
  }

  if (m_fcd)
  {
    if (!m_fcd->mri_norm)
    {
      cerr << "Did not find norm volume for FCD data" << endl;
      return false;
    }
    ::SetProgressCallback(ProgressCallback, 50, 70);
    MakeAllLayers();
    ::SetProgressCallback(ProgressCallback, 70, 100);
    DoCompute(false);
  }

  return (m_fcd != NULL);
}

void LayerFCD::MakeAllLayers()
{
  if (true)
  {
    LayerMRI* mri = m_mri_norm;
    if (m_layerSource)
    {
      mri->SetRefVolume(m_layerSource->GetSourceVolume());
    }
    mri->SetName(GetName() + "_norm");
    mri->SetFileName(m_fcd->mri_norm->fname);
    if ( mri->CreateFromMRIData((void*)m_fcd->mri_norm) )
    {
      if (!m_layerSource)
      {
        m_layerSource = mri;
        InitializeData();
      }
    }
    else
    {
      cerr << "Failed to create norm layer" << endl;
      delete mri;
      m_mri_norm = NULL;
    }
  }

  if (m_fcd->mri_flair)
  {
    LayerMRI* mri = m_mri_flair;
    if (m_layerSource)
    {
      mri->SetRefVolume(m_layerSource->GetSourceVolume());
    }
    mri->SetName(GetName() + "_flair");
    mri->SetFileName(m_fcd->mri_flair->fname);
    if ( !mri->CreateFromMRIData((void*)m_fcd->mri_flair) )
    {
      delete m_mri_flair;
      m_mri_flair = NULL;
    }
  }
  else
  {
    delete m_mri_flair;
    m_mri_flair = NULL;
  }

  if (m_fcd->mri_t2)
  {
    LayerMRI* mri = m_mri_t2;
    if (m_layerSource)
    {
      mri->SetRefVolume(m_layerSource->GetSourceVolume());
    }
    mri->SetName(GetName() + "_t2");
    mri->SetFileName(m_fcd->mri_t2->fname);
    if ( !mri->CreateFromMRIData((void*)m_fcd->mri_t2) )
    {
      delete m_mri_t2;
      m_mri_t2 = NULL;
    }
  }
  else
  {
    delete m_mri_t2;
    m_mri_t2 = NULL;
  }

  if (m_fcd->mri_aseg)
  {
    LayerMRI* mri = m_mri_aseg;
    if (m_layerSource)
    {
      mri->SetRefVolume(m_layerSource->GetSourceVolume());
    }
    mri->SetName(GetName() + "_aseg");
    mri->SetFileName(m_fcd->mri_aseg->fname);
    if ( mri->CreateFromMRIData((void*)m_fcd->mri_aseg) )
    {
      mri->GetProperty()->SetColorMap(LayerPropertyMRI::LUT);
      mri->SetVisible(false);
    }
    else
    {
      delete m_mri_aseg;
      m_mri_aseg = NULL;
    }
  }
  else
  {
    delete m_mri_aseg;
    m_mri_aseg = NULL;
  }

  if (m_fcd->mri_thickness_difference)
  {
    LayerMRI* mri = m_mri_difference;
    if (m_layerSource)
    {
      mri->SetRefVolume(m_layerSource->GetSourceVolume());
    }
    mri->SetName(GetName() + "_thickness_difference");
//    mri->SetFileName(m_fcd->mri_thickness_increase->fname);
    if ( mri->CreateFromMRIData((void*)m_fcd->mri_thickness_difference) )
    {
      mri->GetProperty()->SetColorMap(LayerPropertyMRI::Heat);
    }
    else
    {
      delete m_mri_difference;
      m_mri_difference = NULL;
    }
  }
  else
  {
    delete m_mri_difference;
    m_mri_difference = NULL;
  }

  if (m_fcd->mris_lh)
  {
    LayerSurface* surf = m_surf_lh;
    if (m_layerSource)
    {
      surf->SetRefVolume(m_layerSource);
    }
    surf->SetName(GetName() + "_lh");
    if (!surf->CreateFromMRIS((void*)m_fcd->mris_lh))
    {
      delete m_surf_lh;
      m_surf_lh = NULL;
    }
  }
  else
  {
    delete m_surf_lh;
    m_surf_lh = NULL;
  }

  if (m_fcd->mris_lh_pial)
  {
    LayerSurface* surf = m_surf_lh_pial;
    if (m_layerSource)
    {
      surf->SetRefVolume(m_layerSource);
    }
    surf->SetName(GetName() + "_lh.pial");
    if (!surf->CreateFromMRIS((void*)m_fcd->mris_lh_pial))
    {
      delete m_surf_lh_pial;
      m_surf_lh_pial = NULL;
    }
    else
    {
      surf->GetProperty()->SetEdgeColor(Qt::green);
    }
  }
  else
  {
    delete m_surf_lh_pial;
    m_surf_lh_pial = NULL;
  }

  if (m_fcd->mris_lh_sphere_d1)
  {
    LayerSurface* surf = m_surf_lh_sphere_d1;
    if (m_layerSource)
    {
      surf->SetRefVolume(m_layerSource);
    }
    surf->SetName(GetName() + "_lh.sphere.d1");
    if (!surf->CreateFromMRIS((void*)m_fcd->mris_lh_sphere_d1))
    {
      delete m_surf_lh_sphere_d1;
      m_surf_lh_sphere_d1 = NULL;
    }
    else
    {
      surf->GetProperty()->SetEdgeColor(Qt::red);
      surf->SetVisible(false);
    }
  }
  else
  {
    delete m_surf_lh_sphere_d1;
    m_surf_lh_sphere_d1 = NULL;
  }

  if (m_fcd->mris_rh)
  {
    LayerSurface* surf = m_surf_rh;
    if (m_layerSource)
    {
      surf->SetRefVolume(m_layerSource);
    }
    surf->SetName(GetName() + "_rh");
    if (!surf->CreateFromMRIS((void*)m_fcd->mris_rh))
    {
      delete m_surf_rh;
      m_surf_rh = NULL;
    }
  }
  else
  {
    delete m_surf_rh;
    m_surf_rh = NULL;
  }

  if (m_fcd->mris_rh_pial)
  {
    LayerSurface* surf = m_surf_rh_pial;
    if (m_layerSource)
    {
      surf->SetRefVolume(m_layerSource);
    }
    surf->SetName(GetName() + "_rh.pial");
    if (!surf->CreateFromMRIS((void*)m_fcd->mris_rh_pial))
    {
      delete m_surf_rh_pial;
      m_surf_rh_pial = NULL;
    }
    else
    {
      surf->GetProperty()->SetEdgeColor(Qt::green);
    }
  }
  else
  {
    delete m_surf_rh_pial;
    m_surf_rh_pial = NULL;
  }

  if (m_fcd->mris_rh_sphere_d1)
  {
    LayerSurface* surf = m_surf_rh_sphere_d1;
    if (m_layerSource)
    {
      surf->SetRefVolume(m_layerSource);
    }
    surf->SetName(GetName() + "_rh.sphere.d1");
    if (!surf->CreateFromMRIS((void*)m_fcd->mris_rh_sphere_d1))
    {
      delete m_surf_rh_sphere_d1;
      m_surf_rh_sphere_d1 = NULL;
    }
    else
    {
      surf->GetProperty()->SetEdgeColor(Qt::red);
      surf->SetVisible(false);
    }
  }
  else
  {
    delete m_surf_rh_sphere_d1;
    m_surf_rh_sphere_d1 = NULL;
  }
}


void LayerFCD::Recompute()
{
  if (m_worker->isRunning())
  {
    m_worker->wait();
  }

  m_worker->start(QThread::HighPriority);
}

void LayerFCD::DoCompute(bool resetProgress)
{
  if (m_fcd)
  {
    if (resetProgress)
    {
      ::SetProgressCallback(ProgressCallback, 0, 50);
    }
    try
    {
      ::FCDcomputeThicknessLabels(m_fcd, 
                                  GetProperty()->GetThicknessThreshold(),
                                  GetProperty()->GetSigma(),
                                  GetProperty()->GetMinArea());
    }
    catch (int ret)
    {
      return;
    }

    m_labelVisibility.clear();
    for (int i = 0; i < m_fcd->nlabels; i++)
    {
      m_labelVisibility << true;
    }
    UpdateRASImage(m_imageData);
    for (int i = 0; i < 3; i++)
    {
      mReslice[i]->Modified();
    }

    if (resetProgress)
    {
      ::SetProgressCallback(ProgressCallback, 50, 60);
    }
    if (m_mri_difference)
    {
      m_mri_difference->UpdateMRIToImage();
    }
    if (resetProgress)
    {
      ::SetProgressCallback(ProgressCallback, 60, 100);
    }
//    if (m_mri_decrease)
//    {
//      m_mri_decrease->UpdateMRIToImage();
//    }

    emit LabelsChanged();
    emit ActorUpdated();
  }
}

void LayerFCD::UpdateRASImage( vtkImageData* rasImage)
{
  if (!m_fcd)
  {
    return;
  }

  int n[3];
  double pos[3];
  int* dim = rasImage->GetDimensions();
  memset( rasImage->GetScalarPointer(),
          0,
          dim[0] * dim[1] * dim[2] * rasImage->GetScalarSize() );
  if ( m_fcd->nlabels == 0 )
  {
    cout << "No labels found\n";
    return;
  }

  for (int cnt = 0; cnt < m_fcd->nlabels; cnt++)
  {
    LABEL* label = m_fcd->labels[cnt];
    FSVolume* ref_vol = m_layerSource->GetSourceVolume();
    for ( int i = 0; i < label->n_points; i++ )
    {
      pos[0] = label->lv[i].x;
      pos[1] = label->lv[i].y;
      pos[2] = label->lv[i].z;
      if ( label->coords == LABEL_COORDS_VOXEL )
      {
        MRIvoxelToWorld(ref_vol->GetMRI(),
                        pos[0], pos[1], pos[2], 
                        pos, pos+1, pos+2);
      }
      else if (label->coords == LABEL_COORDS_TKREG_RAS)
      {
        ref_vol->TkRegToNativeRAS( pos, pos );
      }
      ref_vol->NativeRASToRAS( pos, pos );
      ref_vol->RASToTargetIndex( pos, n );
      if ( n[0] >= 0 && n[0] < dim[0] && n[1] >= 0 && n[1] < dim[1] &&
           n[2] >= 0 && n[2] < dim[2] )
      {
        rasImage->SetScalarComponentFromFloat
          ( n[0], n[1], n[2], 0, label->lv[i].vno );
      }
    }
  }
}


void LayerFCD::InitializeActors()
{
  if ( m_layerSource == NULL )
  {
    return;
  }

  for ( int i = 0; i < 3; i++ )
  {
    // The reslice object just takes a slice out of the volume.
    //
    mReslice[i] = vtkSmartPointer<vtkImageReslice>::New();
    mReslice[i]->SetInput( m_imageData );
//  mReslice[i]->SetOutputSpacing( sizeX, sizeY, sizeZ );
    mReslice[i]->BorderOff();

    // This sets us to extract slices.
    mReslice[i]->SetOutputDimensionality( 2 );

    // This will change depending what orienation we're in.
    mReslice[i]->SetResliceAxesDirectionCosines( 1, 0, 0,
        0, 1, 0,
        0, 0, 1 );

    // This will change to select a different slice.
    mReslice[i]->SetResliceAxesOrigin( 0, 0, 0 );
    //
    // Image to colors using color table.
    //
    mColorMap[i] = vtkSmartPointer<vtkImageMapToColors>::New();
    mColorMap[i]->SetLookupTable( GetProperty()->GetLookupTable() );
    mColorMap[i]->SetInputConnection( mReslice[i]->GetOutputPort() );
    mColorMap[i]->SetOutputFormatToRGBA();
    mColorMap[i]->PassAlphaToOutputOn();

    //
    // Prop in scene with plane mesh and texture.
    //
    m_sliceActor2D[i]->SetInput( mColorMap[i]->GetOutput() );
    m_sliceActor3D[i]->SetInput( mColorMap[i]->GetOutput() );

    // Set ourselves up.
    this->OnSlicePositionChanged( i );
  }

  this->UpdateOpacity();
  this->UpdateColorMap();
}

void LayerFCD::UpdateOpacity()
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetOpacity( GetProperty()->GetOpacity() );
    m_sliceActor3D[i]->SetOpacity( GetProperty()->GetOpacity() );
  }

  emit ActorUpdated();
}

void LayerFCD::UpdateColorMap ()
{
  for ( int i = 0; i < 3; i++ )
  {
    mColorMap[i]->SetLookupTable( GetProperty()->GetLookupTable() );
  }
  // m_sliceActor2D[i]->GetProperty()->SetColor(1, 0, 0);

  emit ActorUpdated();
}

void LayerFCD::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  renderer->AddViewProp( m_sliceActor2D[nPlane] );
}

void LayerFCD::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( bSliceVisibility == NULL || bSliceVisibility[i] )
    {
      renderer->AddViewProp( m_sliceActor3D[i] );
    }
  }
}

void LayerFCD::OnSlicePositionChanged( int nPlane )
{
  if ( !m_layerSource )
  {
    return;
  }

  vtkSmartPointer<vtkMatrix4x4> matrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  switch ( nPlane )
  {
  case 0:
    m_sliceActor2D[0]->PokeMatrix( matrix );
    m_sliceActor2D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
    m_sliceActor2D[0]->RotateX( 90 );
    m_sliceActor2D[0]->RotateY( -90 );
    m_sliceActor3D[0]->PokeMatrix( matrix );
    m_sliceActor3D[0]->SetPosition( m_dSlicePosition[0], 0, 0 );
    m_sliceActor3D[0]->RotateX( 90 );
    m_sliceActor3D[0]->RotateY( -90 );

    // Putting negatives in the reslice axes cosines will flip the
    // image on that axis.
    mReslice[0]->SetResliceAxesDirectionCosines( 0, -1, 0,
        0, 0, 1,
        1, 0, 0 );
    mReslice[0]->SetResliceAxesOrigin( m_dSlicePosition[0], 0, 0  );
    mReslice[0]->Modified();
    break;
  case 1:
    m_sliceActor2D[1]->PokeMatrix( matrix );
    m_sliceActor2D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
    m_sliceActor2D[1]->RotateX( 90 );
    // m_sliceActor2D[1]->RotateY( 180 );
    m_sliceActor3D[1]->PokeMatrix( matrix );
    m_sliceActor3D[1]->SetPosition( 0, m_dSlicePosition[1], 0 );
    m_sliceActor3D[1]->RotateX( 90 );
    // m_sliceActor3D[1]->RotateY( 180 );

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

    mReslice[2]->SetResliceAxesDirectionCosines( 1, 0, 0,
        0, 1, 0,
        0, 0, 1 );
    mReslice[2]->SetResliceAxesOrigin( 0, 0, m_dSlicePosition[2]  );
    mReslice[2]->Modified();
    break;
  }
}

void LayerFCD::SetVisible( bool bVisible )
{
  for ( int i = 0; i < 3; i++ )
  {
    m_sliceActor2D[i]->SetVisibility( bVisible ? 1 : 0 );
    m_sliceActor3D[i]->SetVisibility( bVisible ? 1 : 0 );
  }
  LayerVolumeBase::SetVisible(bVisible);
}

bool LayerFCD::IsVisible()
{
  return m_sliceActor2D[0]->GetVisibility() > 0;
}

bool LayerFCD::HasProp( vtkProp* prop )
{
  for ( int i = 0; i < 3; i++ )
  {
    if ( prop == m_sliceActor2D[i].GetPointer() || 
         prop == m_sliceActor3D[i].GetPointer() )
    {
      return true;
    }
  }
  return false;
}

void LayerFCD::GetLabelCentroidPosition(int nLabelIndex, double *pos)
{
  if (nLabelIndex >= m_fcd->nlabels)
  {
    return;
  }

  LABEL* label = m_fcd->labels[nLabelIndex];
  if (label->n_points > 0)
  {
    double x = 0, y = 0, z = 0;
    for ( int i = 0; i < label->n_points; i++ )
    {
      x += label->lv[i].x;
      y += label->lv[i].y;
      z += label->lv[i].z;
    }
    pos[0] = x / label->n_points;
    pos[1] = y / label->n_points;
    pos[2] = z / label->n_points;

    FSVolume* ref_vol = m_layerSource->GetSourceVolume();
    if ( label->coords == LABEL_COORDS_VOXEL )
    {
      MRIvoxelToWorld(ref_vol->GetMRI(), 
                      pos[0], pos[1], pos[2], 
                      pos, pos+1, pos+2);
    }
    else if (label->coords == LABEL_COORDS_TKREG_RAS)
    {
      ref_vol->TkRegToNativeRAS( pos, pos );
    }
    ref_vol->NativeRASToRAS( pos, pos );
    m_layerSource->RASToTarget(pos, pos);
  }
}

void LayerFCD::SetLabelVisible(int nIndex, bool visible)
{
  LABEL* label = m_fcd->labels[nIndex];
  FSVolume* ref_vol = m_layerSource->GetSourceVolume();
  int n[3];
  double pos[3];
  int* dim = m_imageData->GetDimensions();
  for ( int i = 0; i < label->n_points; i++ )
  {
    pos[0] = label->lv[i].x;
    pos[1] = label->lv[i].y;
    pos[2] = label->lv[i].z;
    if ( label->coords == LABEL_COORDS_VOXEL )
    {
      MRIvoxelToWorld(ref_vol->GetMRI(), 
                      pos[0], pos[1], pos[2], 
                      pos, pos+1, pos+2);
    }
    else if (label->coords == LABEL_COORDS_TKREG_RAS)
    {
      ref_vol->TkRegToNativeRAS( pos, pos );
    }
    ref_vol->NativeRASToRAS( pos, pos );
    ref_vol->RASToTargetIndex( pos, n );
    if ( n[0] >= 0 && n[0] < dim[0] && n[1] >= 0 && n[1] < dim[1] &&
         n[2] >= 0 && n[2] < dim[2] )
    {
      m_imageData->SetScalarComponentFromFloat
        ( n[0], n[1], n[2], 0, visible?label->lv[i].vno:0 );
    }
  }
  for (int i = 0; i < 3; i++)
  {
    mReslice[i]->Modified();
  }
  emit ActorUpdated();
}

void LayerFCD::SetMRILayerCTAB(COLOR_TABLE *ctab)
{
  QList<LayerMRI*> layers = GetMRILayers();
  foreach (LayerMRI* mri, layers)
  mri->GetProperty()->SetLUTCTAB(ctab);
}

QList<LayerMRI*> LayerFCD::GetMRILayers()
{
  QList<LayerMRI*> layers;
  if (m_mri_norm)
  {
    layers << m_mri_norm;
  }
  if (m_mri_flair)
  {
    layers << m_mri_flair;
  }
  if (m_mri_t2)
  {
    layers << m_mri_t2;
  }
  if (m_mri_aseg)
  {
    layers << m_mri_aseg;
  }
  if (m_mri_difference)
  {
    layers << m_mri_difference;
  }

  return layers;
}

QList<LayerSurface*> LayerFCD::GetSurfaceLayers()
{
  QList<LayerSurface*> layers;
  if (m_surf_lh)
  {
    layers << m_surf_lh;
  }
  if (m_surf_lh_pial)
  {
    layers << m_surf_lh_pial;
  }
  if (m_surf_lh_sphere_d1)
  {
    layers << m_surf_lh_sphere_d1;
  }
  if (m_surf_rh)
  {
    layers << m_surf_rh;
  }
  if (m_surf_rh_pial)
  {
    layers << m_surf_rh_pial;
  }
  if (m_surf_rh_sphere_d1)
  {
    layers << m_surf_rh_sphere_d1;
  }

  return layers;
}

QThread* LayerFCD::GetWorkerThread()
{
  return m_worker;
}

void LayerFCD::OnLayerDestroyed()
{
  Layer* layer = qobject_cast<Layer*>(sender());
  if (layer == m_surf_lh)
  {
    m_surf_lh = NULL;
  }
  else if (layer == m_surf_rh)
  {
    m_surf_rh = NULL;
  }
  else if (layer == m_surf_lh_pial)
  {
    m_surf_lh_pial = NULL;
  }
  else if (layer == m_surf_rh_pial)
  {
    m_surf_rh_pial = NULL;
  }
  else if (layer == m_surf_lh_sphere_d1)
  {
    m_surf_lh_sphere_d1 = NULL;
  }
  else if (layer == m_surf_rh_sphere_d1)
  {
    m_surf_rh_sphere_d1 = NULL;
  }
  else if (layer == m_mri_norm)
  {
    m_mri_norm = NULL;
  }
  else if (layer == m_mri_flair)
  {
    m_mri_flair = NULL;
  }
  else if (layer == m_mri_t2)
  {
    m_mri_t2 = NULL;
  }
  else if (layer == m_mri_aseg)
  {
    m_mri_aseg = NULL;
  }
  else if (layer == m_mri_difference)
  {
    m_mri_difference = NULL;
  }
}
