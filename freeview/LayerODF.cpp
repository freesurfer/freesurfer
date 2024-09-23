#include "LayerODF.h"
#include "LayerPropertyODF.h"
#include "LayerMRI.h"
#include "mri.h"
#include <QFileInfo>
#include <QDebug>
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"
#include <vtkPointData.h>
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkDecimatePro.h"
#include <QFile>
#include "MyVTKUtils.h"
#include "vtkColorTransferFunction.h"
#include <QSettings>
#include "FSVolume.h"
#include <QRegularExpression>
#include "MigrationDefs.h"

#define SUB_ACTOR_COUNT 4

LayerODF::LayerODF( LayerMRI* ref, QObject* parent, int nMainView ) : LayerMRI( ref, parent ),
  m_mask(NULL), m_nMainView(nMainView)
{
  m_strTypeNames.push_back( "ODF" );
  m_sPrimaryType = "ODF";
  mProperty = new LayerPropertyODF(this);

  connect(mProperty, SIGNAL(OdfPropertyChanged()), this, SIGNAL(UpdateActorRequested()));
  connect(this, SIGNAL(UpdateActorRequested(int)), SLOT(UpdateActors(int)), Qt::QueuedConnection);
  connect(mProperty, SIGNAL(ColorCodeChanged()), SLOT(OnColorCodeChanged()), Qt::QueuedConnection);
  connect(mProperty, SIGNAL(ShowInAllChanged()), SLOT(OnShowInAllChanged()), Qt::QueuedConnection);

  // initialize mesh data
  QFile file(":/resource/DSI_vectors_181.dat");
  file.open(QIODevice::ReadOnly);
  m_nVectors = 181;
  m_nMesh = 720;
  for (int i = 0; i < m_nVectors; i++)
  {
    file.read((char*)m_odfVector[i], sizeof(float)*3);
  }
  file.close();
  file.setFileName(":/resource/DSI_mesh_181.dat");
  file.open(QIODevice::ReadOnly);
  for (int i = 0; i < m_nMesh; i++)
  {
    file.read((char*)m_odfMesh[i], sizeof(int)*3);
  }
  file.close();

  m_actor3D = vtkSmartPointer<vtkActor>::New();
  for (int i = 0; i < 3; i++)
    m_glyphActor2D[i]->SetVisibility(i == m_nMainView);

#ifdef USE_ACTOR_LIST
  for (int i = 0; i < SUB_ACTOR_COUNT; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      vtkActor* actor = vtkActor::New();
      actor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
      m_listActor2D[j] << actor;
      actor = vtkActor::New();
      actor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
      m_listActor3D[j] << actor;
    }
  }
#endif

  SetOdfMask(ref, false);
}

LayerODF::~LayerODF()
{
  ClearSliceActors();

  QSettings s;
  QVariantMap root = s.value("OdfSettings").toMap();
  QVariantMap map = root[GetFileName()].toMap();
  if (!property("mask_fn").toString().isEmpty())
  {
    map[property("mask_fn").toString()] = m_odfMaskThreshold[0];
  }
  map["skip"] = GetProperty()->GetOdfSkip();
  map["scale"] = GetProperty()->GetOdfScale();
  map["color_code"] = GetProperty()->GetOdfColorCode();
  map["inversion"] = GetProperty()->GetOdfInversion();
  double th[2];
  GetProperty()->GetMagnitudeThreshold(th);
  map["color_min"] = th[0];
  map["color_max"] = th[1];
  root[GetFileName()] = map;
  s.setValue("OdfSettings",root);
  s.sync();
}

void LayerODF::ClearSliceActors(int nPlane)
{
#ifdef USE_ACTOR_LIST
  int nStart = 0, nEnd = 2;
  if (nPlane >= 0)
    nStart = nEnd = nPlane;
  for (int i = nStart; i <= nEnd; i++)
  {
    foreach (vtkActor* actor, m_listActor2D[i])
    {
      actor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    }
    foreach (vtkActor* actor, m_listActor3D[i])
    {
      actor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    }
  }
#endif
}

void LayerODF::SetOdfMask(LayerMRI *mri, bool bRefresh)
{
  if (m_mask)
    m_mask->setProperty("odf_mask_threshold", m_odfMaskThreshold[0]);
  m_mask = mri;
  if (mri && mri->property("odf_mask_threshold").isNull())
  {
    m_odfMaskThreshold[0] = mri->GetProperty()->GetMinGenericThreshold();
    m_odfMaskThreshold[1] = mri->GetProperty()->GetMaxValue();
  }

  setProperty("mask_fn", m_mask?m_mask->GetFileName():"");

  if (bRefresh)
    emit UpdateActorRequested();
}

void LayerODF::GetOdfMaskThreshold(double *th)
{
  th[0] = m_odfMaskThreshold[0];
  th[1] = m_odfMaskThreshold[1];
}

void LayerODF::SetOdfMaskThreshold(double *th)
{
  m_odfMaskThreshold[0] = th[0];
  m_odfMaskThreshold[1] = th[1];
  emit UpdateActorRequested();
}

bool LayerODF::Load(const QString &fn, const QString& vertex_fn, const QString& face_fn, bool bPermute, bool bHemisphere)
{
  m_bDtkFormat = vertex_fn.isEmpty();
  m_bHemisphere = (m_bDtkFormat || bHemisphere);
  if (!m_bDtkFormat)
  {
    QFile file(vertex_fn);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      cerr << "Could not load vertex file" << endl;
      return false;
    }

    int n = 0;
    while (!file.atEnd())
    {
      QString line = file.readLine().trimmed();
      if (!line.isEmpty())
      {
        QStringList list = line.split(QRegularExpression("\\s+"), MD_SkipEmptyParts);
        if (list.size() == 3)
        {
          m_odfVector[n][0] = list[0].toDouble();
          m_odfVector[n][1] = list[1].toDouble();
          m_odfVector[n][2] = list[2].toDouble();
          n++;
        }
      }
    }
    m_nVectors = n;
    file.close();

    file.setFileName(face_fn);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      cerr << "Could not load face file" << endl;
      return false;
    }

    n = 0;
    while (!file.atEnd())
    {
      QString line = file.readLine().trimmed();
      if (!line.isEmpty())
      {
        QStringList list = line.split(QRegularExpression("\\s+"), MD_SkipEmptyParts);
        if (list.size() == 3)
        {
          m_odfMesh[n][0] = (int)list[0].toDouble();
          m_odfMesh[n][1] = (int)list[1].toDouble();
          m_odfMesh[n][2] = (int)list[2].toDouble();
          n++;
        }
      }
    }
    m_nMesh = n;
    file.close();
  }

  MRI* mri = ::MRIread(qPrintable(fn));
  if (!mri)
  {
    cerr << "Failed to read file " << qPrintable(fn) << endl;
    return false;
  }
  SetName(QFileInfo(fn).completeBaseName());
  SetFileName(fn);

  // un-permute mri if dtk format
  if (m_bDtkFormat && bPermute)
  {
    MRI* mri2 = ::MRIallocSequence(mri->height, mri->depth, mri->nframes, mri->type, mri->width);
    MRI* mri_ref = m_volumeRef->GetMRI();
    ::MRIcopyHeader(mri_ref, mri2);
    for (int i = 0; i < mri->width; i++)
      for (int j = 0; j < mri->height; j++)
        for (int k = 0; k < mri->depth; k++)
          for (int n = 0; n < mri->nframes; n++)
          {
            switch ( mri->type )
            {
            case MRI_UCHAR:
            MRIseq_vox( mri2, j, k, n, i) = MRIseq_vox( mri, i, j, k, n);
            break;
            case MRI_INT:
            MRIIseq_vox( mri2, j, k, n, i ) = MRIIseq_vox( mri, i, j, k, n );
            break;
            case MRI_LONG:
            MRILseq_vox( mri2, j, k, n, i ) = MRILseq_vox( mri, i, j, k, n  );
            break;
            case MRI_FLOAT:
            MRIFseq_vox( mri2, j, k, n, i ) = MRIFseq_vox( mri, i, j, k, n  );
            break;
            case MRI_SHORT:
            MRISseq_vox( mri2, j, k, n, i ) = MRISseq_vox( mri, i, j, k, n );
            break;
            case MRI_USHRT:
            MRIUSseq_vox( mri2, j, k, n, i ) = MRIUSseq_vox( mri, i, j, k, n );
            break;
            default:
            break;
            }
          }
    MRI* temp = mri;
    mri = mri2;
    ::MRIfree(&temp);
  }

  if ( CreateFromMRIData((void*)mri) )
  {
    m_imageData->GetScalarRange(m_dScalarRange);
    GetProperty()->SetMagnitudeThreshold(m_dScalarRange);

    QSettings s;
    QVariantMap root = s.value("OdfSettings").toMap();
    if (root.contains(fn))
    {
      QVariantMap map = root[fn].toMap();
      GetProperty()->SetOdfSkip(map["skip"].toInt());
      if (map.contains("color_min"))
      {
        double th[2];
        th[0] = map["color_min"].toDouble();
        th[1] = map["color_max"].toDouble();
        GetProperty()->SetMagnitudeThreshold(th);
      }
      if (m_mask && map.contains(m_mask->GetFileName()))
      {
        m_odfMaskThreshold[0] = map[m_mask->GetFileName()].toDouble();
      }
      if (map.contains("scale"))
        GetProperty()->SetOdfScale(map["scale"].toDouble());
      if (map.contains("color_code"))
        GetProperty()->SetOdfColorCode(map["color_code"].toInt());
      if (map.contains("inversion"))
        GetProperty()->SetOdfInversion(map["inversion"].toInt());
    }

    ::MRIfree( &mri );
    return true;
  }
  else
  {
    cerr << "Failed to create " << endl;
    ::MRIfree( &mri );
    return false;
  }
}

void LayerODF::OnSlicePositionChanged(int nPlane)
{
  if (!m_imageData)
    return;

  BuildSlice(nPlane);
}

void LayerODF::BuildSlice(int nPlane)
{
  vtkImageData* imagedata = m_imageData;
  double* pos = GetSlicePosition();
  double* orig = imagedata->GetOrigin();
  double* voxel_size = imagedata->GetSpacing();
  int* dim = imagedata->GetDimensions();
  int slice[3] = {0}, slice_range[3][2];
  for (int i = 0; i < 3; i++)
  {
    slice[i] = (int)((m_dSlicePosition[i]-orig[i])/voxel_size[i]);
    slice_range[i][0] = 0;
    slice_range[i][1] = dim[i]-1;
  }
  ClearSliceActors(nPlane);
  if (slice[nPlane] < 0 || slice[nPlane] >= dim[nPlane] ||
      (!GetProperty()->GetShowInAllViews() && m_nMainView != nPlane && m_nMainView < 3))
  {
#ifndef USE_ACTOR_LIST
    m_glyphActor2D[nPlane]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    m_glyphActor3D[nPlane]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
#endif
    return;
  }
  slice_range[nPlane][0] = slice_range[nPlane][1] = slice[nPlane];

  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkUnsignedCharArray> scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  vtkSmartPointer<vtkFloatArray> scalars_2 = vtkSmartPointer<vtkFloatArray>::New();
  scalars->SetNumberOfComponents( 4 );

  int nskip = GetProperty()->GetOdfSkip()+1;
  int x, y, z;
  double min_size = qMin(voxel_size[0], qMin(voxel_size[1], voxel_size[2]));
  double odf_range[2] = {m_dScalarRange[0], m_dScalarRange[1]};
  double fOffset = 0;
  double l_scale = GetProperty()->GetOdfScale()*2;
  double sharpness = 1;
  bool bNormalize = false;
  if (odf_range[0] < 0)
  {
    fOffset = -odf_range[0];
    odf_range[1] += fOffset;
    odf_range[0] = 0;
  }
  fOffset *= 1.2;
  double over_all = l_scale*min_size/odf_range[1]*3/8;

  char* ptr = (char*)imagedata->GetScalarPointer();
  int scalar_type = imagedata->GetScalarType();
  char* mask_ptr = NULL;
  int mask_scalar_type = 0;
  vtkImageData* maskdata = NULL;
  if (m_mask)
  {
    maskdata = m_mask->GetImageData();
    mask_ptr = (char*)maskdata->GetScalarPointer();
    mask_scalar_type = maskdata->GetScalarType();
  }
  int nFrames = GetNumberOfFrames();
  double raws[4096];
  bool bInvert[3] = {false, false, false};
  if (GetProperty()->GetOdfInversion() > 0)
    bInvert[GetProperty()->GetOdfInversion()-1] = true;
  unsigned char c[4] = {0,0,0,255};
  int vox_cnt = 0, total_vox = 0, cnt = 0;
#ifdef USE_ACTOR_LIST
  for (x = slice_range[0][0]; x <= slice_range[0][1]; x+=nskip)
  {
    for (y = slice_range[1][0]; y <= slice_range[1][1]; y+=nskip)
    {
      for (z = slice_range[2][0]; z <= slice_range[2][1]; z+=nskip)
      {
        if (maskdata)
        {
          double val = MyVTKUtils::GetImageDataComponent(mask_ptr, dim, 1, x, y, z, 0, mask_scalar_type);
          if (val < m_odfMaskThreshold[0] || val > m_odfMaskThreshold[1])
            continue;
        }
        total_vox++;
      }
    }
  }

  int nSubMax = total_vox/SUB_ACTOR_COUNT+1;
  if (nSubMax < 100)
    nSubMax = total_vox+1;
#endif

  int nSubId = 0;
  for (x = slice_range[0][0]; x <= slice_range[0][1]; x+=nskip)
  {
    for (y = slice_range[1][0]; y <= slice_range[1][1]; y+=nskip)
    {
      for (z = slice_range[2][0]; z <= slice_range[2][1]; z+=nskip)
      {
        if (maskdata)
        {
          double val = MyVTKUtils::GetImageDataComponent(mask_ptr, dim, 1, x, y, z, 0, mask_scalar_type);
          if (val < m_odfMaskThreshold[0] || val > m_odfMaskThreshold[1])
            continue;
        }

        for (int i = 0; i < nFrames; i++)
          raws[i] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, x, y, z, i, scalar_type);

        double fMax = qAbs(raws[0]+fOffset), fMin = (raws[0]+fOffset);
        double new_offset = fOffset;
        if (bNormalize || sharpness > 1)
        {
          for (int k = 1; k < nFrames; k++)
          {
            if (fMax < (raws[k]+fOffset))
              fMax = (raws[k]+fOffset);
            if (fMin > (raws[k]+fOffset))
              fMin = (raws[k]+fOffset);
          }
          if (fMin < 0)
          {
            new_offset -= fMin;
            fMax -= fMin;
          }
          if (bNormalize)
            over_all = l_scale*min_size/fMax*7/16;
        }

        for (int k = 0; k < nFrames; k++)
        {
          float scale = (raws[k]+new_offset)*over_all;
          double ov[3];
          for (int nv = 0; nv < 3; nv++)
            ov[nv] = bInvert[nv]?(-m_odfVector[k][nv]):m_odfVector[k][nv];
          if (sharpness == -1)
          {
            scale = over_all*fMax/qMax(qMax(qAbs(ov[0]), qAbs(ov[1])), qAbs(ov[2]));
          }
          else if (sharpness == 0)
          {
            scale = over_all*fMax;
          }
          else if (sharpness >= 2)
          {
            for (int i = 0; i < sharpness-1; i++)
              scale = scale*(raws[k]+new_offset)/fMax;
          }
          pts->InsertNextPoint(orig[0] + voxel_size[0]*x + ov[0]*scale,
                               orig[1] + voxel_size[1]*y + ov[1]*scale,
                               orig[2] + voxel_size[2]*z + ov[2]*scale);

          c[0] = (int)(fabs( ov[0] *255 ) );
          c[1] = (int)(fabs( ov[1] *255 ) );
          c[2] = (int)(fabs( ov[2] *255 ) );
#if VTK_MAJOR_VERSION > 5
          scalars->InsertNextTypedTuple( c );
#else
          scalars->InsertNextTupleValue( c );
#endif
          scalars_2->InsertNextValue(raws[k]);
        }

        if (m_bHemisphere)
        {
          for (int k = 0; k < nFrames; k++)
          {
            float scale = (raws[k]+new_offset)*over_all;
            double ov[3];
            for (int nv = 0; nv < 3; nv++)
              ov[nv] = bInvert[nv]?(-m_odfVector[k][nv]):m_odfVector[k][nv];
            if (sharpness == -1)
            {
              scale = over_all*fMax/qMax(qMax(qAbs(ov[0]), qAbs(ov[1])), qAbs(ov[2]));
            }
            else if (sharpness == 0)
            {
              scale = over_all*fMax;
            }
            else if (sharpness == 2)
            {
              scale = scale*(raws[k]+new_offset)/fMax;
            }
            pts->InsertNextPoint(orig[0] + voxel_size[0]*x - ov[0]*scale,
                                 orig[1] + voxel_size[1]*y - ov[1]*scale,
                                 orig[2] + voxel_size[2]*z - ov[2]*scale);

            c[0] = (int)(fabs( ov[0] *255 ) );
            c[1] = (int)(fabs( ov[1] *255 ) );
            c[2] = (int)(fabs( ov[2] *255 ) );
#if VTK_MAJOR_VERSION > 5
            scalars->InsertNextTypedTuple( c );
#else
            scalars->InsertNextTupleValue( c );
#endif
            scalars_2->InsertNextValue(raws[k]);
          }
        }

        if (m_bDtkFormat)
        {
          for (int k = 0; k < m_nMesh; k++)
          {
            polys->InsertNextCell(3);
            for (int kk = 0; kk < 3; kk++)
            {
              if (m_odfMesh[k][kk] > 0)
                polys->InsertCellPoint(cnt+m_odfMesh[k][kk]-1);
              else
                polys->InsertCellPoint(cnt+nFrames-m_odfMesh[k][kk]-1);
            }
          }
          cnt += nFrames*2;
        }
        else
        {
          for (int k = 0; k < m_nMesh; k++)
          {
            polys->InsertNextCell(3);
            for (int kk = 0; kk < 3; kk++)
            {
              polys->InsertCellPoint(cnt+m_odfMesh[k][kk]);
            }
          }
          cnt += nFrames;
          if (m_bHemisphere)
          {
            cnt += nFrames;
          }
        }

        vox_cnt++;
#ifdef USE_ACTOR_LIST
        if (vox_cnt >= nSubMax && nSubId < SUB_ACTOR_COUNT-1)
        {
          BuildActor(m_listActor2D[nPlane][nSubId], m_listActor3D[nPlane][nSubId], pts, polys, scalars, scalars_2);
          cnt = 0;
          vox_cnt = 0;
          nSubId++;
          pts = vtkSmartPointer<vtkPoints>::New();
          polys = vtkSmartPointer<vtkCellArray>::New();
          scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
          scalars_2 = vtkSmartPointer<vtkFloatArray>::New();
          scalars->SetNumberOfComponents( 4 );
        }
#endif
      }
    }
  }
  if (vox_cnt > 0)
  {
#ifdef USE_ACTOR_LIST
    BuildActor(m_listActor2D[nPlane][nSubId], m_listActor3D[nPlane][nSubId], pts, polys, scalars, scalars_2);
#else
    m_glyphActor2D[nPlane]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    m_glyphActor3D[nPlane]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    m_polydata[nPlane] = BuildActor(m_glyphActor2D[nPlane], m_glyphActor3D[nPlane], pts, polys, scalars, scalars_2);
#endif
  }

  double actor_pos[3] = {0,0,0};
  actor_pos[nPlane] = voxel_size[nPlane]*(nPlane==2?-dim[nPlane]:dim[nPlane])/2;
  m_glyphActor2D[nPlane]->SetPosition(actor_pos);

#ifdef USE_ACTOR_LIST
  foreach (vtkActor* actor, m_listActor2D[nPlane])
    actor->SetPosition(actor_pos);
#endif
}

vtkPolyData* LayerODF::BuildActor(vtkActor* actor2D, vtkActor* actor3D, vtkPoints* pts, vtkCellArray* polys, vtkUnsignedCharArray* scalars, vtkFloatArray* scalars_2)
{
  double fdeci = 1;
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(pts);
  polydata->SetPolys(polys);
  scalars->SetName("Directional");
  polydata->GetPointData()->SetScalars(scalars);
  scalars_2->SetName("Magnitude");
  polydata->GetPointData()->AddArray(scalars_2);
  polydata->GetPointData()->SetActiveScalars(qPrintable(GetProperty()->GetOdfColorCode()==0?"Directional":"Magnitude"));

  vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetFeatureAngle(90);
  if (fdeci < 1)
  {
    vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
    deci->SetTargetReduction(fdeci);
#if VTK_MAJOR_VERSION > 5
    deci->SetInputData(polydata);
#else
    deci->SetInput(polydata);
#endif
    normals->SetInputConnection(deci->GetOutputPort());
  }
  else
#if VTK_MAJOR_VERSION > 5
    normals->SetInputData(polydata);
#else
    normals->SetInput(polydata);
#endif

  /*	vtkSmoothPolyDataFilter* smoother = vtkSmoothPolyDataFilter::New();
    smoother->SetInput(normals->GetOutput());
    smoother->SetFeatureAngle(90);
    smoother->SetNumberOfIterations(10);*/

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkPolyDataMapper::SafeDownCast(actor2D->GetMapper());
  mapper->SetInputConnection(normals->GetOutputPort());
  mapper->SetScalarVisibility(1);  
  mapper->SetLookupTable(GetProperty()->GetOdfLut());
  mapper = vtkPolyDataMapper::SafeDownCast( actor3D->GetMapper() );
  mapper->SetInputConnection(normals->GetOutputPort());
  mapper->SetScalarVisibility(1);
  mapper->SetLookupTable(GetProperty()->GetOdfLut());

  return polydata;
}

void LayerODF::Append2DProps( vtkRenderer* renderer, int nPlane )
{
#ifdef USE_ACTOR_LIST
  foreach (vtkActor* actor, m_listActor2D[nPlane])
    renderer->AddViewProp(actor);
#else
  renderer->AddViewProp( m_glyphActor2D[nPlane] );
#endif
}

void LayerODF::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  for (int i = 0; i < 3; i++)
  {
#ifdef USE_ACTOR_LIST
    foreach (vtkActor* actor, m_listActor3D[i])
      renderer->AddViewProp(actor);
#else
    renderer->AddViewProp( m_glyphActor3D[i] );
#endif
  }
}

void LayerODF::SetVisible(bool bVisible)
{
  for (int i = 0; i < 3; i++)
  {
#ifdef USE_ACTOR_LIST
    foreach (vtkActor* actor, m_listActor2D[i])
      actor->SetVisibility((bVisible && GetProperty()->GetShowIn2DView())?1:0);
    foreach (vtkActor* actor, m_listActor3D[i])
      actor->SetVisibility(bVisible?1:0);
#else
    m_glyphActor2D[i]->SetVisibility((bVisible && (GetProperty()->GetShowInAllViews() || m_nMainView == i))?1:0);
    m_glyphActor3D[i]->SetVisibility((bVisible && (GetProperty()->GetShowInAllViews() || m_nMainView == 3))?1:0);
#endif
  }
  m_actor3D->SetVisibility((bVisible && (GetProperty()->GetShowInAllViews() || m_nMainView == 3))?1:0);
  emit ActorUpdated();
}

bool LayerODF::IsVisible()
{
  if (m_nMainView < 3)
    return m_glyphActor2D[m_nMainView]->GetVisibility();
  else
    return m_actor3D->GetVisibility();
}

void LayerODF::UpdateActors(int nPlane)
{
  this->blockSignals( true );
  if (nPlane >= 0)
    BuildSlice(nPlane);
  else
  {
    for ( int i = 0; i < 3; i++ )
    {
      BuildSlice(i);
    }
  }
  if (nPlane < 0)
    OnColorCodeChanged();

  this->blockSignals( false );

  emit ActorUpdated();
}

void LayerODF::OnColorCodeChanged()
{
  int nCode = GetProperty()->GetOdfColorCode();
  for (int i = 0; i < 3; i++)
  {
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[i]->GetMapper() );
    if (m_polydata[i].GetPointer())
    {
      m_polydata[i]->GetPointData()->SetActiveScalars(qPrintable(nCode==0?"Directional":"Magnitude"));
      mapper->SetLookupTable(GetProperty()->GetOdfLut());
      mapper = vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[i]->GetMapper() );
      mapper->SetLookupTable(GetProperty()->GetOdfLut());
    }
  }
  SetVisible(IsVisible());
}

void LayerODF::OnMainViewChanged(int nView)
{
  bool bVisible = IsVisible();
  m_nMainView = nView;
  if (!GetProperty()->GetShowInAllViews())
  {
    UpdateActors(-1);
    SetVisible(bVisible);
  }
}

void LayerODF::OnShowInAllChanged()
{
  if (GetProperty()->GetShowInAllViews())
  {
    UpdateActors(-1);
  }
  SetVisible(IsVisible());
}
