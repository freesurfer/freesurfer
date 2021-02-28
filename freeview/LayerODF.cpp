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

#define ODF_NUM     181

LayerODF::LayerODF( LayerMRI* ref, QObject* parent ) : LayerMRI( ref, parent ),
  m_mask(NULL)
{
  m_strTypeNames.push_back( "ODF" );
  m_sPrimaryType = "ODF";
  mProperty = new LayerPropertyODF(this);
  for (int i = 0; i < 3; i++)
    m_bInvert[i] = false;

  // initialize mesh data
  QFile file(":/resource/DSI_vectors_181.dat");
  file.open(QIODevice::ReadOnly);
  for (int i = 0; i < ODF_NUM; i++)
  {
    file.read((char*)m_odfVector[i], sizeof(float)*3);
    //        if (m_bInvert[0])
    //            m_odfVector[i][0] = -m_odfVector[i][0];
    //        else if (m_bInvert[1])
    //            m_odfVector[i][1] = -m_odfVector[i][1];
    //        else if (m_bInvert[2])
    //            m_odfVector[i][2] = -m_odfVector[i][2];
  }
  file.close();
  file.setFileName(":/resource/DSI_mesh_181.dat");
  file.open(QIODevice::ReadOnly);
  for (int i = 0; i < 720; i++)
  {
    file.read((char*)m_odfMesh[i], sizeof(int)*3);
  }
  file.close();

  m_actor3D = vtkSmartPointer<vtkActor>::New();
}

LayerODF::~LayerODF()
{

}

void LayerODF::SetMask(LayerMRI *mri)
{
  m_mask = mri;
  UpdateActors();
}

bool LayerODF::Load(const QString &fn)
{
  MRI* mri = ::MRIread(qPrintable(fn));
  SetName(QFileInfo(fn).completeBaseName());
  SetFileName(fn);
  if ( CreateFromMRIData((void*)mri) )
  {
    m_imageData->GetScalarRange(m_dScalarRange);
    //    qDebug() << m_dScalarRange[0] << m_dScalarRange[1];
    return true;
  }
  else
  {
    cerr << "Failed to create norm layer" << endl;
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
    qDebug() << "slice " << i << slice[i];
  }
  if (slice[nPlane] < 0 || slice[nPlane] >= dim[nPlane])
  {
    m_glyphActor2D[nPlane]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    m_glyphActor3D[nPlane]->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
    return;
  }
  slice_range[nPlane][0] = slice_range[nPlane][1] = slice[nPlane];

  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIntArray> scalars = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkFloatArray> scalars_2 = vtkSmartPointer<vtkFloatArray>::New();

  int nskip = 1;
  int x, y, z;
  int cnt = 0;
  double min_size = qMin(voxel_size[0], qMin(voxel_size[1], voxel_size[2]));
  double odf_range[2] = {m_dScalarRange[0], m_dScalarRange[1]};
  double fOffset = 0;
  double l_scale = 1;
  double sharpness = 2;
  bool bNormalize = false;
  double fdeci = 1;
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
  int nFrames = 181;
  double raws[181];
  for (x = slice_range[0][0]; x <= slice_range[0][1]; x+=nskip)
  {
    for (y = slice_range[1][0]; y < slice_range[1][1]; y+=nskip)
    {
      for (z = slice_range[2][0]; z < slice_range[2][1]; z+=nskip)
      {
        for (int i = 0; i < nFrames; i++)
          raws[i] = MyVTKUtils::GetImageDataComponent(ptr, dim, nFrames, x, y, z, i, scalar_type);

        bool bGood = true;
        if (bGood)
        {
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
              ov[nv] = m_bInvert[nv]?(-m_odfVector[k][nv]):m_odfVector[k][nv];
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
            pts->InsertNextPoint(voxel_size[0]*(x + 0.5) + ov[0]*scale,
                voxel_size[1]*(y + 0.5) + ov[1]*scale,
                voxel_size[2]*(z + 0.5) + ov[2]*scale);
            scalars->InsertNextValue(k);
            scalars_2->InsertNextValue(scale);
          }
          for (int k = 0; k < nFrames; k++)
          {
            float scale = (raws[k]+new_offset)*over_all;
            double ov[3];
            for (int nv = 0; nv < 3; nv++)
              ov[nv] = m_bInvert[nv]?(-m_odfVector[k][nv]):m_odfVector[k][nv];
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
            pts->InsertNextPoint(voxel_size[0]*(x + 0.5) - ov[0]*scale,
                voxel_size[1]*(y + 0.5) - ov[1]*scale,
                voxel_size[2]*(z + 0.5) - ov[2]*scale);
            scalars->InsertNextValue(k);
            scalars_2->InsertNextValue(scale);
          }
          for (int k = 0; k < 720; k++)
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
      }
    }
  }
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(pts);
  polydata->SetPolys(polys);
  scalars->SetName("Index");
  polydata->GetPointData()->SetScalars(scalars);
  scalars_2->SetName("Magnitude");
  polydata->GetPointData()->AddArray(scalars_2);

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

  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkPolyDataMapper::SafeDownCast( m_glyphActor2D[nPlane]->GetMapper() );
  mapper->SetInputConnection(normals->GetOutputPort());
  mapper->UseLookupTableScalarRangeOn();
  mapper = vtkPolyDataMapper::SafeDownCast( m_glyphActor3D[nPlane]->GetMapper() );
  mapper->SetInputConnection(normals->GetOutputPort());
  mapper->UseLookupTableScalarRangeOn();
  /*
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
*/

}

void LayerODF::Append2DProps( vtkRenderer* renderer, int nPlane )
{
  renderer->AddViewProp( m_glyphActor2D[nPlane] );
}

void LayerODF::Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility )
{
  for (int i = 0; i < 3; i++)
    renderer->AddViewProp( m_glyphActor3D[i] );
}

void LayerODF::SetVisible(bool bVisible)
{
  for (int i = 0; i < 3; i++)
  {
    m_glyphActor2D[i]->SetVisibility(bVisible?1:0);
    m_glyphActor3D[i]->SetVisibility(bVisible?1:0);
  }
  m_actor3D->SetVisibility(bVisible?1:0);
}

bool LayerODF::IsVisible()
{
  return m_actor3D->GetVisibility();
}

void LayerODF::UpdateActors()
{
  this->blockSignals( true );
  for ( int i = 0; i < 3; i++ )
  {
    BuildSlice(i);
  }
  this->blockSignals( false );

  emit ActorChanged();
}
