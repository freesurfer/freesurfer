#include "LayerConnectomeMatrix.h"
#include "LayerMRI.h"
#include "LayerPropertyConnectomeMatrix.h"

LayerConnectomeMatrix::LayerConnectomeMatrix(LayerMRI* ref, QObject *parent) :
  Layer(parent),
  m_mriRef(ref),
  m_mriParcel(NULL),
  m_cmat(NULL),
  m_ctab(NULL)
{
  this->m_strTypeNames << "CMAT";
  mProperty = new LayerPropertyConnectomeMatrix( this );
}

LayerConnectomeMatrix::~LayerConnectomeMatrix()
{
  if (m_cmat)
    ::CMATfree(&m_cmat);
  if (m_ctab)
    ::CTABfree(&m_ctab);
  if (m_mriParcel)
    delete m_mriParcel;
}

bool LayerConnectomeMatrix::LoadFromFile(const QString &fn_cmat, const QString &fn_parcel, const QString &fn_ctab)
{
  if (m_cmat)
    ::CMATfree(&m_cmat);

  m_cmat = ::CMATread(qPrintable(fn_cmat));
  if (!m_cmat)
  {
    cerr << "Could not load CMAT file " << qPrintable(fn_cmat) << endl;
    return false;
  }

  if (m_mriParcel)
    delete m_mriParcel;

  m_mriParcel = new LayerMRI(m_mriRef);
  m_mriParcel->SetFileName(fn_parcel);
  if (!m_mriParcel->LoadVolumeFromFile())
  {
    cerr << "Could not load parcellation file " << qPrintable(fn_parcel) << endl;
    delete m_mriParcel;
    ::CMATfree(&m_cmat);
    m_mriParcel = NULL;
    return false;
  }

  // update label coords to target coords
  double pt[3];
  for (int i = 0; i < m_cmat->nlabels; i++)
  {
    for (int j = 0; j < m_cmat->nlabels; j++)
    {
      LABEL* label = m_cmat->splines[i][j];
      for (int n = 0; n < label->n_points; n++)
      {
        pt[0] = label->lv[n].x;
        pt[1] = label->lv[n].y;
        pt[2] = label->lv[n].z;
        m_mriParcel->RASToTarget(pt, pt);
        label->lv[n].x = pt[0];
        label->lv[n].y = pt[1];
        label->lv[n].z = pt[2];
      }
    }
  }

  return true;
}

bool LayerConnectomeMatrix::LoadFromFile()
{
  return LoadFromFile(m_sFilename, m_sFilenameParcel, m_sFilenameCTAB);
}

void LayerConnectomeMatrix::Append2DProps(vtkRenderer *renderer, int nPlane)
{

}

void LayerConnectomeMatrix::Append3DProps(vtkRenderer *renderer, bool *bPlaneVisibility)
{

}

bool LayerConnectomeMatrix::HasProp(vtkProp *prop)
{

}

bool LayerConnectomeMatrix::IsVisible()
{

}

void LayerConnectomeMatrix::OnSlicePositionChanged(int nPlane)
{

}
