#include "DialogMovePoint.h"
#include "ui_DialogMovePoint.h"
#include "LayerPointSet.h"
#include "RenderView2D.h"
#include "LayerMRI.h"

DialogMovePoint::DialogMovePoint(QWidget *parent) :
  QDialog(parent), m_bUndoable(false), m_layerPointSet(NULL),
  ui(new Ui::DialogMovePoint)
{
  ui->setupUi(this);
  setWindowFlags(Qt::Dialog | Qt::WindowStaysOnTopHint);

  connect(ui->pushButtonGo, SIGNAL(clicked(bool)), SLOT(OnButtonTest()));
  connect(ui->pushButtonRestore, SIGNAL(clicked(bool)), SLOT(OnButtonRestore()));
}

DialogMovePoint::~DialogMovePoint()
{
  delete ui;
}

void DialogMovePoint::SetData(RenderView2D *view, LayerPointSet *wp, int nId)
{
  m_view = view;
  m_layerPointSet = wp;
  m_nIndex = nId;
  wp->GetPoint(nId, m_dPrevPos);
  m_bUndoable = false;
  ui->labelInfo->setText(tr("PointSet: %1   point #%2").arg(wp->GetName()).arg(nId+1));
}

void DialogMovePoint::OnPointSetPicked(LayerPointSet *wp, int n)
{
  SetData(qobject_cast<RenderView2D*>(sender()), wp, n);
}

void DialogMovePoint::OnButtonTest()
{
  if (!m_layerPointSet)
    return;

  double ras[3], v[3];
  LayerPointSet* wp = m_layerPointSet;
  wp->GetNormalAtPoint(m_nIndex, v, m_view->GetViewPlane());
  LayerMRI* mri = wp->GetReferenceVolume();
  wp->GetPoint(m_nIndex, ras);
  mri->LocateLocalMaximumAtRAS(ras, v[0], v[1], v[2], ras, ui->spinBoxSigma->value(), ui->spinBoxSize->value());
  if (!m_bUndoable)
    wp->SaveForUndo();
  wp->UpdatePoint(m_nIndex, ras);
  m_bUndoable = true;
}

void DialogMovePoint::OnButtonRestore()
{
  if (m_bUndoable)
  {
    m_layerPointSet->Undo();
    m_bUndoable = false;
  }
  else
    m_layerPointSet->UpdatePoint(m_nIndex, m_dPrevPos);
}

double DialogMovePoint::GetSigma()
{
  return ui->spinBoxSigma->value();
}

double DialogMovePoint::GetNeighborSize()
{
  return ui->spinBoxSize->value();
}
