#include "WindowLayerInfo.h"
#include "ui_WindowLayerInfo.h"
#include <QLayoutItem>
#include <QLabel>
#include "LayerMRI.h"
#include "FSVolume.h"
#include "LayerSurface.h"
#include "FSSurface.h"
#include <QSettings>

WindowLayerInfo::WindowLayerInfo(QWidget *parent) :
  QWidget(parent),
  ui(new Ui::WindowLayerInfo)
{
  ui->setupUi(this);
  setWindowFlags( Qt::Tool );

  QSettings settings;
  restoreGeometry(settings.value("WindowLayerInfo/Geometry").toByteArray());
}

WindowLayerInfo::~WindowLayerInfo()
{
  QSettings settings;
  settings.setValue("WindowLayerInfo/Geometry", this->saveGeometry());
  delete ui;
}

void WindowLayerInfo::Clear()
{
  QLayoutItem* item;
  while ( ( item = ui->gridLayout->takeAt( 0 ) ) != NULL )
  {
    delete item->widget();
    delete item;
  }
}

void WindowLayerInfo::AddLine(const QString &name, const QString &value, bool word_wrap)
{
  QLabel* l = new QLabel(name);
  l->setAlignment(Qt::AlignRight | (word_wrap?Qt::AlignTop:Qt::AlignVCenter));
  int nRow = ui->gridLayout->rowCount();
  ui->gridLayout->addWidget(l, nRow, 0);
  l = new QLabel(value);
  l->setTextInteractionFlags(Qt::TextSelectableByMouse);
  if (word_wrap)
    l->setWordWrap(true);
  ui->gridLayout->addWidget(l, nRow, 1);
}

void WindowLayerInfo::SetCaption(const QString &text)
{
  QLabel* l = new QLabel(text);
  l->setMinimumWidth(250);
  l->setTextInteractionFlags(Qt::TextSelectableByMouse);
  l->setAlignment(Qt::AlignCenter);
  l->setWordWrap(true);
  l->setStyleSheet("padding-bottom:10px");
  ui->gridLayout->addWidget(l, 0, 0, 1, 2);
}

void WindowLayerInfo::UpdateInfo(Layer* layer)
{
  QString type = layer->GetPrimaryType();
  if (type == "MRI")
  {
    LayerMRI* layer_mri = qobject_cast<LayerMRI*>(layer);
    MRI* mri = layer_mri->GetSourceVolume()->GetMRI();
    Clear();
    setWindowTitle("Volume Information");
//    SetCaption(QString("Volume information for %1").arg(layer_mri->GetFileName()));
    AddLine("file name:", layer_mri->GetFileName(), true);
    QString val;
    val = QString("%1 x %2 x %3").arg(mri->width).arg(mri->height).arg(mri->depth);
    if (mri->nframes > 1)
      val += QString(" x %1").arg(mri->nframes);
    AddLine("dimensions:", val);
    AddLine("voxel sizes:", QString("%1, %2, %3")
            .arg(mri->xsize, 0, 'f', 6)
            .arg(mri->ysize, 0, 'f', 6)
            .arg(mri->zsize, 0, 'f', 6));
    AddLine("number of frames:", QString("%1 ").arg(mri->nframes));
    AddLine("type:", QString("%1 (%2)").arg(
                            mri->type == MRI_UCHAR   ? "UCHAR" :
                            mri->type == MRI_SHORT   ? "SHORT" :
                            mri->type == MRI_USHRT   ? "USHORT" :
                            mri->type == MRI_INT     ? "INT" :
                            mri->type == MRI_LONG    ? "LONG" :
                            mri->type == MRI_BITMAP  ? "BITMAP" :
                            mri->type == MRI_TENSOR  ? "TENSOR" :
                            mri->type == MRI_FLOAT   ? "FLOAT" : "UNKNOWN").arg(mri->type));
    AddLine("TR:", QString("%1 msec").arg(mri->tr));
    AddLine("TE:", QString("%1 msec").arg(mri->te));
    AddLine("TI:", QString("%1 msec").arg(mri->ti));
    AddLine("flip angle:", QString("%1 degrees").arg(mri->flip_angle));
  }
  else if (type == "Surface")
  {
    LayerSurface* surf = qobject_cast<LayerSurface*>(layer);
    MRIS* mris = surf->GetSourceSurface()->GetMRIS();
    Clear();
    setWindowTitle("Surface Information");
//    SetCaption(QString("Surface information for %1").arg(surf->GetFileName()));
    AddLine("file name:", surf->GetFileName(), true);
    AddLine("num vertices:", QString::number(mris->nvertices));
    AddLine("num faces:", QString::number(mris->nfaces));
    AddLine("num strips:", QString::number(mris->nstrips));
    AddLine("surface area:", QString::number(mris->total_area));
  }
}
