#include "WindowEditAnnotation.h"
#include "ui_WindowEditAnnotation.h"
#include "LayerSurface.h"
#include "SurfaceAnnotation.h"
#include "LayerPropertySurface.h"
#include "MainWindow.h"
#include "RenderView.h"
#include <QTreeWidgetItem>
#include <QSettings>

#define ANNOT_COLOR_ROLE  Qt::UserRole
#define ANNOT_INDEX_ROLE  Qt::UserRole+1
#define ANNOT_NAME_ROLE  Qt::UserRole+2

WindowEditAnnotation::WindowEditAnnotation(QWidget *parent) :
  QWidget(parent), m_layerSurface(NULL), m_ctab(NULL),
  ui(new Ui::WindowEditAnnotation)
{
  ui->setupUi(this);
  setWindowFlags( Qt::Tool );
  connect(ui->treeWidgetExistingLabels, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)), SLOT(OnExistingLabelClicked(QTreeWidgetItem*)));
  connect(ui->treeWidgetExistingLabels, SIGNAL(itemChanged(QTreeWidgetItem*,int)), SLOT(OnExistingLabelItemChanged(QTreeWidgetItem*)));
  connect(ui->checkBoxShowAll, SIGNAL(stateChanged(int)), SLOT(OnCheckBoxShowAllLabels(int)));
  connect(ui->lineEditName, SIGNAL(returnPressed()), SLOT(OnEditNameReturnPressed()));
  connect(ui->colorpicker, SIGNAL(colorChanged(QColor)), SLOT(OnColorChanged(QColor)));
  connect(ui->treeWidgetAvailableLabels, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)), SLOT(OnAvailableLabelClicked(QTreeWidgetItem*)));
  connect(ui->pushButtonSet, SIGNAL(clicked(bool)), SLOT(OnButtonSet()));

  QSettings settings;
  QVariant v = settings.value("WindowEditAnnotation/Geometry");
  if (v.isValid())
  {
    this->restoreGeometry(v.toByteArray());
  }
}

WindowEditAnnotation::~WindowEditAnnotation()
{
  QSettings settings;
  settings.setValue("WindowEditAnnotation/Geometry", this->saveGeometry());

  delete ui;
}

void WindowEditAnnotation::showEvent(QShowEvent *e)
{
  MainWindow* wnd = MainWindow::GetMainWindow();
  connect(wnd->GetRenderView(3), SIGNAL(SurfaceVertexClicked(LayerSurface*)),
          SLOT(OnSurfaceVertexClicked(LayerSurface*)), Qt::UniqueConnection);
  UpdateUI();
}

void WindowEditAnnotation::OnActiveSurfaceChanged(Layer* layer)
{
  if (m_layerSurface)
  {
    disconnect(m_layerSurface, 0, this, 0);
    disconnect(this, 0, m_layerSurface, 0);
  }
  m_layerSurface = qobject_cast<LayerSurface*>(layer);
  if (m_layerSurface)
  {
    connect(m_layerSurface, SIGNAL(ActiveAnnotationChanged(int)), SLOT(UpdateUI()), Qt::QueuedConnection);
    connect(m_layerSurface, SIGNAL(ActiveAnnotationChanged(int)), SLOT(PopulateAvailableColorTable()), Qt::QueuedConnection);
    connect(this, SIGNAL(LabelClicked(int)), m_layerSurface, SLOT(SetHighlightedLabelOnAnnotation(int)), Qt::QueuedConnection);
    UpdateUI();
  }
}

void WindowEditAnnotation::UpdateUI()
{
  if (m_layerSurface)
  {
    QList<QWidget*> list = this->findChildren<QWidget*>();
    foreach (QWidget* w, list)
      w->blockSignals(true);

    PopulateColorTable();
    ui->pushButtonSet->setEnabled(ui->treeWidgetAvailableLabels->currentItem() != NULL);

    foreach (QWidget* w, list)
      w->blockSignals(false);
  }
}

void WindowEditAnnotation::PopulateAvailableColorTable()
{
  if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
    return;

  ui->treeWidgetAvailableLabels->blockSignals(true);
  SurfaceAnnotation* annot = m_layerSurface->GetActiveAnnotation();
  COLOR_TABLE* ct = annot->GetColorTable();
  if (ct && ct != m_ctab)
  {
    m_ctab = ct;
    ui->treeWidgetAvailableLabels->clear();
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( ct, &nTotalCount );
    int nValid = 0;
    char name[1000];
    for ( int i = 0; i < nTotalCount; i++ )
    {
      CTABisEntryValid( ct, i, &nValid );
      if ( nValid )
      {
        CTABcopyName( ct, i, name, 1000 );
        int nr, ng, nb;
        CTABrgbAtIndexi( ct, i, &nr, &ng, &nb );
        QTreeWidgetItem* item = new QTreeWidgetItem(ui->treeWidgetAvailableLabels);
        UpdateLabelItem(item, i, name,  QColor( nr, ng, nb ));
      }
    }
  }
  ui->treeWidgetAvailableLabels->blockSignals(false);
}

void WindowEditAnnotation::UpdateLabelItem(QTreeWidgetItem* item, int i,
                                           const QString& name, const QColor& color)
{
  if (i < UNASSIGNED_ANNOT_BASE)
    item->setText(0, QString("%1 %2").arg(i).arg(name));
  else
    item->setText(0, name);
  item->setToolTip( 0, name );
  QPixmap pix(13, 13);
  pix.fill( color );
  item->setIcon(0, QIcon(pix) );
  item->setData(0, ANNOT_COLOR_ROLE, color);
  item->setData(0, ANNOT_INDEX_ROLE, i);
  item->setData(0, ANNOT_NAME_ROLE, name);
}

void WindowEditAnnotation::PopulateColorTable(bool bForce)
{
  if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
    return;

  ui->treeWidgetExistingLabels->blockSignals(true);
  SurfaceAnnotation* annot = m_layerSurface->GetActiveAnnotation();
  COLOR_TABLE* ct = annot->GetColorTable();
  if (ct) // && (bForce || ct != m_ctab))
  {
    QTreeWidgetItem* item = ui->treeWidgetExistingLabels->currentItem();
    int last_index = -1;
    if (item)
      last_index = item->data(0, ANNOT_INDEX_ROLE).toInt();
    ui->treeWidgetExistingLabels->clear();
    int nTotalCount = 0;
    CTABgetNumberOfTotalEntries( ct, &nTotalCount );
    int nValid = 0;
    int nSel = -1;
    QList<int> annotIds = annot->GetExistingAnnotations();
    QList<int> selectedLabels = annot->GetVisibleLabels();
    int nValidCount = 0;
    bool bHasSelected = false, bHasUnselected = false;
    for ( int n = 0; n < annotIds.size(); n++ )
    {
      int i = annotIds[n];
      bool bUnassigned = false;
      if (i >= UNASSIGNED_ANNOT_BASE)
        bUnassigned = true;
      else
        CTABisEntryValid( ct, i, &nValid );
      QString name;
      QColor color;
      if ( bUnassigned || nValid )
      {
        QTreeWidgetItem* item = new QTreeWidgetItem(ui->treeWidgetExistingLabels);
        if (bUnassigned)
        {
          NewAnnotationLabel nl = annot->GetNewLabels().value(i);
          name = nl.name;
          color = nl.color;
        }
        else
        {
          char ch[1000];
          CTABcopyName( ct, i, ch, 1000 );
          name = ch;
          int nr, ng, nb;
          CTABrgbAtIndexi( ct, i, &nr, &ng, &nb );
          color = QColor( nr, ng, nb );
        }
        UpdateLabelItem(item, i, name, color);
        item->setCheckState(0,  selectedLabels.contains(i)?Qt::Checked:Qt::Unchecked);
        if (i > 0)
        {
          if (item->checkState(0) == Qt::Checked)
            bHasSelected = true;
          else
            bHasUnselected = true;
        }
        if ( i == last_index )
        {
          nSel = nValidCount;
        }
        nValidCount++;
      }
    }
//    QMap<int, NewAnnotationLabel> temp_labels = annot->GetNewLabels();
//    QList<int> temp_ids = temp_labels.keys();
//    foreach(int id, temp_ids)
//    {
//      NewAnnotationLabel nl = temp_labels[id];
//      QTreeWidgetItem* item = new QTreeWidgetItem(ui->treeWidgetAvailableLabels);
//      item->setText( 0, QString("%1").arg(nl.name) );
//      item->setToolTip( 0, nl.name );
//      QPixmap pix(13, 13);
//      pix.fill( nl.color );
//      item->setIcon(0, QIcon(pix) );
//      item->setData(0, ANNOT_COLOR_ROLE, nl.color);
//      item->setData(0, ANNOT_INDEX_ROLE, nl.id);
//      item->setData(0, ANNOT_NAME_ROLE, nl.name);
//      item->setCheckState(0,  selectedLabels.contains(nl.id)?Qt::Checked:Qt::Unchecked);
//    }
    if ( nSel >= 0 )
    {
      ui->treeWidgetExistingLabels->setCurrentItem( ui->treeWidgetExistingLabels->topLevelItem( nSel ) );
    }
    if (bHasSelected && !bHasUnselected)
      ui->checkBoxShowAll->setCheckState(Qt::Checked);
    else if (bHasSelected)
      ui->checkBoxShowAll->setCheckState(Qt::PartiallyChecked);
    else
      ui->checkBoxShowAll->setCheckState(Qt::Unchecked);
  }
  ui->treeWidgetExistingLabels->blockSignals(false);
}

void WindowEditAnnotation::OnExistingLabelClicked(QTreeWidgetItem *item)
{
  if (item)
  {
    if (isVisible())
      emit LabelClicked(item->data(0, ANNOT_INDEX_ROLE).toInt());

    ui->lineEditName->blockSignals(true);
    ui->lineEditName->setText(item->data(0, ANNOT_NAME_ROLE).toString());
    ui->lineEditName->blockSignals(false);
    ui->colorpicker->blockSignals(true);
    ui->colorpicker->setCurrentColor(item->data(0, ANNOT_COLOR_ROLE).value<QColor>());
    ui->colorpicker->blockSignals(false);

    ui->treeWidgetAvailableLabels->setCurrentItem(NULL);
  }
}

void WindowEditAnnotation::OnExistingLabelItemChanged(QTreeWidgetItem *item)
{
  ui->checkBoxShowAll->blockSignals(true);
  ui->checkBoxShowAll->setCheckState(Qt::PartiallyChecked);
  if ( m_layerSurface && m_layerSurface->GetActiveAnnotation() )
  {
    int nVal = item->data(0, ANNOT_INDEX_ROLE).toInt();
    m_layerSurface->GetActiveAnnotation()->SetSelectLabel(nVal, item->checkState(0) == Qt::Checked);
    ui->checkBoxShowAll->setCheckState(m_layerSurface->GetActiveAnnotation()->GetVisibleLabels().isEmpty()?Qt::Unchecked:Qt::PartiallyChecked);
  }

  ui->checkBoxShowAll->blockSignals(false);
}

void WindowEditAnnotation::OnCheckBoxShowAllLabels(int nState)
{
  ui->treeWidgetExistingLabels->blockSignals(true);
  if (nState == Qt::PartiallyChecked)
  {
    ui->checkBoxShowAll->blockSignals(true);
    ui->checkBoxShowAll->setCheckState(Qt::Checked);
    ui->checkBoxShowAll->blockSignals(false);
  }
  for ( int i = 0; i < ui->treeWidgetExistingLabels->topLevelItemCount(); i++ )
  {
    QTreeWidgetItem* item = ui->treeWidgetExistingLabels->topLevelItem( i );
    item->setCheckState(0, nState == Qt::Unchecked ? Qt::Unchecked : Qt::Checked);
  }
  ui->treeWidgetExistingLabels->blockSignals(false);

  if ( m_layerSurface && m_layerSurface->GetActiveAnnotation() )
  {
    if (nState == Qt::Unchecked)
      m_layerSurface->GetActiveAnnotation()->SetUnselectAllLabels();
    else
      m_layerSurface->GetActiveAnnotation()->SetSelectAllLabels();
  }
}

void WindowEditAnnotation::OnSurfaceVertexClicked(LayerSurface *surf)
{
  if (surf == m_layerSurface)
  {
    if (surf->GetActiveAnnotation() && m_layerSurface->GetCurrentVertex() >= 0)
    {
      int n = surf->GetActiveAnnotation()->GetIndexAtVertex(m_layerSurface->GetCurrentVertex());
      for ( int i = 0; i < ui->treeWidgetExistingLabels->topLevelItemCount(); i++ )
      {
        QTreeWidgetItem* item = ui->treeWidgetExistingLabels->topLevelItem( i );
        if (item->data(0, ANNOT_INDEX_ROLE).toInt() == n)
        {
          ui->treeWidgetExistingLabels->setCurrentItem(item);
        }
      }
    }
  }
}

int WindowEditAnnotation::GetCurrentIndex()
{
  QTreeWidgetItem* item = ui->treeWidgetExistingLabels->currentItem();
  if (item)
    return item->data(0, ANNOT_INDEX_ROLE).toInt();
  else
    return -1;
}

void WindowEditAnnotation::OnEditNameReturnPressed()
{
  QTreeWidgetItem* item = ui->treeWidgetExistingLabels->currentItem();
  QString name = ui->lineEditName->text().trimmed();
  if (!name.isEmpty() && item)
  {
    int i = item->data(0, ANNOT_INDEX_ROLE).toInt();
    item->setText(0, QString("%1 %2").arg(i).arg(name));
    item->setData(0, ANNOT_NAME_ROLE, name);
    if (m_layerSurface && m_layerSurface->GetActiveAnnotation())
      m_layerSurface->GetActiveAnnotation()->UpdateLabelInfo(i, name);
  }
}

void WindowEditAnnotation::OnColorChanged(const QColor &color)
{
  QTreeWidgetItem* item = ui->treeWidgetExistingLabels->currentItem();
  if (item)
  {
    int i = item->data(0, ANNOT_INDEX_ROLE).toInt();
    item->setData(0, ANNOT_COLOR_ROLE, color);
    QPixmap pix(13, 13);
    pix.fill( color );
    item->setIcon(0, QIcon(pix));
    if (m_layerSurface && m_layerSurface->GetActiveAnnotation())
      m_layerSurface->GetActiveAnnotation()->UpdateLabelInfo(i, "", color);
  }
}

void WindowEditAnnotation::OnAvailableLabelClicked(QTreeWidgetItem *item)
{
  if (item)
  {
    ui->lineEditName->blockSignals(true);
    ui->lineEditName->setText(item->data(0, ANNOT_NAME_ROLE).toString());
    ui->lineEditName->blockSignals(false);
    ui->colorpicker->blockSignals(true);
    ui->colorpicker->setCurrentColor(item->data(0, ANNOT_COLOR_ROLE).value<QColor>());
    ui->colorpicker->blockSignals(false);
    ui->pushButtonSet->setEnabled(true);
  }
}

void WindowEditAnnotation::OnButtonSet()
{
  QTreeWidgetItem* itemTo = ui->treeWidgetExistingLabels->currentItem();
  QTreeWidgetItem* itemFrom = ui->treeWidgetAvailableLabels->currentItem();
  if (itemTo && itemFrom)
  {
    int i = itemFrom->data(0, ANNOT_INDEX_ROLE).toInt();
    int old_i = itemTo->data(0, ANNOT_INDEX_ROLE).toInt();
    QString name = ui->lineEditName->text().trimmed();
    QColor color = ui->colorpicker->currentColor();
    UpdateLabelItem(itemTo, i, name, color);
    if (m_layerSurface && m_layerSurface->GetActiveAnnotation())
      m_layerSurface->GetActiveAnnotation()->ReassignNewLabel(old_i, i);
  }
}
