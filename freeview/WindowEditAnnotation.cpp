#include "WindowEditAnnotation.h"
#include "ui_WindowEditAnnotation.h"
#include <cstddef>
#include "LayerSurface.h"
#include "SurfaceAnnotation.h"
#include "LayerPropertySurface.h"
#include "MainWindow.h"
#include "RenderView.h"
#include <QTreeWidgetItem>
#include <QSettings>
#include <QFileDialog>
#include "MigrationDefs.h"
#ifdef Q_OS_MAC
#include "MacHelper.h"
#endif

#define ANNOT_COLOR_ROLE  Qt::UserRole
#define ANNOT_INDEX_ROLE  Qt::UserRole+1
#define ANNOT_NAME_ROLE  Qt::UserRole+2

WindowEditAnnotation::WindowEditAnnotation(QWidget *parent) :
  QWidget(parent), m_layerSurface(NULL), m_ctab(NULL),
  ui(new Ui::WindowEditAnnotation)
{
  ui->setupUi(this);
  setWindowFlags( Qt::Tool );
#ifdef Q_OS_MAC
  if (MacHelper::IsDarkMode())
  {
      QIcon icn(":/resource/icons/surface_path_dm.png");
      icn.addPixmap(QPixmap(":/resource/icons/surface_path_dm.png"), QIcon::Normal, QIcon::On);
      ui->actionMakePath->setIcon(QIcon(":/resource/icons/surface_path_make_dm.png"));
      ui->actionMakeClosedPath->setIcon(QIcon(":/resource/icons/surface_path_make_closed_dm.png"));
      ui->actionDeletePath->setIcon(QIcon(":/resource/icons/surface_path_delete_dm.png"));
  }
#endif

  connect(ui->treeWidgetExistingLabels, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)),
          SLOT(OnExistingLabelClicked(QTreeWidgetItem*)), Qt::QueuedConnection);
  connect(ui->treeWidgetExistingLabels, SIGNAL(itemChanged(QTreeWidgetItem*,int)), SLOT(OnExistingLabelItemChanged(QTreeWidgetItem*)));
  connect(ui->checkBoxShowAll, SIGNAL(stateChanged(int)), SLOT(OnCheckBoxShowAllLabels(int)));
  connect(ui->lineEditName, SIGNAL(returnPressed()), SLOT(OnButtonSet()));
  connect(ui->lineEditIndex, SIGNAL(returnPressed()), SLOT(OnButtonSet()));
  connect(ui->lineEditColor, SIGNAL(returnPressed()), SLOT(OnButtonSet()));
  connect(ui->lineEditColor, SIGNAL(textChanged(QString)), SLOT(OnEditColorTextChanged()));
  connect(ui->colorpicker, SIGNAL(colorChanged(QColor)), SLOT(OnColorChanged(QColor)), Qt::QueuedConnection);
  connect(ui->treeWidgetAvailableLabels, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)),
          SLOT(UpdateActions()), Qt::QueuedConnection);
  connect(ui->actionRedo, SIGNAL(triggered(bool)), SLOT(OnButtonRedo()));
  connect(ui->actionUndo, SIGNAL(triggered(bool)), SLOT(OnButtonUndo()));
  connect(ui->pushButtonSet, SIGNAL(clicked(bool)), SLOT(OnButtonSet()));
  connect(ui->toolButtonFromCTab, SIGNAL(clicked(bool)), SLOT(OnButtonFromCTab()));
  connect(ui->pushButtonDelete, SIGNAL(clicked(bool)), SLOT(OnButtonDelete()));
  connect(ui->pushButtonCleanUp, SIGNAL(clicked(bool)), SLOT(OnButtonCleanUp()));
  connect(ui->pushButtonLoadSegmentation, SIGNAL(clicked(bool)), SLOT(OnButtonLoadSegmentation()));
  connect(ui->pushButtonLoadLUT, SIGNAL(clicked(bool)), SLOT(OnButtonLoadColorTable()));

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

void WindowEditAnnotation::UpdateUI(int nIndex)
{
  if (m_layerSurface)
  {
    QList<QWidget*> list = this->findChildren<QWidget*>();
    foreach (QWidget* w, list)
      w->blockSignals(true);

    PopulateColorTable(nIndex);
    ui->toolButtonFromCTab->setEnabled(ui->treeWidgetAvailableLabels->currentItem() != NULL);

    foreach (QWidget* w, list)
      w->blockSignals(false);

    UpdateActions();
  }
}

void WindowEditAnnotation::UpdateActions()
{
  SurfaceAnnotation* annot = m_layerSurface->GetActiveAnnotation();
  ui->actionRedo->setEnabled(annot && annot->HasRedo());
  ui->actionUndo->setEnabled(annot && annot->HasUndo());
  ui->toolButtonFromCTab->setEnabled(!ui->treeWidgetAvailableLabels->selectedItems().isEmpty() &&
                                !ui->treeWidgetExistingLabels->selectedItems().isEmpty());
}

void WindowEditAnnotation::PopulateAvailableColorTable(bool bForce)
{
  if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
    return;

  ui->treeWidgetAvailableLabels->blockSignals(true);
  SurfaceAnnotation* annot = m_layerSurface->GetActiveAnnotation();
  COLOR_TABLE* ct = annot->GetColorTable();
  if (ct && (bForce || ct != m_ctab))
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

void WindowEditAnnotation::PopulateColorTable(int nIndex)
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
    if (nIndex >= 0)
      last_index = nIndex;
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
    if ( nSel >= 0 )
    {
      ui->treeWidgetExistingLabels->setCurrentItem( ui->treeWidgetExistingLabels->topLevelItem( nSel ) );
      UpdateInfoFromItem(ui->treeWidgetExistingLabels->currentItem());
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

    UpdateInfoFromItem(item);
  }
  UpdateActions();
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

void WindowEditAnnotation::OnEditColorTextChanged()
{
  QStringList list = ui->lineEditColor->text().split(" ", MD_SkipEmptyParts);
  if (list.size() == 3)
  {
    int rgb[3];
    for (int i = 0; i < 3; i++)
    {
      bool ok;
      rgb[i] = list[i].toInt(&ok);
      if (!ok || rgb[i] < 0 || rgb[i] > 255)
        return;
    }
    QColor c(rgb[0], rgb[1], rgb[2]);
    ui->colorpicker->blockSignals(true);
    ui->colorpicker->setCurrentColor(c);
    ui->colorpicker->blockSignals(false);
  }
}

void WindowEditAnnotation::OnButtonSet()
{
  QTreeWidgetItem* item = ui->treeWidgetExistingLabels->currentItem();
  if (item)
  {
    ui->treeWidgetExistingLabels->blockSignals(true);
    int index_old = item->data(0, ANNOT_INDEX_ROLE).toInt();
    bool ok;
    int index_new = ui->lineEditIndex->text().toInt(&ok);
    if (!ok || index_new < 0)
      index_new = index_old;
    QColor color = ui->colorpicker->currentColor();
    QString name = ui->lineEditName->text().trimmed();
    if (name.isEmpty())
    {
      QMessageBox::warning(this, "Error", "Please provide an unique name for annotation");
      return;
    }
    for ( int n = 0; n < ui->treeWidgetExistingLabels->topLevelItemCount(); n++ )
    {
      QTreeWidgetItem* aitem = ui->treeWidgetExistingLabels->topLevelItem( n );
      if (aitem != item && aitem->data(0, ANNOT_INDEX_ROLE).toInt() == index_new)
      {
        if (QMessageBox::question(this, "Confirm", "This label already exists in the annotation. Do you want to combine them?")
            != QMessageBox::Yes)
          return;
        else
          delete ui->treeWidgetExistingLabels->takeTopLevelItem(n);
      }
    }

    UpdateLabelItem(item, index_new, name, color);

    ui->treeWidgetExistingLabels->blockSignals(false);
    if (m_layerSurface && m_layerSurface->GetActiveAnnotation())
    {
      if (index_new != index_old)
        m_layerSurface->GetActiveAnnotation()->ReassignNewLabel(index_old, index_new, name, color);
      else
        m_layerSurface->GetActiveAnnotation()->UpdateLabelInfo(index_old, name, color);
      PopulateAvailableColorTable(true);
    }
    UpdateActions();
  }
}

void WindowEditAnnotation::OnColorChanged(const QColor &color)
{
  QTreeWidgetItem* item = ui->treeWidgetExistingLabels->currentItem();
  if (item)
  {
    for ( int i = 0; i < ui->treeWidgetExistingLabels->topLevelItemCount(); i++ )
    {
      QTreeWidgetItem* other_item = ui->treeWidgetExistingLabels->topLevelItem( i );
      if (item != other_item && other_item->data(0, ANNOT_COLOR_ROLE).value<QColor>() == color)
      {
        QMessageBox::information(this, "Existing Color", "This color already exists in the color table. Please choose another one.");
        ui->colorpicker->blockSignals(true);
        ui->colorpicker->setCurrentColor(item->data(0, ANNOT_COLOR_ROLE).value<QColor>());
        ui->colorpicker->blockSignals(false);
        return;
      }
    }
//    ui->treeWidgetExistingLabels->blockSignals(true);
//    int i = item->data(0, ANNOT_INDEX_ROLE).toInt();
//    item->setData(0, ANNOT_COLOR_ROLE, color);
//    QPixmap pix(13, 13);
//    pix.fill( color );
//    item->setIcon(0, QIcon(pix));
//    ui->treeWidgetExistingLabels->blockSignals(false);
//    if (m_layerSurface && m_layerSurface->GetActiveAnnotation())
//      m_layerSurface->GetActiveAnnotation()->UpdateLabelInfo(i, "", color);

    ui->lineEditColor->blockSignals(true);
    ui->lineEditColor->setText(QString("%1 %2 %3").arg(color.red()).arg(color.green()).arg(color.blue()));
    ui->lineEditColor->blockSignals(false);
  }
}

void WindowEditAnnotation::UpdateInfoFromItem(QTreeWidgetItem *item)
{
  ui->lineEditName->blockSignals(true);
  ui->lineEditName->setText(item->data(0, ANNOT_NAME_ROLE).toString());
  ui->lineEditName->blockSignals(false);

  ui->colorpicker->blockSignals(true);
  QColor c = item->data(0, ANNOT_COLOR_ROLE).value<QColor>();
  ui->colorpicker->setCurrentColor(c);
  ui->colorpicker->blockSignals(false);
  ui->lineEditColor->blockSignals(true);
  ui->lineEditColor->setText(QString("%1 %2 %3").arg(c.red()).arg(c.green()).arg(c.blue()));
  ui->lineEditColor->blockSignals(false);

  ui->lineEditIndex->blockSignals(true);
  int nIndex = item->data(0, ANNOT_INDEX_ROLE).toInt();
  if (nIndex < UNASSIGNED_ANNOT_BASE)
    ui->lineEditIndex->setText(QString::number(nIndex));
  else
    ui->lineEditIndex->clear();
  ui->lineEditIndex->blockSignals(false);
}

void WindowEditAnnotation::OnButtonFromCTab()
{
  QTreeWidgetItem* itemTo = ui->treeWidgetExistingLabels->currentItem();
  QTreeWidgetItem* itemFrom = ui->treeWidgetAvailableLabels->currentItem();
  if (itemTo && itemFrom)
  {
    UpdateInfoFromItem(itemFrom);
    OnButtonSet();
  }
}

void WindowEditAnnotation::OnButtonUndo()
{
  if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
      return;

  m_layerSurface->GetActiveAnnotation()->Undo();
  UpdateUI();
  PopulateAvailableColorTable(true);
}

void WindowEditAnnotation::OnButtonRedo()
{
  if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
      return;

  m_layerSurface->GetActiveAnnotation()->Redo();
  UpdateUI();
  PopulateAvailableColorTable(true);
}

void WindowEditAnnotation::OnButtonDelete()
{
  if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
      return;

  QTreeWidgetItem* item = ui->treeWidgetExistingLabels->currentItem();
  if (item)
  {
    m_layerSurface->GetActiveAnnotation()->DeleteLabel(item->data(0, ANNOT_INDEX_ROLE).toInt());
    int n = ui->treeWidgetExistingLabels->indexOfTopLevelItem(item);
    ui->treeWidgetExistingLabels->takeTopLevelItem(n);
    delete item;
    if (n < ui->treeWidgetExistingLabels->topLevelItemCount())
      ui->treeWidgetExistingLabels->setCurrentItem(ui->treeWidgetExistingLabels->topLevelItem(n));
    else if (ui->treeWidgetExistingLabels->topLevelItemCount() > 0)
      ui->treeWidgetExistingLabels->setCurrentItem(ui->treeWidgetExistingLabels->topLevelItem(0));

    UpdateActions();
  }
}

void WindowEditAnnotation::OnButtonCleanUp()
{
  if (m_layerSurface && m_layerSurface->GetActiveAnnotation())
  {
    m_layerSurface->GetActiveAnnotation()->CleanUpColorTable();
    UpdateActions();
    PopulateAvailableColorTable(true);
  }
}

void WindowEditAnnotation::OnButtonLoadSegmentation()
{
    if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
        return;

    if (ui->treeWidgetAvailableLabels->topLevelItemCount() == 0)
    {
        QMessageBox::information(this, "Load Annotation", "Please load color table first");
        return;
    }
    QString fn = QFileDialog::getOpenFileName(this, "Load Segmentation", MainWindow::GetMainWindow()->AutoSelectLastDir("mri" ),
                                 "Segmentation files (*.mgz *.mgh *.nii *.nii.gz)");
    if (!fn.isEmpty())
    {
        if (!m_layerSurface->GetActiveAnnotation()->LoadFromSegmentation(fn))
            QMessageBox::warning(this, "Annotation", "Failed to load annotation from segmentation file");
        else
            UpdateUI();
    }
}

void WindowEditAnnotation::OnButtonLoadColorTable()
{
    if (!m_layerSurface || !m_layerSurface->GetActiveAnnotation())
        return;

    QString fn = QFileDialog::getOpenFileName(this, "Load Color Table", MainWindow::GetMainWindow()->AutoSelectLastDir("label" ),
                                 "Look Up Table files (*)");
    if (!fn.isEmpty())
    {
        if (!m_layerSurface->GetActiveAnnotation()->LoadColorTable(fn))
            QMessageBox::warning(this, "Color Table", "Failed to load color table");
        else
            PopulateAvailableColorTable();
    }
}
