#include "ToolWindowLesionPopup.h"
#include "ui_ToolWindowLesionPopup.h"
#include "LayerPointSet.h"
#include "LayerPropertyPointSet.h"

ToolWindowLesionPopup::ToolWindowLesionPopup(QWidget *parent) :
    QWidget(parent), UIUpdateHelper(), m_curLayer(NULL),
    ui(new Ui::ToolWindowLesionPopup)
{
  ui->setupUi(this);
  this->setWindowFlags( Qt::Tool | Qt::WindowTitleHint | Qt::CustomizeWindowHint );

  connect(ui->spinBoxGoToPoint, SIGNAL(valueChanged(int)), SLOT(OnSpinBoxGoToPoint(int)));
  connect(ui->pushButtonGoToPoint, SIGNAL(clicked()), SIGNAL(GoToPointTriggered()));

  connect(ui->colorpickerPointColor, SIGNAL(colorChanged(QColor)), SIGNAL(PointColorChanged(QColor)));
  connect(ui->lineEditRadius, SIGNAL(textChanged(QString)), SIGNAL(RadiusTextChanged(QString)));

  connect(ui->checkBoxShowLesionOnly, SIGNAL(toggled(bool)), SLOT(OnCheckBoxShowRegion(bool)), Qt::QueuedConnection);
  connect(ui->checkBoxShowNonLesionOnly, SIGNAL(toggled(bool)), SLOT(OnCheckBoxShowRegion(bool)), Qt::QueuedConnection);
  connect(ui->checkBoxShowVRSpaceOnly, SIGNAL(toggled(bool)), SLOT(OnCheckBoxShowRegion(bool)), Qt::QueuedConnection);

  connect(ui->checkBoxRegionLesion, SIGNAL(toggled(bool)), SLOT(OnCheckBoxSetRegion()), Qt::QueuedConnection);
  connect(ui->checkBoxRegionNonLesion, SIGNAL(toggled(bool)), SLOT(OnCheckBoxSetRegion()), Qt::QueuedConnection);
  connect(ui->checkBoxRegionVRSpace, SIGNAL(toggled(bool)), SLOT(OnCheckBoxSetRegion()), Qt::QueuedConnection);

  connect(ui->checkBoxNonLesionCloseToCortex, SIGNAL(toggled(bool)), SLOT(OnCheckBoxNonLesionReason(bool)), Qt::QueuedConnection);
  connect(ui->checkBoxNonLesionCloseToGM, SIGNAL(toggled(bool)), SLOT(OnCheckBoxNonLesionReason(bool)), Qt::QueuedConnection);
  connect(ui->checkBoxNonLesionOutsideOfBrain, SIGNAL(toggled(bool)), SLOT(OnCheckBoxNonLesionReason(bool)), Qt::QueuedConnection);
  connect(ui->checkBoxNonLesionHighB0Distortion, SIGNAL(toggled(bool)), SLOT(OnCheckBoxNonLesionReason(bool)), Qt::QueuedConnection);
  connect(ui->checkBoxNonLesionPoorImageQuality, SIGNAL(toggled(bool)), SLOT(OnCheckBoxNonLesionReason(bool)), Qt::QueuedConnection);

  connect(ui->sliderConfidenceScore, SIGNAL(valueChanged(int)), SLOT(OnConfidenceScoreChanged(int)));
  connect(ui->spinBoxConfidenceScore, SIGNAL(valueChanged(int)), SLOT(OnConfidenceScoreChanged(int)));

  connect(ui->textEditAdditionalComment, SIGNAL(textChanged()), SLOT(OnTextAdditionalComment()));

  connect(ui->doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), SIGNAL(OpacityChanged(double)));
  connect(ui->sliderOpacity, SIGNAL(valueChanged(int)), SLOT(OnSliderOpacity(int)));
}

ToolWindowLesionPopup::~ToolWindowLesionPopup()
{
  delete ui;
}

void ToolWindowLesionPopup::UpdateUI(LayerPointSet *layer)
{
  m_curLayer = layer;
  if (!layer)
  {
 //   hide();
    return;
  }

  BlockAllSignals(true);

  ui->lineEditFileName->clear();
  ui->sliderOpacity->setValue( (int)( layer->GetProperty()->GetOpacity() * 100 ) );
  ChangeDoubleSpinBoxValue( ui->doubleSpinBoxOpacity, layer->GetProperty()->GetOpacity() );
  double* rgb = layer->GetProperty()->GetColor();
  ui->colorpickerPointColor->setCurrentColor( QColor( (int)(rgb[0]*255), (int)(rgb[1]*255), (int)(rgb[2]*255) ) );
  ui->lineEditFileName->setText(layer->GetFileName());
  ui->lineEditFileName->setCursorPosition( ui->lineEditFileName->text().size() );
  ChangeLineEditNumber( ui->lineEditRadius, layer->GetProperty()->GetRadius() );
  ui->labelNumberTotal->setText(QString("%1").arg(layer->GetNumberOfPoints()));
  ui->spinBoxGoToPoint->setRange(1, layer->GetNumberOfPoints());

  ui->checkBoxShowLesionOnly->setChecked(layer->GetProperty()->GetShowLesionsOnly());
  ui->checkBoxShowNonLesionOnly->setChecked(layer->GetProperty()->GetShowNonLesionsOnly());
  ui->checkBoxShowVRSpaceOnly->setChecked(layer->GetProperty()->GetShowVRSpaceOnly());

  BlockAllSignals(false);
}

void ToolWindowLesionPopup::BlockAllSignals(bool bBlock)
{
  QList<QWidget*> allWidgets = this->findChildren<QWidget*>();
  for ( int i = 0; i < allWidgets.size(); i++ )
  {
    allWidgets[i]->blockSignals( bBlock );
  }
}

void ToolWindowLesionPopup::UpdatePointInfo(int nIndex, ControlPoint* p)
{
  SetCurrentPoint(nIndex);
  BlockAllSignals(true);
  QString region = p->info.value("pathological_region").toString().toLower();
  ui->checkBoxRegionLesion->setChecked(region == "lesion");
  ui->checkBoxRegionNonLesion->setChecked(region == "non_lesion");
  ui->checkBoxRegionVRSpace->setChecked(region == "vr_space");
  ui->sliderConfidenceScore->setValue(p->info.value("confidence_score").toInt());
  ui->spinBoxConfidenceScore->setValue(p->info.value("confidence_score").toInt());
  QVariantMap non_lesion_reason = p->info.value("non_lession_reason").toMap();
  ui->checkBoxNonLesionCloseToCortex->setChecked(non_lesion_reason.value("close_to_cortex").toBool());
  ui->checkBoxNonLesionCloseToGM->setChecked(non_lesion_reason.value("close_to_gm").toBool());
  ui->checkBoxNonLesionOutsideOfBrain->setChecked(non_lesion_reason.value("outside_of_brain").toBool());
  ui->checkBoxNonLesionHighB0Distortion->setChecked(non_lesion_reason.value("high_b0_distortion").toBool());
  ui->checkBoxNonLesionPoorImageQuality->setChecked(non_lesion_reason.value("poor_image_quality").toBool());
  ui->textEditAdditionalComment->setPlainText(p->info.value("additional_comment").toString());
  BlockAllSignals(false);
}

void ToolWindowLesionPopup::SetCurrentPoint(int nIndex)
{
  ui->spinBoxGoToPoint->blockSignals(true);
  ui->spinBoxGoToPoint->setValue(nIndex+1);
  ui->spinBoxGoToPoint->blockSignals(false);
}

void ToolWindowLesionPopup::OnCheckBoxShowRegion(bool bShow)
{
  if (!m_curLayer)
    return;

  if (sender() == ui->checkBoxShowLesionOnly)
    m_curLayer->GetProperty()->SetShowLesionsOnly(bShow);
  else if (sender() == ui->checkBoxShowNonLesionOnly)
    m_curLayer->GetProperty()->SetShowNonLesionsOnly(bShow);
  else if (sender() == ui->checkBoxShowVRSpaceOnly)
    m_curLayer->GetProperty()->SetShowVRSpaceOnly(bShow);
}

void ToolWindowLesionPopup::OnCheckBoxSetRegion()
{
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (m_curLayer && nIndex < m_curLayer->GetNumberOfPoints())
  {
    if (sender() == ui->checkBoxRegionLesion && ui->checkBoxRegionLesion->isChecked())
      m_curLayer->UpdatePoint(nIndex, "pathological_region", "lesion");
    if (sender() == ui->checkBoxRegionNonLesion && ui->checkBoxRegionNonLesion->isChecked())
      m_curLayer->UpdatePoint(nIndex, "pathological_region", "non_lesion");
    if (sender() == ui->checkBoxRegionVRSpace && ui->checkBoxRegionVRSpace->isChecked())
      m_curLayer->UpdatePoint(nIndex, "pathological_region", "vr_space");

    ControlPoint p = m_curLayer->GetPoint(nIndex);
    UpdatePointInfo(nIndex, &p);
  }
}

void ToolWindowLesionPopup::OnSpinBoxGoToPoint(int n)
{
  emit GoToPointChanged(n);
}

void ToolWindowLesionPopup::OnTextAdditionalComment()
{
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (m_curLayer && nIndex < m_curLayer->GetNumberOfPoints())
  {
    m_curLayer->UpdatePoint(nIndex, "additional_comment", ui->textEditAdditionalComment->toPlainText());
  }
}

void ToolWindowLesionPopup::OnCheckBoxNonLesionReason(bool b)
{
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (m_curLayer && nIndex < m_curLayer->GetNumberOfPoints())
  {
    QVariantMap non_lesion_reason = m_curLayer->GetPoint(nIndex).info.value("non_lession_reason").toMap();
    if (sender() == ui->checkBoxNonLesionCloseToCortex)
      non_lesion_reason["close_to_cortex"] = b;
    else if (sender() == ui->checkBoxNonLesionCloseToGM)
      non_lesion_reason["close_to_gm"] = b;
    else if (sender() == ui->checkBoxNonLesionOutsideOfBrain)
      non_lesion_reason["outside_of_brain"] = b;
    else if (sender() == ui->checkBoxNonLesionHighB0Distortion)
      non_lesion_reason["high_b0_distortion"] = b;
    else if (sender() == ui->checkBoxNonLesionPoorImageQuality)
      non_lesion_reason["poor_image_quality"] = b;

    m_curLayer->UpdatePoint(nIndex, "non_lession_reason", non_lesion_reason);
  }
}

void ToolWindowLesionPopup::OnConfidenceScoreChanged(int val)
{
  int nIndex = ui->spinBoxGoToPoint->value()-1;
  if (m_curLayer && nIndex < m_curLayer->GetNumberOfPoints())
  {
    m_curLayer->UpdatePoint(nIndex, "confidence_score", val);
    ControlPoint p = m_curLayer->GetPoint(nIndex);
    UpdatePointInfo(nIndex, &p);
  }
}

void ToolWindowLesionPopup::OnSliderOpacity(int n)
{
  emit OpacityChanged(n/100.0);
}
