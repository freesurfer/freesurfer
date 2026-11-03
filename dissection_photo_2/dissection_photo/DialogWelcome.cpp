#include "DialogWelcome.h"
#include <QLabel>
#include <QScreen>
#include <QGuiApplication>
#include "ui_DialogWelcome.h"
#include "AppTheme.h"

#include <QBuffer>
#include <QHBoxLayout>
#include <QPushButton>
#include <QVBoxLayout>

DialogWelcome::DialogWelcome(QWidget *parent) :
  QDialog(parent),
  ui(new Ui::DialogWelcome)
{
  ui->setupUi(this);

  // When shown, this precedes the Hub, so it carries the appearance override
  // too (it can be skipped permanently). Inserted above the Designer content in
  // code, matching DialogProjectsHub.
  m_btnDarkMode = new QPushButton(this);
  m_btnDarkMode->setCheckable(true);
  m_btnDarkMode->setChecked(AppTheme::EffectiveDarkMode());
  m_btnDarkMode->setCursor(Qt::PointingHandCursor);
  m_btnDarkMode->setText(m_btnDarkMode->isChecked() ? "☀  Light Mode" : "☽  Dark Mode");
  connect(m_btnDarkMode, &QPushButton::toggled, this, [this](bool bDark) {
    m_btnDarkMode->setText(bDark ? "☀  Light Mode" : "☽  Dark Mode");
    OnToggleDarkMode(bDark);
  });
  QHBoxLayout* header = new QHBoxLayout();
  header->addStretch();
  header->addWidget(m_btnDarkMode);
  if (QVBoxLayout* vLayout = qobject_cast<QVBoxLayout*>(layout()))
    vLayout->insertLayout(0, header);

  // The link glyph goes inside the anchor so glyph and text are one clickable
  // unit. Rich text renders data: URLs, so the tinted glyph is embedded as a
  // PNG; the .ui's text-only anchor is the fallback. 32px so it draws 1:1 in a
  // 16px box at DPR 2.
  const QPixmap glyph = AppTheme::TintedPixmap("link", QColor("#c95659"));
  if (!glyph.isNull())
  {
    QByteArray png;
    QBuffer buf(&png);
    buf.open(QIODevice::WriteOnly);
    glyph.scaled(32, 32, Qt::KeepAspectRatio, Qt::SmoothTransformation).save(&buf, "PNG");
    ui->label_10->setText(QString("<a style='color:#c95659;font-size:16px;text-decoration:none' "
                                  "href='https://surfer.nmr.mgh.harvard.edu/fswiki/PhotoTools'>"
                                  "<img src='data:image/png;base64,%1' width='16' height='16' "
                                  "style='vertical-align:middle'/>&nbsp;Photo Tools Wiki Page &gt;</a>")
                          .arg(QString::fromLatin1(png.toBase64())));
  }

  Retheme(AppTheme::EffectiveDarkMode());
}

void DialogWelcome::Retheme(bool bDark)
{
  // Widget-level sheets (this and the divider's in the .ui) must be re-set by
  // hand on a mode change.
  m_btnDarkMode->setStyleSheet(AppTheme::ToggleStyleSheet(bDark));
  // The .ui's light-grey hairline is too bright on dark.
  ui->frameWelcomeDivider->setStyleSheet(QString("background: %1;").arg(bDark ? "#4a4a4a" : "#d0d0d0"));
}

DialogWelcome::~DialogWelcome()
{
  delete ui;
}

bool DialogWelcome::GetDoNotShow()
{
  return ui->checkBoxDoNotShow->isChecked();
}

void DialogWelcome::OnToggleDarkMode(bool bDark)
{
  AppTheme::SetMode(bDark ? "dark" : "light");
  // Apply(), not setPalette(): the app-wide sheet must be swapped too (see
  // AppTheme::BaseStyleSheet()).
  AppTheme::Apply(bDark);
  Retheme(bDark);
}

void DialogWelcome::UpdateBodyHeights()
{
  // The four wrapped description paragraphs (see resizeEvent() in the header).
  QLabel* bodies[] = { ui->label_3, ui->label_7, ui->labelWelcomeCCBody, ui->label_5 };
  for (QLabel* l : bodies)
  {
    if (!l || l->width() <= 0)
      continue;
    const int nNeed = l->heightForWidth(l->width());
    if (nNeed > 0 && nNeed != l->minimumHeight())
      l->setMinimumHeight(nNeed);
  }
}

void DialogWelcome::resizeEvent(QResizeEvent* e)
{
  QDialog::resizeEvent(e);
  UpdateBodyHeights();
}

void DialogWelcome::showEvent(QShowEvent* e)
{
  QDialog::showEvent(e);
  UpdateBodyHeights();
  // Resize to the rows' new height, clamped to the screen.
  adjustSize();
  // QWidget::screen() is Qt 5.14; FreeSurfer still pins 5.12 on several
  // platforms, so go through the application instead.
  QScreen* pScreen = QGuiApplication::primaryScreen();
  if (pScreen)
  {
    const QRect rcAvail = pScreen->availableGeometry();
    const int nMax = int(rcAvail.height() * 0.92);
    if (height() > nMax)
      resize(width(), nMax);
  }
}
