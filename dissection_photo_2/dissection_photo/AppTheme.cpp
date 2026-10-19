#include "AppTheme.h"

#include <QApplication>
#include <QFile>
#include <QGuiApplication>
#include <QPainter>
#include <QSettings>
#include <QStyleHints>

bool AppTheme::IsSystemDarkMode()
{
#if QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
  Qt::ColorScheme scheme = QGuiApplication::styleHints()->colorScheme();
  if (scheme != Qt::ColorScheme::Unknown)
    return scheme == Qt::ColorScheme::Dark;
#endif
  // Fallback: the platform palette follows the OS (Fusion's standardPalette()
  // is always light). Only valid before BuildPalette() is applied.
  return qApp && qApp->palette().color(QPalette::Window).lightness() < 128;
}

QString AppTheme::Mode()
{
  QSettings s;
  QString mode = s.value("Appearance/Mode", "system").toString();
  if (mode != "light" && mode != "dark")
    mode = "system";
  return mode;
}

void AppTheme::SetMode(const QString& mode)
{
  QSettings s;
  s.setValue("Appearance/Mode", mode);
}

bool AppTheme::EffectiveDarkMode()
{
  QString mode = Mode();
  if (mode == "dark")
    return true;
  if (mode == "light")
    return false;
  return IsSystemDarkMode();
}

QPalette AppTheme::BuildPalette(bool bDark)
{
  QPalette pal;
  if (bDark)
  {
    pal.setColor(QPalette::Window, QColor(45, 45, 45));
    pal.setColor(QPalette::WindowText, Qt::white);
    pal.setColor(QPalette::Base, QColor(35, 35, 35));
    pal.setColor(QPalette::AlternateBase, QColor(50, 50, 50));
    pal.setColor(QPalette::ToolTipBase, QColor(45, 45, 45));
    pal.setColor(QPalette::ToolTipText, Qt::white);
    pal.setColor(QPalette::Text, Qt::white);
    pal.setColor(QPalette::Button, QColor(55, 55, 55));
    pal.setColor(QPalette::ButtonText, Qt::white);
    pal.setColor(QPalette::BrightText, Qt::red);
    pal.setColor(QPalette::Link, QColor(0x6f, 0xb3, 0xc9));
    pal.setColor(QPalette::Highlight, QColor(0x45, 0x7c, 0x92));
    pal.setColor(QPalette::HighlightedText, Qt::white);
    pal.setColor(QPalette::Disabled, QPalette::Text, QColor(120, 120, 120));
    pal.setColor(QPalette::Disabled, QPalette::WindowText, QColor(120, 120, 120));
    pal.setColor(QPalette::Disabled, QPalette::ButtonText, QColor(120, 120, 120));
  }
  else
  {
    // Soft near-black; pure black is too harsh for this muted palette.
    pal.setColor(QPalette::Window, QColor(240, 240, 240));
    pal.setColor(QPalette::WindowText, QColor(0x26, 0x26, 0x26));
    pal.setColor(QPalette::Base, Qt::white);
    pal.setColor(QPalette::AlternateBase, QColor(245, 245, 245));
    pal.setColor(QPalette::ToolTipBase, Qt::white);
    pal.setColor(QPalette::ToolTipText, QColor(0x26, 0x26, 0x26));
    pal.setColor(QPalette::Text, QColor(0x26, 0x26, 0x26));
    pal.setColor(QPalette::Button, QColor(240, 240, 240));
    pal.setColor(QPalette::ButtonText, QColor(0x26, 0x26, 0x26));
    pal.setColor(QPalette::BrightText, Qt::red);
    pal.setColor(QPalette::Link, QColor(0x45, 0x7c, 0x92));
    pal.setColor(QPalette::Highlight, QColor(0x45, 0x7c, 0x92));
    pal.setColor(QPalette::HighlightedText, Qt::white);
    pal.setColor(QPalette::Disabled, QPalette::Text, QColor(160, 160, 160));
    pal.setColor(QPalette::Disabled, QPalette::WindowText, QColor(160, 160, 160));
    pal.setColor(QPalette::Disabled, QPalette::ButtonText, QColor(160, 160, 160));
  }
  return pal;
}

QString AppTheme::BaseStyleSheet(bool bDark)
{
  QString strBase;
  QFile styleFile(":/resource/style.qss");
  if (styleFile.open(QFile::ReadOnly | QFile::Text))
    strBase = QString::fromUtf8(styleFile.readAll());
  if (!bDark)
    return strBase;

  // Same tones as BuildPalette(true).
  static const char* kDarkCss =
      "\n/* ---- dark-mode overrides, appended by AppTheme::BaseStyleSheet() ---- */\n"
      "QWidget { background-color: #2d2d2d; color: #f0f0f0; }\n"
      // The bare QWidget rule above overrides every bare type rule in
      // style.qss (equal specificity, appended later), so those are re-declared.
      // Without this, dialogs' default buttons turned invisible.
      "QPushButton { background: #457c92; color: white;"
      " border: 1px solid transparent; }\n"
      "QPushButton:hover { background: #3c6d80; }\n"
      "QPushButton:pressed { background: #34606f; }\n"
      "QPushButton:focus { outline: none; border-color: #6fb3c9; }\n"
      // style.qss's light grey would look more prominent than an enabled button.
      "QPushButton:disabled { background: #3a3a3a; color: #7a7a7a; }\n"
      // Bare type selectors again, as for QPushButton.
      "QToolTip { background-color: #1e1e1e; color: #f0f0f0;"
      " border: 1px solid #555555; }\n"
      "QScrollBar:vertical, QScrollBar:horizontal { background: transparent; }\n"
      "QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {\n"
      "  background: #333333; color: #f0f0f0;\n"
      "  border-top: 1px solid #262626; border-left: 1px solid #262626;\n"
      "  border-bottom: 1px solid #444444; border-right: 1px solid #444444;\n"
      "  selection-background-color: #457c92; }\n"
      "QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {\n"
      "  border: 1px solid #6fb3c9; }\n"
      "QComboBox QAbstractItemView {\n"
      "  background: #383838; color: #f0f0f0; border: 1px solid #555555;\n"
      "  selection-background-color: #457c92; selection-color: #ffffff; }\n"
      // Otherwise QWidget background with palette text: unreadable in one theme.
      "QMenu { background-color: #383838; color: #f0f0f0; border: 1px solid #555555; }\n"
      "QMenu::item:selected { background-color: #457c92; color: #ffffff; }\n"
      "QMenu::item:disabled { color: #787878; }\n"
      "QMenu::separator { background: #4f4f4f; }\n"
      // Same rules as style.qss's ::indicator (see there), dark colours.
      "QCheckBox::indicator, QRadioButton::indicator {\n"
      "  width: 14px; height: 14px;\n"
      "  background-color: #333333; border: 1px solid #5f5f5f; }\n"
      "QCheckBox::indicator { border-radius: 3px; }\n"
      "QRadioButton::indicator { border-radius: 8px; }\n"
      "QCheckBox::indicator:unchecked:hover, QRadioButton::indicator:unchecked:hover {\n"
      "  border: 1px solid #8a8a8a; }\n"
      "QCheckBox::indicator:checked, QRadioButton::indicator:checked {\n"
      "  background-color: #c8c8c8; border: 1px solid #c8c8c8; }\n"
      "QCheckBox::indicator:checked:hover, QRadioButton::indicator:checked:hover {\n"
      "  background-color: #e2e2e2; border: 1px solid #e2e2e2; }\n"
      "QCheckBox::indicator:disabled, QRadioButton::indicator:disabled {\n"
      "  background-color: #2f2f2f; border: 1px solid #474747; }\n"
      "QCheckBox::indicator:checked:disabled, QRadioButton::indicator:checked:disabled {\n"
      "  background-color: #6a6a6a; border: 1px solid #6a6a6a; }\n"
      "QScrollBar::handle:vertical, QScrollBar::handle:horizontal { background: #575757; }\n"
      "QScrollBar::handle:vertical:hover, QScrollBar::handle:horizontal:hover { background: #6f6f6f; }\n"
      ;
  return strBase + QString::fromLatin1(kDarkCss);
}

void AppTheme::Apply(bool bDark)
{
  if (!qApp)
    return;
  qApp->setPalette(BuildPalette(bDark));
  qApp->setStyleSheet(BaseStyleSheet(bDark));
}

QString AppTheme::ToggleStyleSheet(bool bDark)
{
  // Every state is explicit: an unstyled :checked button falls back to
  // Fusion's Highlight fill (a blue box). The border is the input-box bevel.
  if (bDark)
  {
    return "QPushButton {"
           "  background: #3a3a3a; color: #d8d8d8;"
           "  font-size: 13px; font-weight: bold; padding: 5px 16px;"
           "  border-radius: 13px;"
           "  border-top: 1px solid #272727; border-left: 1px solid #272727;"
           "  border-bottom: 1px solid #4e4e4e; border-right: 1px solid #4e4e4e;}"
           "QPushButton:hover {background: #454545; color: #ffffff;}"
           "QPushButton:pressed {background: #303030;}"
           "QPushButton:checked {"
           "  background: #457c92; color: #ffffff;"
           "  border-top: 1px solid #2c4c58; border-left: 1px solid #2c4c58;"
           "  border-bottom: 1px solid #6f9fb2; border-right: 1px solid #6f9fb2;}"
           "QPushButton:checked:hover {background: #4d8aa2;}"
           "QPushButton:checked:pressed {background: #3c6d80;}";
  }
  return "QPushButton {"
         "  background: #e4e4e4; color: #4a4a4a;"
         "  font-size: 13px; font-weight: bold; padding: 5px 16px;"
         "  border-radius: 13px;"
         "  border-top: 1px solid #c9c9c9; border-left: 1px solid #c9c9c9;"
         "  border-bottom: 1px solid #fbfbfb; border-right: 1px solid #fbfbfb;}"
         "QPushButton:hover {background: #dcdcdc; color: #262626;}"
         "QPushButton:pressed {background: #d2d2d2;}"
         "QPushButton:checked {"
         "  background: #457c92; color: #ffffff;"
         "  border-top: 1px solid #2c4c58; border-left: 1px solid #2c4c58;"
         "  border-bottom: 1px solid #6f9fb2; border-right: 1px solid #6f9fb2;}"
         "QPushButton:checked:hover {background: #3c6d80;}"
         "QPushButton:checked:pressed {background: #34606f;}";
}

QPixmap AppTheme::TintedPixmap(const QString& strName, const QColor& color, bool bMirror)
{
  QPixmap pix(QString(":/resource/icons/%1@2x.png").arg(strName));
  if (pix.isNull())
    return QPixmap();   // callers degrade to text-only

  // Qt tags an "@2x" file with DPR 2; reset so the composite is 1:1.
  pix.setDevicePixelRatio(1.0);
  if (bMirror)
    pix = pix.transformed(QTransform().scale(-1, 1), Qt::SmoothTransformation);

  QPixmap out(pix.size());
  out.fill(Qt::transparent);
  QPainter p(&out);
  p.drawPixmap(0, 0, pix);
  p.setCompositionMode(QPainter::CompositionMode_SourceIn);
  p.fillRect(out.rect(), color);
  p.end();
  out.setDevicePixelRatio(2.0);
  return out;
}

QIcon AppTheme::TintedIcon(const QString& strName, const QColor& color, bool bMirror)
{
  const QPixmap pix = TintedPixmap(strName, color, bMirror);
  return pix.isNull() ? QIcon() : QIcon(pix);
}
