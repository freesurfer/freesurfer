#ifndef APPTHEME_H
#define APPTHEME_H

#include <QColor>
#include <QIcon>
#include <QPalette>
#include <QPixmap>
#include <QString>

// Shared light/dark palettes and the active-mode decision, so every window
// agrees, including the dialogs that exist before MainWindow.
namespace AppTheme
{
// Queried each call (no cross-platform change signal before Qt 6.5).
bool IsSystemDarkMode();

// Persisted user override: "system" (default), "light" or "dark".
QString Mode();
void SetMode(const QString& mode);

// Mode() resolved against IsSystemDarkMode().
bool EffectiveDarkMode();

QPalette BuildPalette(bool bDark);

// resource/style.qss (light) plus, in dark mode, dark counterparts. A
// stylesheet beats a QPalette, and more specific selectors beat an appended
// QWidget rule, so every hardcoded light colour needs a dark counterpart.
QString BaseStyleSheet(bool bDark);

// Applies BuildPalette() + BaseStyleSheet() to qApp. Doesn't call SetMode(),
// so a "system" preference is never persisted as "light"/"dark".
void Apply(bool bDark);

// The light/dark pill toggle on the Projects Hub and Welcome dialog.
QString ToggleStyleSheet(bool bDark);

// Recolours a single-colour resource/icons glyph through its alpha. Loads the
// 48px @2x file tagged DPR 2 (automatic @2x lookup doesn't apply to a built
// pixmap). Null if the resource is missing. bMirror flips horizontally.
QPixmap TintedPixmap(const QString& strName, const QColor& color, bool bMirror = false);
QIcon TintedIcon(const QString& strName, const QColor& color, bool bMirror = false);
}

#endif // APPTHEME_H
