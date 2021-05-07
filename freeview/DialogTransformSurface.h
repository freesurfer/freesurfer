#ifndef DIALOGTRANSFORMSURFACE_H
#define DIALOGTRANSFORMSURFACE_H

#include <QDialog>
#include <QList>
#include "UIUpdateHelper.h"

class QLineEdit;
class QScrollBar;

namespace Ui {
class DialogTransformSurface;
}

class DialogTransformSurface : public QDialog, public UIUpdateHelper
{
  Q_OBJECT

public:
  explicit DialogTransformSurface(QWidget *parent = nullptr);
  ~DialogTransformSurface();

  void showEvent(QShowEvent *);

public slots:
  void UpdateUI();
  void OnScrollBarValueChanged(int n);
  void OnTextEdit(const QString& text);
  void OnButtonRestore();
  void OnButtonUndo();
  void OnButtonSaveTransform();

protected slots:
  void ResetUI(const QList<QWidget*>& excluded = QList<QWidget*>());

private:
  Ui::DialogTransformSurface *ui;
  QScrollBar*    m_scrollRotate[3];
  QLineEdit*     m_textRotate[3];
  QScrollBar*    m_scrollTranslate[3];
  QLineEdit*     m_textTranslate[3];
  QScrollBar*    m_scrollScale[3];
  QLineEdit*     m_textScale[3];
};

#endif // DIALOGTRANSFORMSURFACE_H
