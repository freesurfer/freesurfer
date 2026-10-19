#ifndef DIALOGWELCOME_H
#define DIALOGWELCOME_H

#include <QDialog>

namespace Ui {
class DialogWelcome;
}

class QPushButton;

class DialogWelcome : public QDialog
{
  Q_OBJECT

public:
  explicit DialogWelcome(QWidget *parent = nullptr);
  ~DialogWelcome();

  bool GetDoNotShow();

private slots:
  void OnToggleDarkMode(bool bDark);

protected:
  // Wrapped QLabels report sizeHint for a default width, and a plain QWidget
  // doesn't propagate heightForWidth, so rows clipped their text. Recompute
  // each label's minimum height from its real width on resize.
  void resizeEvent(QResizeEvent* e) override;
  void showEvent(QShowEvent* e) override;

private:
  void UpdateBodyHeights();

  // Re-applies this dialog's widget-level colours for the mode.
  void Retheme(bool bDark);

  Ui::DialogWelcome *ui;
  QPushButton* m_btnDarkMode;
};

#endif // DIALOGWELCOME_H
