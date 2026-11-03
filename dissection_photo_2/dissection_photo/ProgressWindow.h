#ifndef PROGRESSWINDOW_H
#define PROGRESSWINDOW_H

#include <QFrame>
#include <QString>

namespace Ui {
class ProgressWindow;
}

class ProgressWindow : public QFrame
{
  Q_OBJECT

public:
  explicit ProgressWindow(QWidget *parent = nullptr);
  ~ProgressWindow();

  // Shows msg over a busy (indeterminate) progress bar.
  void SetIndeterminate(const QString& msg);

private:
  Ui::ProgressWindow *ui;
};

#endif // PROGRESSWINDOW_H
