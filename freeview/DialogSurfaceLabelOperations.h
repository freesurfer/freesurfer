#ifndef DIALOGSURFACELABELOPERATIONS_H
#define DIALOGSURFACELABELOPERATIONS_H

#include <QDialog>

namespace Ui {
class DialogSurfaceLabelOperations;
}

class DialogSurfaceLabelOperations : public QDialog
{
  Q_OBJECT

public:
  explicit DialogSurfaceLabelOperations(QWidget *parent = 0);
  ~DialogSurfaceLabelOperations();

signals:
  void OperationTriggered(const QVariantMap& op);

public slots:
  void OnButtonClicked();

private:
  Ui::DialogSurfaceLabelOperations *ui;
};

#endif // DIALOGSURFACELABELOPERATIONS_H
