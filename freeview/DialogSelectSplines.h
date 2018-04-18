#ifndef DIALOGSELECTSPLINES_H
#define DIALOGSELECTSPLINES_H

#include <QDialog>
#include <QList>

namespace Ui {
class DialogSelectSplines;
}

class LayerPointSet;
class Layer;

class DialogSelectSplines : public QDialog
{
  Q_OBJECT

public:
  explicit DialogSelectSplines(QWidget *parent = 0);
  ~DialogSelectSplines();

  void SetPointSets(const QList<Layer*>& list);

  QList<LayerPointSet*> GetSelectedPointSets();

public slots:
  void OnButtonUp();
  void OnButtonDown();
  void OnButtonOK();

private:
  Ui::DialogSelectSplines *ui;
};

#endif // DIALOGSELECTSPLINES_H
