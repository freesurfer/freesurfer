#ifndef DIALOGMOVEPOINT_H
#define DIALOGMOVEPOINT_H

#include <QDialog>

namespace Ui {
class DialogMovePoint;
}

class LayerPointSet;
class RenderView2D;

class DialogMovePoint : public QDialog
{
  Q_OBJECT

public:
  explicit DialogMovePoint(QWidget *parent = nullptr);
  ~DialogMovePoint();

  void SetData(RenderView2D* view, LayerPointSet* wp, int nId);

  double GetSigma();
  double GetNeighborSize();

public slots:
  void OnButtonTest();
  void OnButtonRestore();
  void OnPointSetPicked(LayerPointSet* wp, int n);

private:
  Ui::DialogMovePoint *ui;

  RenderView2D*   m_view;
  LayerPointSet*  m_layerPointSet;
  int     m_nIndex;
  double  m_dPrevPos[3];
  int     m_nPlane;
  bool    m_bUndoable;
};

#endif // DIALOGMOVEPOINT_H
