#ifndef BINARYTREEVIEW_H
#define BINARYTREEVIEW_H

#include <QGraphicsView>

#include <fstream>
#include <iostream>
#include <QVariantMap>
#include <QThread>

class BinaryTreeNode;
class BinaryTreeView;

class TreeDataLoader : public QObject
{
  Q_OBJECT
public:
  TreeDataLoader(QObject* parent = 0);

signals:
  void LoadRequested(const QString& dirName, BinaryTreeView* view);
  void DataLoaded(const QVariantMap& data);

public slots:
  void Load(const QString& dirName, BinaryTreeView* view)
  {
    emit LoadRequested(dirName, view);
  }

protected slots:
  void DoLoad(const QString& dirName, BinaryTreeView* view);
};

class BinaryTreeView : public QGraphicsView
{
  friend class TreeDataLoader;

  Q_OBJECT
public:
  BinaryTreeView(QWidget* parent = 0);
  ~BinaryTreeView();

  QVariantMap GetTractMapData()
  {
    return m_mapData;
  }

signals:
  void TreeNodeActivated(const QStringList& trk_files);
  void TreeDataLoaded(const QVariantMap& data);

public slots:
  void Load(const QString& dirName);
  void ZoomIn();
  void ZoomOut();
  void OnDataLoaded(const QVariantMap& data);
  void SetData(const QVariantMap& data);

protected:
  void keyPressEvent(QKeyEvent* e);
  void mousePressEvent(QMouseEvent* e);
  void wheelEvent(QWheelEvent* e);

  void drawBackground(QPainter *painter, const QRectF &rect);
  void ScaleView(qreal scaleFactor);

private:
  QVariantMap m_mapData;
  QMap<QString, BinaryTreeNode*> m_mapNode;
  BinaryTreeNode* m_selectedNode;
  BinaryTreeNode* m_nodeStart;

  QGraphicsScene* m_scene;
  TreeDataLoader* m_dataLoader;
  QThread   m_threadLoad;
  QString   m_strDataDir;
};

#endif // BINARYTREEVIEW_H
