#include "BinaryTreeView.h"
#include "BinaryTreeEdge.h"
#include "BinaryTreeNode.h"
#include <QMouseEvent>
#include <QDebug>
#include <QDir>
#include "math.h"
#include <QRegularExpression>

TreeDataLoader::TreeDataLoader(QObject *parent) : QObject(parent)
{
  connect(this, SIGNAL(LoadRequested(QString,BinaryTreeView*)), SLOT(DoLoad(QString,BinaryTreeView*)));
}

void TreeDataLoader::DoLoad(const QString &dirPath, BinaryTreeView *view)
{
  int max_level;
  max_level = 0;
  //1) Make a dictionary of tract values (set in header)
  //find tract directory
  QDir mDir(dirPath);
  mDir.setNameFilters(QStringList()<<"*.trk");
  mDir.setSorting(QDir::Name);

  QFileInfoList info_list = mDir.entryInfoList();
  QStringList tract_list;
  foreach (QFileInfo fn, info_list)
  {
    if (!fn.baseName().contains(QRegularExpression("[^0-1]")))
      tract_list << fn.fileName();
  }

  //read tract files
  QVariantMap tract_map;
  for (int i = 0; i <tract_list.size();i++)
  {
    QString current_key;
    current_key = tract_list.at(i);
    //remove ".trk"
    current_key.remove(-4,4);
    //determine the max tree level
    if (current_key.size() > max_level){
      max_level = current_key.size();
    }

    //read each character of tract
    //populate dictionary with tract lists for each key
    for (int j = current_key.size()+1; (j--)>1;)
    {
      QString current_tract;
      current_tract = current_key.left(j);
      if (!tract_map.contains(current_tract))
      {
        //create list of associated tracts (keys) for each tract
        QStringList list_of_children;

        //find all children and insert into tract map
        for (int k = 0; k <tract_list.size();k++)
        {
          QString possible_child = tract_list.at(k);
          possible_child.remove(-4,4);
          if (current_tract.size()<=possible_child.size())
          {
            if (possible_child.startsWith(current_tract))
            {
              //qDebug() << "current key = "+ current_key;
              list_of_children.append(possible_child);
            }
          }
        }
        tract_map.insert(current_tract,list_of_children);
      }
    }
    //qDebug() << "Filename " + QString::number(i) + " = " + list.at(i);
  }

  QVariantMap data;
  if (!tract_map.isEmpty())
  {
    data["tree_map"] = tract_map;
    data["max_level"] = max_level;
    data["path"] = dirPath;
  }
  emit DataLoaded(data);
}

BinaryTreeView::BinaryTreeView(QWidget* parent) : QGraphicsView(parent), m_selectedNode(NULL)
{
  setMouseTracking(true);
  m_scene = new QGraphicsScene(this);
  m_scene->setItemIndexMethod(QGraphicsScene::NoIndex);
  setScene(m_scene);
  setCacheMode(CacheBackground);
  setViewportUpdateMode(BoundingRectViewportUpdate);
  setRenderHint(QPainter::Antialiasing);
  setTransformationAnchor(AnchorUnderMouse);
  scale(qreal(0.5), qreal(0.5));

  m_dataLoader = new TreeDataLoader;
  m_dataLoader->moveToThread(&m_threadLoad);
  m_threadLoad.start();
  connect(m_dataLoader, SIGNAL(DataLoaded(QVariantMap)), SLOT(OnDataLoaded(QVariantMap)), Qt::QueuedConnection);
}

BinaryTreeView::~BinaryTreeView()
{
  m_threadLoad.quit();
  m_threadLoad.wait();
  m_dataLoader->deleteLater();
}

void BinaryTreeView::Load(const QString &dirName)
{
  m_dataLoader->Load(dirName, this);
}

void BinaryTreeView::mousePressEvent(QMouseEvent *event)
{
  QGraphicsItem* gitem = itemAt(event->pos());
  if (!gitem || gitem->type() != (QGraphicsItem::UserType+1))
    return;
  BinaryTreeNode *item = (BinaryTreeNode*)gitem; //Get the node at the position
  if (item) //if there is a node at that position
  {
    QString tract_name = m_mapNode.key(item);
    QStringList list = m_mapData.value("tree_map").toMap().value(tract_name).toStringList();
    qDebug() << "\nNode" << tract_name; //print the tract name for the node
    qDebug() << "Leaves" << list; //print list of all leaves tracts to show for this node
    if (!list.isEmpty())
    {
      QStringList filenames;
      foreach (QString fn, list)
      {
        filenames << QFileInfo(QDir(m_strDataDir), fn + ".trk").absoluteFilePath();
      }
      emit TreeNodeActivated(filenames);
    }
    if (m_selectedNode)
      m_selectedNode->SetHighlighted(false, true);
    m_selectedNode = item;
    item->SetHighlighted(true, true);
  }
  QGraphicsView::mousePressEvent(event);
}

void BinaryTreeView::keyPressEvent(QKeyEvent *event)
{
  switch (event->key()) {
  case Qt::Key_Up:
    m_nodeStart->moveBy(0, -20);
    break;
  case Qt::Key_Down:
    m_nodeStart->moveBy(0, 20);
    break;
  case Qt::Key_Left:
    m_nodeStart->moveBy(-20, 0);
    break;
  case Qt::Key_Right:
    m_nodeStart->moveBy(20, 0);
    break;
  case Qt::Key_Plus:
    ZoomIn();
    break;
  case Qt::Key_Minus:
    ZoomOut();
    break;
  case Qt::Key_Space:
  case Qt::Key_Enter:
    break;
  default:
    QGraphicsView::keyPressEvent(event);
  }
}

void BinaryTreeView::wheelEvent(QWheelEvent *event)
{
  ScaleView(pow((double)2, -event->angleDelta().y() / 240.0));
}

void BinaryTreeView::drawBackground(QPainter *painter, const QRectF &rect)
{
  Q_UNUSED(rect);

//  // Shadow
//  QRectF sceneRect = this->sceneRect();
//  QRectF rightShadow(sceneRect.right(), sceneRect.top() + 5, 5, sceneRect.height());
//  QRectF bottomShadow(sceneRect.left() + 5, sceneRect.bottom(), sceneRect.width(), 5);
//  if (rightShadow.intersects(rect) || rightShadow.contains(rect))
//    painter->fillRect(rightShadow, Qt::darkGray);
//  if (bottomShadow.intersects(rect) || bottomShadow.contains(rect))
//    painter->fillRect(bottomShadow, Qt::darkGray);

//  // Fill
//  QLinearGradient gradient(sceneRect.topLeft(), sceneRect.bottomRight());
//  gradient.setColorAt(0, Qt::white);
//  gradient.setColorAt(1, Qt::lightGray);
//  painter->fillRect(rect.intersected(sceneRect), gradient);
//  painter->setBrush(Qt::NoBrush);
//  painter->drawRect(sceneRect);

//  // Text
//  QRectF textRect(sceneRect.left() + 4, sceneRect.top() + 4,
//                  sceneRect.width() - 4, sceneRect.height() - 4);
//  QString message(tr(""));

//  QFont font = painter->font();
//  font.setBold(true);
//  font.setPointSize(14);
//  painter->setFont(font);
//  painter->setPen(Qt::lightGray);
//  painter->drawText(textRect.translated(2, 2), message);
//  painter->setPen(Qt::black);
//  painter->drawText(textRect, message);
}

void BinaryTreeView::ScaleView(qreal scaleFactor)
{
  qreal factor = transform().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
  if (factor < 0.07 || factor > 100)
    return;

  scale(scaleFactor, scaleFactor);
}

void BinaryTreeView::ZoomIn()
{
  ScaleView(qreal(1.2));
}

void BinaryTreeView::ZoomOut()
{
  ScaleView(1 / qreal(1.2));
}

void BinaryTreeView::OnDataLoaded(const QVariantMap& data)
{
  if (!data.isEmpty())
    SetData(data);
  QStringList list = data.value("tree_map").toMap().value("1").toStringList();
  QStringList filenames;
  if (!list.isEmpty())
  {
    foreach (QString fn, list)
    {
      filenames << QFileInfo(QDir(m_strDataDir), fn + ".trk").absoluteFilePath();
    }
  }
  QVariantMap map = m_mapData;
  map["filenames"] = filenames;
  emit TreeDataLoaded(map);
}

void BinaryTreeView::SetData(const QVariantMap &data)
{
  m_selectedNode = NULL;
  m_mapData = data;
  if (m_mapData.contains("path"))
    m_strDataDir = m_mapData["path"].toString();
  QVariantMap mapTract = data["tree_map"].toMap();
  int max_level = data["max_level"].toInt();

  double max_dist;
  //summation 2^n is 2^(n)-1
  max_dist = 6*(pow(2,max_level-1)-1);

  m_scene->clear();
  //node map
  //create all nodes based on tract_map
  QStringList keys = mapTract.keys();
  for (int i=0; i< keys.size(); i++)
  {
    QString key;
    key = keys[i];
    //assign new node to node map
    m_mapNode[key] = new BinaryTreeNode(this);
    //parent node
    if (key=="1")
    {
      m_scene->addItem(m_mapNode[key]);
      m_mapNode[key]->setPos(0,-800);
    }
    else
    {
      double dist;
      dist = 0;
      m_scene->addItem(m_mapNode[key]);

      //set distance of node based on label
      for (int j = 2; j<=key.size(); j++)
      {
        if (key.left(j).endsWith("1"))
        {
          dist -=6*pow(2,(max_level-1-key.left(j).size()));
        }
        else
        {
          dist +=6*pow(2,(max_level-1-key.left(j).size()));
        }
      }
      m_mapNode[key] -> setPos(dist,(-800+(key.size()-1)*100));

      //assign to node pos map
      //node_pos_map[key] = QPoint(dist,(-800+(key.size()-1)*100));

      //add edges to nodes
      if (key.size() > 1)
      {
        m_scene->addItem(new BinaryTreeEdge(m_mapNode[key.left(key.size()-1)], m_mapNode[key]));
      }
    }
  }
}

