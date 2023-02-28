/*=========================================================================

  Copyright 2004 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/

/*========================================================================
 For general information about using VTK and Qt, see:
 http://www.trolltech.com/products/3rdparty/vtksupport.html
=========================================================================*/

/*========================================================================
 !!! WARNING for those who want to contribute code to this file.
 !!! If you use a commercial edition of Qt, you can modify this code.
 !!! If you use an open source version of Qt, you are free to modify
 !!! and use this code within the guidelines of the GPL license.
 !!! Unfortunately, you cannot contribute the changes back into this
 !!! file.  Doing so creates a conflict between the GPL and BSD-like VTK
 !!! license.
=========================================================================*/

// .NAME QVTKPaintEngine - directs QPainter calls to a VTK window


#include "QVTKPaintEngine.h"
#include "QVTKWidget.h"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkImageData.h"
#include "vtkSmartPointer.h"

#include <QCache>
#include <QPainterPath>

class QVTKPaintEngineInternal
{
public:
  // cache of pixmaps
  QCache<qint64, vtkSmartPointer<vtkImageData> > mImageCache;
};

QVTKPaintEngine::QVTKPaintEngine()
  : QPaintEngine(QPaintEngine::PaintOutsidePaintEvent |
                 QPaintEngine::AlphaBlend)
{
  this->Internal = new QVTKPaintEngineInternal;
}

QVTKPaintEngine::~QVTKPaintEngine()
{
  delete this->Internal;
}

bool QVTKPaintEngine::begin(QPaintDevice* dev)
{
  Widget = static_cast<QVTKWidget*>(dev);
  return true;
}

bool QVTKPaintEngine::end()
{
  Widget = NULL;
  return true;
}

QPaintEngine::Type QVTKPaintEngine::type() const
{
  return QPaintEngine::User;
}

void QVTKPaintEngine::updateState(const QPaintEngineState&)
{
}

// at a minimum, we only need to re-implement this drawPixmap function.
// Qt can do all other drawing to create a pixmap and then we draw it here.
void QVTKPaintEngine::drawPixmap(const QRectF& r, const QPixmap& pm, const QRectF& sr)
{
  if(!this->Widget)
  {
    return;
  }
  QRect ri = r.toRect();
  QRect sri = sr.toRect();

  QPixmap pix = pm.copy(sri);
  if(sri.size() != ri.size())
  {
    pix = pix.scaled(ri.size());
  }

  QImage img = pix.toImage().mirrored().rgbSwapped();

  // blend the pixels from QImage into the vtkRenderWindow's buffer
  vtkRenderWindow* renWin = this->Widget->GetRenderWindow();
  if (renWin)
    renWin->SetRGBACharPixelData(ri.left(), this->Widget->height() - ri.top() - ri.height(),
                               ri.left()+img.width() - 1,
                               this->Widget->height() - ri.top() - 1,
                               img.bits(),
                               renWin->GetDoubleBuffer() ? 0 : 1,
                               1);

  // NOTE: this would perform much better if textures were used and caching of those
  // textures was done (probably vtkActor2D and vtkImageMapper)
}

void QVTKPaintEngine::drawPath(const QPainterPath& path)
{
  // drawPath in base class does nothing so here we make it do something
  QRectF box = path.boundingRect();
  QPixmap img((int)box.width(), (int)box.height());
  img.fill(Qt::transparent);
  QPainter p(&img);
  p.translate(-box.left(), -box.right());
  p.drawPath(path);
  this->drawPixmap(QRectF(QPoint(0,0), img.size()), img, box);
}

