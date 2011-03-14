/**
 * @file  LayerPropertyROI.h
 * @brief The common properties available to MRI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Kevin Teich
 * Reimplemented by: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2007-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef LayerPropertyROI_h
#define LayerPropertyROI_h

#include "vtkSmartPointer.h"
#include "LayerProperty.h"
#include <QColor>

extern "C"
{
#include "colortab.h"
}

class vtkRGBAColorTransferFunction;

class LayerPropertyROI : public LayerProperty
{
  Q_OBJECT
public:
  LayerPropertyROI ( QObject* parent = NULL );
  ~LayerPropertyROI ();

  double GetOpacity() const;

  double* GetColor()
  {
    return mRGB;
  }
  void SetColor( double r, double g, double b );

  vtkRGBAColorTransferFunction* GetLookupTable() const;

public slots:
  void SetOpacity( double opacity );
  void SetColor( const QColor& c )
  {
      SetColor( c.redF(), c.greenF(), c.blueF() );
  }

signals:
  void OpacityChanged( double opacity );
  void ColorMapChanged();

private:
  void SetColorMapChanged();

  //BTX

  // Color tables --------------------------------------------------------
  vtkSmartPointer<vtkRGBAColorTransferFunction> mLUTTable;
  // ---------------------------------------------------------------------

  double mOpacity;
  double mRGB[3];
};

#endif
