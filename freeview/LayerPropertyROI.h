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
 *    $Author: nicks $
 *    $Date: 2011/03/14 23:44:47 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
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
