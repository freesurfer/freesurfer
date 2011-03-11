/**
 * @file  LayerPropertiesROI.h
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
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
 *    $Revision: 1.1 $
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
 */

#ifndef LayerPropertiesROI_h
#define LayerPropertiesROI_h

#include "vtkSmartPointer.h"
#include "LayerProperties.h"

extern "C"
{
#include "colortab.h"
}

class vtkRGBAColorTransferFunction;

class LayerPropertiesROI : public LayerProperties
{

public:
  LayerPropertiesROI ();
  ~LayerPropertiesROI ();

  double GetOpacity() const;
  void SetOpacity( double opacity );

  double* GetColor()
  {
    return mRGB;
  }
  void SetColor( double r, double g, double b );

  vtkRGBAColorTransferFunction* GetLookupTable() const;

private:
  void ColorMapChanged ();

  //BTX

  // Color tables --------------------------------------------------------
  vtkSmartPointer<vtkRGBAColorTransferFunction> mLUTTable;
  // ---------------------------------------------------------------------

  double mOpacity;
  double mRGB[3];
};

#endif
