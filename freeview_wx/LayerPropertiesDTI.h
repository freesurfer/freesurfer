/**
 * @file  LayerPropertiesDTI.h
 * @brief Layer properties available to DTI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
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

#ifndef LayerPropertiesDTI_h
#define LayerPropertiesDTI_h

#include "vtkSmartPointer.h"
#include "LayerPropertiesMRI.h"

class vtkLookupTable;

class LayerPropertiesDTI : public LayerPropertiesMRI
{
public:
  LayerPropertiesDTI ();
  ~LayerPropertiesDTI ();

  enum DirectionCode { RGB = 0, RBG, GRB, GBR, BRG, BGR };

  vtkLookupTable* GetDirectionCodedTable () const;

  DirectionCode GetDirectionCode()
  {
    return m_nDirectionCode;
  }

  void SetDirectionCode( DirectionCode nCode );

  void ColorMapChanged();

private:
  DirectionCode m_nDirectionCode;
  vtkSmartPointer<vtkLookupTable> mDirectionCodedTable;
};

#endif
