/**
 * @brief Layer properties available to DTI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef LayerPropertyDTI_h
#define LayerPropertyDTI_h

#include "vtkSmartPointer.h"
#include "LayerPropertyMRI.h"

class vtkLookupTable;

class LayerPropertyDTI : public LayerPropertyMRI
{
  Q_OBJECT
public:
  LayerPropertyDTI ( QObject* parent = NULL );
  ~LayerPropertyDTI ();

  enum DirectionCode { RGB = 0, RBG, GRB, GBR, BRG, BGR };

  vtkLookupTable* GetDirectionCodedTable () const;

  int GetDirectionCode()
  {
    return m_nDirectionCode;
  }

  void OnColorMapChanged();

public slots:
  void SetDirectionCode( int nCode );

private:
  int m_nDirectionCode;
  vtkSmartPointer<vtkLookupTable> mDirectionCodedTable;
};

#endif
