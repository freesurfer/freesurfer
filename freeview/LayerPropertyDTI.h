/**
 * @file  LayerPropertyDTI.h
 * @brief Layer properties available to DTI layers
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
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
