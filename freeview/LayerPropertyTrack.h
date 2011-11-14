/**
 * @file  LayerPropertyTrack.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/11/14 16:30:24 $
 *    $Revision: 1.5 $
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
#ifndef LAYERPROPERTYTRACK_H
#define LAYERPROPERTYTRACK_H

#include "LayerProperty.h"
#include <QColor>

class LayerPropertyTrack : public LayerProperty
{
  Q_OBJECT
public:
  LayerPropertyTrack(QObject* parent = 0);

  QColor  m_color;
};

#endif // LAYERPROPERTYTRACK_H
