/**
 * @file  LayerVolumeTracks.h
 * @brief Layer class for tracks saved in a multi-frame volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/21 21:27:40 $
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
 *
 */

#ifndef LAYERVOLUMETRACKS_H
#define LAYERVOLUMETRACKS_H

#include "LayerMRI.h"

class LayerVolumeTracks : public LayerMRI
{
public:
    LayerVolumeTracks( LayerMRI* ref, QObject* parent = NULL );
};

#endif // LAYERVOLUMETRACKS_H
