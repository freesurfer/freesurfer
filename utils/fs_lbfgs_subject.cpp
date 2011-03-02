/**
 * @file  fs_lbfgs_subject.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author$
 *    $Date$
 *    $Revision$
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


#include "fs_vnl/fs_lbfgs_subject.h"

fs_lbfgs_subject::fs_lbfgs_subject()
{
  mObserver = NULL;
}

fs_lbfgs_subject::~fs_lbfgs_subject()
{}

void fs_lbfgs_subject::setObserver( fs_lbfgs_observer* observer )
{
  mObserver = observer;
}

void fs_lbfgs_subject::notify( double bestF, vnl_vector< double >* bestX )
{
  // we might want a number of observers, but right now we only have one
  if ( mObserver != NULL )
  {
    mObserver->update( bestF, bestX );
  }
}
