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
 * Copyright (C) 2002-2007,
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
