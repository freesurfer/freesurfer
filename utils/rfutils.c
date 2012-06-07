/**
 * @file  rfutils.c
 * @brief Utilitie functions for Random Forests on surfaces and in volumes.
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/06/07 12:57:46 $
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "diag.h"
#include "error.h"
#include "rfa.h"
#include "cma.h"
#include "transform.h"
#include "talairachex.h"
#include "gcamorph.h"
#include "utils.h"
#include "rfutils.h"
#include "macros.h"

