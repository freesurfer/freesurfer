/**
 * @brief Utilitie functions for Random Forests on surfaces and in volumes.
 *
 */
/*
 * Original Author: Bruce Fischl
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
 */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cma.h"
#include "diag.h"
#include "error.h"
#include "gcamorph.h"
#include "macros.h"
#include "rfa.h"
#include "rfutils.h"
#include "talairachex.h"
#include "transform.h"
#include "utils.h"
