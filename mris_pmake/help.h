/**
 * @file  help.h
 * @brief help text utils
 *
 */
/*
 * Original Author: Rudolph Pienaar / Christian Haselgrove
 * CVS Revision Info:
 *    $Author: rudolph $
 *    $Date: 2012/10/19 19:01:35 $
 *    $Revision: 1.9 $
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


#ifndef __HELP_H__
#define __HELP_H__

#include <getopt.h>


#include "general.h"
#include "env.h"

#ifdef __cplusplus
extern  "C" {
#endif

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"

#ifdef __cplusplus
}
#endif
#include <unistd.h>
#include <string>
using namespace std;

extern string   G_VERSION;

static struct option const longopts[] = {
    {"optionsFile",     required_argument,      NULL, 'o'},
    {"dir",             required_argument,      NULL, 'D'},
    {"version",         no_argument,            NULL, 'v'},
    {"subject",         required_argument,      NULL, 'S'},
    {"hemi",            required_argument,      NULL, 'h'},
    {"surface",         required_argument,      NULL, 's'},
    {"surface1",        required_argument,      NULL, 't'},
    {"curv",            required_argument,      NULL, 'c'},
    {"curv1",           required_argument,      NULL, 'd'},
    {"mpmProg",         required_argument,      NULL, 'm'},
    {"mpmArgs",         required_argument,      NULL, 'M'},
    {"port",            required_argument,      NULL, 'p'},
    {"useAbsCurvs",     no_argument,            NULL, 'a'},
    {"mpmOverlay",	required_argument,	NULL, 'O'},
    {"mpmOverlayArgs",	required_argument,	NULL, 'V'},
    {NULL, 0, NULL, 0}
};

string
commandLineOptions_process(
    int         argc,
    char**      ppch_argv,
    s_env&      st_env                           
);

void synopsis_show(void);

void asynchEvent_processHELP(
  s_env&   st_env,
  string   str_event
);

#endif //__HELP_H__


