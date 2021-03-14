/**
 * @brief help text utils
 *
 */
/*
 * Original Author: Rudolph Pienaar / Christian Haselgrove
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

#ifndef __HELP_H__
#define __HELP_H__

#include <getopt.h>

#include "env.h"
#include "general.h"

#include "error.h"
#include "label.h"
#include "mri.h"
#include "mrisurf.h"

#include <string>
#include <unistd.h>

extern std::string G_VERSION;

static struct option const longopts[] = {
    {"optionsFile", required_argument, nullptr, 'o'},
    {"dir", required_argument, nullptr, 'D'},
    {"version", no_argument, nullptr, 'v'},
    {"subject", required_argument, nullptr, 'S'},
    {"hemi", required_argument, nullptr, 'h'},
    {"surface", required_argument, nullptr, 's'},
    {"surface1", required_argument, nullptr, 't'},
    {"curv", required_argument, nullptr, 'c'},
    {"curv1", required_argument, nullptr, 'd'},
    {"mpmProg", required_argument, nullptr, 'm'},
    {"mpmArgs", required_argument, nullptr, 'M'},
    {"port", required_argument, nullptr, 'p'},
    {"useAbsCurvs", no_argument, nullptr, 'a'},
    {"mpmOverlay", required_argument, nullptr, 'O'},
    {"mpmOverlayArgs", required_argument, nullptr, 'V'},
    {nullptr, 0, nullptr, 0}};

std::string commandLineOptions_process(int argc, char **ppch_argv,
                                       s_env &st_env);

void synopsis_show();

void asynchEvent_processHELP(s_env &st_env, std::string str_event);

#endif //__HELP_H__
