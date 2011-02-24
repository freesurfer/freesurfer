/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis       *
 *    Massachusetts General Hospital        *
 *    Building 149, 13th St.         *
 *    Charlestown, MA 02129         *
 *    {ch|rudolph}@nmr.mgh.harvard.edu      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/// \file help.h
///
/// \brief Brief description
/// The help menus / data.
///
/// \b DESCRIPTION
/// Internal help data.
///
/// \b HISTORY
/// 15 March 2005 - Initial consolidation from several other sources.
/// $Id: help.h,v 1.5 2011/02/24 21:14:30 rudolph Exp $
///
///

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

#include <string>
using namespace std;

extern string   G_VERSION;

static struct option const longopts[] = {
    {"optionsFile",     required_argument,      NULL, 'o'},
    {"dir",             required_argument,      NULL, 'D'},
    {"version",         no_argument,            NULL, 'v'},
    {"subject",         required_argument,      NULL, 'S'},
    {"hemi",            required_argument,      NULL, 'h'},
    {"surface0",        required_argument,      NULL, 's'},
    {"surface1",        required_argument,      NULL, 't'},
    {"curv0",           required_argument,      NULL, 'c'},
    {"curv1",           required_argument,      NULL, 'd'},
    {"mpmProg",         required_argument,      NULL, 'm'},
    {"mpmArgs",         required_argument,      NULL, 'M'},
    {"useAbsCurvs",     no_argument,            NULL, 'a'},
    {"mpmOverlay",	required_argument,	NULL, 'O'},
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


