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
/// \file general.h
///
/// \brief Brief description
/// Some functions shared between the various dijkstra_p1 components.
///
/// \b DESCRIPTION
/// This file defines some "general" functions that are used by
/// several different dijkstra_p1 components.
///
///

#ifndef __GENERAL_H__
#define __GENERAL_H__

using namespace std;

#include <string.h>
#include "env.h"
#include "c_SMessage.h"
#include "pathconvert.h"

#define PENV  str_env.

#define ULOUT( msg )  if(st_env.pcsm_userlog) st_env.pcsm_userlog->dump(st_env.b_syslogPrepend, ( msg ) );
#define nULOUT( msg ) if(st_env.pcsm_userlog) st_env.pcsm_userlog->dump(false, ( msg ) );
#define SLOUT( msg )  if(st_env.pcsm_syslog)  st_env.pcsm_syslog->dump(st_env.b_syslogPrepend, ( msg ) );
#define nSLOUT( msg ) if(st_env.pcsm_syslog)  st_env.pcsm_syslog->dump(false, ( msg ) );
#define RLOUT( msg )  if(st_env.pcsm_resultlog)  st_env.pcsm_resultlog->dump(st_env.b_syslogPrepend, ( msg ) );
#define nRLOUT( msg ) if(st_env.pcsm_resultlog)  st_env.pcsm_resultlog->dump(false, ( msg ) );

#define Cn( msg ) cout << ( msg ) << endl;
#define CW(x, msg) cout.width((x))  ; cout << (msg);
#define CWn(x, msg) cout.width((x))  ; cout << (msg); cout << endl;
#define CCn(w, msg) for(unsigned __co=0; __co<(w-strlen(msg))/2; __co++)\
       cout << " ";    \
   cout << (msg); cout << endl;

typedef enum {
  e_gaussian, e_mean
} e_CURVATURE;

typedef enum _FILEACCESS {
        e_UNSPECIFIED           = -10,
        e_WRONGMAGICNUMBER      = -1,
        e_OK                    =  0,
        e_READACCESSERROR       =  1,
        e_WRITEACCESSERROR      =  2
} e_FILEACCESS;

short CURV_arrayProgress_print(
  int   asize,
  int   acurrent,
  char* apch_message
);

void    lprintf(int lw, const char* format, ...);
void    colprintf(int lw, int rw,
                    const char* pch_lstr, const char* format, ...);
char*   colsprintf(int lw, int rw, char* pch_buffer,
                    const char* pch_lstr, const char* format, ...);

/// A *very* simple structure to handle some recurring 3D vector operations.
/// It would make more sense to use a CMatrix class, but that is perhaps
/// too much of an overkill right now. Same argument for creating a
/// custom vector3D class.
typedef struct _vector3D {
  float  f_x;
  float  f_y;
  float  f_z;
}
st_V3D;

/// \fn V3D_normalizedDirection_find(st_V3D* pV_A, st_V3D* pV_B, st_V3D* pV_C)
/// \brief Calculate the normalized direction vector between A and B
/// \param  pV_A  The "start" point
/// \param  pV_B The "end" point
/// \param  pV_C The normalized direction vector from A to B
/// \return    (float) magnitude of the distance
float  V3D_normalizedDirection_find(
  st_V3D&  V_A,
  st_V3D&  V_B,
  st_V3D*  pV_C
);

/// \fn float  V3D_dot(st_V3D& V_A, st_v3D& V_B)
/// \brief Calculate the dot product between A and B
/// \param  pV_A  The first vector
/// \param  pV_B The second vector
/// \return    (float) dot product
float  V3D_dot( st_V3D&  V_A,
                st_V3D&  V_B);

/// \fn float  V3D_distance(st_V3D& V_A, st_v3D& V_B)
/// \brief Calculate the distance between A and B
/// \param  pV_A  The first vector
/// \param  pV_B The second vector
/// \return    (float) dot product
float  V3D_distance( st_V3D&  V_A,
                     st_V3D&  V_B);

/// \fn bool str_3parse(string& astr_input, string& astr_first, string& astr_second, string& astr_third)
/// \brief Break (on whitespace) an input string into three substrings. Any remainder is returned in the input.
/// \param astr_input  Source string to parse. Destructive parsing.
/// \param astr_first  The first token component.
/// \param astr_second  The second token component.
/// \param astr_third  The third token component.
/// \return int   Number of substrings successfully extracted.
int
str_3parse(
  string&  astr_input,
  string&  astr_first,
  string&  astr_second,
  string&  astr_third
);

/// \fn bool str_leadingWSremove(string& astr_input)
/// \brief Remove (destructively) any leading whitespace in passed string.
/// \param astr_input  Source string to parse. Destructive parsing.
/// \return bool   OK
bool
str_leadingWSremove(
  string&  astr_input
);

/// \fn bool relDirSpec_test(string& astr_dirSpec)
/// \brief Checks if the current directory specification is relative or absolute.
/// \param astr_dirSpec  An input string containing a directory (or file) spec
/// \return bool   TRUE - if astr_dirSpec is a relative specification
bool
relDirSpec_test(
  string&  astr_dirSpec
);

/// \fn bool str_rel2absDirSpec_change(string& astr_rel, string& astr_abs)
/// \brief Converts the relative dir spec in astr_rel to an absolute directory in astr_abs.
/// \param astr_rel  An input string containing a relative directory
/// \param astr_abs  An output string containing the input as an absolute dir spec
/// \return bool   OK
bool
str_rel2absDirSpec_change(
  string&  astr_rel,
  string&  astr_abs
);


void warn(
  string          str_action,
  string          str_errorMsg,
  int  errorCode
);

void
error_exit(
  string          str_action,
  string          str_errorMsg,
  int  errorCode
);

void
command_line_error(
  char   *fmt, ...
);


#endif //__GENERAL_H__

