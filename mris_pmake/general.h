/**
 * @brief Some functions shared between the various dijkstra_p1 components.
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


#ifndef __GENERAL_H__
#define __GENERAL_H__

#include <string.h>
#include "env.h"
#include "c_SMessage.h"
#include "pathconvert.h"

using namespace std;

#define PENV  str_env.

#define ULOUT( msg )  if(st_env.pcsm_userlog)        st_env.pcsm_userlog->dump(st_env.b_syslogPrepend, ( msg ) );
#define nULOUT( msg ) if(st_env.pcsm_userlog)        st_env.pcsm_userlog->dump(false, ( msg ) );
#define SLOUT( msg )  if(st_env.pcsm_syslog)         st_env.pcsm_syslog->dump(st_env.b_syslogPrepend, ( msg ) );
#define nSLOUT( msg ) if(st_env.pcsm_syslog)         st_env.pcsm_syslog->dump(false, ( msg ) );
#define RLOUT( msg )  if(st_env.pcsm_resultlog)	     st_env.pcsm_resultlog->dump(st_env.b_syslogPrepend, ( msg ) );
#define nRLOUT( msg ) if(st_env.pcsm_resultlog)      st_env.pcsm_resultlog->dump(false, ( msg ) );

#define pULOUT( msg )  if(mps_env->pcsm_userlog)   mps_env->pcsm_userlog->dump(mps_env->b_syslogPrepend, ( msg ) );
#define pnULOUT( msg ) if(mps_env->pcsm_userlog)   mps_env->pcsm_userlog->dump(false, ( msg ) );
#define pSLOUT( msg )  if(mps_env->pcsm_syslog)    mps_env->pcsm_syslog->dump(mps_env->b_syslogPrepend, ( msg ) );
#define pnSLOUT( msg ) if(mps_env->pcsm_syslog)    mps_env->pcsm_syslog->dump(false, ( msg ) );
#define pRLOUT( msg )  if(mps_env->pcsm_resultlog) mps_env->pcsm_resultlog->dump(mps_env->b_syslogPrepend, ( msg ) );
#define pnRLOUT( msg ) if(mps_env->pcsm_resultlog) mps_env->pcsm_resultlog->dump(false, ( msg ) );

#define Cn( msg ) cout << ( msg ) << endl;
#define CW(x, msg) cout.width((x))  ; cout << (msg);
#define CWn(x, msg) cout.width((x))  ; cout << (msg); cout << endl;
#define CCn(w, msg) for(unsigned __co=0; __co<(w-strlen(msg))/2; __co++)\
       cout << " ";    \
   cout << (msg); cout << endl;

// A simple structure to hold some statistical information
typedef struct _e_stats {
    float       f_max;
    int         indexMax;
    float       f_min;
    int         indexMin;
} e_stats;


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
char* 	lsprintf(int lw, char* pch_bufferOut, const char* format, ...);
void    colprintf(int lw, int rw,
                    const char* pch_lstr, const char* format, ...);
char*   colsprintf(int lw, int rw, char* pch_buffer,
                    const char* pch_lstr, const char* format, ...);

int
arr_stats(e_stats& a_stats, float* arr, int asize);

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

/// \fn void str_findAndReplace(string& astr_source, const string astr_find, string astr_replace)
/// \brief Find and replace.
/// \param astr_source 	The source string
/// \param astr_find  	The string to find in the source
/// \param astr_replace	The replacement string
/// \return bool	TRUE: Performed replacement; FALSE: No replacement
bool
str_findAndReplace(
    string& 		astr_source,
    const string 	astr_find,
    string		astr_replace
);


void warn(
  string          str_action,
  string          str_errorMsg,
  int             errorCode
);

void
error_exit(
  string          str_action,
  string          str_errorMsg,
  int             errorCode
);

void
command_line_error(
  char   *fmt, ...
);

#endif //__GENERAL_H__

