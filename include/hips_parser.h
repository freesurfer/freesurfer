/**
 * @file  hips_parser.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.2 $
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


/*
 * Copyright (c) 1991 Michael Landy
 *
 * Disclaimer:  No guarantees of performance accompany this software,
 * nor is any responsibility assumed on the part of the authors.  All the
 * software has been tested extensively and every effort has been made to
 * insure its reliability.
 */

/*
 * hips_parser.h - definitions related to the HIPS argument parser
 *
 * David Wilson - 30/1/91
 */

/*
 * Definitions for the command line parser lib/sources/parseargs.c
 */

/*
 * Make the library visible.
 */

extern int  parseargs ();

/*
 * Define an enumeration for the accepted trailing filename formats. Note
 * that these should keep in sequence since they are used as index values
 * in data structures internal to the parser.
 */

#define FFNONE 0
#define FFONE 1
#define FFTWO 2
#define FFLIST  3

typedef int Filename_Format;

/*
 * Define a flag as a character string of arbitrary length.
 */

typedef char *Flag;

#define LASTFLAG (Flag) 0 /* Terminator for flag list */

/*
 * Define an enumeration for the accepted parameter types. Note that
 * these values should keep in sequence since they are used as index
 * values in data structures internal to the parser.
 */

#define PTNULL    -1  /* Null type used as error/terminator */
#define PTBOOLEAN   0
#define PTCHAR     1
#define PTSTRING    2
#define PTINT     3
#define PTDOUBLE    4
#define PTFILENAME  5
#define PTLIST      6

typedef int  Parameter_Type;

#define LASTPARAMETER PTNULL /* Terminator for parameter type list */

/*
 * Define a structure for a parameter to a flag option. This holds the
 * type of the parameter and its default value. Note that the default
 * value is given as a string.
 */

struct parameter
{
  Parameter_Type type;
  char *pdefault;
  char *par_usage;
};

typedef struct parameter Parameter;

/*
 * Define a structure for the format of the flag options accepted by a HIPS
 * filter. This consists of the flag option itself, a list of all the
 * mutually exclusive flag options, a count of the minimum  number of
 * parameters which must be present and a list of the types and default
 * values of all parameters associated with the option.
 *
 * Note that this is the structure which is filled in by the user for
 * each flag recognised by a filter and passed to the parser.
 */

#define MAX_MUTEX_FLAGS  20 /* Arbitrary limits - may be extended */
#define MAX_PARAMETERS   20

struct flag_format
{
  Flag  value;
  Flag  mutex_flags [MAX_MUTEX_FLAGS];
  int  min_parameters;
  Parameter parameters [MAX_PARAMETERS];
};

typedef struct flag_format Flag_Format;

/*
 * Define a structure to hold a list of trailing filenames.
 */

struct file_list
{
  int   *count;
  char  ***list;
};

typedef struct file_list  File_List;

/*
 * Define a pointer to the filenames associated with each of the
 * accepted filename formats.
 */

union filename_ptr
{
  char       **filename;
  char       **filepair [2];
  File_List  filenames;
};

typedef union filename_ptr Filename_Ptr;

/*
 * Define a structure to hold a PTLIST parameter
 */

struct listarg
{
  int     argcount;
  char    **args;
};

typedef struct listarg Listarg;

/*
 * Define a generic pointer to a parameter of arbitrary type.
 */

union generic_ptr
{
  int         *boolean_ptr;
  char        *char_ptr;
  char        **string_ptr;
  int         *int_ptr;
  double      *double_ptr;
  char        **filename_ptr;
  Listarg     *listarg_ptr;
};

typedef union generic_ptr Generic_Ptr;

/*
 * Define a structure to hold a list of pointers to each of the
 * parameters assoociated with a flag. The variables pointed to by these
 * pointers are set to contain the required operating mode for the
 * filter. The structure also holds a signal to indicate whether a flag
 * is mutually excluded.
 */

struct flag_key
{
  Flag_Format    *format;
  Generic_Ptr    parameter_ptrs[MAX_PARAMETERS];
  int     locked;  /* Really h_boolean*/
  h_boolean     specified;
};

typedef struct flag_key Flag_Key;
