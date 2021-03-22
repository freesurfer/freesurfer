#ifndef  MINC_FILES_H
#define  MINC_FILES_H

/*
 * Original Author: David MacDonald, modified to compile within freesurfer/utils by Bevin Brett
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

/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
@VERSION    : $Header: /private-cvsroot/minc/volume_io/Include/volume_io/files.h,v 1.8.2.2 2005/03/31 17:39:49 bert Exp $
---------------------------------------------------------------------------- */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : files.h
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Types for use with the general file io routines of the library.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#include  <stdio.h>

#include "minc_volume_io.h"

typedef int VIO_BOOL;	// Note - this is written to and read from files, so can't be changed to bool or other types
#define VIO_TRUE  1
#define VIO_FALSE 0

typedef  enum  { ASCII_FORMAT, BINARY_FORMAT }          VIO_File_formats;
typedef  enum  { READ_FILE, WRITE_FILE, APPEND_FILE }   VIO_IO_types;

VIO_Status  input_character(
    FILE *file,
    char *ch );

VIO_Status  input_nonwhite_character(
    FILE   *file,
    char   *ch );

VIO_Status  input_string(
    FILE    *file,
    char*   *str,
    char    termination_char );
   
VIO_Status  input_double(
    FILE   *file,
    double *d );

VIO_Status  input_int(
    FILE  *file,
    int   *i );
    
VIO_Status  input_newline(
    FILE *file );

VIO_Status  io_int(
    FILE            *file,
    VIO_IO_types    io_flag,
    VIO_File_formats format,
    int             *i );

VIO_Status  unget_character(
    FILE  *file,
    char  ch );
    
VIO_Status  io_binary_data(
    FILE            *file,
    VIO_IO_types     io_flag,
    void            *data,
    size_t           element_size,
    int              n );

char*  expand_filename(
    const char*  filename );

char*  extract_directory(
    const char*    filename );

char*  get_absolute_filename(
    const char*    filename,		// expands and prepends path if necessary
    const char*    directory );
    
VIO_Status  open_file(
    const char*        filename,
    VIO_IO_types       io_type,
    VIO_File_formats   file_format,
    FILE             **file );
    
VIO_Status open_file_with_default_suffix(
    const char*        filename,
    const char*        default_suffix,
    VIO_IO_types       io_type,
    VIO_File_formats   file_format,
    FILE             **file );

VIO_Status  set_file_position(
    FILE     *file,
    long     byte_position );
    
VIO_Status  close_file(
    FILE     *file );

VIO_Status  copy_file(
    const char*  src,
    const char*  dest );

void  remove_file(
    const char*  filename );

#endif
