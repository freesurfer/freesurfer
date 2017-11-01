#if defined(BEVIN_EXCLUDE_MINC)

/*
 * Original Author: David MacDonald, modified to compile within freesurfer/utils by Bevin Brett
 * CVS Revision Info:
 *    $Author: Brett $
 *    $Date: 2017/11/01 12:46:00 $
 *    $Revision: 1.00 $
 *
 * Copyright © 2017 The General Hospital Corporation (Boston, MA) "MGH"
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
---------------------------------------------------------------------------- */

#include  "files.h"

#if HAVE_PWD_H
#include  <pwd.h>
#endif /* HAVE_PWD_H */

#include  <stdlib.h>

#if HAVE_UNISTD_H
#include  <unistd.h>
#else
#error "not HAVE_UNISTD_H"
#endif /* HAVE_UNISTD_H */
#include  <errno.h>

#include  <string.h>

#if TIME_WITH_SYS_TIME
#include  <sys/time.h>
#include  <time.h>
#else
#if HAVE_SYS_TIME_H
#include  <sys/time.h>
#else
#include  <time.h>
#endif
#endif

#if HAVE_UNISTD_H
#include  <unistd.h>
#endif

static void delete_string(const char* s) {
    free((void*)s);
}

static int string_length(const char* s) {
    return s ? (int)strlen(s) : 0;
}

static void concat_to_string(
    const char* *lhs,
    const char*  rhs) {
    if (!rhs) return;
    size_t rhs_strlen = strlen(rhs);
    if (!*lhs) *lhs = (const char*)malloc(rhs_strlen + 1);
    size_t lhs_strlen = strlen(*lhs);
    *lhs = (const char*)realloc((void*)*lhs, lhs_strlen + rhs_strlen + 1);
    memcpy((void*)(*lhs + lhs_strlen), (void*)rhs, rhs_strlen + 1);
}

const char* create_string(const char* s) {
    const char* p = NULL;
    concat_to_string(&p, s);
    return p;
}

const char* get_date() {

    time_t currentTime;
    time( &currentTime );

    struct tm  localtimeBuffer;
    struct tm* localTime = localtime_r(&currentTime, &localtimeBuffer);

    char   asctimeBuffer[32];
    char * ascTime = asctime_r(localTime, asctimeBuffer);

    return create_string(ascTime);
}

static  bool          has_no_extension( const char* );
static  const char*   compressed_endings[] = { ".z", ".Z", ".gz" };


/* ----------------------------- MNI Header -----------------------------------
@NAME       : print_system_error
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Prints the most recent system error.
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    :        , 1996    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  print_system_error( void )
{
    const char   *error = strerror( errno );
    fprintf(stderr, "\nSystem message: %s\n", error );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : file_exists
@INPUT      : filename
@OUTPUT     : 
@RETURNS    : VIO_TRUE or VIO_FALSE if file exists
@DESCRIPTION: Checks if the file of the given name exists
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_BOOL  file_exists(
    const char*        filename )
{
    VIO_BOOL  exists;
    FILE     *file;
    const char*   expanded;

    expanded = expand_filename( filename );

    file = fopen( expanded, "r" );

    if( file != NULL )
    {
        (void) fclose( file );
        exists = VIO_TRUE;
    }
    else
        exists = VIO_FALSE;

    delete_string( expanded );

    return( exists );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : file_directory_exists
@INPUT      : filename
@OUTPUT     : 
@RETURNS    : VIO_TRUE if directory containing file exists.
@DESCRIPTION: Checks if the directory contained in the path name exists.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov. 2, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_BOOL  file_directory_exists(
    const char*        filename )
{
    VIO_BOOL  exists;
    const char*   dir;

    dir = extract_directory( filename );

    if( string_length( dir ) != 0 )
        exists = file_exists( dir );
    else
        exists = VIO_TRUE;

    delete_string( dir );

    return( exists );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_clobber_file
@INPUT      : filename
@OUTPUT     : 
@RETURNS    : VIO_TRUE if can write file
@DESCRIPTION: Checks if the file exists.  If so, asks the user for permission
              to overwrite the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_BOOL  check_clobber_file(
    const char*   filename )
{
    char     ch;
    VIO_BOOL  okay;
    const char*   expanded;

    okay = VIO_TRUE;

    if( file_exists( filename ) )
    {
        expanded = expand_filename( filename );

        printf( "File <%s> exists, do you wish to overwrite (y or n): ",
               expanded );

        delete_string( expanded );

        while( input_character( stdin, &ch ) == OK && ch != 'y' && ch != 'n' &&
               ch != 'N' && ch != 'Y' )
        {
            if( ch == '\n' )
                printf( "  Please type y or n: " );
        }

        (void) input_newline( stdin );

        okay = (ch == 'y' || ch == 'Y');
    }

    return( okay );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_clobber_file_default_suffix
@INPUT      : filename
              default_suffix
@OUTPUT     : 
@RETURNS    : VIO_TRUE if can write file
@DESCRIPTION: Checks if the file exists (adding the default suffix if
              necessary).  If the file exists, asks the user for permission
              to overwrite the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_BOOL  check_clobber_file_default_suffix(
    const char*   filename,
    const char*   default_suffix )
{
    const char*   expanded;
    VIO_BOOL  can_write;

    expanded = expand_filename( filename );

    if( has_no_extension( expanded ) )
    {
        concat_to_string( &expanded, "." );
        concat_to_string( &expanded, default_suffix );
    }

    can_write = check_clobber_file( expanded );

    delete_string( expanded );

    return( can_write );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_backup_filename
@INPUT      : filename
@OUTPUT     : 
@RETURNS    : const char* - a backup filename
@DESCRIPTION: Creates a backup filename that is filename.{date}.bkp
              If this already exists (not very likely), then it tries
              appending _1, _2, ...
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    : Feb.  3, 1997    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  const char*  create_backup_filename(
    const char*   filename )
{
    const char * const expanded = expand_filename( filename );
    const char * const date     = get_date();

    char* const backup_filename = (char*)malloc( string_length( expanded ) + string_length( date ) + 100 );

    int count = 0;
    do
    {
        if( count == 0 )
        {
            (void) sprintf( backup_filename, "%s.%s.bkp",
                            expanded, date );
        }
        else
        {
            (void) sprintf( backup_filename, "%s.%s.bkp_%d",
                            expanded, date, count );
        }

        int len = string_length( backup_filename );
        while( len > 0 && 
               (backup_filename[len-1] == ' ' ||
                backup_filename[len-1] == '\t' ||
                backup_filename[len-1] == '\n') )
        {
            --len;
        }
        backup_filename[len] = (char) 0;

        int i;
        for (i = 0; i < len; i++)
        {
            if( backup_filename[i] == ' ' || 
	        backup_filename[i] == '\t' ||
                backup_filename[i] == '\n' )
                backup_filename[i] = '_';

            /* remove ':' for windows */
            if( backup_filename[i] == ':')
               backup_filename[i] = '-';
        }

        ++count;
    }
    while( file_exists( backup_filename ) );

    delete_string( expanded );
    delete_string( date );

    return( backup_filename );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_backup_file
@INPUT      : filename
@OUTPUT     : backup_filename
@RETURNS    : OK or ERROR
@DESCRIPTION: If the file exists, creates a backup of the file, and passes
              back the name of the backup file, which must be passed to
              cleanup_backup_file after the write of filename is performed.
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    : Feb.  3, 1997    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  make_backup_file(
    const char*   filename,
    const char*   *backup_filename )
{
    VIO_Status   status;

    status = OK;

    if( file_exists( filename ) )
    {
        *backup_filename = create_backup_filename( filename );

        status = copy_file( filename, *backup_filename );

        if( status != OK )
        {
            fprintf(stderr, "Error making backup file for: %s\n", filename );
            *backup_filename = NULL;
        }
    }
    else
        *backup_filename = NULL;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : cleanup_backup_file
@INPUT      : filename
              backup_filename
              status_of_write
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: This function is called after writing a file.  If a backup file
              was made before the write, then it is deleted, if the write was
              successful, or copied to the original file, otherwise.
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    : Feb.  3, 1997    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  cleanup_backup_file(
    const char*   filename,
    const char*   backup_filename,
    VIO_Status   status_of_write )
{
    VIO_BOOL  can_remove;

    if( backup_filename != NULL )
    {
        can_remove = VIO_TRUE;
        if( status_of_write != OK )
        {
            if( copy_file( backup_filename, filename ) != OK )
            {
                fprintf(stderr, "File %s was corrupted during a failed write,\n",
                             filename );
                fprintf(stderr,
                   "File %s contains the state prior to the write attempt.\n",
                  backup_filename );
                can_remove = VIO_FALSE;
            }
        }

        if( can_remove )
            remove_file( backup_filename );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : remove_file
@INPUT      : filename
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the given file.
@METHOD     : Makes a system call to unlink().
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  remove_file(
    const char*  filename )
{
    const char*   expanded;

    expanded = expand_filename( filename );

    if( unlink( expanded ) != 0 )
    {
        fprintf(stderr, "Error removing %s.  ", expanded );
        print_system_error();
    }

    delete_string( expanded );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : copy_file
@INPUT      : src
              dest
@OUTPUT     : 
@RETURNS    : OK or ERROR
@DESCRIPTION: Copies the src file to the dest file.
@METHOD     : Makes a UNIX system call, using /bin/cp
@GLOBALS    : 
@CALLS      :  
@CREATED    : Feb.  3, 1997    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  copy_file(
    const char*  src,
    const char*  dest )
{
    VIO_Status   status;
    const char*   src_expanded, dest_expanded, command;

    src_expanded = expand_filename( src );
    dest_expanded = expand_filename( dest );

    command = concat_strings( "/bin/cp ", src_expanded );
    concat_to_string( &command, " " );
    concat_to_string( &command, dest_expanded );

    if( system( command ) != 0 )
    {
        fprintf(stderr, "Error copying file %s to %s: ",
                     src_expanded, dest_expanded );
        print_system_error();
        status = ERROR;
    }
    else
        status = OK;

    delete_string( src_expanded );
    delete_string( dest_expanded );
    delete_string( command );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : move_file
@INPUT      : src
              dest
@OUTPUT     : 
@RETURNS    : OK or ERROR
@DESCRIPTION: Move the src file to the dest file.
@METHOD     : Makes a UNIX system call, using /bin/mv
@GLOBALS    : 
@CALLS      :  
@CREATED    : Feb.  3, 1997    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  move_file(
    const char*  src,
    const char*  dest )
{
    VIO_Status   status;
    const char*   src_expanded, dest_expanded, command;

    src_expanded = expand_filename( src );
    dest_expanded = expand_filename( dest );

    command = concat_strings( "/bin/cp -f ", src_expanded );
    concat_to_string( &command, " " );
    concat_to_string( &command, dest_expanded );

    if( system( command ) != 0 )
    {
        fprintf(stderr, "Error moving file %s to %s: ",
                     src_expanded, dest_expanded );
        print_system_error();
        status = ERROR;
    }
    else
        status = OK;

    delete_string( src_expanded );
    delete_string( dest_expanded );
    delete_string( command );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_user_home_directory
@INPUT      : user_name
@OUTPUT     : 
@RETURNS    : Pointer to home directory string.
@DESCRIPTION: Returns the home directory of the specified user.
@METHOD     : UNIX password file utilities
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  const char*  get_user_home_directory(
    const char*   user_name )
{
#if HAVE_GETPWNAM
    struct   passwd  *p;

    p = getpwnam( user_name );

    if( p == NULL )
        return( NULL );
    else
        return( p->pw_dir );
#else
    return (".");
#endif /* HAVE_GETPWNAM */
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : expand_filename
@INPUT      : filename
@OUTPUT     : expanded_filename
@RETURNS    : 
@DESCRIPTION: Expands certain strings in the filename, if present:

                  environment variables, e.g.   "$DATA_DIR/filename.txt"
                  ~                      e.g.   "~david/filename.txt"

              If a dollar sign or backslash is desired, it must be preceded
              by a backslash.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

const char*  expand_filename(
    const char*  filename )
{
    int      i, new_i, dest, len, env_index;
    VIO_BOOL  tilde_found, prev_was_backslash;
    char     *expand_value;
    int      n_alloced, n_env_alloced;
    const char*   env, expanded;

    /* --- copy from filename to expanded_filename, changing environment
           variables and home directories */

    len = string_length( filename );

    prev_was_backslash = VIO_FALSE;
    i = 0;
    dest = 0;

    n_alloced = 0;

    n_env_alloced = 0;
    env = NULL;
    expanded = NULL;

    while( i < len+1 )
    {
        /* --- if not escaped by backslash, and is either a '~' at the
               beginning or a '$' anywhere, expand it */

        if( !prev_was_backslash &&
            ((i == 0 && filename[i] == '~') || filename[i] == '$') )
        {
            /* --- pick up the environment variable name or user name, by
                   searching until the next '/' or a '.' or end of string */

            new_i = i;
            tilde_found = (filename[new_i] == '~');
            ++new_i;

            env_index = 0;
            while( filename[new_i] != '/' &&
                   filename[new_i] != '.' &&
                   filename[new_i] != END_OF_STRING )
            {
                ADD_ELEMENT_TO_ARRAY_WITH_SIZE( env, n_env_alloced, env_index,
                                                filename[new_i],
                                                DEFAULT_CHUNK_SIZE );
                ++new_i;
            }

            ADD_ELEMENT_TO_ARRAY_WITH_SIZE( env, n_env_alloced, env_index,
                                            END_OF_STRING, DEFAULT_CHUNK_SIZE );

            /* --- if expanding a '~', find the corresponding home directory */

            if( tilde_found )
            {
                if( string_length( env ) == 0 )
                    expand_value = getenv( "HOME" );
                else
                    expand_value = get_user_home_directory( env );
            }
            else               /* --- get the environment variable value */
                expand_value = getenv( env );

            /* --- if an expansion is found, copy it, otherwise just copy char*/

            if( expand_value != NULL )
            {
                SET_ARRAY_SIZE( expanded, n_alloced,
                                n_alloced + string_length(expand_value),
                                DEFAULT_CHUNK_SIZE );
                n_alloced += string_length(expand_value);
                (void) strcpy( &expanded[dest], expand_value );
                dest += string_length( expand_value );
                i = new_i;
            }
            else
            {
                SET_ARRAY_SIZE( expanded, n_alloced, n_alloced + 1,
                                DEFAULT_CHUNK_SIZE );
                ++n_alloced;
                expanded[dest] = filename[i];
                ++dest;
                ++i;
            }

            prev_was_backslash = VIO_FALSE;
        }
        else
        {
            /* --- if not a backslash or if it is escaped, add character */

            if( filename[i] != '\\' || prev_was_backslash )
            {
                SET_ARRAY_SIZE( expanded, n_alloced, n_alloced + 1,
                                DEFAULT_CHUNK_SIZE );
                ++n_alloced;
                expanded[dest] = filename[i];
                ++dest;
                prev_was_backslash = VIO_FALSE;
            }
            else
                prev_was_backslash = VIO_TRUE;
            ++i;
        }
    }

    if( n_env_alloced > 0 )
        delete_string( env );

    return( expanded );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : filename_extension_matches
@INPUT      : filename
              extension
@OUTPUT     : 
@RETURNS    : VIO_TRUE if filename extension matches
@DESCRIPTION: Checks if the filename ends in a period, then the given
              extension.  Note that the filename first undergoes expansion
              for home directories and environment variables, and any
              ending of ".z", ".Z", or ".gz" is first removed.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_BOOL  filename_extension_matches(
    const char*   filename,
    const char*   extension )
{
    int       len, i;
    const char*    filename_no_z, ending;
    VIO_BOOL   matches;

    filename_no_z = expand_filename( filename );

    len = string_length( filename_no_z );

    for_less( i, 0, SIZEOF_STATIC_ARRAY(compressed_endings) )
    {
        if( string_ends_in( filename_no_z, compressed_endings[i] ) )
        {
            filename_no_z[len-string_length(compressed_endings[i])] =
                                   END_OF_STRING;
        }
    }

    ending = concat_strings( ".", extension );

    matches = string_ends_in( filename_no_z, ending );

    delete_string( filename_no_z );
    delete_string( ending );

    return( matches );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : remove_directories_from_filename
@INPUT      : filename
@OUTPUT     : filename_no_directories
@RETURNS    : 
@DESCRIPTION: Creates a new filename with no directories in it.
              E.G.  if filename equals  "/usr/people/david/test.c"
              filename_no_directories will be set to "test.c"
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

const char*  remove_directories_from_filename(
    const char*  filename )
{
    const char*   expanded, no_directories;
    int      i;

    expanded = expand_filename( filename );

    i = string_length( expanded );

    while( i >= 0 && expanded[i] != '/' )
        --i;

    ++i;

    no_directories = create_string( &expanded[i] );

    delete_string( expanded );

    return( no_directories );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : file_exists_as_compressed
@INPUT      : filename
@OUTPUT     : compressed_filename
@RETURNS    : VIO_TRUE if a compressed file exists
@DESCRIPTION: Checks to see if a compressed version of the file exists.  If so,
              passes back the name of the compressed file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 21, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_BOOL  file_exists_as_compressed(
    const char*       filename,
    const char*       *compressed_filename )
{
    int      i;
    const char*   compressed, expanded;
    VIO_BOOL  gzipped;

    gzipped = VIO_FALSE;

    expanded = expand_filename( filename );

    /* --- check to see if file.z or file.Z, etc, exists */

    for_less( i, 0, SIZEOF_STATIC_ARRAY( compressed_endings ) )
    {
        compressed = concat_strings( expanded, compressed_endings[i] );

        if( file_exists( compressed ) )
        {
            if( *compressed_filename == filename )
                delete_string( filename );

            *compressed_filename = compressed;
            gzipped = VIO_TRUE;
            break;
        }

        delete_string( compressed );
    }

    delete_string( expanded );

    return( gzipped );
}

const char*  get_temporary_filename( void )
{
    return micreate_tempfile();
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : open_file
@INPUT      : filename
            : io_type        READ_FILE or WRITE_FILE
            : file_format    ASCII_FORMAT or BINARY_FORMAT
@OUTPUT     : file
@RETURNS    : 
@DESCRIPTION: Opens the given filename for ascii or binary input or output.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  open_file(
    const char*             filename,
    IO_types           io_type,
    File_formats       file_format,
    FILE               **file )
{
    VIO_Status   status;
    int      i;
    char     *tmp_name;
    char     command[EXTREMELY_LARGE_STRING_SIZE];
    const char*   access_str, expanded;
    VIO_BOOL  gzipped;
    int      command_status;

    /* --- determine what mode of file access */

    switch( io_type )
    {
    case APPEND_FILE:   access_str = create_string( "a" );  break;

    case WRITE_FILE:    access_str = create_string( "w" );  break;

    case READ_FILE:
    default:            access_str = create_string( "r" );  break;
    }

    /* --- check if ascii or binary */

    if( file_format == BINARY_FORMAT )
        concat_to_string( &access_str, "b" );

    /* --- expand ~ and $ in filename */

    expanded = expand_filename( filename );

    gzipped = VIO_FALSE;

    /* --- if reading the file, check if it is in compressed format */

    if( io_type == READ_FILE )
    {
        /* --- check if the filename ends in one of the compressed suffixes */

        for_less( i, 0, SIZEOF_STATIC_ARRAY( compressed_endings ) )
        {
            if( string_ends_in( expanded, compressed_endings[i] ) )
            {
                gzipped = VIO_TRUE;
                break;
            }
        }

        /* --- if the filename does not have a compressed suffix and
               the file does not exist, check to see if file.z or file.Z, etc,
               exists */

        if( !gzipped && !file_exists( expanded ) )
            gzipped = file_exists_as_compressed( expanded, &expanded );
    }

    /* --- if reading from a compressed file, decompress it to a temp file */

    status = OK;

    if( gzipped )
    {
        /* --- uncompress to a temporary file */

        tmp_name = get_temporary_filename();

        (void) sprintf( command, "gunzip -c %s > %s", expanded, tmp_name );
        command_status = system( command );

        /* Try again, using bzip2 */
        if( command_status != 0 )
        {
            (void) sprintf( command, "bunzip2 -c %s > %s", expanded, tmp_name );
            command_status = system( command );
        }

        /* Check for failure */
        if( command_status != 0 )
        {
            fprintf(stderr, "Error uncompressing %s into %s using gunzip and bunzip2\n",
                        expanded, tmp_name );
            status = ERROR;
        }
        else
            replace_string( &expanded, create_string(tmp_name) );

        free(tmp_name);
    }

    /* --- finally, open the file */

    if( status == OK )
    {
        *file = fopen( expanded, access_str );

        if( *file == NULL )          /* --- print error message if needed */
        {
            fprintf(stderr, "Error:  could not open file \"%s\".  ", expanded );
            print_system_error();
            status = ERROR;
        }
        else if( gzipped )           /* if reading a decompressed temp file, */
            remove_file( expanded ); /* unlink it, so that when the program  */
                                     /* closes the file or dies, the file is */
                                     /* removed                              */
    }

    delete_string( access_str );
    delete_string( expanded );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : open_file_with_default_suffix
@INPUT      : filename
            : default_suffix  - e.g. ".obj"
            : io_type        READ_FILE or WRITE_FILE
            : file_format    ASCII_FORMAT or BINARY_FORMAT
@OUTPUT     : file
@RETURNS    : 
@DESCRIPTION: Opens the given filename for ascii or binary input or output.
            : On output, if the file has no suffix, it adds the default suffix.
            : On input, if the file does not exist as given, then it tries to
            : find the file with the default_suffix.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  open_file_with_default_suffix(
    const char*             filename,
    const char*             default_suffix,
    IO_types           io_type,
    File_formats       file_format,
    FILE               **file )
{
    VIO_Status   status;
    VIO_BOOL  suffix_added;
    const char*   used_filename, expanded;

    expanded = expand_filename( filename );

    if( io_type == READ_FILE )
    {
        suffix_added = VIO_FALSE;

        if( !file_exists(expanded) && has_no_extension( expanded ) )
        {
            used_filename = concat_strings( expanded, "." );
            concat_to_string( &used_filename, default_suffix );
            if( file_exists( used_filename ) )
                suffix_added = VIO_TRUE;
            else
                delete_string( used_filename );
        }

        if( !suffix_added )
            used_filename = create_string( expanded );
    }
    else if( has_no_extension( expanded ) )
    {
        used_filename = concat_strings( expanded, "." );
        concat_to_string( &used_filename, default_suffix );
    }
    else
    {
        used_filename = create_string( expanded );
    }

    status = open_file( used_filename, io_type, file_format, file );

    delete_string( expanded );
    delete_string( used_filename );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : has_no_extension
@INPUT      : filename
@OUTPUT     : 
@RETURNS    : VIO_TRUE if there is no . extension
@DESCRIPTION: Checks if there is an extension on the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  VIO_BOOL  has_no_extension(
    const char*   filename )
{
    const char*   base_name;
    VIO_BOOL  dot_found;

    base_name = remove_directories_from_filename( filename );

    dot_found = (find_character( base_name, '.' ) >= 0);

    delete_string( base_name  );

    return( !dot_found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_file_position
@INPUT      : file
            : byte_position
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the file position to the given offset from the start.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  set_file_position(
    FILE     *file,
    long     byte_position )
{
    VIO_Status   status;

    if( fseek( file, byte_position, 0 ) == 0 )
    {
        status = OK;
    }
    else
    {
        fprintf(stderr, "Error setting the file position.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : close_file
@INPUT      : file
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Closes the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  close_file(
    FILE     *file )
{
    if( file != NULL )
    {
        (void) fclose( file );
        return( OK );
    }
    else
        return( ERROR );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : extract_directory
@INPUT      : filename
@OUTPUT     : directory
@RETURNS    : 
@DESCRIPTION: Extracts the directory from the filename by copying the string
            : from the beginning up to the last '/'.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

const char*  extract_directory(
    const char*    filename )
{
    int     i, slash_index;
    const char*  expanded, directory;

    expanded = expand_filename( filename );

    slash_index = string_length(expanded) - 1;

    while( slash_index >= 0 && expanded[slash_index] != '/' )
        --slash_index;

    if( slash_index < 0 )
        directory = create_string( "." );
    else
    {
        ++slash_index;

        directory = alloc_string( slash_index );

        for_less( i, 0, slash_index )
            directory[i] = expanded[i];

        directory[slash_index] = END_OF_STRING;
    }

    delete_string( expanded );

    return( directory );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_absolute_filename
@INPUT      : filename
            : directory
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Given a filename and a default directory, determines the correct
            : filename by checking if the filename is a relative or absolute
            : pathname, and prepending the directory, if the former.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

const char*  get_absolute_filename(
    const char*    filename,
    const char*    directory )
{
    const char*  abs_filename, expanded;

    /* if the directory is non-null and the filename is not already
       absolute (begins with '/'), then prefix the directory to the filename */

    expanded = expand_filename( filename );

    if( string_length( directory ) > 0 && expanded[0] != '/' )
    {
        if( directory[string_length(directory)-1] == '/' )
            abs_filename = create_string( directory );
        else
            abs_filename = concat_strings( directory, "/" );
    }
    else
    {
        abs_filename = create_string( NULL );
    }

    concat_to_string( &abs_filename, expanded );

    delete_string( expanded );

    return( abs_filename );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : flush_file
@INPUT      : file
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Flushes the output buffer for the given file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  flush_file(
    FILE     *file )
{
    VIO_Status   status;

    if( fflush( file ) == 0 )
    {
        status = OK;
    }
    else
    {
        fprintf(stderr, "Error flushing file.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_character
@INPUT      : file
@OUTPUT     : ch
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs one character from the file, returning ERROR if eof.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_character(
    FILE  *file,
    char   *ch )
{
    VIO_Status   status;
    int      c;

    c = fgetc( file );

    if( c == EOF )
    {
        status = ERROR;
    }
    else
    {
        *ch = (char) c;
        status = OK;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : unget_character
@INPUT      : file
@OUTPUT     : ch
@RETURNS    : VIO_Status
@DESCRIPTION: Ungets one character back to the file, returning status.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  unget_character(
    FILE  *file,
    char  ch )
{
    VIO_Status   status;
    int      c;

    c = ungetc( (int) ch, file );

    if( c == EOF )
        status = ERROR;
    else
        status = OK;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_nonwhite_character
@INPUT      : file
@OUTPUT     : ch
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs the next nonwhite (tab, space, newline) character.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_nonwhite_character(
    FILE   *file,
    char   *ch )
{
    VIO_Status   status;

    do
    {
        status = input_character( file, ch );
    }
    while( status == OK && (*ch == ' ' || *ch == '\t' || *ch == '\n') );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_character
@INPUT      : file
            : ch
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs the character to the file, returning the status.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_character(
    FILE   *file,
    char   ch )
{
    VIO_Status   status;

    if( fputc( (int) ch, file ) != ch )
    {
        status = ERROR;
    }
    else
    {
        status = OK;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : skip_input_until
@INPUT      : file
            : search_char
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Skips characters in the file, up to and including the first match
            : of the search_char;
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status   skip_input_until(
    FILE   *file,
    char   search_char )
{
    VIO_Status   status;
    char     ch;

    status = OK;

    do
    {
        status = input_character( file, &ch );
    }
    while( status == OK && ch != search_char );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_string
@INPUT      : file
            : str
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs the string to the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_string(
    FILE    *file,
    const char*  str )
{
    VIO_Status   status;

    if( fprintf( file, "%s", str ) == string_length(str) )
        status = OK;
    else
    {
        fprintf(stderr, "Error outputting string.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_string
@INPUT      : file
            : termination_char
@OUTPUT     : str 
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs a string from the file.  First it skips white space, then
            : inputs all characters until the termination_char is found.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_string(
    FILE    *file,
    const char*  *str,
    char    termination_char )
{
    char    ch;
    VIO_Status  status;

    status = input_nonwhite_character( file, &ch );

    *str = create_string( NULL );

    while( status == OK && ch != termination_char && ch != '\n' )
    {
        concat_char_to_string( str, ch );

        status = input_character( file, &ch );
    }

    if( termination_char != '\n' && ch == '\n' )
        (void) unget_character( file, ch );

    if( status != OK )
    {
        delete_string( *str );
        *str = NULL;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_quoted_string
@INPUT      : file
@OUTPUT     : str
@RETURNS    : VIO_Status
@DESCRIPTION: Skips to the next nonwhitespace character, checks if it is a
            : quotation mark ( ", ', or ` ), then reads characters into the
            : string until the : next quotation mark.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_quoted_string(
    FILE            *file,
    const char*          *str )
{
    char     ch, quote;
    VIO_Status   status;

    status = input_nonwhite_character( file, &quote );

    if( status == OK && quote != '"' && quote != '\'' && quote != '`' )
        status = ERROR;

    if( status == OK )
        status = input_character( file, &ch );

    *str = create_string( NULL );

    while( status == OK && ch != quote )
    {
        concat_char_to_string( str, ch );

        status = input_character( file, &ch );
    }

    if( status != OK )
    {
        delete_string( *str );
        *str = NULL;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_possibly_quoted_string
@INPUT      : file
            : str
            : str_length    - size of string storage
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Skips to the next nonwhitespace character, checks if it is a
            : quotation mark, then reads characters into the string until the
            : next quotation mark.  If it is not a quotation mark, reads to
            : the next white space.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_possibly_quoted_string(
    FILE            *file,
    const char*          *str )
{
    VIO_BOOL  quoted;
    char     ch, quote;
    VIO_Status   status;

    status = input_nonwhite_character( file, &quote );

    if( status == OK )
    {
        if( quote == '"' || quote == '\'' || quote == '`' )
        {
            quoted = VIO_TRUE;
            status = input_character( file, &ch );
        }
        else
        {
            quoted = VIO_FALSE;
            ch = quote;
        }
    }

    *str = create_string( NULL );

    while( status == OK &&
           (quoted && ch != quote ||
            !quoted && ch != ' ' && ch != '\t' && ch != '\n') )
    {
        concat_char_to_string( str, ch );

        status = input_character( file, &ch );
    }

    if( !quoted )
        (void) unget_character( file, ch );

    if( status != OK )
    {
        delete_string( *str );
        *str = NULL;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_quoted_string
@INPUT      : file
            : str
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs the given string, with quotation marks around it.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_quoted_string(
    FILE            *file,
    const char*          str )
{
    VIO_Status   status;

    if( fprintf( file, " \"%s\"", str ) > 0 )
        status = OK;
    else
        status = ERROR;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_binary_data
@INPUT      : file
            : element_size       size of each element
            : n                  number of elements
@OUTPUT     : data               array of elements to input
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs the data in binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_binary_data(
    FILE            *file,
    void            *data,
    size_t          element_size,
    int             n )
{
    VIO_Status   status;
    int      n_done;

    status = OK;

    n_done = (int) fread( data, element_size, (size_t) n, file );
    if( n_done != n )
    {
        fprintf(stderr, "Error inputting binary data.\n" );
        fprintf(stderr, "     (%d out of %d items of size %ld).  ", n_done, n,
                     element_size );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_binary_data
@INPUT      : file
            : data               array of elements to output
            : element_size       size of each element
            : n                  number of elements
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs the data in binary format.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_binary_data(
    FILE            *file,
    void            *data,
    size_t          element_size,
    int             n )
{
    VIO_Status   status;
    int      n_done;

    status = OK;

    n_done = (int) fwrite( data, element_size, (size_t) n, file );
    if( n_done != n )
    {
        fprintf(stderr, "Error outputting binary data.\n" );
        fprintf(stderr, "     (%d out of %d items of size %ld).  ", n_done, n,
                     element_size );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_newline
@INPUT      : file
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Skips to after the next newline in the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_newline(
    FILE            *file )
{
    VIO_Status   status;

    status = skip_input_until( file, '\n' );

    if( status != OK )
    {
        fprintf(stderr, "Error inputting newline.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_newline
@INPUT      : file
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs a newline to the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_newline(
    FILE            *file )
{
    VIO_Status   status;

    if( fprintf( file, "\n" ) > 0 )
        status = OK;
    else
    {
        fprintf(stderr, "Error outputting newline.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_line
@INPUT      : line         - string to input to
            : str_length   - storage allocated to the string
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs all characters upto the next newline.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_line(
    FILE    *file,
    const char*  *line )
{
    VIO_Status   status;
    char     ch;

    *line = create_string( NULL );

    status = input_character( file, &ch );

    while( status == OK && ch != '\n' )
    {
        concat_char_to_string( line, ch );

        status = input_character( file, &ch );
    }

    if( status != OK )
    {
        delete_string( *line );
        *line = NULL;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_boolean
@INPUT      : file
@OUTPUT     : b
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs a VIO_BOOL value from a file, by looking for an 'f' or 't'.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_boolean(
    FILE            *file,
    VIO_BOOL         *b )
{
    VIO_Status   status;
    char     ch;

    status = input_nonwhite_character( file, &ch );

    if( status == OK )
    {
        if( ch == 'f' || ch == 'F' )
            *b = VIO_FALSE;
        else if( ch == 't' || ch == 'T' )
            *b = VIO_TRUE;
        else
            status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_boolean
@INPUT      : file
            : b
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs a T or F to the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_boolean(
    FILE            *file,
    VIO_BOOL         b )
{
    VIO_Status   status;
    const char*   str;

    status = OK;

    if( b )
        str = "T";
    else
        str = "F";

    if( fprintf( file, " %s", str ) <= 0 )
    {
        fprintf(stderr, "Error outputting VIO_BOOL.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_short
@INPUT      : file
@OUTPUT     : s
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs an ascii short.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_short(
    FILE            *file,
    short           *s )
{
    VIO_Status   status;

    if( fscanf( file, "%hd", s ) == 1 )
        status = OK;
    else
        status = ERROR;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_short
@INPUT      : file
            : s
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs an ascii short.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_short(
    FILE            *file,
    short           s )
{
    VIO_Status   status;

    if( fprintf( file, " %d", s ) > 0 )
        status = OK;
    else
    {
        fprintf(stderr, "Error outputting short.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_unsigned_short
@INPUT      : file
@OUTPUT     : s
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs an ascii unsigned short.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_unsigned_short(
    FILE            *file,
    unsigned short  *s )
{
    int      i;
    VIO_Status   status;

    if( fscanf( file, "%d", &i ) == 1 )
    {
        *s = (unsigned short) i;
        status = OK;
    }
    else
        status = ERROR;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_unsigned_short
@INPUT      : file
            : s
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs an ascii unsigned short.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_unsigned_short(
    FILE            *file,
    unsigned short  s )
{
    VIO_Status   status;

    if( fprintf( file, " %d", (int) s ) > 0 )
        status = OK;
    else
    {
        fprintf(stderr, "Error outputting unsigned short.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_int
@INPUT      : file
@OUTPUT     : i
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs an ascii integer.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_int(
    FILE  *file,
    int   *i )
{
    VIO_Status   status;

    if( fscanf( file, "%d", i ) == 1 )
        status = OK;
    else
        status = ERROR;

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_int
@INPUT      : file
            : i
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs an ascii integer.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_int(
    FILE            *file,
    int             i )
{
    VIO_Status   status;

    if( fprintf( file, " %d", i ) > 0 )
        status = OK;
    else
    {
        fprintf(stderr, "Error outputting int.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_float
@INPUT      : file
@OUTPUT     : f
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs an ascii float.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_float(
    FILE            *file,
    float           *f )
{
    VIO_Status   status;

    if( fscanf( file, "%f", f ) == 1 )
        status = OK;
    else
    {
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_float
@INPUT      : file
            : f
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs an ascii float value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_float(
    FILE            *file,
    float           f )
{
    VIO_Status   status;

    if( fprintf( file, " %g", f ) > 0 )
        status = OK;
    else
    {
        fprintf(stderr, "Error outputting float.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : input_double
@INPUT      : file
@OUTPUT     : d
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs an ascii double.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  input_double(
    FILE            *file,
    double          *d )
{
    VIO_Status   status;

    if( fscanf( file, "%lf", d ) == 1 )
        status = OK;
    else
    {
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : output_double
@INPUT      : file
            : d
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Outputs an ascii double value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  output_double(
    FILE            *file,
    double          d )
{
    VIO_Status   status;

    if( fprintf( file, " %g", d ) > 0 )
        status = OK;
    else
    {
        fprintf(stderr, "Error outputting double.  " );
        print_system_error();
        status = ERROR;
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_binary_data
@INPUT      : file
            : io_flag
            : data
            : element_size
            : n
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs binary data, depending on io_flag.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_binary_data(
    FILE            *file,
    IO_types        io_flag,
    void            *data,
    size_t          element_size,
    int             n )
{
    VIO_Status   status;

    if( io_flag == READ_FILE )
        status = input_binary_data( file, data, element_size, n );
    else
        status = output_binary_data( file, data, element_size, n );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_newline
@INPUT      : file
            : io_flag
            : data
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary newline char, as appropriate.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_newline(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format )
{
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
            status = OK;
        else
            status = output_newline( file );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_quoted_string
@INPUT      : file
            : io_flag
            : format
            : str
            : str_length
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary quoted string.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_quoted_string(
    FILE           *file,
    IO_types        io_flag,
    File_formats    format,
    const char*    *str )
{
    int          length;
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
            status = input_quoted_string( file, str );
        else
            status = output_quoted_string( file, *str );
    }
    else
    {
        if( io_flag == WRITE_FILE )
            length = string_length( *str );

        status = io_int( file, io_flag, format, &length );

        if( io_flag == READ_FILE )
            *str = alloc_string( length );

        if( status == OK )
        {
            status = io_binary_data( file, io_flag, (void *) (*str),
                                     sizeof((*str)[0]), length );
        }

        str[length] = END_OF_STRING;
    }

    if( status != OK )
        fprintf(stderr, "Error in quoted string in file.\n" );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_boolean
@INPUT      : file
            : io_flag
            : format
            : b              boolean value
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary boolean value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_boolean(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    VIO_BOOL         *b )
{
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
            status = input_boolean( file, b );
        else
            status = output_boolean( file, *b );
    }
    else
        status = io_binary_data( file, io_flag, (void *) b, sizeof(*b), 1 );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_short
@INPUT      : file
            : io_flag
            : format
            : short_int              short value
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary short value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_short(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    short           *short_int )
{
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
            status = input_short( file, short_int );
        else
            status = output_short( file, *short_int );
    }
    else
        status = io_binary_data( file, io_flag, (void *) short_int,
                                 sizeof(*short_int), 1 );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_unsigned_short
@INPUT      : file
            : io_flag
            : format
            : unsigned_short              short value
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary unsigned short value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_unsigned_short(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    unsigned short  *unsigned_short )
{
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
            status = input_unsigned_short( file, unsigned_short );
        else
            status = output_unsigned_short( file, *unsigned_short );
    }
    else
        status = io_binary_data( file, io_flag, (void *) unsigned_short,
                                 sizeof(*unsigned_short), 1 );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_unsigned_char
@INPUT      : file
            : io_flag
            : format
            : c              unsigned char value
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary unsigned char.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_unsigned_char(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    unsigned  char  *c )
{
    int      i;
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
        {
            if( fscanf( file, "%d", &i ) == 1 )
                *c = (unsigned char) i;
            else
            {
                fprintf(stderr, "Error inputting unsigned char.  " );
                print_system_error();
                status = ERROR;
            }
        }
        else
        {
            if( fprintf( file, "%d", (int) *c ) != 1 )
            {
                fprintf(stderr, "Error outputting unsigned char.  " );
                print_system_error();
                status = ERROR;
            }
        }
    }
    else
        status = io_binary_data( file, io_flag, (void *) c, sizeof(*c), 1 );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_int
@INPUT      : file
            : io_flag
            : format
            : i              integer value
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary integer value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_int(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             *i )
{
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
            status = input_int( file, i );
        else
            status = output_int( file, *i );
    }
    else
        status = io_binary_data( file, io_flag, (void *) i, sizeof(*i), 1 );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_double
@INPUT      : file
            : io_flag
            : format
            : d              double value
@OUTPUT     :
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs an ascii or binary double value.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_double(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    double          *d )
{
    VIO_Status   status;

    status = OK;

    if( format == ASCII_FORMAT )
    {
        if( io_flag == READ_FILE )
            status = input_double( file, d );
        else
            status = output_double( file, *d );
    }
    else
        status = io_binary_data( file, io_flag, (void *) d, sizeof(*d), 1 );

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_ints
@INPUT      : file
            : io_flag
            : format
            : n               number of ints
            : ints            array of ints
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs a list of ascii or binary integers.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_ints(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    int             *ints[] )
{
    VIO_Status   status;
    int      i;
#define      INTS_PER_LINE   8

    status = OK;

    if( io_flag == READ_FILE )
    {
        ALLOC( *ints, n );
    }

    if( format == ASCII_FORMAT )
    {
        for_less( i, 0, n )
        {
            status = io_int( file, io_flag, format, &(*ints)[i] );

            if( status == OK )
            {
                if( i == n - 1 || (i+1) % INTS_PER_LINE == 0 )
                    status = io_newline( file, io_flag, format );
            }

            if( status == ERROR )
                break;
        }
    }
    else
    {
        status = io_binary_data( file, io_flag, (void *) *ints,
                                 sizeof((*ints)[0]), n );
    }

    return( status );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : io_unsigned_chars
@INPUT      : file
            : io_flag
            : format
            : n               number of unsigned chars
            : unsigned_chars  array of unsigned chars
@OUTPUT     : 
@RETURNS    : VIO_Status
@DESCRIPTION: Inputs or outputs a list of ascii or binary unsigned chars.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :                      David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Status  io_unsigned_chars(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    unsigned char   *unsigned_chars[] )
{
    VIO_Status   status;
    int      i;

    status = OK;

    if( io_flag == READ_FILE )
        ALLOC( *unsigned_chars, n );

    if( format == ASCII_FORMAT )
    {
        for_less( i, 0, n )
        {
            status = io_unsigned_char( file, io_flag, format,
                                       &(*unsigned_chars)[i] );

            if( status == OK )
            {
                if( i == n - 1 || (i+1) % INTS_PER_LINE == 0 )
                    status = io_newline( file, io_flag, format );
            }

            if( status == ERROR )
                break;
        }
    }
    else
    {
        status = io_binary_data( file, io_flag, (void *) (*unsigned_chars),
                                 sizeof((*unsigned_chars)[0]), n );
    }

    return( status );
}

#endif // BEVIN_EXCLUDE_MINC
