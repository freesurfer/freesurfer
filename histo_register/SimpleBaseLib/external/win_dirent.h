/*
 * dirent.h - operating system independent dirent implementation
 * 
 * Copyright (C) 1998-2002  Toni Ronkko
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * ``Software''), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED ``AS IS'', WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL TONI RONKKO BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 * 
 * 
 * Aug 20, 2006, Toni Ronkko
 * Removed all remarks about MSVC 1.0, which is antiqued now.  Simplified
 * comments by removing SGML tags.
 *
 * May 14 2002, Toni Ronkko
 * Embedded the function definitions directly to the header so that no source
 * modules must be included in the MS Visual C project for using the
 * interface.  Removed all the dependencies to other projects so
 * that this very header can be used independently.
 *
 * May 28 1998, Toni Ronkko
 * First version.
 */
#ifndef DIRENT_H
#define DIRENT_H
#define DIRENT_H_INCLUDED

/* find out platform */
#if defined(MSDOS)                             /* MS-DOS */
#elif defined(__MSDOS__)                       /* Turbo C/Borland */
# define MSDOS
#elif defined(__DOS__)                         /* Watcom */
# define MSDOS
#endif

#if defined(WIN32)                             /* MS-Windows */
#elif defined(__NT__)                          /* Watcom */
# define WIN32
#elif defined(_WIN32)                          /* Microsoft */
# define WIN32
#elif defined(__WIN32__)                       /* Borland */
# define WIN32
#endif

/*
 * See what kind of dirent interface we have unless autoconf has already
 * determinated that.
 */
#if !defined(HAVE_DIRENT_H) && !defined(HAVE_DIRECT_H) && !defined(HAVE_SYS_DIR_H) && !defined(HAVE_NDIR_H) && !defined(HAVE_SYS_NDIR_H) && !defined(HAVE_DIR_H)
# if defined(_MSC_VER)                         /* Microsoft C/C++ */
    /* no dirent.h */
# elif defined(__BORLANDC__)                   /* Borland C/C++ */
#   define HAVE_DIRENT_H
#   define VOID_CLOSEDIR
# elif defined(__TURBOC__)                     /* Borland Turbo C */
    /* no dirent.h */
# elif defined(__WATCOMC__)                    /* Watcom C/C++ */
#   define HAVE_DIRECT_H
# elif defined(__apollo)                       /* Apollo */
#   define HAVE_SYS_DIR_H
# elif defined(__hpux)                         /* HP-UX */
#   define HAVE_DIRENT_H
# elif defined(__alpha) || defined(__alpha__)  /* Alpha OSF1 */
#   error "not implemented"
# elif defined(__sgi)                          /* Silicon Graphics */
#   define HAVE_DIRENT_H
# elif defined(sun) || defined(_sun)           /* Sun Solaris */
#   define HAVE_DIRENT_H
# elif defined(__FreeBSD__)                    /* FreeBSD */
#   define HAVE_DIRENT_H
# elif defined(__linux__)                      /* Linux */
#   define HAVE_DIRENT_H
# elif defined(__GNUC__)                       /* GNU C/C++ */
#   define HAVE_DIRENT_H
# else
#   error "not implemented"
# endif
#endif

/* include proper interface headers */
#if defined(HAVE_DIRENT_H)
# include <dirent.h>
# ifdef FREEBSD
#   define NAMLEN(dp) ((int)((dp)->d_namlen))
# else
#   define NAMLEN(dp) ((int)(strlen((dp)->d_name)))
# endif

#elif defined(HAVE_NDIR_H)
# include <ndir.h>
# define NAMLEN(dp) ((int)((dp)->d_namlen))

#elif defined(HAVE_SYS_NDIR_H)
# include <sys/ndir.h>
# define NAMLEN(dp) ((int)((dp)->d_namlen))

#elif defined(HAVE_DIRECT_H)
# include <direct.h>
# define NAMLEN(dp) ((int)((dp)->d_namlen))

#elif defined(HAVE_DIR_H)
# include <dir.h>
# define NAMLEN(dp) ((int)((dp)->d_namlen))

#elif defined(HAVE_SYS_DIR_H)
# include <sys/types.h>
# include <sys/dir.h>
# ifndef dirent
#   define dirent direct
# endif
# define NAMLEN(dp) ((int)((dp)->d_namlen))

#elif defined(MSDOS) || defined(WIN32)

  /* figure out type of underlaying directory interface to be used */
# if defined(WIN32)
#   define DIRENT_WIN32_INTERFACE
# elif defined(MSDOS)
#   define DIRENT_MSDOS_INTERFACE
# else
#   error "missing native dirent interface"
# endif

  /*** WIN32 specifics ***/
# if defined(DIRENT_WIN32_INTERFACE)
#   include <windows.h>
#   if !defined(DIRENT_MAXNAMLEN)
#     define DIRENT_MAXNAMLEN (MAX_PATH)
#   endif


  /*** MS-DOS specifics ***/
# elif defined(DIRENT_MSDOS_INTERFACE)
#   include <dos.h>

    /* Borland defines file length macros in dir.h */
#   if defined(__BORLANDC__)
#     include <dir.h>
#     if !defined(DIRENT_MAXNAMLEN)
#       define DIRENT_MAXNAMLEN ((MAXFILE)+(MAXEXT))
#     endif
#     if !defined(_find_t)
#       define _find_t find_t
#     endif

    /* Turbo C defines ffblk structure in dir.h */
#   elif defined(__TURBOC__)
#     include <dir.h>
#     if !defined(DIRENT_MAXNAMLEN)
#       define DIRENT_MAXNAMLEN ((MAXFILE)+(MAXEXT))
#     endif
#     define DIRENT_USE_FFBLK

    /* MSVC */
#   elif defined(_MSC_VER)
#     if !defined(DIRENT_MAXNAMLEN)
#       define DIRENT_MAXNAMLEN (12)
#     endif

    /* Watcom */
#   elif defined(__WATCOMC__)
#     if !defined(DIRENT_MAXNAMLEN)
#       if defined(__OS2__) || defined(__NT__)
#         define DIRENT_MAXNAMLEN (255)
#       else
#         define DIRENT_MAXNAMLEN (12)
#       endif
#     endif

#   endif
# endif

  /*** generic MS-DOS and MS-Windows stuff ***/
# if !defined(NAME_MAX) && defined(DIRENT_MAXNAMLEN)
#   define NAME_MAX DIRENT_MAXNAMLEN
# endif
# if NAME_MAX < DIRENT_MAXNAMLEN
#   error "assertion failed: NAME_MAX >= DIRENT_MAXNAMLEN"
# endif


  /*
   * Substitute for real dirent structure.  Note that d_name field is a
   * true character array although we have it copied in the implementation
   * dependent data.  We could save some memory if we had declared d_name
   * as a pointer refering the name within implementation dependent data.
   * We have not done that since some code may rely on sizeof(d_name) to be
   * something other than four.  Besides, directory entries are typically so
   * small that it takes virtually no time to copy them from place to place.
   */
  typedef struct dirent {
    char d_name[NAME_MAX + 1];

    /*** Operating system specific part ***/
# if defined(DIRENT_WIN32_INTERFACE)       /*WIN32*/
    WIN32_FIND_DATA data;
# elif defined(DIRENT_MSDOS_INTERFACE)     /*MSDOS*/
#   if defined(DIRENT_USE_FFBLK)
    struct ffblk data;
#   else
    struct _find_t data;
#   endif
# endif
  } dirent;

  /*
   * DIR substitute structure containing directory name.  The name is
   * essential for the operation of rewinndir() function.
   */
  typedef struct DIR {
    char          *dirname;                    /* directory being scanned */
    dirent        current;                     /* current entry */
    int           dirent_filled;               /* is current un-processed? */

  /*** Operating system specific part ***/
#  if defined(DIRENT_WIN32_INTERFACE)
    HANDLE        search_handle;
#  elif defined(DIRENT_MSDOS_INTERFACE)
#  endif
  } DIR;

# ifdef __cplusplus
extern "C" {
# endif

/* prototypes of public dirent functions */
static DIR *opendir (const char *dirname);
static struct dirent *readdir (DIR *dirp);
static int closedir (DIR *dirp);
static void rewinddir (DIR *dirp);

/*
 * Implement dirent interface as static functions so that the user does not
 * need to change his project in any way to use dirent function.  With this
 * it is sufficient to include this very header from source modules using
 * dirent functions and the functions will be pulled in automatically.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

/* use ffblk instead of _find_t if requested */
#if defined(DIRENT_USE_FFBLK)
# define _A_ARCH   (FA_ARCH)
# define _A_HIDDEN (FA_HIDDEN)
# define _A_NORMAL (0)
# define _A_RDONLY (FA_RDONLY)
# define _A_SUBDIR (FA_DIREC)
# define _A_SYSTEM (FA_SYSTEM)
# define _A_VOLID  (FA_LABEL)
# define _dos_findnext(dest) findnext(dest)
# define _dos_findfirst(name,flags,dest) findfirst(name,dest,flags)
#endif

/* internal prototypes */
static int _initdir (DIR *p);
static const char *_getdirname (const struct dirent *dp);
static void _setdirname (struct DIR *dirp);


/*
 * Open directory stream DIRNAME for read and return pointer to the
 * internal working area that is used for retrieving individual directory
 * entries.  The internal working area has no fields of your interest.
 */
static DIR *
opendir(
    const char *dirname)
{
  DIR *dirp;
  assert (dirname != NULL);
  
  dirp = (DIR*)malloc (sizeof (struct DIR));
  if (dirp != NULL) {
    char *p;
    
    /* allocate room for directory name */
    dirp->dirname = (char*) malloc (strlen (dirname) + 1 + strlen ("\\*.*"));
    if (dirp->dirname == NULL) {
      /* failed to duplicate directory name.  errno set by malloc() */
      free (dirp);
      return NULL;
    }
    
    /*
     * Copy directory name while appending directory separator and "*.*".
     * Directory separator is not appended if the name already ends with
     * drive or directory separator.  Directory separator is assumed to be
     * '/' or '\' and drive separator is assumed to be ':'.
     */
    strcpy (dirp->dirname, dirname);
    p = strchr (dirp->dirname, '\0');
    if (dirp->dirname < p  &&
        *(p - 1) != '\\'  &&  *(p - 1) != '/'  &&  *(p - 1) != ':')
    {
      strcpy (p++, "\\");
    }
# ifdef DIRENT_WIN32_INTERFACE
    strcpy (p, "*"); /*scan files with and without extension in win32*/
# else
    strcpy (p, "*.*"); /*scan files with and without extension in DOS*/
# endif

    /* open stream */
    if (_initdir (dirp) == 0) {
      /* initialization failed */
      free (dirp->dirname);
      free (dirp);
      return NULL;
    }
  }
  return dirp;
}

  

/*
 * Read a directory entry and return pointer to a structure containing the
 * name of the entry.  Individual directory entries returned by this very
 * function include normal files, sub-directories, pseudo-directories "."
 * and ".." and also volume labels, hidden files and system files in MS-DOS
 * and MS-Windows.   You might want to use stat(2) function to determinate
 * which one you are dealing with.
 *
 * The dirent structure contains several system dependent fields that
 * generally have no interest to you.  The only interesting one is char
 * d_name[] that is also portable across different systems.  The d_name
 * field contains the name of the directory entry without leading path.
 * While d_name is portable across different systems, the actual storage
 * capacity of d_name varies from system to system and there is no portable
 * way to find out it at compile time.  Different systems define the
 * capacity of d_name with different macros and some systems do not define
 * capacity at all (besides actual declaration of the field).  If you really
 * need to find out storage capacity of d_name then you might want to try
 * NAME_MAX macro.  The NAME_MAX is defined in POSIX standard althought
 * there are many MS-DOS and MS-Windows implementations those do not define
 * it.  There are also systems that declare d_name as "char d_name[1]" and
 * then allocate suitable amount of memory at run-time.  Thanks to Alain
 * Decamps (Alain dot Decamps atsign advalvas dot be) for pointing it out to
 * me.
 * 
 * This all leads to the fact that it is difficult to allocate space
 * for the directory names when the very same program is being compiled on
 * number of operating systems.  Therefore I suggest that you always
 * allocate space for directory names dynamically.
 */
static struct dirent *
readdir (DIR *dirp)
{
  assert (dirp != NULL);
  if (dirp == NULL) {
    errno = EBADF;
    return NULL;
  }

#if defined(DIRENT_WIN32_INTERFACE)
  if (dirp->search_handle == INVALID_HANDLE_VALUE) {
    /* directory stream was opened/rewound incorrectly or it ended normally */
    errno = EBADF;
    return NULL;
  }
#endif

  if (dirp->dirent_filled != 0) {
    /*
     * Directory entry has already been retrieved and there is no need to
     * retrieve a new one.  Directory entry will be retrieved in advance
     * when the user calls readdir function for the first time.  This is so
     * because real dirent has separate functions for opening and reading
     * the stream whereas Win32 and DOS dirents open the stream
     * automatically when we retrieve the first file.  Therefore, we have to
     * save the first file when opening the stream and later we have to
     * return the saved entry when the user tries to read the first entry.
     */
    dirp->dirent_filled = 0;
  } else {
    /* fill in entry and return that */
#if defined(DIRENT_WIN32_INTERFACE)
    if (FindNextFile (dirp->search_handle, &dirp->current.data) == FALSE) {
      /* Last file has been processed or an error occured */
      FindClose (dirp->search_handle);
      dirp->search_handle = INVALID_HANDLE_VALUE;
      errno = ENOENT;
      return NULL;
    }

# elif defined(DIRENT_MSDOS_INTERFACE)
    if (_dos_findnext (&dirp->current.data) != 0) {
      /* _dos_findnext and findnext will set errno to ENOENT when no
       * more entries could be retrieved. */
      return NULL;
    }
# endif

    _setdirname (dirp);
    assert (dirp->dirent_filled == 0);
  }
  return &dirp->current;
}



/*
 * Close directory stream opened by opendir() function.  Close of the
 * directory stream invalidates the DIR structure as well as any previously
 * read directory entry.
 *
 * Usually, the function returns 0 on success and -1 on failure.  However,
 * on some systems, the function has been declared to return void.  If you
 * want to keep your code portable, it may be better to ignore the return
 * value altogether.
 */
static int
closedir (DIR *dirp)
{   
  int retcode = 0;

  /* make sure that dirp points to legal structure */
  assert (dirp != NULL);
  if (dirp == NULL) {
    errno = EBADF;
    return -1;
  }
 
  /* free directory name */
  if (dirp->dirname != NULL) {
    free (dirp->dirname);
  }

  /* release search handle */
#if defined(DIRENT_WIN32_INTERFACE)
  if (dirp->search_handle != INVALID_HANDLE_VALUE) {
    if (FindClose (dirp->search_handle) == FALSE) {
      /* Unknown error */
      retcode = -1;
      errno = EBADF;
    }
  }
#endif                     

  /* clear dirp structure to make sure that it cannot be used anymore*/
  memset (dirp, 0, sizeof (*dirp));
# if defined(DIRENT_WIN32_INTERFACE)
  dirp->search_handle = INVALID_HANDLE_VALUE;
# endif

  free (dirp);
  return retcode;
}



/*
 * Rewind the directory stream to re-read the directory entries.  If new
 * files or directories have been created since opendir(), then also these
 * new files will be available with readdir().
 *
 * Almost every dirent implementation ensures that rewinddir() updates its
 * internal caches to reflect the current directory structure on disk.
 * However, if you really depend on this feature, it may be better to close
 * and re-open the directory stream.  I recall one particual dirent
 * interface, where rewinddir() did not update its internal cache.
 *
 * BUGS: If the current working directory is changed between opendir() and
 * rewinddir(), then the rewound directory stream may refer a non-existent
 * directory, or at worst case, to another directory.
 */
static void
rewinddir (DIR *dirp)
{   
  /* make sure that dirp is legal */
  assert (dirp != NULL);
  if (dirp == NULL) {
    errno = EBADF;
    return;
  }
  assert (dirp->dirname != NULL);
  
  /* close previous stream */
#if defined(DIRENT_WIN32_INTERFACE)
  if (dirp->search_handle != INVALID_HANDLE_VALUE) {
    if (FindClose (dirp->search_handle) == FALSE) {
      /* Unknown error */
      errno = EBADF;
    }
  }
#endif

  /* re-open previous stream */
  if (_initdir (dirp) == 0) {
    /* initialization failed but we cannot deal with error.  User will notice
     * error later when she tries to retrieve first directory enty. */
    /*EMPTY*/;
  }
}



/*
 * Open native directory stream object and retrieve first file.
 * Be sure to close previous stream before opening new one.
 */
static int
_initdir (DIR *dirp)
{ 
  assert (dirp != NULL);
  assert (dirp->dirname != NULL);
  dirp->dirent_filled = 0;

# if defined(DIRENT_WIN32_INTERFACE)
  /* Open stream and retrieve first file */
  dirp->search_handle = FindFirstFile (dirp->dirname, &dirp->current.data);
  if (dirp->search_handle == INVALID_HANDLE_VALUE) {
    /* something went wrong but we don't know what.  GetLastError() could
     * give us more information about the error, but then we should map
     * the error code into errno. */
    errno = ENOENT;
    return 0;
  }

# elif defined(DIRENT_MSDOS_INTERFACE)
  if (_dos_findfirst (dirp->dirname,
          _A_SUBDIR | _A_RDONLY | _A_ARCH | _A_SYSTEM | _A_HIDDEN,
          &dirp->current.data) != 0)
  {
    /* _dos_findfirst and findfirst will set errno to ENOENT when no 
     * more entries could be retrieved. */
    return 0;
  }
# endif

  /* initialize DIR and it's first entry */
  _setdirname (dirp);
  dirp->dirent_filled = 1;
  return 1;
}



/*
 * Return implementation dependent name of the current directory entry.
 */
static const char *
_getdirname (const struct dirent *dp)
{
#if defined(DIRENT_WIN32_INTERFACE)
  return dp->data.cFileName;
  
#elif defined(DIRENT_USE_FFBLK)
  return dp->data.ff_name;
  
#else
  return dp->data.name;
#endif  
}



/*
 * Copy name of implementation dependent directory entry to the d_name field.
 */
static void
_setdirname (struct DIR *dirp) {
  /* make sure that d_name is long enough */
  assert (strlen (_getdirname (&dirp->current)) <= NAME_MAX);
  
  strncpy (dirp->current.d_name,
      _getdirname (&dirp->current),
      NAME_MAX);
  dirp->current.d_name[NAME_MAX] = '\0'; /*char d_name[NAME_MAX+1]*/
}
  
# ifdef __cplusplus
}
# endif
# define NAMLEN(dp) ((int)(strlen((dp)->d_name)))

#else
# error "missing dirent interface"
#endif


#endif /*DIRENT_H*/
