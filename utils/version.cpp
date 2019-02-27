/**
 * @file  version.c
 * @brief freesurfer version functions defined here
 *
 * Looks for the --version, -version, --all-info, or -all-info tag in the
 * argv and if found, prints out version information, namely this:
 * ProgramVersion, TimeStamp, CVS, User, Machine, Platform, PlatformVersion
 * CompilerName, and CompilerVersion.
 */
/*
 * Original Author: Kevin Teich
 * Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 8c0544cb51c408c315be1f6cae68f38fc534b7da $
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

#include "version.h"
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>
#include "const.h"
#include "error.h"
#include "utils.h"

/* Set our compiler name */
#if defined(__INTEL_COMPILER)
#undef __GNUC__
#define COMPILER_NAME "INTEL_COMPILER"
#elif defined(__GNUC__)
#define COMPILER_NAME "GCC"
#else
#define COMPILER_NAME "Non-GCC"
#endif

/* If GCC (probably is) get the version number */
#if defined(__GNUC__)
#if defined(__GNU_PATCHLEVEL__)
#define COMPILER_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#else
#define COMPILER_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100)
#endif
#else
#define COMPILER_VERSION 0
#endif

/* Figure out the platform. */
#if defined(Linux) || defined(linux) || defined(__linux)
#define PLATFORM "Linux"
#endif
#if defined(Darwin) || defined(__MACOSX__) || defined(__APPLE__)
#define PLATFORM "Darwin"
#endif
#if defined(IRIX) || defined(sgi) || defined(mips) || defined(_SGI_SOURCE)
#define PLATFORM "IRIX"
#endif
#if defined(sun) || defined(__sun)
#if defined(__SVR4) || defined(__svr4__)
#define PLATFORM "Solaris"
#else
#define PLATFORM "SunOS"
#endif
#endif
#if defined(Windows_NT)
#define PLATFORM "Windows"
#endif

#ifndef PLATFORM
#error "PLATFORM not defined!"
#endif

/* This function looks for the --version, or -version tag in the
   argv and if found, prints out version information. This can be used
   in any binary. It will return the number of options processed and
   copy the remaining items in argv back, so the caller can shorten
   their own argc variable, like this:

   nargs = handle_version_option (argc, argv, "dollarIDdollar");
   argc -= nargs ;

   (This is done so that it can easily be used with older programs
   that make assumptions about argument counts instead of scanning for
   options.)

   It will print out the id string passed in, which should just be
   dollar sign + ID + dollar sign, which CVS will expand to include
   the CVS information including CVS file, revision number, date,
   peson who checked it in, and tag, as well as the OS and GCC
   information.

   The binary may also want to exit if there are no other options to
   handle, i.e.

   nargs = handle_version_option (argc, argv, "dollarIDdollar");
   if (nargs && argc - nargs == 1)
   exit (0);
   argc -= nargs;
*/
int make_cmd_version_string(
    int argc, char **argv, const char *id_string, const char *version_string, char *return_string)
{
  int nnarg = 0;
  char stripped_version_string[1024];
  int length;
  time_t seconds;
  struct tm broken_time;
  struct utsname kernel_info;
  int result;
  char *begin;
  char program_name[1024];
  char arguments[1024];
  char current_time_stamp[1024];
  char user[1024];
  char machine[1024];
  char platform_version[1024];
  struct passwd *pw;

  if (strlen(version_string) > 7) {
    strcpy(stripped_version_string, &(version_string[7]));
    length = strlen(stripped_version_string);
    if (length > 2) {
      stripped_version_string[length - 2] = '\0';
    }
  }
  else {
    strcpy(stripped_version_string, version_string);
  }

  begin = argv[0];
  strcpy(program_name, begin);

  /* Copy the arguments to the arguments string. */
  strcpy(arguments, "");
  if (argc > 1) {
    strncpy(arguments, argv[1], 1023);
    for (nnarg = 2; nnarg < argc; nnarg++) {
      // on Slackware Linux, libc does not support having the same source and
      // destination, like this:
      // sprintf (arguments, "%s %s", arguments, argv[nnarg]);
      // the correct way to do this is:
      strcat(arguments, " ");
      if (strlen(arguments) + strlen(argv[nnarg]) >= 1023) break;
      strcat(arguments, argv[nnarg]);
    }
  }

  /* Find the time string. */
  seconds = time(NULL);
  gmtime_r(&seconds, &broken_time);
  sprintf(current_time_stamp,
          "20%02d/%02d/%02d-%02d:%02d:%02d-GMT",
          broken_time.tm_year % 100, /* mod here to change 103 to 03 */
          broken_time.tm_mon + 1,    /* +1 here because tm_mon is 0-11 */
          broken_time.tm_mday,
          broken_time.tm_hour,
          broken_time.tm_min,
          broken_time.tm_sec);

  /* As suggested by the getlogin() manpage, use getpwuid(geteuid())
     to get the user controlling this process. don't use getlogin()
     as that returns the name of the user logged in on the controlling
     terminal of the process (ie the person sitting at the terminal),
     and don't use cuserid() because the manpage says not to.
  */
  pw = getpwuid(geteuid());
  if ((pw != NULL) && (pw->pw_name != NULL))
    strcpy(user, pw->pw_name);
  else
    strcpy(user, "UNKNOWN");

  /* Call uname to get the machine. */
  result = uname(&kernel_info);
  if (0 != result) {
    // fprintf (stderr, "uname() returned %d\n", result);
    strcpy(machine, "UNKNOWN");
    strcpy(platform_version, "UNKNOWN");
  }
  else {
    strcpy(machine, kernel_info.nodename);
    strcpy(platform_version, kernel_info.release);
  }

  // TODO: build_timestamp really ought to be passed-in as a parameter
  // from the calling binary, to get a more accurate build time, but for
  // now, the build time of this version.c (libutils) is better than nothing
  char build_timestamp[] = __DATE__ " " __TIME__;
  char argstr[CMD_LINE_LEN];

  if (strlen(arguments) > CMD_LINE_LEN / 2)
    strncpy(argstr, arguments, CMD_LINE_LEN / 2);
  else
    strcpy(argstr, arguments);

  /* Build the info string. */
  sprintf(return_string,
          "%s %s "
          "ProgramVersion: %s  TimeStamp: %s  "
          "BuildTimeStamp: %s  Id: %s  User: %s  "
          "Machine: %s  Platform: %s  PlatformVersion: %s  "
          "CompilerName: %s  CompilerVersion: %d  ",
          program_name,
          argstr,
          version_string,
          current_time_stamp,
          build_timestamp,
          id_string,
          user,
          machine,
          PLATFORM,
          platform_version,
          COMPILER_NAME,
          COMPILER_VERSION);

  return (NO_ERROR);
}

int handle_version_option(int argc, char **argv, const char *id_string, const char *version_string)
{
  int narg = 0;
  int nnarg = 0;
  int num_processed_args = 0;
  char stripped_version_string[1024];
  int length;
  char *option = NULL;
  time_t seconds;
  struct tm broken_time;
  struct utsname kernel_info;
  int result;
  char *begin;
  char program_name[1024];
  char arguments[1024];
  char current_time_stamp[1024];
  char user[1024];
  char machine[1024];
  char platform_version[1024];
  struct passwd *pw;
  char *myarg;

  /* Go through each option looking for --version, -version,
     --all-info, or -all-info */
  for (narg = 1; narg < argc; narg++) {
    option = argv[narg];

    if (!strncmp(option, "--version", 9) || !strncmp(option, "-version", 8)) {
#if 0
      /* Since BIRN is now using --all-info to get version stuff,
         I want to keep this as simple as possible. So I'm
         commenting out some of this stuff, let's see if anybody
         complains. */

      /* Print out the entire command line. */
      for (nnarg = 0; nnarg < argc; nnarg++)
        fprintf (stdout, "%s ", argv[nnarg]);
      fprintf (stdout, "\n");

      fprintf (stdout, "%s Platform: %s C lib: %d\n",
               id_string, PLATFORM, COMPILER_VERSION);
#else
      /* Strip the "{Name: " and "}" from the id string. */
      if (strlen(version_string) > 7) {
        strcpy(stripped_version_string, &(version_string[7]));
        length = strlen(stripped_version_string);
        if (length > 2) {
          stripped_version_string[length - 2] = '\0';
        }
      }
      else {
        strcpy(stripped_version_string, version_string);
      }

      if (strcmp(" $", stripped_version_string) == 0) {
        // on the dev build, where a sticky tag does not exist,
        // just a dollar sign is printed, which isnt very helpful,
        // so print something...
        strcpy(stripped_version_string, "dev build (use --all-info flag for full version info)");
      }
      fprintf(stdout, "%s\n", stripped_version_string);

#endif
      num_processed_args++;

      /* Copy later args one step back. */
      for (nnarg = narg; nnarg < argc - num_processed_args; nnarg++) {
        myarg = (char *)malloc(strlen(argv[nnarg + 1]) + 1);
        strcpy(myarg, argv[nnarg + 1]);
        argv[nnarg] = myarg;
      }
    }

    if (!strncmp(option, "--all-info", 11) || !strncmp(option, "-all-info", 10)) {
      /* Copy argv[0] without the path into program_name. */
      begin = strrchr(argv[0], (int)'/');
      if (NULL == begin) {
        begin = argv[0];
      }
      else {
        begin = begin + 1;
      }
      strcpy(program_name, begin);

      /* Copy the arguments to the arguments string. */
      strcpy(arguments, "");
      if (argc > 1) {
        strcpy(arguments, argv[1]);
        for (nnarg = 2; nnarg < argc; nnarg++) {
          // on Slackware Linux, libc does not support having the same
          // source and destination, like this:
          // sprintf (arguments, "%s %s", arguments, argv[nnarg]);
          // the correct way to do this is:
          strcat(arguments, " ");
          strcat(arguments, argv[nnarg]);
        }
      }

      /* Find the time string. */
      seconds = time(NULL);
      gmtime_r(&seconds, &broken_time);
      sprintf(current_time_stamp,
              "20%02d/%02d/%02d-%02d:%02d:%02d-GMT",
              broken_time.tm_year % 100, /* mod here to change 103 to 03 */
              broken_time.tm_mon + 1,    /* +1 here because tm_mon is 0-11 */
              broken_time.tm_mday,
              broken_time.tm_hour,
              broken_time.tm_min,
              broken_time.tm_sec);

      /* As suggested by the getlogin() manpage, use getpwuid(geteuid())
         to get the user controlling this process. don't use getlogin()
         as that returns the name of the user logged in on the controlling
         terminal of the process (ie the person sitting at the terminal),
         and don't use cuserid() because the manpage says not to.
      */
      pw = getpwuid(geteuid());
      if ((pw != NULL) && (pw->pw_name != NULL))
        strcpy(user, pw->pw_name);
      else
        strcpy(user, "UNKNOWN");

      /* Call uname to get the machine. */
      result = uname(&kernel_info);
      if (0 != result) {
        // fprintf (stderr, "uname() returned %d\n", result);
        strcpy(machine, "UNKNOWN");
        strcpy(platform_version, "UNKNOWN");
      }
      else {
        strcpy(machine, kernel_info.nodename);
        strcpy(platform_version, kernel_info.release);
      }

      // TODO: build_timestamp really ought to be passed-in as a parameter
      // from the calling binary, to get a more accurate build time, but for
      // now, the build time of version.c (libutils) is better than nothing
      char build_timestamp[] = __DATE__ " " __TIME__;

      /* Build the info string. */
      fprintf(stdout,
              "ProgramName: %s  ProgramArguments: %s  "
              "ProgramVersion: %s  TimeStamp: %s  "
              "BuildTimeStamp: %s  Id: %s  User: %s  "
              "Machine: %s  Platform: %s  PlatformVersion: %s  "
              "CompilerName: %s  CompilerVersion: %d \n",
              program_name,
              arguments,
              version_string,
              current_time_stamp,
              build_timestamp,
              id_string,
              user,
              machine,
              PLATFORM,
              platform_version,
              COMPILER_NAME,
              COMPILER_VERSION);

      num_processed_args++;

      /* Copy later args one step back. */
      for (nnarg = narg; nnarg < argc - num_processed_args; nnarg++) {
        myarg = (char *)malloc(strlen(argv[nnarg + 1]) + 1);
        strcpy(myarg, argv[nnarg + 1]);
        argv[nnarg] = myarg;
      }
    }
  }

  /* Return the number of arguments processed. */
  return num_processed_args;
}

/*------------------------------------------------------------------------
  argv2cmdline() - converts argv into a single string.
  ------------------------------------------------------------------------*/
char *argv2cmdline(int argc, char *argv[])
{
  int len, n;
  char *cmdline;

  len = 0;
  for (n = 0; n < argc; n++) len += strlen(argv[n]);

  cmdline = (char *)calloc(sizeof(char), len + n + 1);
  for (n = 0; n < argc; n++) {
    strcat(cmdline, argv[n]);
    strcat(cmdline, " ");
  }

  return (cmdline);
}

/*----------------------------------------------------------
  VERuser(void) - returns the user id (or UNKNOWN). Allocates
  char, so caller must free.
  *----------------------------------------------------------*/
char *VERuser(void)
{
  char *user;
  struct passwd *pw;

  /* As suggested by the getlogin() manpage, use getpwuid(geteuid())
     to get the user controlling this process. don't use getlogin()
     as that returns the name of the user logged in on the controlling
     terminal of the process (ie the person sitting at the terminal),
     and don't use cuserid() because the manpage says not to.
  */
  pw = getpwuid(geteuid());
  if ((pw != NULL) && (pw->pw_name != NULL))
    user = strcpyalloc(pw->pw_name);
  else
    user = strcpyalloc("UNKNOWN");

  return (user);
}

/*----------------------------------------------------------
  VERfileTimeStamp(fname)
  *----------------------------------------------------------*/
char *VERfileTimeStamp(char *fname)
{
  struct stat buf;
  struct tm *lt = NULL;
  int err;
  char *timestamp, tmpstr[1000];

  err = stat(fname, &buf);
  if (err) {
    printf("ERROR: stating file %s\n", fname);
    return (NULL);
  }
  lt = localtime(&buf.st_mtime);
  sprintf(tmpstr,
          "%04d/%02d/%02d %02d:%02d:%02d",
          lt->tm_year + 1900,
          lt->tm_mon + 1, /* +1 here because tm_mon is 0-11 */
          lt->tm_mday,
          lt->tm_hour,
          lt->tm_min,
          lt->tm_sec);
  // free(lt); // Dies here
  timestamp = strcpyalloc(tmpstr);
  return (timestamp);
}

/*----------------------------------------------------------
  VERcurTimeStamp() - time stamp at the time this function
  is called.
  *----------------------------------------------------------*/
char *VERcurTimeStamp(void)
{
  char tmpstr[2000], *current_time_stamp;
  time_t seconds;
  struct tm broken_time;
  seconds = time(NULL);
  gmtime_r(&seconds, &broken_time);
  sprintf(tmpstr,
          "20%02d/%02d/%02d-%02d:%02d:%02d-GMT",
          broken_time.tm_year % 100, /* mod here to change 103 to 03 */
          broken_time.tm_mon + 1,    /* +1 here because tm_mon is 0-11 */
          broken_time.tm_mday,
          broken_time.tm_hour,
          broken_time.tm_min,
          broken_time.tm_sec);
  // Note: NOT Y3K compliant!
  current_time_stamp = strcpyalloc(tmpstr);
  return (current_time_stamp);
}
