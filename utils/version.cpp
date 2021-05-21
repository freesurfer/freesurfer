#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>

#include "version.h"
#include "version_info.h"
#include "const.h"
#include "error.h"
#include "utils.h"


// set the platform
#ifdef __linux__
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

// set the compiler name
#if defined(__INTEL_COMPILER)
  #undef __GNUC__
  #define COMPILER_NAME "Intel"
#elif defined(__clang__)
  #define COMPILER_NAME "Clang"
#elif defined(__GNUC__)
  #define COMPILER_NAME "GCC"
#else
  #define COMPILER_NAME "Unknown"
#endif

static std::string ver_string(int a, int b, int c) {
  std::ostringstream ss;
  ss << a << '.' << b << '.' << c;
  return ss.str();
}

// set the compiler version
#if defined(__clang__)
  #define COMPILER_VERSION ver_string(__clang_major__, __clang_minor__, __clang_patchlevel__)
#elif defined(__GNUC__)
  #if defined(__GNU_PATCHLEVEL__)
    #define COMPILER_VERSION ver_string(__GNUC__, __GNUC_MINOR__)
  #else
    #define COMPILER_VERSION ver_string(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)
  #endif
#else
  #define COMPILER_VERSION 0
#endif


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



/*
  Returns the freesurfer version string.
*/
std::string getVersion()
{
  return FS_VERSION;
}


/*
  Returns the freesurfer build stamp ID.
*/
std::string getBuildStamp()
{
  return FS_BUILD_STAMP;
}


/*
  This function looks for the --version, or -version tag in the
  argv and if found, prints out version information. This can be used
  in any binary. It will return the number of options processed and
  copy the remaining items in argv back, so the caller can shorten
  their own argc variable, like this:

  nargs = handleVersionOption(argc, argv, "dollarIDdollar");
  argc -= nargs;

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

  nargs = handleVersionOption(argc, argv, "dollarIDdollar");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
*/
int handleVersionOption(int argc, char **argv, const char *progname)
{
  char *myarg;
  int num_processed_args = 0;

  for (int narg = 1; narg < argc; narg++) {
    char *option = argv[narg];

    if (!strncmp(option, "--version", 9) || !strncmp(option, "-version", 8)) {
      std::cout << progname << " freesurfer " << FS_VERSION << std::endl;
    }
    else if (!strncmp(option, "--all-info", 11) || !strncmp(option, "-all-info", 10)) {
      std::cout << getAllInfo(argc, argv, progname) << std::endl;
    }
    else {
      continue;
    }

    // copy later args one step back if processed the flag
    num_processed_args++;
    for (int nnarg = narg; nnarg < argc - num_processed_args; nnarg++) {
      myarg = (char *)malloc(strlen(argv[nnarg + 1]) + 1);
      strcpy(myarg, argv[nnarg + 1]);
      argv[nnarg] = myarg;
    }
  }

  // return the number of arguments processed
  return num_processed_args;
}


std::string getAllInfo(int argc, char **argv, const std::string& progname)
{
  // get current time
  char current_time_stamp[1024];
  struct tm broken_time;
  time_t seconds = time(NULL);
  gmtime_r(&seconds, &broken_time);
  sprintf(current_time_stamp,
    "20%02d/%02d/%02d-%02d:%02d:%02d-GMT",
    broken_time.tm_year % 100,
    broken_time.tm_mon + 1,
    broken_time.tm_mday,
    broken_time.tm_hour,
    broken_time.tm_min,
    broken_time.tm_sec);

  // TODO this timestamp really ought to be passed-in as a parameter
  // from the calling binary, to get a more accurate build time, but for
  // now, the build time of libutils is better than nothing
  char build_time[] = __DATE__ " " __TIME__;

  // get username
  char username[1024];
  struct passwd *pw = getpwuid(geteuid());
  if ((pw != NULL) && (pw->pw_name != NULL)) {
    strcpy(username, pw->pw_name);
  } else {
    strcpy(username, "UNKNOWN");
  }

  // get current machine name and platform version
  char machine[1024];
  char platform_version[1024];
  struct utsname kernel_info;
  if (uname(&kernel_info) != 0) {
    strcpy(machine, "UNKNOWN");
    strcpy(platform_version, "UNKNOWN");
  } else {
    strcpy(machine, kernel_info.nodename);
    strcpy(platform_version, kernel_info.release);
  }

  // get args as string
  std::ostringstream argvss;
  argvss << argv[0];
  for (int i = 1; i < argc ; i++) argvss << " " << argv[i];

  std::ostringstream infoss;
  infoss << "ProgramName: " << progname
         << "  ProgramArguments: " << argvss.str()
         << "  ProgramVersion: " << getVersion()
         << "  TimeStamp: " << current_time_stamp
         << "  BuildTime: " << build_time
         << "  BuildStamp: " << getBuildStamp()
         << "  User: " << username
         << "  Machine: " << machine
         << "  Platform: " << PLATFORM
         << "  PlatformVersion: " << platform_version
         << "  CompilerName: " << COMPILER_NAME
         << "  CompilerVersion: " << COMPILER_VERSION;

  return infoss.str();
}
