/**
 * @file  version.h
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
 *    $Date: 2011/03/02 00:04:10 $
 *    $Id$
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


#ifndef VERSION_H
#define VERSION_H

/* This function looks for the --version, or -version tag in the
   argv and if found, prints out version information. This can be used
   in any binary. It will return the number of options processed and
   copy the remaining items in argv back, so the caller can shorten
   their own argc variable, like this:

    nargs = handle_version_option (argc, argv,
                                   "dollarIDdollar", "dollarNamedollar");
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

    nargs = handle_version_option (argc, argv,
                                   "dollarIDdollar", "dollarNamedollar");
    if (nargs && argc - nargs == 1)
      exit (0);
    argc -= nargs;

   This will also handle the --all-info option which prints out
   more information in the format necessary for the BIRN provenance
   spec, including the version_string.

*/

int handle_version_option(int argc, char** argv, const char* id_string, const char* version_string);
int make_cmd_version_string(int argc, char** argv,  const char* id_string, const char* version_string, char *return_string);
char *argv2cmdline(int argc, char *argv[]);
char *VERuser(void);
char *VERfileTimeStamp(char *fname);
char *VERcurTimeStamp(void);

#endif
