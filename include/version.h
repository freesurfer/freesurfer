// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author $
// Revision Date  : $Date $
// Revision       : $Revision $


#ifndef VERSION_H
#define VERSION_H

/* This function looks for the -v, --version, or -version tag in the
   argv and if found, prints out version information. This can be used
   in any binary. It will return the number of argc/argv options
   processed, so the caller can shorten their own argc/argv variables,
   like this:

    nargs = handle_version_option (argc, argv, "$Id: version.h,v 1.1 2003/03/19 18:04:51 kteich Exp $");
    argc -= nargs ;
    argv += nargs ;

   It will print out the id string passed in, which should just be
   "$Id", which CVS will expand to include the CVS information
   including CVS file, revision number, date, peson who checked it in,
   and tag, as well as the OS and GCC information.
*/

int handle_version_option (int argc, char** argv, char* id_string);


#endif
