// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author $
// Revision Date  : $Date $
// Revision       : $Revision $


#ifndef VERSION_H
#define VERSION_H

/* This function looks for the -v, --version, or -version tag in the
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

    if (1 == argc)
      exit (0);
*/

int handle_version_option (int argc, char** argv, char* id_string);


#endif
