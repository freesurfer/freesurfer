#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "version.h"


// char version[] = "version 1.0, Wed Apr 26 11:11:26 EDT 2000" ;


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
int
handle_version_option (int argc, char** argv, char* id_string) 
{
  
#if defined(__GNUC__)
# if defined(__GNU_PATCHLEVEL__)
#  define __GNUC_VERSION__ (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100 \
                            + __GNUC_PATCHLEVEL__)
# else
#  define __GNUC_VERSION__ (__GNUC__ * 10000 \
                            + __GNUC_MINOR__ * 100)
# endif
#endif

#ifdef Linux
#  define __PLATFORM__ "Linux"
#endif
#ifdef Darwin
#  define __PLATFORM__ "Darwin"
#endif
#ifdef IRIX
#  define __PLATFORM__ "IRIX"
#endif
#ifdef Solaris
#  define __PLATFORM__ "Solaris"
#endif

  int narg = 0;
  int nnarg = 0;
  int num_processed_args = 0;
  char *option = NULL;

  /* Go through each option looking for --version, or -version */
  for (narg = 1; narg < argc; narg++) 
    {
      option = argv[narg];
      
      if (!strncmp(option,"--version",9) ||
	  !strncmp(option,"-version",8))
	{

	  /* Print out the entire command line. */
	  for (nnarg = 0; nnarg < argc; nnarg++)
	    fprintf (stdout, "%s ", argv[nnarg]);
	  fprintf (stdout, "\n");
	  
	  fprintf (stdout, "%s Platform: %s C lib: %d\n",
		   id_string, __PLATFORM__, __GNUC_VERSION__);
	  
	  num_processed_args++;

	  /* Copy later args one step back. */
	  for (nnarg = narg; nnarg < argc - num_processed_args; nnarg++)
	    {
	      strcpy (argv[nnarg], argv[nnarg+1] );
	    }
	}
    }
  
  /* Return the number of arguments processed. */
  return num_processed_args;
}
