#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "version.h"


// char version[] = "version 1.0, Wed Apr 26 11:11:26 EDT 2000" ;


/* This function looks for the -v, --version, or -version tag in the
   argv and if found, prints out version information. This can be used
   in any binary. It will return the number of argc/argv options
   processed, so the caller can shorten their own argc/argv variables,
   like this:

    nargs = handle_version_option (argc, argv, "$Id: version.c,v 1.2 2003/03/19 18:04:26 kteich Exp $");
    argc -= nargs ;
    argv += nargs ;

   It will print out the id string passed in, which should just be
   "$Id", which CVS will expand to include the CVS information
   including CVS file, revision number, date, peson who checked it in,
   and tag, as well as the OS and GCC information.

   The binary may also want to exit if there are no other options to
   handle, i.e.

    if (1 == argc)
      exit (0);
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

  int narg = 0;
  int nnarg = 0;
  int num_processed_args = 0;
  char *option = NULL;
  char *os = NULL;

  /* Go through each option looking for -v, --version, or -version */
  for (narg = 1; narg < argc; narg++) 
    {
      option = argv[narg];
      
      if (!strncmp(option,"--version",9) ||
	  !strncmp(option,"-version",8) ||
	  !strncmp(option,"-v",2)) 
	{

	  /* Print out the entire command line. */
	  for (nnarg = 0; nnarg < argc; nnarg++)
	    fprintf (stdout, "%s ", argv[nnarg]);
	  fprintf (stdout, "\n");
	  
	  /* See if we can get the OS string. Print out the version
	     information with or without it. */
	  os = getenv("OS");
	  if (NULL != os)
	    {
	      fprintf (stdout, "%s Platform: %s GCC: %d\n",
		       id_string, os, __GNUC_VERSION__);
	    }
	  else
	    {
	      fprintf (stdout, "%s Platform: UNKNOWN GCC: %d\n",
		       id_string, __GNUC_VERSION__);
	    }
	  
	  num_processed_args++;
	}
    }
  
  /* Return the number of arguments processed. */
  return num_processed_args;
}
