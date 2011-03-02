/**
 * @file  dummy.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:01 $
 *    $Revision: 1.10 $
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



/*!
\file dummy.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/


// $Id: dummy.c,v 1.10 2011/03/02 00:04:01 nicks Exp $

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: dummy.c,v 1.10 2011/03/02 00:04:01 nicks Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *TempVolFile=NULL;
char *subject, *hemi, *SUBJECTS_DIR;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  return 0;
}
/* ------ Doxygen markup starts on the line below ---- */
/*!
\fn int parse_commandline(int argc, char **argv)
\brief Parses the command-line arguments
\param argc - number of command line arguments
\param argv - pointer to a character pointer
*/
/* ------ Doxygen markup ends on the line above ---- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--temp-vol")) {
      if (nargc < 1) CMDargNErr(option,1);
      TempVolFile = pargv[0];
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void usage_exit(void)
\brief Prints usage and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_usage(void)
\brief Prints usage and returns (does not exit)
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --temp-vol volfile : template volume \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_help(void)
\brief Prints help and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_version(void)
\brief Prints version and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void check_options(void)
\brief Checks command-line options
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void check_options(void) {
  return;
}

/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void dump_options(FILE *fp)
\brief Prints command-line options to the given file pointer
\param FILE *fp - file pointer
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  return;
}
