/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

// cmdargs.c - utilities for handling command-line arguments

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cmdargs.h"

/* --------------------------------------------------------------------------
   CMDargNErr() - prints a msg as to the number of arguments needed by a flag
   and then exits.
   ------------------------------------------------------------------------ */
void CMDargNErr(char *option, int n)
{
  if (n == 1)
    fprintf(stderr, "ERROR: %s flag needs %d argument\n", option, n);
  else
    fprintf(stderr, "ERROR: %s flag needs %d arguments\n", option, n);
  exit(-1);
}
/*---------------------------------------------------------------
  CMDisFlag() - returns 1 if the arg string starts with a double-dash.
  The double-dash, ie, "--" indicates that the arg is a flag.
  ---------------------------------------------------------------*/
int CMDisFlag(char *arg)
{
  int len;
  len = strlen(arg);
  if (len < 2) return (0);

  if (arg[0] == '-' && arg[1] == '-') return (1);
  return (0);
}
/*-------------------------------------------------------------------------
  CMDsingleDash() - returns 1 if the arg string starts with a single-dash
  "-".  This often happens when the user accidentally specifies a flag
  with a single dash instead of a  double-dash.
  -----------------------------------------------------------------------*/
int CMDsingleDash(char *arg)
{
  int len;
  len = strlen(arg);
  if (len < 2) return (0);

  if (arg[0] == '-' && arg[1] != '-') return (1);
  return (0);
}
/*---------------------------------------------------------------
  CMDnthIsArg() - Checks that nth arg exists and is NOT a flag (ie,
  it does NOT start with a double-dash). This is good for reading in
  a variable number of arguments for a flag. nth is 0-based.
  ---------------------------------------------------------------*/
int CMDnthIsArg(int nargc, char **argv, int nth)
{
  /* check that there are enough args for nth to exist */
  if (nargc <= nth) return (0);
  /* check whether the nth arg is a flag */
  if (CMDisFlag(argv[nth])) return (0);
  return (1);
}
/*------------------------------------------------------------
  CMDstringMatch() - returns 1 if the two strings match
  ------------------------------------------------------------*/
int CMDstringMatch(const char *str1, const char *str2)
{
  if (str1 == NULL && str2 != NULL) return (0);
  if (str2 == NULL && str1 != NULL) return (0);
  if (!strcmp(str1, str2)) return (1);
  return (0);
}

/*------------------------------------------------------------
  CMDprintUsage() - eventually will find a text file with
  usage info and print it out, but currently does nothing.
  ------------------------------------------------------------*/
int CMDprintUsage(FILE *fp, char *ProgName) { return (0); }
/*------------------------------------------------------------
  CMDusageExit() -  same as CMDprintUsage() but exits.
  ------------------------------------------------------------*/
int CMDusageExit(char *ProgName)
{
  CMDprintUsage(stdout, ProgName);
  exit(1);
}
/*------------------------------------------------------------
  CMDprintHelp() - eventually will find a text file with
  help info and print it out, but currently does nothing.
  ------------------------------------------------------------*/
int CMDprintHelp(FILE *fp, char *ProgName) { return (0); }
