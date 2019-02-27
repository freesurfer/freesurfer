/**
 * @file  cmdargs.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.5 $
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


// cmdargs.h - include file for cmdargs.c (used for
//   handling command-line arguments)
// $Id: cmdargs.h,v 1.5 2011/03/02 00:04:09 nicks Exp $


#ifndef CMDARGS_H
#define CMDARGS_H

const char * CMDSrcVersion(void);
void CMDargNErr(char *option, int n);
int CMDsingleDash(char *flag);
int CMDisFlag(char *flag);
int CMDnthIsArg(int nargc, char **argv, int nth);
int CMDstringMatch(char *str1, char *str2);

int CMDprintUsage(FILE *fp, char *ProgName);
int CMDusageExit(char *ProgName);
int CMDprintHelp(FILE *fp, char *ProgName);

#endif
