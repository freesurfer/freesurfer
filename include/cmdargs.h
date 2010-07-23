/**
 * @file  cmdargs.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2010/07/23 21:07:44 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


// cmdargs.h - include file for cmdargs.c (used for
//   handling command-line arguments)
// $Id: cmdargs.h,v 1.4 2010/07/23 21:07:44 ayendiki Exp $


#ifndef CMDARGS_H
#define CMDARGS_H

#if defined(__cplusplus)
extern "C" {
#endif

const char * CMDSrcVersion(void);
void CMDargNErr(char *option, int n);
int CMDsingleDash(char *flag);
int CMDisFlag(char *flag);
int CMDnthIsArg(int nargc, char **argv, int nth);
int CMDstringMatch(char *str1, char *str2);

int CMDprintUsage(FILE *fp, char *ProgName);
int CMDusageExit(char *ProgName);
int CMDprintHelp(FILE *fp, char *ProgName);

#if defined(__cplusplus)
};
#endif

#endif



