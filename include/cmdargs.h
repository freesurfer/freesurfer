// cmdargs.h - include file for cmdargs.c (used for 
//   handling command-line arguments)
// $Id: cmdargs.h,v 1.2 2005/09/14 18:20:51 greve Exp $


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



