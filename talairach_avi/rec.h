/*$Header: /space/repo/1/dev/dev/talairach_avi/rec.h,v 1.1 2007/05/04 22:34:03 nicks Exp $*/
/*$Log: rec.h,v $
/*Revision 1.1  2007/05/04 22:34:03  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.5  2006/09/23  20:53:17  avi
 * startrecle()  startrece () final argument now control
 *
 * Revision 1.4  2006/09/23  06:29:44  avi
 * extern int startrece()  extern int startrecle()
 *
 * Revision 1.3  2004/03/03  02:33:23  avi
 **/
/**************************************************************************************/
/* Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006                           */
/* Washington University, Mallinckrodt Institute of Radiology.                        */
/* All Rights Reserved.                                                               */
/* This software may not be reproduced, copied, or distributed without written        */
/* permission of Washington University. For further information contact A. Z. Snyder. */
/**************************************************************************************/

extern int printrec	(char *string);
extern int catrec	(char *file);
extern int startrec	(char *outfile, int argc, char *argv[], char *rcsid);
extern int startrecl	(char *outfile, int argc, char *argv[], char *rcsid);
extern int startrece	(char *outfile, int argc, char *argv[], char *rcsid, char control);
extern int startrecle	(char *outfile, int argc, char *argv[], char *rcsid, char control);
extern int endrec	(void);
