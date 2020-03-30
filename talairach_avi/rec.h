/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
 *
 * Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
 * Washington University, Mallinckrodt Institute of Radiology.
 * All Rights Reserved.
 *
 * This software may not be reproduced, copied, or distributed without 
 * written permission of Washington University. For further information 
 * contact A. Z. Snyder.
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

extern int printrec	(char *string);
extern int catrec	(char *file);
extern int startrec	(char *outfile, int argc, char *argv[], char *rcsid);
extern int startrecl	(char *outfile, int argc, char *argv[], char *rcsid);
extern int startrece	(char *outfile, int argc, char *argv[], char *rcsid, char control);
extern int startrecle	(char *outfile, int argc, char *argv[], char *rcsid, char control);
extern int endrec	(void);

const char* current_date_time();
