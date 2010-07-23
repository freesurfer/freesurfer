/**
 * @file  pdf.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: ayendiki $
 *    $Date: 2010/07/23 21:07:45 $
 *    $Revision: 1.6 $
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



#ifndef PDF_INC
#define PDF_INC

#if defined(__cplusplus)
extern "C" {
#endif

unsigned long PDFtodSeed(void);
double PDFgaussian(void);
double PDFerlang(int order);
int PDFloadCDF(char *fname, double **xcdf, double **cdf, int *ncdf);
double PDFsampleCDF(double *xcdf, double *cdf, int ncdf);
int PDFsearchOrderedTable(double u, double *y, int ny);


#if defined(__cplusplus)
};
#endif

#endif //#ifndef PDF_INC
