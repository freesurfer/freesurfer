/**
 * @file  hips_error.h
 * @brief definitions related to the HIPS error handler
 *
 */
/*
 * Original Author: Michael Landy
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/11/18 03:03:32 $
 *    $Revision: 1.3 $
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef HIPS_ERROR_H
#define HIPS_ERROR_H

/*
 * Copyright (c) 1991 Michael Landy
 *
 * Disclaimer:  No guarantees of performance accompany this software,
 * nor is any responsibility assumed on the part of the authors.  All the
 * software has been tested extensively and every effort has been made to
 * insure its reliability.
 */

/*
 * hips_error.h - definitions related to the HIPS error handler
 *
 * Michael Landy - 12/28/90
 */

/* error-related variables */

extern int hipserrno;
extern int hipserrlev; /* if errlev <= errno - then print&die */
extern int hipserrprt; /* if errprt <= errno < errlev - then print&return */
extern char hipserr[];

struct h_errstruct
{  /* the error structure defined in herrs.c */
  char *h_errstr;  /* sprintf string for error code */
  int h_errfmt;  /* list of variables for sprintf */
  int h_errsev;  /* error severity */
};

#define HIPS_ERROR -1 /* error-return from hips routines */
#define HIPS_OK  0 /* one possible nonerror-return value */

/* Error code severities */
#define HEL_INFORM 1 /* informational */
#define HEL_WARN 2 /* warning */
#define HEL_ERROR 3 /* error */
#define HEL_SEVERE 4 /* severe error */

/*
 * Standard error codes, the strings are defined in lib/herrs.c
 *
 * Note that these codes must be allocated sequentially, starting from 1.
 * perr.c depends on this fact.
 */

#define HE_ALLOC 1 /* "can't allocate memory" */
#define HE_ALLOCSUBR 2 /* "%s: can't allocate memory" */
#define HE_FREE 3 /* "can't free memory" */
#define HE_READFR 4 /* "error reading frame %d" */
#define HE_READFRFILE 5 /* "error reading frame %d from file %s" */
#define HE_READFILE 6 /* "error reading from file %s" */
#define HE_READ 7 /* "error during read" */
#define HE_WRITEFR 8 /* "error writing frame %d" */
#define HE_WRITEFRFILE 9 /* "error writing frame %d to file %s" */
#define HE_OPEN 10 /* "can't open file: `%s'" */
#define HE_SEEK 11 /* "can't perform seek" */
#define HE_FMTMATCH 12 /* "pixel format mismatch" */
#define HE_FMTMATCHFILE 13 /* "pixel format mismatch, file %s" */
#define HE_FMT 14 /* "can't handle pixel format %s" */
#define HE_FMTFILE 15 /* "can't handle pixel format %s, file: %s" */
#define HE_POW2 16 /* "image dimensions must be powers of 2" */
#define HE_SNEG 17 /* "number of frames must be zero or positive" */
#define HE_SMATCH 18 /* "size mismatch" */
#define HE_FRMATCH 19 /* "mismatch of number of frames" */
#define HE_LARGE 20 /* "frame dimensions are too large" */
#define HE_ZNEG 21 /* "zero or negative dimension" */
#define HE_BUF 22 /* "%s: buffer limit exceeded" */
#define HE_CODE 23 /* "%s: unknown op-code" */
#define HE_CUT1 24 /* "cut_frame: 1 intersection" */
#define HE_CUTI 25 /* "cut_frame: intersections not between?" */
#define HE_FFTI 26 /* "%s: strange index=%d" */
#define HE_PYRTLZ 27 /* "pyrnumpix: toplev less than zero?" */
#define HE_PYRTL 28 /* "pyrnumpix: toplev too large?" */
#define HE_REFL 29 /* "%s: Invalid reflection type %d" */
#define HE_FRMEND 30 /* "read_frame: did not find frame-end" */
#define HE_HDRREAD 31 /* "error while reading header, file %s" */
#define HE_HDRPREAD 32 /* "error reading extended parameters, file %s" */
#define HE_HDRBREAD 33 /* "error reading binary parameters, file %s" */
#define HE_HDRPTYPE 34 /* "invalid extended parameter format %d, file %s" */
#define HE_HDRXOV 35
/* "header parameters overflow, file: %s, size should be %d and was %d" */
#define HE_HDRWRT 36 /* "error while writing header, file %s" */
#define HE_HDRPWRT 37 /* "error while writing header parameters, file %s" */
#define HE_HDRBWRT 38 /* "error while writing header binary area, file %s" */
#define HE_HDWPTYPE 39
/* "invalid extended parameter format %d during write, file %s" */
#define HE_REQ 40 /* "%s: reqested %d, got %d" */
#define HE_BADFMT 41 /* "invalid format code %d" */
#define HE_HDPTYPE 42 /* "%s: invalid extended parameter format %d" */
#define HE_MISSPAR 43 /* "%s: can't find extended parameter %s" */
#define HE_C_ROW 44 /* "mismatched number of rows, file: %s" */
#define HE_C_COL 45 /* "mismatched number of columns, file: %s" */
#define HE_C_FRM 46 /* "mismatched number of frames, file: %s" */
#define HE_C_FMT 47 /* "mismatched pixel format (%s), file: %s" */
#define HE_C_NCL 48 /* "mismatched number of colors, file: %s" */
#define HE_C_NLV 49 /* "mismatched number of pyramid levels, file: %s" */
#define HE_FMTSUBR 50 /* "%s: can't handle pixel format %s" */
#define HE_CTORTP 51 /* "%s: unknown complex-to-real conversion: %d" */
#define HE_METH 52 /* "%s: unknown method (%d), file: %s" */
#define HE_FMTSUBRFILE 53 /* "%s: can't handle pixel format %s, file: %s" */
#define HE_SETFP 54 /* "setformat: can't handle pyramid formats" */
#define HE_SETPF 55 /* "setpyrformat: can only handle pyramid formats" */
#define HE_RTOCTP 56 /* "%s: unknown real-to-complex conversion: %d" */
#define HE_MSG 57 /* "%s" */
#define HE_ROI 58
/* "setroi: ROI out of bounds, first=(%d,%d), size=(%d,%d)" */
#define HE_CONV 59 /* "converting from %s to %s, file: %s" */
#define HE_CONVI 60 /* "converting from %s to %s via integer, file: %s" */
#define HE_SMALL 61 /* "frame dimensions are too small" */
#define HE_RNG 62 /* "%s: invalid range" */
#define HE_XINC 63
/* "header parameters inconsistency for (%s %d %d %d) offset is %d, file: %s" */
#define HE_ROI8 64 /* "%s: packed image ROI - columns not multiple of 8" */
#define HE_ROI8F 65
/* "%s: packed image ROI - columns not multiple of 8, file: %s" */
#define HE_ROI8C 66
/* "%s: packed image ROI - columns not multiple of 8, clearing ROI" */
#define HE_IMSG 67 /* "%s" */
#define HE_UNKFLG 68 /* "unrecognised flag option %s\n%s" */
#define HE_MUTEX 69 /* "flags are mutually exclusive\n%s" */
#define HE_MISSFPAR 70 /* "missing parameter for flag %s\n%s" */
#define HE_INVFPAR 71 /* "invalid parameter %s for flag %s\n%s" */
#define HE_SYNTAX 72 /* "invalid syntax\n%s" */
#define HE_FILECNT 73 /* "too many image filenames\n%s" */
#define HE_MISSFILE 74 /* "missing image filenames\n%s" */
#define HE_INVFILE 75 /* "invalid image filename %s" */
#define HE_STDIN 76 /* "can't open the standard input twice" */
#define HE_FMT3SUBR 77
/* "%s: can't handle pixel format combination %s/%s/%s" */
#define HE_FMT2SUBR 78
/* "%s: can't handle pixel format combination %s/%s" */
#define HE_RCSELN 79 /* "%s: row/column selection out of range" */
#define HE_C_OROW 80 /* "mismatched total number of rows, file: %s" */
#define HE_C_OCOL 81 /* "mismatched total number of columns, file: %s" */
#define HE_C_FRMC 82
/* "mismatched number of frames and numcolor>1 or numdepth > 1, file: %s" */
#define HE_MULT2 83 /* "image dimensions must be multiples of 2" */
#define HE_MSKFUNFILE 84/* "bad mask function number %d, file: %s" */
#define HE_MSKFUNSUBR 85/* "%s: bad mask function number %d" */
#define HE_MSKCNT 86 /* "%s: mask function %d, bad mask count %d" */
#define HE_SIGMA 87 /* "%s: sigma less than or equal to zero" */
#define HE_PTWICE 88 /* "flag `-%s' specified more than once" */
#define HE_WINDSZ 89 /* "%s: invalid window size (%d)" */
#define HE_COMPOV 90 /* "component table overflow" */
#define HE_PNTOV 91 /* "point table overflow" */
#define HE_FILTPAR 92 /* "bad filter parameter" */
#define HE_HDRPTYPES 93 /* "invalid extended parameter format %s, file %s" */
#define HE_COLSPEC 94 /* "%s: invalid color specification: `%s'" */
#define HE_PNAME 95 /* "%s: parameter name contains white space: `%s'" */
#define HE_PCOUNT 96 /* "%s: supplied count doesn't match that of `%s'" */
#define HE_COLOVF 97 /* "%s: colormap table overflow, file: %s */
#define HE_FMTSLEN 98 /* "length of formats array mismatch, file: %s" */
#define HE_COL1 99 /* "input 3-color image has numcolor>1, file %s" */
#define HE_COL3 100
/* "can't convert 1-color image with numcolor!=3 to 3-color, file %s" */
#define HE_C_DEPTH 101 /* "mismatched number of depths, file: %s" */
#define HE_C_DPTHC 102
/* "mismatched number of depths and numcolor>1 or numdepth > 1, file: %s" */
#define HE_BADDEPTH 103 /* "bad value of depth, should be > 0" */
#define MAXERR 103

/* error print formats */

#define HEP_N 0  /* no arguments */
#define HEP_D 1  /* %d */
#define HEP_S 2  /* %s */
#define HEP_SD 3 /* %s %d */
#define HEP_DS 4 /* %d %s */
#define HEP_SS 5 /* %s %s */
#define HEP_SDD 6 /* %s %d %d */
#define HEP_SDS 7 /* %s %d %s */
#define HEP_SSS 8 /* %s %s %s */
#define HEP_DDDD 9 /* %d %d %d %d */
#define HEP_SDDDDS 10 /* %s %d %d %d %d %s */
#define HEP_SSSS 11 /* %s %s %s %s */

#endif
