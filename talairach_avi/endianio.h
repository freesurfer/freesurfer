/**************************************************************************************/
/* Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006                           */
/* Washington University, Mallinckrodt Institute of Radiology.                        */
/* All Rights Reserved.                                                               */
/* This software may not be reproduced, copied, or distributed without written        */
/* permission of Washington University. For further information contact A. Z. Snyder. */
/**************************************************************************************/
/*$Header: /space/repo/1/dev/dev/talairach_avi/endianio.h,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: endianio.h,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.1  2006/09/23  05:29:38  avi
 * Initial revision
 **/

#include <ANALYZE.h>

extern void 	swab2 			(char *a);
extern void 	swab4 			(char *a);
extern void 	swab_hdr 		(struct dsr *phdr);
extern int 	CPU_is_bigendian 	(void);
extern void 	errf 			(char* program);
extern void 	errm 			(char* program);
extern void 	errr 			(char* program, char* filespc);
extern void 	errw 			(char* program, char* filespc);
extern void 	getroot 		(char *filespc, char *imgroot);
extern int 	eread 			(float *imgt, int n, int isbig, FILE *fp);
extern int 	ewrite 			(float *imgt, int n, char control, FILE *fp);
extern int 	gread 			(char *imgt, size_t bytes, int n, FILE *fp, int isbig);
extern int 	gwrite 			(char *imgt, size_t bytes, int n, FILE *fp, char control);
extern void 	load_4dfp_frame 	(char *fileroot, int *imgdim, int frame, int isbig, float *fimg);
extern int 	get_4dfp_dimoe 		(char *fileroot, int *imgdim, float *voxsiz, int *orient, int *isbig);
extern int 	get_4dfp_dimoe_quiet 	(char *fileroot, int *imgdim, float *voxsiz, int *orient, int *isbig);

/***************************************************************************************************************
swab2	 		byte swaps.
swab4 			byte swaps 4 bytes.
swab_hdr 		byte swaps the Analyze header (see ANALYZE.h for struct dsr).
CPU_is_bigendian 	returns 1 if the CPU on which it is executed is BIG ENDIAN architecture.
errf 			prints illegal byte/word message to stderr.
errm 			prints memory allocation error message to stderr.
errr 			prints read error message to stderr.
errw 			prints write error message to stderr.
getroot 		peels off the extensions from filespc and populates imgroot.
eread 			reads n*sizeof(float) bytes from fp and populates imgt. eread also performs swaps
		        if necessary depending on isbig and CPU_is_bigendian.
ewrite 			writes n*sizeof(float) bytes to fp from imgt. ewrite also performs swaps
			if necessary depending on control and CPU_is_bigendian.
gread 			reads n*bytes from fp and populates imgt. gread also performs swaps
			if necessary depending on isbig and CPU_is_bigendian.
gwrite 			writes n*bytes to fp from imgt. gwrite also performs swaps
			if necessary depending on control and CPU_is_bigendian.
load_4dfp_frame 	loads the volume identified by frame into fimg. It performs the necessary swaps.
get_4dfp_dimoe 		reads the IFH file and returns the image dimensions into imgdim,voxel size into voxsiz,
			orientation into orient and returns the isbig status.
get_4dfp_dimoe_quiet 	is the non-verbose version of get_4dfp_dimoe
***************************************************************************************************************/
