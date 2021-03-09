/**
 * @brief Identify MRI volume format based on filename extension
 */
/*
 * Original Author: Christian Haselgrove
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


#ifndef MRI_IDENTIFY_H
#define MRI_IDENTIFY_H

#define GE_MAGIC  (0x494d4746)  /* GE magic number
(first four bytes, = "IMGF") */

#define NIFTI1_MAGIC  "ni1\0"
#define NII_MAGIC  "n+1\0"
#define NRRD_MAGIC  "NRRD"

/* ge compression codes */
#define GE_COMPRESSION_ASIS                  0
#define GE_COMPRESSION_RECTANGULAR           1
#define GE_COMPRESSION_PACKED                2
#define GE_COMPRESSION_COMPRESSED            3
#define GE_COMPRESSION_COMPRESSED_AND_PACKED 4

int string_to_type(const char *string);
char *type_to_string(int type);

int mri_identify(const char *fname);

int   IDtypeFromStem(const char *stem);
char *IDnameFromStem(const char *stem);
char *IDstemFromName(const char *name);
char *IDextensionFromName(const char *name);

int is_cor(const char *fname);
int is_genesis(const char *fname);
int is_ge_lx(const char *fname);
int is_mgh(const char *fname);
int is_mnc(const char *fname);
int is_analyze(const char *fname);
int is_siemens(const char *fname);
int is_brik(const char *fname);
int is_bhdr(const char *fname);
int is_bshort(const char *fname);
int is_bfloat(const char *fname);
int is_sdt(const char *fname);
int is_gdf(const char *fname);
int is_otl(const char *fname);
int is_ximg(const char *fname);
int is_nifti1(const char *fname);
int is_nii(const char *fname);
int is_nrrd(const char *fname);
int IDisCurv(const char *curvfile);
char * bhdr_stem(const char *fname);
char * bhdr_precisionstring(const char *fname);
int bhdr_precision(const char *fname);
char * bhdr_firstslicefname(const char *fname);

/* EOF */
#endif
