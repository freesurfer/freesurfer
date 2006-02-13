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

int string_to_type(char *string);
char *type_to_string(int type);

int mri_identify(char *fname);

int is_cor(char *fname);
int is_genesis(char *fname);
int is_ge_lx(char *fname);
int is_mgh(char *fname);
int is_mnc(char *fname);
int is_analyze(char *fname);
int is_siemens(char *fname);
int is_brik(char *fname);
int is_bhdr(char *fname);
int is_bshort(char *fname);
int is_bfloat(char *fname);
int is_sdt(char *fname);
int is_gdf(char *fname);
int is_otl(char *fname);
int is_ximg(char *fname);
int is_nifti1(char *fname);
int is_nii(char *fname);
int is_nrrd(char *fname);
int IDisCurv(char *curvfile);
char * bhdr_stem(char *fname);
char * bhdr_precisionstring(char *fname);
int bhdr_precision(char *fname);
char * bhdr_firstslicefname(char *fname);

/* EOF */
#endif
