#define GE_MAGIC  (0x494d4746)  /* GE magic number
                                   (first four bytes, = "IMGF") */

/* ge compression codes */
#define GE_COMPRESSION_ASIS                  0
#define GE_COMPRESSION_RECTANGULAR           1
#define GE_COMPRESSION_PACKED                2
#define GE_COMPRESSION_COMPRESSED            3
#define GE_COMPRESSION_COMPRESSED_AND_PACKED 4

int mri_identify(char *fname);

int is_cor(char *fname);
int is_genesis(char *fname);
int is_ge_lx(char *fname);
int is_mgh(char *fname);
int is_mnc(char *fname);
int is_analyze(char *fname);
int is_siemens(char *fname);
int is_brik(char *fname);
int is_bshort(char *fname);

/* EOF */
