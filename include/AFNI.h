//
// AFNI.h
//
// define afni header structure

typedef struct 
{
  // mandatory attributes
  int dataset_rank[2];
  int dataset_dimensions[3];
  char typestring[16];
  int scene_data[3];
  int orient_specific[3];
  float origin[3];
  float delta[3];
  // almost mandatory attributes
  char *idcode_string; int numchars;
  char byteorder_string[10];
  float *brick_stats; int numstats;
  int *brick_types; int numtypes;
  float *brick_float_facs; int numfacs;
} AFNI_HEADER, AF;

MRI *afniRead(char *fname, int read_volume);
int afniWrite(MRI *mri, char *fname);
int readAFNIHeader(FILE *fp, AF *paf);
void AFinit(AF *pAF);
void AFclean(AF *pAF);
void printAFNIHeader(AF *pAF);
