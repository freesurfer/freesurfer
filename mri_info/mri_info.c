////////////////////////////////////////////////////////////////////
// mri_info.c
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: tosa $
// Revision Date  : $Date: 2004/08/25 19:05:20 $
// Revision       : $Revision: 1.29 $
//
////////////////////////////////////////////////////////////////////
char *MRI_INFO_VERSION = "$Revision: 1.29 $";
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "const.h"
#include "machine.h"
#include "fio.h"
#include "utils.h"
#include "mri.h"
#include "volume_io.h"
#include "analyze.h"
#include "mri_identify.h"
#include "error.h"
#include "diag.h"
#include "version.h"
#include "mghendian.h"
#include "fio.h"

struct ge_header {
  int magic;
  int width;
  int height;
  int depth;
  int exam_header_offset;
  int series_header_offset;
  int image_header_offset;
  int compression;
};

static void do_file(char *fname);
#if 0
static void read_ge_5x_file(char *fname, struct stat stat_buf);
static void read_ge_8x_file(char *fname, struct stat stat_buf);
static void read_siemens_file(char *fname, struct stat stat_buf);
static void read_minc_file(char *fname, struct stat stat_buf);
static void read_mgh_file(char *fname, struct stat stat_buf);
static void read_analyze_file(char *fname, struct stat stat_buf);
static void read_sdt_file(char *fname, struct stat stat_buf);
static void read_cor(char *fname);
static void read_brik_file(char *fname, struct stat stat_buf);
static void read_short(FILE *fp, int offset, short *val);
static void read_int(FILE *fp, int offset, int *val);
static void read_string(FILE *fp, int offset, char *val, int length);
static void read_float(FILE *fp, int offset, float *val);
static void read_double(FILE *fp, int offset, double *val);
static void read_analyze_header(char *fname, dsr *bufptr);
static void read_ge_header(FILE *fp, int offset, struct ge_header *header);
static void flip_analyze_header(dsr *header);


static char *month[] = {
  "Month Zero",
  "January", "February", "March", "April", "May", "June", 
  "July", "August", "September", "October", "November", "December"
};

static char *type_text[] = { "unsigned char", "int", "long", "float", "short", "bitmap" };
#endif

char *Progname ;

static void usage(char *prog_name, int exit_val)
{

  fprintf(stderr, "usage: %s filename\n", prog_name);
  exit(exit_val);

} /* end usage() */



int main(int argc, char *argv[])
{
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_info.c,v 1.29 2004/08/25 19:05:20 tosa Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;


  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc < 2)
    {
    usage(argv[0], 1);
    }

  for(argv++;*argv;argv++)
    do_file(*argv);

  exit(0);

} /* end main() */

int PrettyMatrixPrint(MATRIX *mat)
{
  int row;

  if (mat == NULL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat = NULL!")) ;

  if (mat->type != MATRIX_REAL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not Real type")) ;
 
  if (mat->rows != 4 || mat->cols != 4)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "mat is not of 4 x 4")) ;
    
  for (row=1; row < 5; ++row)
    printf("              %8.4f %8.4f %8.4f %10.4f\n",
	   mat->rptr[row][1], mat->rptr[row][2], mat->rptr[row][3], mat->rptr[row][4]);
  return (NO_ERROR);
}

static void do_file(char *fname)
{
#if 0
  struct stat stat_buf;
  int type;
  char *at;
  char fname2[STRLEN];
  int i;
#endif
  MRI *mri ;
  MATRIX *m ;

  mri = MRIreadHeader(fname, MRI_VOLUME_TYPE_UNKNOWN) ;
  if (!mri)
    return;

  printf("Volume information for %s\n", fname);
  // mri_identify has been called but the result is not stored and thus I have to call it again
  printf("          type: %s\n", type_to_string(mri_identify(fname)));
  if (mri->nframes > 1)
    printf("    dimensions: %d x %d x %d x %d\n", mri->width, mri->height, mri->depth, mri->nframes) ;
  else
    printf("    dimensions: %d x %d x %d\n", mri->width, mri->height, mri->depth) ;
  printf("   voxel sizes: %6.4f, %6.4f, %6.4f\n", mri->xsize, mri->ysize, mri->zsize) ;
  printf("          type: %s (%d)\n",
	 mri->type == MRI_UCHAR   ? "UCHAR" :
	 mri->type == MRI_SHORT   ? "SHORT" :
	 mri->type == MRI_INT     ? "INT" :
	 mri->type == MRI_LONG    ? "LONG" :
	 mri->type == MRI_BITMAP  ? "BITMAP" :
	 mri->type == MRI_TENSOR  ? "TENSOR" :
	 mri->type == MRI_FLOAT   ? "FLOAT" : "UNKNOWN", mri->type) ;
  printf("           fov: %2.3f\n", mri->fov) ;
  printf("        xstart: %2.1f, xend: %2.1f\n", mri->xstart*mri->xsize, mri->xend*mri->xsize) ;
  printf("        ystart: %2.1f, yend: %2.1f\n", mri->ystart*mri->ysize, mri->yend*mri->ysize) ;
  printf("        zstart: %2.1f, zend: %2.1f\n", mri->zstart*mri->zsize, mri->zend*mri->zsize) ;
  printf("            TR: %2.2f msec, TE: %2.2f msec, TI: %2.2f msec, flip angle: %2.2f degrees\n",
	 mri->tr, mri->te, mri->ti, DEGREES(mri->flip_angle)) ;
  printf("       nframes: %d\n", mri->nframes) ;
  printf("ras xform %spresent\n", mri->ras_good_flag ? "" : "not ") ;
  printf("    xform info: x_r = %8.4f, y_r = %8.4f, z_r = %8.4f, c_r = %10.4f\n",
	 mri->x_r, mri->y_r, mri->z_r, mri->c_r);
  printf("              : x_a = %8.4f, y_a = %8.4f, z_a = %8.4f, c_a = %10.4f\n",
	 mri->x_a, mri->y_a, mri->z_a, mri->c_a);
  printf("              : x_s = %8.4f, y_s = %8.4f, z_s = %8.4f, c_s = %10.4f\n",
	 mri->x_s, mri->y_s, mri->z_s, mri->c_s);

  if (fio_IsDirectory(fname))
    printf("\ntalairach xfm : %s\n", mri->transform_fname);
  else
  {
    char *ext = 0;
    ext = fio_extension(fname);
    if (ext)
    {
      if (strcmp(ext, "mgz") == 0 || strcmp(ext, "mgh")==0)
	printf("\ntalairach xfm : %s\n", mri->transform_fname);
      free(ext);
    }
  }
  m = MRIgetVoxelToRasXform(mri) ; // extract_i_to_r(mri) (just macto)
  printf("\nvoxel to ras transform:\n") ; PrettyMatrixPrint(m) ;
  MatrixFree(&m) ;
  m = extract_r_to_i(mri);
  printf("\nras to voxel transform:\n"); PrettyMatrixPrint(m);
  MatrixFree(&m);
  MRIfree(&mri);
  
  return;

#if 0 
  strcpy(fname2, fname) ;
  at = strrchr(fname2, '@');

  if(at)
  {
    *at = '\0';
    at++;
  }

  if (at)
  {
#if 1
    type  = string_to_type(at) ;
#else
    for(i = 0;at[i] != '\0';i++)
      at[i] = (at[i] >= 'a' && at[i] <= 'z' ? at[i] += 'A' - 'a' : at[i]);
    if (!strcmp(at, "MNC"))
      type = MRI_MINC_FILE ;
    else if (!strcmp(at, "MINC"))
      type = MRI_MINC_FILE ;
    else if (!strcmp(at, "BRIK"))
      type = BRIK_FILE ;
    else if(!strcmp(at, "SIEMENS"))
      type = SIEMENS_FILE ;
    else if (!strcmp(at, "MGH"))
      type = MRI_MGH_FILE ;
    else if (!strcmp(at, "MR"))
      type = GENESIS_FILE ;
    else if (!strcmp(at, "GE"))
      type = GE_LX_FILE ;
    else if (!strcmp(at, "IMG"))
      type = MRI_ANALYZE_FILE ;
    else if (!strcmp(at, "SDT"))
      type = SDT_FILE;
    else if(!strcmp(at, "COR"))
      type = MRI_CORONAL_SLICE_DIRECTORY ;
    else
    {
      fprintf(stderr, "unknown file type %s\n", at);
      return;
    }
#endif
  }
  else  /* no '@' found */
  {

    type = -1;

    if(stat(fname2, &stat_buf) < 0)
    {
      fprintf(stderr, "can't stat %s\n", fname2);
      return;
    }

    if(S_ISDIR(stat_buf.st_mode))
      type = MRI_CORONAL_SLICE_DIRECTORY;
    else
    {

      if(is_genesis(fname2))
	type = GENESIS_FILE;
      else if(is_ge_lx(fname2))
	type = GE_LX_FILE;
      else if(is_brik(fname2))
	type = BRIK_FILE;
      else if(is_siemens(fname2))
	type = SIEMENS_FILE;
      else if(is_sdt(fname2))
	type = SDT_FILE;
      else if(is_analyze(fname2))
	type = MRI_ANALYZE_FILE;
      else if(is_brik(fname2))
	type = BRIK_FILE;
      else if(is_mgh(fname2))
	type = MRI_MGH_FILE;
      else if(is_mnc(fname2))
	type = MRI_MINC_FILE;
    }

    if(type == -1)
    {
      fprintf(stderr, "unrecognized file type for file %s\n", fname2);
      return;
    }

  }

  switch(type)
  {
    case GENESIS_FILE:
      read_ge_5x_file(fname2, stat_buf);
      break;

    case GE_LX_FILE:
      read_ge_8x_file(fname2, stat_buf);
      break;

    case SIEMENS_FILE:
      read_siemens_file(fname2, stat_buf);
      break;

    case MRI_ANALYZE_FILE:
      read_analyze_file(fname2, stat_buf);
      break;

    case MRI_MGH_FILE:
      read_mgh_file(fname2, stat_buf);
      break;

    case MRI_MINC_FILE:
      read_minc_file(fname2, stat_buf);
      break;

    case BRIK_FILE:
      read_brik_file(fname2, stat_buf);
      break;

    case SDT_FILE:
      read_sdt_file(fname2, stat_buf);
      break;

    case MRI_CORONAL_SLICE_DIRECTORY:
      read_cor(fname2);
      break;

    default:
      break;
  }
#endif
} /* end do_file */

#if 0 
static void read_short(FILE *fp, int offset, short *val)
{

  short buf;

  fseek(fp, offset, 0);
  fread(&buf, 2, 1, fp);
  *val = orderShortBytes(buf);

} /* end read_short() */

static void read_int(FILE *fp, int offset, int *val)
{

  int buf;

  fseek(fp, offset, 0);
  fread(&buf, 4, 1, fp);
  *val = orderIntBytes(buf);

} /* end read_int() */

static void read_string(FILE *fp, int offset, char *val, int length)
{

  fseek(fp, offset, 0);
  fread(val, 1, length, fp);
  val[length] = '\0';

} /* end read_string() */

static void read_float(FILE *fp, int offset, float *val)
{

  float buf;

  fseek(fp, offset, 0);
  fread(&buf, 4, 1, fp);
  *val = orderFloatBytes(buf);

} /* end read_float() */

static void read_double(FILE *fp, int offset, double *val)
{

  double buf;

  fseek(fp, offset, 0);
  fread(&buf, 8, 1, fp);
  *val = orderDoubleBytes(buf);

} /* end read_double() */

static void read_brik_file(char *fname, struct stat stat_buf)
{

  FILE *fp;
  char type_line[STRLEN], name_line[STRLEN], type[STRLEN], name[STRLEN], count_line[STRLEN];
  int count;
  char buf[10000];
  int vals[50];
  int i;
  char fname2[STRLEN];
  char *ext;

  strcpy(fname2, fname);
  if((ext = strstr(fname2, "BRIK")) != NULL)
    sprintf(ext, "HEAD");

  if((fp = fopen(fname2, "r")) == NULL)
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: BRIK\n");

  while(!feof(fp))
  {
    type_line[0] = '\0';
    while(!feof(fp) && strncmp(type_line, "type", 4))
      fgets(type_line, STRLEN, fp);
    if(!feof(fp))
    {
      name_line[0] = '\0';
      fgets(name_line, STRLEN, fp);
      fgets(count_line, STRLEN, fp);
      sscanf(type_line, "type = %s", type);
      sscanf(name_line, "name = %s", name);
      sscanf(count_line, "count = %d", &count);

      if(strcmp(type, "integer-attribute") == 0)
        for(i = 0;i < count;i++)
          fscanf(fp, "%d", &vals[i]);
      if(strcmp(type, "string-attribute") == 0)
      {
        if(count < 10000)
          fread(buf, 1, count, fp);
        else
          count = 0;
        buf[count] = '\0';
        for(i = 0;i < count;i++)
          if(buf[i] == '~')
            buf[i] = '\0';
      }

      if(strcmp(name, "BYTEORDER_STRING") == 0)
      {
        if(strcmp(buf, "`MSB_FIRST") == 0)
          printf("byteorder: big-endian\n");
        if(strcmp(buf, "`LSB_FIRST") == 0)
          printf("byteorder: little-endian\n");
      }
      if(strcmp(name, "DATASET_DIMENSIONS") == 0)
      {
        printf("width: %d\n", vals[0]);
        printf("height: %d\n", vals[1]);
      }
      if(strcmp(name, "IDCODE_DATE") == 0)
        printf("idcode date: %s\n", &buf[1]);
      if(strcmp(name, "IDCODE_STRING") == 0)
        printf("idcode string: %s\n", &buf[1]);
      if(strcmp(name, "TYPESTRING") == 0)
        printf("typestring: %s\n", &buf[1]);
      if(strcmp(name, "DATASET_NAME") == 0)
        printf("dataset name: %s\n", &buf[1]);
      if(strcmp(name, "ANATOMY_PARENTNAME") == 0)
        printf("anatomy parentname: %s\n", &buf[1]);
      if(strncmp(name, "LABEL_", 6) == 0)
        printf("label %s: %s\n", &name[6], &buf[1]);

    }
  }

  fclose(fp);

}  /*  end read_brik_file()  */

static void read_ge_5x_file(char *fname, struct stat stat_buf)
{

  FILE *fp;
  struct ge_header header;
  char string[50];
  short s;
  float f;
  int i;

  if((fp = fopen(fname, "r")) == NULL)
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: GE Signa 5.x\n");

  read_ge_header(fp, 0, &header);

  read_string(fp, 170, string, 20);
  printf("exam number/series number/first image number: %s\n", string);

  read_string(fp, header.exam_header_offset + 0, string, 4);
  printf("suite ID: %s\n", string);
  read_string(fp, header.exam_header_offset + 97, string, 25);
  printf("patient name: %s\n", string);
  read_short(fp, header.exam_header_offset + 122, &s);
  printf("patient age: %hd\n", s);
  read_short(fp, header.exam_header_offset + 126, &s);
  printf("patient sex: %hd\n", s);
  read_string(fp, header.exam_header_offset + 305, string, 3);
  printf("exam type: %s\n", string);
  read_int(fp, header.exam_header_offset + 208, &i);
  strcpy(string, ctime((time_t *)&i));
  string[strlen(string)-1] = '\0';
  printf("exam date: %s\n", string);

  read_string(fp, header.exam_header_offset + 10, string, 33);
  printf("hospital name: %s\n", string);
  read_int(fp, header.exam_header_offset + 80, &i);
  printf("magnet strength in gauss: %d\n", i);
  read_string(fp, header.exam_header_offset + 86, string, 61);
  printf("patient history: %s\n", string);
  read_string(fp, header.exam_header_offset + 212, string, 33);
  printf("referring physician: %s\n", string);
  read_string(fp, header.exam_header_offset + 245, string, 33);
  printf("diagnostician/radiologist: %s\n", string);
  read_string(fp, header.exam_header_offset + 282, string, 23);
  printf("exam description: %s\n", string);
  read_string(fp, header.exam_header_offset + 318, string, 9);
  printf("creator suite and host: %s\n", string);

  read_short(fp, header.series_header_offset + 10, &s);
  printf("series number: %hd\n",s);
  read_string(fp, header.series_header_offset + 84, string, 3);
  printf("anatomical reference: %s\n", string);
  read_short(fp, header.series_header_offset + 72, &s);
  printf("most like plane: %hd\n",s);
  read_string(fp, header.series_header_offset + 20, string, 30);
  printf("series description: %s\n", string);
  read_short(fp, header.series_header_offset + 140, &s);
  printf("number of acquisitions: %hd\n", s);
  read_int(fp, header.series_header_offset + 16, &i);
  strcpy(string, ctime((time_t *)&i));
  string[strlen(string)-1] = '\0';
  printf("series date: %s\n", string);

  read_short(fp, header.image_header_offset + 12, &s);
  printf("image number: %hd\n", s);
  read_float(fp, header.image_header_offset + 26, &f);
  printf("slice thickness: %g mm\n", f);
  read_float(fp, header.image_header_offset + 34, &f);
  printf("display field of view x: %g mm\n", f);
  read_float(fp, header.image_header_offset + 38, &f);
  printf("display field of view y: %g mm\n", f);
  read_float(fp, header.image_header_offset + 42, &f);
  printf("image dimension x: %g\n", f);
  read_float(fp, header.image_header_offset + 46, &f);
  printf("image dimension y: %g\n", f);
  read_float(fp, header.image_header_offset + 50, &f);
  printf("pixel size x: %g mm\n", f);
  read_float(fp, header.image_header_offset + 54, &f);
  printf("pixel size y: %g mm\n", f);
  read_string(fp, header.image_header_offset + 72, string, 17);
  printf("iv contrast agent: %s\n",string);
  read_string(fp, header.image_header_offset + 89, string, 17);
  printf("oral contrast agent: %s\n",string);
  read_string(fp, header.image_header_offset + 308, string, 33);
  printf("pulse sequence name: %s\n",string);
  read_string(fp, header.image_header_offset + 362, string, 17);
  printf("coil name: %s\n",string);
  read_int(fp, header.image_header_offset + 18, &i);
  strcpy(string, ctime((time_t *)&i));
  string[strlen(string)-1] = '\0';
  printf("image date: %s\n", string);

  read_float(fp, header.image_header_offset + 130, &f);
  printf("center R: %g\n", f);
  read_float(fp, header.image_header_offset + 134, &f);
  printf("center A: %g\n", f);
  read_float(fp, header.image_header_offset + 138, &f);
  printf("center S: %g\n", f);
  read_float(fp, header.image_header_offset + 142, &f);
  printf("normal R: %g\n", f);
  read_float(fp, header.image_header_offset + 146, &f);
  printf("normal A: %g\n", f);
  read_float(fp, header.image_header_offset + 150, &f);
  printf("normal S: %g\n", f);
  read_float(fp, header.image_header_offset + 154, &f);
  printf("top left R: %g\n", f);
  read_float(fp, header.image_header_offset + 158, &f);
  printf("top left A: %g\n", f);
  read_float(fp, header.image_header_offset + 162, &f);
  printf("top left S: %g\n", f);
  read_float(fp, header.image_header_offset + 166, &f);
  printf("top right R: %g\n", f);
  read_float(fp, header.image_header_offset + 170, &f);
  printf("top right A: %g\n", f);
  read_float(fp, header.image_header_offset + 174, &f);
  printf("top right S: %g\n", f);
  read_float(fp, header.image_header_offset + 178, &f);
  printf("bottom right R: %g\n", f);
  read_float(fp, header.image_header_offset + 182, &f);
  printf("bottom right A: %g\n", f);
  read_float(fp, header.image_header_offset + 186, &f);
  printf("bottom right S: %g\n", f);

  read_short(fp, header.image_header_offset + 398, &s);
  printf("number of slices in scan group: %hd\n", s);

  fclose(fp);

} /* end read_ge_5x_file() */

static void read_ge_8x_file(char *fname, struct stat stat_buf)
{

  FILE *fp;
  struct ge_header header;
  short s;
  int i;
  float f;
  char string[100];

  if((fp = fopen(fname, "r")) == NULL)
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: GE Signa 8.x\n");

  read_ge_header(fp, 3228, &header);

  read_string(fp, 116 + 0, string, 4);
  printf("suite id: %s\n", string);
  read_short(fp, 116 + 8, &s);
  printf("exam number: %hd\n", s);
  read_int(fp, 116 + 84, &i);
  printf("magnet strength in gauss: %d\n",i);
  read_string(fp, 116 + 10, string, 33);
  printf("hospital name: %s\n", string);
  read_string(fp, 116 + 88, string, 13);
  printf("patient id: %s\n", string);
  read_string(fp, 116 + 101, string, 25);
  printf("patient name: %s\n", string);
  read_short(fp, 116 + 126, &s);
  printf("patient age: %hd\n", s);
  read_short(fp, 116 + 128, &s);
  printf("patient age notation: %hd\n", s);
  read_short(fp, 116 + 130, &s);
  printf("patient sex: %hd\n", s);
  read_int(fp, 116 + 132, &i);
  printf("patient weight: %d\n",i);
  read_string(fp, 116 + 138, string, 61);
  printf("patient history: %s\n", string);
  read_string(fp, 116 + 216, string, 33);
  printf("referring physician: %s\n", string);
  read_string(fp, 116 + 249, string, 33);
  printf("diagnostician/radiologist: %s\n", string);
  read_string(fp, 116 + 286, string, 23);
  printf("exam description: %s\n", string);
  read_string(fp, 116 + 324, string, 9);
  printf("creator suite and host: %s\n", string);

  read_short(fp, 1156 + 8, &s);
  printf("series number: %hd\n", s);
  read_string(fp, 1156 + 92, string, 25);
  printf("scan protocol: %s\n", string);
  read_string(fp, 1156 + 92, string, 25);
  printf("scan protocol: %s\n", string);

  read_short(fp, 2184 + 12, &s);
  printf("image number: %hd\n", s);
  read_float(fp, 2184 + 28, &f);
  printf("slice thickness: %g\n", f);
  read_short(fp, 2184 + 32, &s);
  printf("image matrix size (x): %hd\n",s);
  read_short(fp, 2184 + 34, &s);
  printf("image matrix size (y): %hd\n",s);
  read_float(fp, 2184 + 36, &f);
  printf("display field of view (x): %g\n",f);
  read_float(fp, 2184 + 40, &f);
  printf("display field of view (y, if different from x): %g\n",f);
  read_float(fp, 2184 + 44, &f);
  printf("image x dimension: %g\n", f);
  read_float(fp, 2184 + 48, &f);
  printf("image y dimension: %g\n", f);
  read_float(fp, 2184 + 52, &f);
  printf("image x pixel size: %g\n", f);
  read_float(fp, 2184 + 56, &f);
  printf("image y pixel size: %g\n", f);
  read_string(fp, 2184 + 74, string, 17);
  printf("iv contrast agent: %s\n", string);
  read_string(fp, 2184 + 91, string, 17);
  printf("oral contrast agent: %s\n", string);
  read_short(fp, 2184 + 116, &s);
  printf("plane type: %hd\n", s);
  read_float(fp, 2184 + 120, &f);
  printf("scan spacing: %g\n", f);
  read_string(fp, 2184 + 320, string, 33);
  printf("pulse sequence name: %s\n", string);
  read_string(fp, 2184 + 376, string, 17);
  printf("coil name: %s\n", string);
  read_int(fp, 2184 + 404, &i);
  printf("calibrated field strength (x10 uGauss): %d\n", i);
  read_short(fp, 2184 + 416, &s);
  printf("number of slices in this scan group: %hd\n", s);
  fclose(fp);

} /* end read_ge_8x_file() */

static void read_siemens_file(char *fname, struct stat stat_buf)
{

  int i1, i2, i3;
  FILE *fp;
  char string[100];
  double d;
  short s;

  if((fp = fopen(fname, "r")) == NULL)
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: Siemens\n");

  read_int(fp, 0, &i1);
  read_int(fp, 4, &i2);
  read_int(fp, 8, &i3);
  printf("study date: %s %d, %d\n", month[i2], i3, i1);
  read_string(fp, 105, string, 27);
  printf("institution name: %s\n", string);
  read_string(fp, 768, string, 25);
  printf("patient name: %s\n", string);
  read_string(fp, 795, string, 12);
  printf("patient id: %s\n", string);
  read_int(fp, 808, &i1);
  read_int(fp, 812, &i2);
  read_int(fp, 816, &i3);
  printf("patient date of birth: %s %d, %d\n", month[i2], i3, i1);
  read_string(fp, 851, string, 4);
  string[4] = '\0';
  printf("patient age: %s\n", string);
  read_double(fp, 1544, &d);
  printf("slice thickness: %g\n", d);
  read_double(fp, 1560, &d);
  printf("repetition time: %g\n", d);
  read_double(fp, 1568, &d);
  printf("echo time: %g\n", d);
  read_string(fp, 1767, string, 16);
  printf("coil: %s\n",string);
  read_double(fp, 2560, &d);
  printf("field strength: %g\n", d);
  read_string(fp, 2944, string, 65);
  printf("parameter file name: %s\n",string);
  read_short(fp, 4994, &s);
  printf("image matrix size (rows): %hd\n", s);
  read_short(fp, 4996, &s);
  printf("image matrix size (columns): %hd\n", s);
  read_int(fp, 4004, &i1);
  printf("nominal number of slices: %d\n", i1);
  read_int(fp, 2864, &i1);
  printf("base raw matrix size: %d\n", i1);
  read_double(fp, 4136, &d);
  printf("slice distance factor: %g\n", d);
  read_int(fp, 1584, &i1);
  printf("number of averages: %d\n", i1);
  read_double(fp, 5000, &d);
  printf("pixel size row: %g\n", d);
  read_double(fp, 5008, &d);
  printf("pixel size column: %g\n", d);
  read_string(fp, 3009, string, 65);
  printf("sequence name: %s\n",string);
  read_string(fp, 5814, string, 7);
  printf("slice direction: %s\n",string);
  read_int(fp, 1052, &i1);
  read_int(fp, 1056, &i2);
  read_int(fp, 1060, &i3);
  printf("registration date: %04d%02d%02d\n", i1, i2, i3);
  read_int(fp, 1064, &i1);
  read_int(fp, 1068, &i2);
  read_int(fp, 1072, &i3);
  printf("registration time: %02d%02d%02d\n", i1, i2, i3);
  read_string(fp, 3904, string, 32);
  printf("experiment name: %s\n", string);
  read_string(fp, 358, string, 25);
  printf("experimenter: %s\n", string);
  read_string(fp, 281, string, 27);
  printf("manufacturer model: %s\n", string);
  read_string(fp, 1612, string, 27);
  printf("device serial number: %s\n", string);

  fclose(fp);

} /* end read_siemens_file() */

static void read_minc_file(char *fname, struct stat stat_buf)
{

  Volume volume;
  char *dim_names[4];
  volume_input_struct input_info;
  static char MIfspace[] = "frame space";
  int sizes[4];
  int ndim;

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: MINC\n");

  dim_names[0] = MIxspace;
  dim_names[1] = MIyspace;
  dim_names[2] = MIzspace;
  dim_names[3] = MIfspace;

  if(start_volume_input(fname, 0, dim_names, NC_UNSPECIFIED, TRUE, 0, 0, TRUE, &volume, NULL, &input_info))
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  get_volume_sizes(volume, sizes);
  ndim = get_volume_n_dimensions(volume);

  printf("width: %d\n", sizes[0]);
  printf("height: %d\n", sizes[1]);
  printf("depth: %d\n", sizes[2]);
  printf("number of frames: %d\n", (ndim < 4 ? 1 : sizes[3]));
  printf("data type: ");
  switch(volume->nc_data_type)
  {
    case NC_BYTE:
      printf("byte\n");
      break;

    case NC_CHAR:
      printf("char\n");
      break;

    case NC_SHORT:
      printf("short\n");
      break;

    case NC_LONG:
      printf("long\n");
      break;

    case NC_FLOAT:
      printf("float\n");
      break;

    case NC_DOUBLE:
      printf("double\n");
      break;

    default:
      printf("unknown\n");
      break;

  }

} /* end read_minc_file() */

static void read_mgh_file(char *fname, struct stat stat_buf)
{

  FILE *fp;
  int width, height, depth, nframes, type, dof, version;

  if((fp = fopen(fname, "r")) == NULL)
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: MGH\n");

  version = freadInt(fp) ;
  width = freadInt(fp) ;
  height = freadInt(fp) ;
  depth =  freadInt(fp) ;
  nframes = freadInt(fp) ;
  type = freadInt(fp) ;
  dof = freadInt(fp) ;

  printf("version: %d\n", version);
  printf("width: %d\n", width);
  printf("height: %d\n", height);
  printf("depth: %d\n", depth);
  printf("number of frames: %d\n", nframes);
  printf("type: %s\n", type_text[type]);
  printf("dof: %d\n", dof);

  fclose(fp);

} /* end read_mgh_file() */

static void read_analyze_file(char *fname, struct stat stat_buf)
{

  FILE *fp;
  dsr header;

  if((fp = fopen(fname, "r")) == NULL)
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: analyze\n");

  read_analyze_header(fname, &header);

  printf("byte order: ");
  if(header.hk.sizeof_hdr != sizeof(header))
  {
    flip_analyze_header(&header);
    if(header.hk.sizeof_hdr != sizeof(header))
    {
      printf("unknown\n");
      return;
    }
#if (BYTE_ORDER == LITTLE_ENDIAN)
    printf("big-endian\n");
#else
    printf("little-endian\n");
#endif
  }
  else
  {
#if (BYTE_ORDER == LITTLE_ENDIAN)
    printf("little-endian\n");
#else
    printf("big-endian\n");
#endif
  }

  printf("width: %d\n", header.dime.dim[1]);
  printf("height: %d\n", header.dime.dim[2]);
  printf("depth: %d\n", header.dime.dim[3]);
  printf("voxel size: %g, %g, %g, %g, %g, %g, %g, %g\n", header.dime.pixdim[0], header.dime.pixdim[1], header.dime.pixdim[2], header.dime.pixdim[3], header.dime.pixdim[4], header.dime.pixdim[5], header.dime.pixdim[6], header.dime.pixdim[7]);
  printf("originator: %hd, %hd, %hd\n", 
*(short *)&header.hist.originator[0], *(short *)&header.hist.originator[2], *(short *)&header.hist.originator[4]);
  printf("data type: ");
  switch(header.dime.datatype)
  {
    case DT_BINARY:
      printf("binary\n");
      break;

    case DT_UNSIGNED_CHAR:
      printf("unsigned char\n");
      break;

    case DT_SIGNED_SHORT:
      printf("short\n");
      break;

    case DT_SIGNED_INT:
      printf("int\n");
      break;

    case DT_FLOAT:
      printf("float\n");
      break;

    case DT_DOUBLE:
      printf("double\n");
      break;

     default:
       printf("unknown\n");
  }

  fclose(fp);

} /* end read_analyze_file() */

static void read_sdt_file(char *fname, struct stat stat_buf)
{

  char hfname[STR_LEN];
  char *dot;
  FILE *fp;
  char line[STR_LEN], *colon, *tc;
  float vs;

  printf("file name: %s\n", fname);
  printf("file size: %d\n", (int)(stat_buf.st_size));
  printf("file type: SDT\n");

  strcpy(hfname, fname);

  if((dot = strrchr(hfname, '.')))
    {
    dot++;
    sprintf(dot, "spr");
    }
  else
    strcat(hfname, ".spr");

  if((fp = fopen(hfname, "r")) == NULL)
  {
    fprintf(stderr, "error opening file %s\n", fname);
    return;
  }

  while(!feof(fp))
  {
    fgets(line, STR_LEN, fp);
    if((colon = strchr(line, ':')))
    {

      if(line[strlen(line)-1] == '\n')
        line[strlen(line)-1] = '\0';

      *colon = '\0';

      do
        {
        colon++;
        } while(isspace((int)*colon));

      if(strcmp(line, "numDim") == 0)
        printf("number of dimensions: %s\n", colon);
      else if(strcmp(line, "dim") == 0)
        printf("dimensions: %s\n", colon);
      else if(strcmp(line, "dataType") == 0)
      {
        if(strncmp(colon, "BYTE", 4) == 0)
          printf("data type: character\n");
        else if(strncmp(colon, "WORD", 4) == 0)
          printf("data type: short\n");
        else if(strncmp(colon, "LWORD", 5) == 0)
          printf("data type: integer\n");
        else if(strncmp(colon, "REAL", 4) == 0)
          printf("data type: float\n");
        else if(strncmp(colon, "COMPLEX", 7) == 0)
          printf("data type: complex\n");
        else
          printf("unknown data type: %s\n", colon);
      }
      else if(strcmp(line, "interval") == 0)
      {
        printf("voxel size:");
        line[strlen(line)+1] = '\0';
        while(*colon != '\0')
          {
          tc = colon;
          while(isdigit((int)*colon) || *colon == '-' || 
                *colon == '+' || *colon == '.')
            colon++;
          while(isspace((int)*colon))
            colon++;
          *(colon-1) = '\0';
          sscanf(tc, "%f", &vs);
          printf(" %g", vs*10.0);
          }
        printf("\n");
      }
      else if(strcmp(line, "sdtOrient") == 0)
      {
        if(strncmp(colon, "sag", 3) == 0)
          printf("orientation: sagittal\n");
        else if(strncmp(colon, "ax", 2) == 0)
          printf("orientation: horizontal\n");
        else if(strncmp(colon, "cor", 3) == 0)
          printf("orientation: coronal\n");
        else
        printf("orientation: unknown (%s)\n", colon);
      }
      else
      {
      }
    }
  }

  fclose(fp);

}

static void read_cor(char *fname)
{

  char header_fname[STRLEN];
  struct stat stat_buf;
  FILE *fp;
  int i, i2;
  float f;

  strcpy(header_fname, fname);

  if(header_fname[strlen(header_fname) - 1] != '/')
  {
    header_fname[strlen(header_fname) + 1] = '\0';
    header_fname[strlen(header_fname)] = '/';
  }

  if(stat(header_fname, &stat_buf) < 0)
  {
    fprintf(stderr, "can't stat %s\n", header_fname);
    return;
  }

  if(!S_ISDIR(stat_buf.st_mode))
  {
    fprintf(stderr, "%s isn't a directory of COR- files\n", header_fname);
    return;
  }

  sprintf(header_fname, "%sCOR-.info", header_fname);

  if((fp = fopen(header_fname, "r")) == NULL)
  {
    fprintf(stderr, "can't open file %s\n", header_fname);
    return;
  }

  printf("file name: %s\n", fname);
  printf("file type: coronal slice directory\n");

  fscanf(fp, "%*s %d", &i);
  fscanf(fp, "%*s %d", &i2);
  printf("number of images: %d\n", i2 - i + 1);
  fscanf(fp, "%*s %d", &i);
  printf("data type: %s\n", type_text[i]);
  fscanf(fp, "%*s %d", &i);
  printf("width: %d\n", i);
  fscanf(fp, "%*s %d", &i);
  printf("height: %d\n", i);
  fscanf(fp, "%*s %f", &f);
  printf("field of view: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("slice thickness: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("pixel size: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("location: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("start x: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("end x: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("start y: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("end y: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("start z: %g\n", f);
  fscanf(fp, "%*s %f", &f);
  printf("end z: %g\n", f);

  fclose(fp);

} /* end read_cor() */


static void read_ge_header(FILE *fp, int offset, struct ge_header *header)
{

  char buf[156];

  fseek(fp, offset, 0);
  fread(buf, 1, 156, fp);

  header->magic = orderIntBytes(*(int *)(&buf[0]));
  header->width = orderIntBytes(*(int *)(&buf[8]));
  header->height = orderIntBytes(*(int *)(&buf[12]));
  header->depth = orderIntBytes(*(int *)(&buf[16]));
  header->exam_header_offset = orderIntBytes(*(int *)(&buf[132]));
  header->series_header_offset = orderIntBytes(*(int *)(&buf[144]));
  header->image_header_offset = orderIntBytes(*(int *)(&buf[148]));
  header->compression = orderIntBytes(*(int *)(&buf[20]));

  printf("width: %d\n",header->width);
  printf("height: %d\n",header->height);
  printf("data point length (in bits): %d\n",header->depth);
  printf("compression: ");
  switch(header->compression)
  {
    case GE_COMPRESSION_ASIS:
      printf("as is (uncompressed)\n");
      break;
    case GE_COMPRESSION_RECTANGULAR:
      printf("rectangular\n");
      break;
    case GE_COMPRESSION_PACKED:
      printf("packed\n");
      break;
    case GE_COMPRESSION_COMPRESSED:
      printf("compressed\n");
      break;
    case GE_COMPRESSION_COMPRESSED_AND_PACKED:
      printf("compressed and packed\n");
      break;
    default:
      printf("unknown\n");
      break;
  }

} /* end read_ge_header() */

/* lifted from mriio.c (and altered) */
static void
read_analyze_header(char *fname, dsr *bufptr)
{
  FILE *fp;
  char hdr_fname[STRLEN], *dot;
  int  nread;

  strcpy(hdr_fname, fname) ;
  dot = strrchr(hdr_fname, '.') ;
  if (dot)
    *dot = 0 ;
  strcat(hdr_fname, ".hdr") ;
  fp = fopen(hdr_fname,"rb");

  nread = fread(bufptr,sizeof(char), sizeof(dsr),fp);
  fclose(fp);

}

static void flip_analyze_header(dsr *header)
{

  int i;
  char c;

  header->hk.sizeof_hdr = swapInt(header->hk.sizeof_hdr);
  header->hk.extents = swapInt(header->hk.extents);
  header->hk.session_error = swapShort(header->hk.session_error);

  for(i = 0;i < 8;i++)
  {
    header->dime.dim[i] = swapShort(header->dime.dim[i]);
    header->dime.pixdim[i] = swapFloat(header->dime.pixdim[i]);
  }

  header->dime.unused1 = swapShort(header->dime.unused1);
  header->dime.datatype = swapShort(header->dime.datatype);
  header->dime.bitpix = swapShort(header->dime.bitpix);
  header->dime.dim_un0 = swapShort(header->dime.dim_un0);
  header->dime.vox_offset = swapFloat(header->dime.vox_offset);
  header->dime.roi_scale = swapFloat(header->dime.roi_scale);
  header->dime.funused1 = swapFloat(header->dime.funused1);
  header->dime.funused2 = swapFloat(header->dime.funused2);
  header->dime.cal_max = swapFloat(header->dime.cal_max);
  header->dime.cal_min = swapFloat(header->dime.cal_min);
  header->dime.compressed = swapInt(header->dime.compressed);
  header->dime.verified = swapInt(header->dime.verified);
  header->dime.glmax = swapInt(header->dime.glmax);
  header->dime.glmin = swapInt(header->dime.glmin);

  header->hist.views = swapInt(header->hist.views);
  header->hist.vols_added = swapInt(header->hist.vols_added);
  header->hist.start_field = swapInt(header->hist.start_field);
  header->hist.field_skip = swapInt(header->hist.field_skip);
  header->hist.omax = swapInt(header->hist.omax);
  header->hist.omin = swapInt(header->hist.omin);
  header->hist.smax = swapInt(header->hist.smax);
  header->hist.smin = swapInt(header->hist.smin);

  c = header->hist.originator[0]; header->hist.originator[0] = header->hist.originator[1] ; header->hist.originator[1] = c;
  c = header->hist.originator[2]; header->hist.originator[2] = header->hist.originator[3] ; header->hist.originator[3] = c;
  c = header->hist.originator[4]; header->hist.originator[4] = header->hist.originator[5] ; header->hist.originator[5] = c;
  c = header->hist.originator[6]; header->hist.originator[6] = header->hist.originator[7] ; header->hist.originator[7] = c;
  c = header->hist.originator[8]; header->hist.originator[8] = header->hist.originator[9] ; header->hist.originator[9] = c;

}  /*  end flip_analyze_header()  */

/* EOF */
#endif
