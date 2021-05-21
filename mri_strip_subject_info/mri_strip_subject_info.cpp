/*
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


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include "mri.h"
#include "error.h"
#include "mri_identify.h"
#include "machine.h"
#include "version.h"

char *get_base_name(char *fullpath);
int fix_genesis(char *fname, char *dname);
int fix_siemens(char *fname, char *dname);

extern int errno;
const char *Progname;

void usage(void) {

  fprintf(stderr, "usage: %s <input file> ... <output directory>\n", Progname);

} /* end usage() */

int main(int argc, char *argv[]) {

  struct stat stat_buf;
  int i;
  int nargs;

  nargs = handleVersionOption(argc, argv, "mri_strip_subject_info");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = get_base_name(argv[0]);

  if (argc < 3) {
    usage();
    exit(1);
  }

  if (stat(argv[argc-1], &stat_buf) != 0) {
    ErrorExit(ERROR_BADFILE, "%s: error on stat(%s)", Progname, argv[argc-1]);
  }

  if (!S_ISDIR(stat_buf.st_mode)) {
    errno = 0;
    ErrorExit(ERROR_BADFILE, "%s: %s is not a directory", Progname, argv[argc-1]);
    exit(1);
  }

  for (i = 1;i < argc-1;i++) {
    printf("%s\n", argv[i]);
    if (stat(argv[i], &stat_buf) != 0) {
      fprintf(stderr, "%s: error accessing %s\n", Progname, argv[i]);
    } else if (S_ISDIR(stat_buf.st_mode)) {
      fprintf(stderr, "%s: %s is a directory\n", Progname, argv[i]);
    } else if (is_genesis(argv[i])) {
      fix_genesis(argv[i], argv[argc-1]);
    } else if (is_siemens(argv[i])) {
      fix_siemens(argv[i], argv[argc-1]);
    } else
      fprintf(stderr, "%s: file %s is not GE or Siemens\n", Progname, argv[i]);
  }

  exit(0);

} /* end main() */

char *get_base_name(char *fullpath) {

  char *bn;

  bn = strrchr(fullpath, '/');
  bn = (bn == NULL ? fullpath : bn+1);

  return(bn);

} /* end get_base_name() */

int fix_genesis(char *fname, char *dname) {

  FILE *fp;
  char *basename;
  char out_fname[STRLEN];
  char *buf;
  int file_length;
  int exam_header_offset;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error opening file %s", fname));
  }

  fseek(fp, 0, SEEK_END);
  file_length = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  buf = (char *)malloc(file_length);
  if (buf == NULL) {
    ErrorReturn(ERROR_NOMEMORY, (ERROR_NOMEMORY, "error: no memory to read file %s", fname));
  }

  if (fread(buf, 1, file_length, fp) != file_length) {
    free(buf);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error reading from file %s", fname));
  }

  fclose(fp);

  memcpy(&exam_header_offset, &(buf[132]), 4);
  exam_header_offset = orderIntBytes(exam_header_offset);

  if (exam_header_offset + 282 + 23 >= file_length) {
    free(buf);
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "file %s is too short", fname));
  }

  memset(&(buf[exam_header_offset+10]), 0x00, 33);     /*  hospital name              */
  memset(&(buf[exam_header_offset+97]), 0x00, 25);     /*  patient name               */
  memset(&(buf[exam_header_offset+122]), 0x00, 2);     /*  patient age                */
  memset(&(buf[exam_header_offset+126]), 0x00, 2);     /*  patient sex                */
  memset(&(buf[exam_header_offset+134]), 0x00, 61);    /*  patient history            */
  memset(&(buf[exam_header_offset+212]), 0x00, 33);    /*  referring physician        */
  memset(&(buf[exam_header_offset+245]), 0x00, 33);    /*  diagnostician/radiologist  */
  memset(&(buf[exam_header_offset+282]), 0x00, 23);    /*  exam description           */

  basename = get_base_name(fname);
  sprintf(out_fname, "%s/%s", dname, basename);

  fp = fopen(out_fname, "w");
  if (fp == NULL) {
    free(buf);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error opening file %s", fname));
  }

  if (fwrite(buf, 1, file_length, fp) != file_length) {
    free(buf);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error writing to file %s", out_fname));
  }

  fclose(fp);

  free(buf);

  return(NO_ERROR);

} /* end fix_genesis() */

int fix_siemens(char *fname, char *dname) {

  FILE *fp;
  char *basename;
  char out_fname[STRLEN];
  char *buf;
  int file_length;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error opening file %s", fname));
  }

  fseek(fp, 0, SEEK_END);
  file_length = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  buf = (char *)malloc(file_length);
  if (buf == NULL) {
    ErrorReturn(ERROR_NOMEMORY, (ERROR_NOMEMORY, "error: no memory to read file %s", fname));
  }

  if (fread(buf, 1, file_length, fp) != file_length) {
    free(buf);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error reading from file %s", fname));
  }

  fclose(fp);

  if (6058 + 27 >= file_length) {
    free(buf);
    errno = 0;
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "file %s is too short", fname));
  }

  memset(&(buf[105]), 0x00, 27);     /*  G08.Ide.InstitutionID              */
  memset(&(buf[132]), 0x00, 27);     /*  G08.Ide.ReferringPhysician         */
  memset(&(buf[308]), 0x00, 25);     /*  G08.Ide.AttendingPhysician         */
  memset(&(buf[333]), 0x00, 25);     /*  G08.Ide.Radiologist                */
  memset(&(buf[358]), 0x00, 25);     /*  G08.Ide.OperatorIdentification     */
  memset(&(buf[638]), 0x00, 25);     /*  G09.Ide.OperatorIdentification     */
  memset(&(buf[768]), 0x00, 27);     /*  G10.Pat.PatientName                */
  memset(&(buf[795]), 0x00, 13);     /*  G10.Pat.PatientId                  */
  memset(&(buf[808]), 0x00, 4);      /*  G10.Pat.PatientBirthdate.Year      */
  memset(&(buf[812]), 0x00, 4);      /*  G10.Pat.PatientBirthdate.Month     */
  memset(&(buf[816]), 0x00, 4);      /*  G10.Pat.PatientBirthdate.Day       */
  memset(&(buf[820]), 0x00, 4);      /*  G10.Pat.PatientSex                 */
  memset(&(buf[824]), 0x00, 27);     /*  G10.Pat.PatientMaidenName          */
  memset(&(buf[851]), 0x00, 5);      /*  G10.Pat.PatientAge                 */
  memset(&(buf[856]), 0x00, 8);      /*  G10.Pat.PatientSize                */
  memset(&(buf[864]), 0x00, 4);      /*  G10.Pat.PatientWeight              */
  memset(&(buf[1080]), 0x00, 4);     /*  G11.Pat.UsedPatientWeight          */
  memset(&(buf[1088]), 0x00, 27);    /*  G11.Pat.CaseIdentification         */
  memset(&(buf[1115]), 0x00, 27);    /*  G11.Pat.RequestIdentification      */
  memset(&(buf[1152]), 0x00, 27);    /*  G13.PatMod.ModifyingPhysician      */
  memset(&(buf[1208]), 0x00, 65);    /*  G13.PatMod.PatientName_DICOM       */
  memset(&(buf[1273]), 0x00, 65);    /*  G13.PatMod.PatientId_DICOM         */
  memset(&(buf[1340]), 0x00, 4);     /*  G13.PatMod.PatientBirthdate.Year   */
  memset(&(buf[1344]), 0x00, 4);     /*  G13.PatMod.PatientBirthdate.Month  */
  memset(&(buf[1348]), 0x00, 4);     /*  G13.PatMod.PatientBirthdate.Day    */
  memset(&(buf[1352]), 0x00, 4);     /*  G13.PatMod.PatientWeight           */
  memset(&(buf[1356]), 0x00, 4);     /*  G13.PatMod.PatientSex              */
  memset(&(buf[1360]), 0x00, 27);    /*  G13.PatMod.Comment1                */
  memset(&(buf[1387]), 0x00, 27);    /*  G13.PatMod.Comment2                */
  memset(&(buf[1612]), 0x00, 27);    /*  G18.Acq.DeviceSerialNumber         */
  memset(&(buf[3904]), 0x00, 27);    /*  G21.Rel1.CM.StudyName              */
  memset(&(buf[5056]), 0x00, 65);    /*  G28.Pre.PatientName_DICOM          */
  memset(&(buf[5121]), 0x00, 65);    /*  G28.Pre.PatientId_DICOM            */
  memset(&(buf[5504]), 0x00, 13);    /*  G51.Txt.PatientNumber              */
  memset(&(buf[5517]), 0x00, 12);    /*  G51.Txt.PatientSexAndAge           */
  memset(&(buf[5655]), 0x00, 27);    /*  G51.Txt.InstallationName           */
  memset(&(buf[5938]), 0x00, 12);    /*  G51.Txt.StudyNumber                */
  memset(&(buf[5956]), 0x00, 12);    /*  G51.Txt.PatientBirthdate           */
  memset(&(buf[6058]), 0x00, 27);    /*  G51.Txt.PatientName                */

  basename = get_base_name(fname);
  sprintf(out_fname, "%s/%s", dname, basename);

  fp = fopen(out_fname, "w");
  if (fp == NULL) {
    free(buf);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error opening file %s", fname));
  }

  if (fwrite(buf, 1, file_length, fp) != file_length) {
    free(buf);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "error writing to file %s", out_fname));
  }

  fclose(fp);

  free(buf);

  return(NO_ERROR);

} /* end fix_siemens() */

/* eof */
