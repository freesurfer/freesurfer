#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "diag.h"
#include "error.h"
#include "mri.h"
#include "mri2.h"
#include "fmriutils.h"
#include "mri_identify.h"
#include "gcamorph.h"
#include "DICOMRead.h"
#include "version.h"
#include "utils.h"
#include "macros.h"
#include "fmriutils.h"

extern int errno;
char *Progname;

void get_string(int argc, char *argv[], int *pos, char *val);
void usage_message(FILE *stream);

int main(int argc, char *argv[]) {

  char in_name[STRLEN];  
  in_name[0] = '\0';
  int i;
  FILE* fp;

  if (argc <= 1)
    printf("Input filename is required! Use -i or --input_volume and the name of the file.\n");

  for(i = 1;i < argc;i++) {
    if(strcmp(argv[i], "-i") == 0 ||
       strcmp(argv[i], "--input_volume") == 0)
      get_string(argc, argv, &i, in_name);
    else
      printf("Input filename is required! Use -i or --input_volume and the name of the file.\n");
  }
  if ((fp = fopen(in_name, "r")) != NULL)
    {
      list_labels_in_otl_file(fp);
    }

  return 0;

}


void get_string(int argc, char *argv[], int *pos, char *val) {

  if (*pos + 1 >= argc) {
    fprintf(stderr, "\n%s: argument %s expects an extra argument; "
            "none found\n", Progname, argv[*pos]);
    usage_message(stdout);
    exit(1);
  }

  strcpy(val, argv[*pos+1]);

  (*pos)++;

} /* end get_string() */

void usage_message(FILE *stream) {

  fprintf(stream, "\n");
  fprintf(stream, "type %s -u for usage\n", Progname);
  fprintf(stream, "\n");

} /* end usage_message() */
