#ifndef ANNOTATION_H
#define ANNOTATION_H

int   read_annotation_table(void) ;
char  *annotation_to_name(int annotation, int *pindex) ;
int   annotation_to_index(int annotation) ;
int   print_annotation_table(FILE *fp);

#endif
