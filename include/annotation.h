#ifndef ANNOTATION_H
#define ANNOTATION_H

int   read_annotation_table(void) ;
int   read_named_annotation_table(char *fname) ;
char  *annotation_to_name(int annotation, int *pindex) ;
int   annotation_to_index(int annotation) ;
int   print_annotation_table(FILE *fp);
int   index_to_annotation(int index) ;

#endif
