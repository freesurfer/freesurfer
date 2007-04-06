/**
 * @file  annotation.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/04/06 06:11:51 $
 *    $Revision: 1.10 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef ANNOTATION_H
#define ANNOTATION_H


#ifdef ANNOTATION_SRC
char *annotation_table_file = NULL;
#else
extern char *annotation_table_file;
#endif

int   read_annotation_table(void) ;
int   read_named_annotation_table(char *fname) ;
char *index_to_name(int index);
char  *annotation_to_name(int annotation, int *pindex) ;
int   annotation_to_index(int annotation) ;
int   print_annotation_table(FILE *fp);
int print_annotation_colortable(FILE *fp);
int   index_to_annotation(int index) ;
LABEL *annotation2label(int annotid, MRIS *Surf);
int set_atable_from_ctable(COLOR_TABLE *pct);
int  MRISdivideAnnotation(MRI_SURFACE *mris, int *nunits) ;
int  MRISdivideAnnotationUnit(MRI_SURFACE *mris, int annot, int nunits) ;

#endif
