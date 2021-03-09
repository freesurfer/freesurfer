/**
 * @brief compare two parcellated (annotated) surfaces and computes accuracy.
 *
 * Compares two parcellated (annotated) surfaces
 * and computes an overall Dice coefficient, and computes mean min distance
 * for each label (indicating in mm the error).
 *
 * Usage:
 *   mris_compute_parc_overlap --s subject --hemi hemi \
 *   --annot1 annotfile --annot2 annotfile
 *
 * Required:
 *   --s subject          subject to check
 *   --hemi hemi          hemisphere: rh or lh
 *   --annot1 annotfile   first .annot file
 *   --annot2 annotfile   second .annot file
 *
 * Optional:
 *   --sd subj_dir        set SUBJECTS_DIR
 *   --use-labels file    file containing labels to check, one per line
 *   --version            version info
 *   --help               this usage info
 *
 * Example:
 *   mris_compute_parc_overlap --s bert --hemi lh \
 *     --annot1 aparc --annot2 aparc.ernie
 *
 *   mris_compute_parc_overlap --s bert --hemi lh \
 *     --label1 precentral --label2 precentral
 */
/*
 * Original Author: Nick Schmansky
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

#include "mrisurf.h"
#include "annotation.h"
#include "version.h"
#include "error.h"
#include "utils.h"

#include "compilerdefs.h"

#define MAX_VNOS 200000
typedef struct _labelInfo
{
  int annotation;      /* its annotation identifier */
  char vnos[MAX_VNOS]; /* if array element is 1, then that array index, which
                          is a vertex number, is a border vertex of this 
                          label */
  float meanMinDist;   /* average of the minimum distance of every boundary
                          vertex to the other surfaces boundary vertices for
                          this same label */
}
LABEL_INFO;
#define MAX_LABELS 100
static LABEL_INFO surf1BoundaryLabels[MAX_LABELS];
static LABEL_INFO surf2BoundaryLabels[MAX_LABELS];

#define MAX_SKIPPED_LABELS 100000
static char skippedLabels[MAX_SKIPPED_LABELS];

#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
#endif
// this mini colortable is used when two labels are being compared
static const COLOR_TABLE_ENTRY unknown = 
{"unknown", 0,0,0,255, 0,0,0,255};
static COLOR_TABLE_ENTRY userLabel = 
{ "user label name gets copied here                   ", 
  220,20,20,255, 0.8,0.08,0.08,1};
static const CTE *entries[2] = {&unknown, &userLabel};
static const COLOR_TABLE miniColorTable = 
{(CTE**)entries, 2, "miniColorTable", 2};
#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic pop
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic pop
#endif

static void addToExcludedLabelsList(COLOR_TABLE *ct, const char *labelToExclude);
static int isExcludedLabel(int colortabIndex);
static void addToIncludedLabelsList(COLOR_TABLE *ct, const char *labelToInclude);
static int isIncludedLabel(int colortabIndex);
static void calcMeanMinLabelDistances(void);
static void usage(int exit_val);
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void dump_options(FILE *fp);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static int  singledash(char *flag);
static void padWhite(char* str, int maxLen);

const char *Progname;
static char *FREESURFER_HOME = NULL;
static char *SUBJECTS_DIR = NULL;
static char *subject = NULL;
static char *hemi = NULL;
static char *annot1 = NULL;
static char *annot2 = NULL;
static int unknown_annot=0;
static char *label1name = NULL;
static char *label2name = NULL;
static LABEL *label1 = NULL;
static LABEL *label2 = NULL;
static MRIS *surface1 = NULL;
static MRIS *surface2 = NULL;
static int debug_overlap = 0;
static int debug_boundaries = 0;
static int debug_labels = 0;
static int use_label1_xyz = 0;
static int use_label2_xyz = 0;
static int check_label1_xyz = 1;
static int check_label2_xyz = 1;
static char tmpstr[2000];
static char *labelsfile=NULL; // optional file containing lists of label to chk
static FILE *logfile = NULL ;

int main(int argc, char *argv[])
{
  char filename[2000];

  // clear our label info arrays
  memset(surf1BoundaryLabels,0,sizeof(surf1BoundaryLabels));
  memset(surf2BoundaryLabels,0,sizeof(surf2BoundaryLabels));
  memset(skippedLabels,0,sizeof(skippedLabels));

  // command-line processing...
  int nargs;
  nargs = handleVersionOption(argc, argv, "mris_compute_parc_overlap");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc--;
  argv++;
  if (argc == 0) usage(1);
  if (argc < 8) usage(1);

  FREESURFER_HOME = getenv("FREESURFER_HOME");
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  if (SUBJECTS_DIR==NULL)
  {
    printf
    ("ERROR: SUBJECTS_DIR not defined in environment or on command-line\n");
    exit(1);
  }

  /* ------ Load subject's inflated surface ------ */
  sprintf(tmpstr,"%s/%s/surf/%s.inflated",SUBJECTS_DIR,subject,hemi);
  printf("\nReading %s inflated surface \n %s\n",hemi,tmpstr);
  surface1 = MRISread(tmpstr);
  if (surface1 == NULL)
  {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  surface2 = MRISclone(surface1);
  MRIS *overlapSurf = MRISclone(surface1); // a debugging surface

  /* --- Load subject's 'reference' annotation or label onto surface1 ---- */
  if (annot1)
  {
    sprintf(filename,"%s/%s/label/%s.%s.annot",
            SUBJECTS_DIR,subject,hemi,annot1);
    printf("\nLoading %s annotations from %s\n",hemi,filename);
    fflush(stdout);
    int err = MRISreadAnnotation(surface1, filename);
    if (err)
    {
      printf("ERROR: MRISreadAnnotation() failed %s\n",filename);
      exit(1);
    }
  }
  else if (label1name) 
  {
    sprintf(filename,"%s.%s",hemi,label1name);
    printf("Loading label file %s...",filename);
    fflush(stdout);
    label1 = LabelRead(subject,filename);
    if (label1==NULL)
    {
      printf("ERROR: LabelRead() failed to read %s\n",filename);
      exit(1);
    }
    printf(" %d vertices.\n",label1->n_points);
    surface1->ct = (CT*)&miniColorTable;
    // put this label name in our color table
    strcpy(userLabel.name,label1name);
  }


  /* ------- Load the 'test' annotation or label onto surface2 --------- */
  if (annot2)
  {
    sprintf(filename,"%s/%s/label/%s.%s.annot",
            SUBJECTS_DIR,subject,hemi,annot2);
    printf("\nLoading %s annotations from %s\n",hemi,filename);
    fflush(stdout);
    int err = MRISreadAnnotation(surface2, filename);
    if (err)
    {
      printf("ERROR: MRISreadAnnotation() failed %s\n",filename);
      exit(1);
    }
    MRISreadAnnotation(overlapSurf, filename);
  }
  else if (label2name) 
  {
    sprintf(filename,"%s.%s",hemi,label2name);
    printf("Loading label file %s...",filename);
    fflush(stdout);
    label2 = LabelRead(subject,filename);
    if (label2==NULL)
    {
      printf("ERROR: LabelRead() failed to read %s\n",filename);
      exit(1);
    }
    printf(" %d vertices.\n",label2->n_points);
    surface2->ct = (CT*)&miniColorTable;
    overlapSurf->ct = (CT*)&miniColorTable;
    strcat(userLabel.name,"->");
    strcat(userLabel.name,label2name);
  }

  /* make sure there are a valid color look-up tables */
  if (surface1->ct == NULL)
  {
    printf("ERROR: %s does not contain color table!\n",annot1);
    exit(1);
  }
  if (surface2->ct == NULL)
  {
    printf("ERROR: %s does not contain color table!\n",annot2);
    exit(1);
  }
  char *ctabname1 = surface1->ct->fname;
  char *ctabname2 = surface2->ct->fname;
  if (strcmp(ctabname1,ctabname2) != 0)
  {
    // just a warning, as the table could nonetheless be identical
    printf("Warning: annotation files based on different colortables:\n");
    printf("\tannot1: %s\n",ctabname1);
    printf("\tannot2: %s\n",ctabname2);
  }
  int ctab_nentries1 = surface1->ct->nentries;
  int ctab_nentries2 = surface2->ct->nentries;
  if (ctab_nentries1 != ctab_nentries2)
  {
    printf("ERROR: annotation files have unequal number of "
           "colortable entries:\n");
    printf("\tannot1: %d entries\n",ctab_nentries1);
    printf("\tannot2: %d entries\n",ctab_nentries2);
    exit(1);
  }

  // stupid sanity check:
  if (surface1->nvertices != surface2->nvertices)
  {
    printf("ERROR: (surface1->nvertices=%d != surface2->nvertices=%d\n)",
           surface1->nvertices,surface2->nvertices);
    exit(1);
  }

  // don't want to overrun our global array of label info.
  if (surface1->nvertices >= MAX_VNOS)
  {
    printf("ERROR: surface1->nvertices=%d >= MAX_VNOS=%d\n",
           surface1->nvertices,MAX_VNOS);
    exit(1); // need to increase MAX_VNOS
  }

  /*
   * if a file containing a list of labels to use was specified, 
   * then initialize list of labels to include in accuracy measurements
   */
  if (labelsfile)
  {
    char *cp, line[STRLEN], label[STRLEN];
    FILE* fp = fopen(labelsfile,"r");
    if (NULL == fp)
    {
      printf("ERROR: could not open file '%s'\n",labelsfile);
      exit(1);
    }
    while ((cp = fgetl(line, STRLEN, fp)) != NULL)
    {
      COLOR_TABLE *ct=surface1->ct;
      if (!sscanf(cp, "%s", label))    
      {
        printf("ERROR: could not scan # of lines from %s",labelsfile);
        exit(1);
      }
      addToIncludedLabelsList(ct, label);
    }
  }

  /*
     get annot info on the label 'unknown', so that we can use it to
     create a .annot file containing vertices where perfect overlap occurs
     which are labeled as unknown, making it easy to see where
     the failure in overlap occurs.
  */
  int unknown_index=0;
  unknown_annot=0;
  CTABfindName(surface1->ct, "unknown", &unknown_index);
  if (unknown_index == -1)
  {
    printf
      ("ERROR: could not retrieve index for label 'unknown'\n");
    exit(1);
  }// else printf("label 'unknown' has index %d\n",unknown_index);
  int err=CTABannotationAtIndex(surface1->ct,unknown_index,&unknown_annot);
  if (err != NO_ERROR)
  {
    printf
      ("ERROR: could not retrieve annotation for label 'unknown'\n");
    exit(1);
  }// else printf("label 'unknown' has annotation 0x%8.8X\n",unknown_annot);

  /*
   * initialize list of labels to exclude from accuracy measurements
   */
  COLOR_TABLE *ct=surface1->ct;
  addToExcludedLabelsList(ct, "unknown"); //Desikan atlas
  addToExcludedLabelsList(ct, "corpuscallosum"); //Desikan atlas
  addToExcludedLabelsList(ct, "Unknown"); //Christophe atlas
  addToExcludedLabelsList(ct, "Corpus_callosum"); //Christophe atlas
  addToExcludedLabelsList(ct, "Medial_wall"); //Christophe atlas

  /*
   * if we are performing a label comparison, then get the annotation info
   * for our'userLabel', from our local mini colortable
   */
  int label_index=0;
  int label_annot=0;
  if (label1 && label2)
  {
    CTABfindName((CT*)&miniColorTable, userLabel.name, &label_index);
    if (label_index == -1)
    {
      printf
        ("ERROR: could not retrieve index for label '%s'\n",label1name);
      exit(1);
    }// else printf("label '%s' has index %d\n",userLabel.name,label_index);
    if (label_index >= MAX_LABELS)
    {
      printf
        ("ERROR: label_index=%d >= MAX_LABELS=%d\n",label_index,MAX_LABELS);
      exit(1);
    }
    int err=CTABannotationAtIndex
      ((CT*)&miniColorTable,label_index,&label_annot);
    if (err != NO_ERROR)
    {
      printf
        ("ERROR: could not retrieve annotation for label '%s'\n",
         userLabel.name);
      exit(1);
    }// else printf("label '%s' has annotation 0x%8.8X\n",
    //userLabel.name,label_annot);
  }

  /*
   * if performing a label comparison, we need to give the surface vertices
   * cooresponding to the labled vertices an annotation; also check that the
   * x,y,z is valid for this surface label
   */
  if (label1 && label2)
  {
    int vno;
    for (vno=0; vno < surface1->nvertices; vno++)
    {
      VERTEX *v1 = &surface1->vertices[vno];
      VERTEX *v2 = &surface2->vertices[vno];
      VERTEX *vo = &overlapSurf->vertices[vno];
      v1->annotation = unknown_annot; // default annotation
      v2->annotation = unknown_annot;
      vo->annotation = unknown_annot;
      int n;
      for (n=0; n < label1->n_points; n++)
      {
        if (label1->lv[n].vno == vno)
        {
          if (use_label1_xyz)
          {
            // replace surface x,y,z coordinates with those found in label file
            MRISsetXYZ(surface1,vno,
              label1->lv[n].x,
              label1->lv[n].y,
              label1->lv[n].z);
          }
          if (check_label1_xyz)
          {
            // this vertex is a label vertex, so check its coordinates
            if (label1->lv[n].x != v1->x)
            {
              printf("ERROR: label1 vno=%d has x=%f, while surface x=%f\n",
                     vno, label1->lv[n].x, v1->x);
            }
            if (label1->lv[n].y != v1->y)
            {
              printf("ERROR: label1 vno=%d has y=%f, while surface y=%f\n",
                     vno, label1->lv[n].y, v1->y);
            }
            if (label1->lv[n].z != v1->z)
            {
              printf("ERROR: label1 vno=%d has z=%f, while surface z=%f\n",
                     vno, label1->lv[n].z, v1->z);
            }
          }
          // looks good, so give it an annotation
          v1->annotation = label_annot;
          vo->annotation = label_annot; // the overlap (debugging) surface
        }
      }

      // repeat for label2
      for (n=0; n < label2->n_points; n++)
      {
        if (label2->lv[n].vno == vno)
        {
          if (use_label2_xyz)
          {
            // replace surface x,y,z coordinates with those found in label file
            MRISsetXYZ(surface2,vno,
              label2->lv[n].x,
              label2->lv[n].y,
              label2->lv[n].z);
          }
          if (check_label2_xyz)
          {
            // this vertex is a label vertex, so check its coordinates
            if (label2->lv[n].x != v2->x)
            {
              printf("ERROR: label2 vno=%d has x=%f, while surface x=%f\n",
                     vno, label2->lv[n].x, v2->x);
            }
            if (label2->lv[n].y != v2->y)
            {
              printf("ERROR: label2 vno=%d has y=%f, while surface y=%f\n",
                     vno, label2->lv[n].y, v2->y);
            }
            if (label2->lv[n].z != v2->z)
            {
              printf("ERROR: label2 vno=%d has z=%f, while surface z=%f\n",
                     vno, label2->lv[n].z, v2->z);
            }
          }
          // looks good, so give it an annotation
          v2->annotation = label_annot;
          vo->annotation = label_annot; // the overlap (debugging) surface
        }
      }
    }
    if (debug_labels)
    {
      sprintf(tmpstr,"%s/%s/label/%s.label1.annot",SUBJECTS_DIR,subject,hemi);
      printf("Writing %s...",tmpstr);
      MRISwriteAnnotation(surface1,tmpstr);
      printf("done.\n");
      sprintf(tmpstr,"%s/%s/label/%s.label2.annot",SUBJECTS_DIR,subject,hemi);
      printf("Writing %s...",tmpstr);
      MRISwriteAnnotation(surface2,tmpstr);
      printf("done.\n");
    }
  }


  /* 
   * ------------------------------ 
   * Do vertex-by-vertex comparison 
   * ------------------------------
   */
  printf("Checking %d surface vertices...\n",surface1->nvertices);
  fflush(stdout);
  fflush(stderr);
  int mismatchCount=0;
  int dice_overlap=0; // for Dice calc
  int dice_surf1=0; // for Dice calc
  int dice_surf2=0; // for Dice calc
  int dice_union=0;
  int n;
  for (n = 0; n < surface1->nvertices; n++)
  {
    /*
     * surface 1: get the colortable index.for this vertex's annotation
     */
    int colorTabIndex1 = 0;
    VERTEX_TOPOLOGY const * const v1t = &surface1->vertices_topology[n];
    VERTEX                * const v1  = &surface1->vertices         [n];
    if ( v1->annotation )
    {
      //printf("v1->annotation=0x%8.8.8X\n",v1->annotation);
      CTABfindAnnotation (surface1->ct, v1->annotation, &colorTabIndex1);
    }
    else
    {
      // for some reason, this vertex has a zero annotation field,
      // so mark it as 'unknown'
      v1->annotation = unknown_annot;
      colorTabIndex1 = unknown_index;
    }

    /*
     * surface 2: get the colortable index.for this vertex's annotation
     */
    int colorTabIndex2 = 0;
    VERTEX_TOPOLOGY const * const v2t = &surface2->vertices_topology[n];
    VERTEX                * const v2  = &surface2->vertices         [n];
    if ( v2->annotation )
    {
      CTABfindAnnotation (surface2->ct, v2->annotation, &colorTabIndex2);
    }
    else
    {
      // for some reason, this vertex has a zero annotation field,
      // so mark it as 'unknown'
      v2->annotation = unknown_annot;
      colorTabIndex2 = unknown_index;
    }

    /*
     * compare: gather Dice info (skipping overlapping excluded labels)
     */
    //printf("colorTabIndex1=%d, colorTabIndex2=%d\n",
    //colorTabIndex1, colorTabIndex2);
    if (( isExcludedLabel(colorTabIndex1)) &&
        ( isExcludedLabel(colorTabIndex2)))
    {
      // skip overlapping excluded labels
    }
    else if ((NULL != labelsfile) &&
             (( ! isIncludedLabel(colorTabIndex1)) ||
             ( ! isIncludedLabel(colorTabIndex2))) )
    {
      // if a label-list file was specified on command-line, 
      // then skip labels not found in that list
    }
    else
    {
      if (colorTabIndex1 != colorTabIndex2)
      {
        //printf("colorTabIndex1=%d != colorTabIndex2=%d\n",
        //     colorTabIndex1, colorTabIndex2);
        mismatchCount++; // hey!  a mismatch!  used later for Dice calc
      }
      else if (colorTabIndex1 == colorTabIndex2) 
      {
        dice_overlap++; // intersection
      }
      dice_union++;
      if (colorTabIndex1) dice_surf1++;
      if (colorTabIndex2) dice_surf2++;
    }

    /*
     * while we are looping in this Dice calc, we need to gather some
     * info for the other way of assessing accuracy (the mean-min-distances): 
     * determine if this vertex is a label boundary (outline) by checking
     * if any of its neighbors have annotations that are different.
     * if so, indicate by setting the 'border' field (we are stealing this
     * field for our own use).
     */
    v1->border=0;
    int neighbor_vno;
    for (neighbor_vno = 0; neighbor_vno < v1t->vnum; neighbor_vno++ )
    {
      VERTEX *v1Neighbor = &surface1->vertices[v1t->v[neighbor_vno]];
      if (v1Neighbor->annotation != v1->annotation)
      {
        v1->border=1; // this neighbor is a foreigner! (a different label)
        //printf("border, surf1: vno=%d\n",n);
        break;
      }
    }
    v2->border=0;
    for (neighbor_vno = 0; neighbor_vno < v2t->vnum; neighbor_vno++ )
    {
      VERTEX *v2Neighbor = &surface2->vertices[v2t->v[neighbor_vno]];
      if (v2Neighbor->annotation != v2->annotation)
      {
        v2->border=1;
        //printf("border, surf2: vno=%d\n",n);
        break;
      }
    }
    /*
     * while we are looping in this Dice calc, we also want to build a 
     * table of vertex numbers composing each label, for use during the
     * function calcMeanMinLabelDistances()
     */
    if (colorTabIndex1 >= MAX_LABELS)
    {
      if (colorTabIndex1 >= MAX_SKIPPED_LABELS)
      {
        printf("ERROR: Colortable1 index (%d) exceeds "
               "MAX_SKIPPED_LABELS (%d)\n",
               colorTabIndex1,MAX_SKIPPED_LABELS);
        exit(1);
      }
      // mark this label as skipped...
      skippedLabels[colorTabIndex1]=1;
    }
    else
    {
      if (surf1BoundaryLabels[colorTabIndex1].annotation == 0)
      {
        // store-away the annotation identifier
        surf1BoundaryLabels[colorTabIndex1].annotation = v1->annotation;
        //if (v1->annotation) printf(
        //"surf1BoundaryLabels[%d].annotation = 0x%8.8X\n",
        //colorTabIndex1,v1->annotation);
      }
      else if (v1->annotation != 0)
      {
        // sanity check of our data structures
        if (surf1BoundaryLabels[colorTabIndex1].annotation != v1->annotation)
        {
          printf("ERROR: "
                 "surf1BoundaryLabels[%d].annotation=0x%8.8X != "
                 "v1->annotation=0x%8.8X\n",
                 colorTabIndex1,
                 surf1BoundaryLabels[colorTabIndex1].annotation, 
                 v1->annotation);
          exit(1);
        }
      }
      // flag this vertex number as occupied, if it is a boundary vertex
      if (v1->border) surf1BoundaryLabels[colorTabIndex1].vnos[n] = 1;
    }

    // repeat for the other surface
    if (colorTabIndex2 >= MAX_LABELS)
    {
      if (colorTabIndex2 >= MAX_SKIPPED_LABELS)
      {
        printf("ERROR: Colortable2 index (%d) exceeds "
               "MAX_SKIPPED_LABELS (%d)\n",
               colorTabIndex2,MAX_SKIPPED_LABELS);
        exit(1);
      }
      // mark this label as skipped...
      skippedLabels[colorTabIndex2]=1;
    }
    else
    {
      if (surf2BoundaryLabels[colorTabIndex2].annotation == 0)
      {
        // store-away the annotation identifier
        surf2BoundaryLabels[colorTabIndex2].annotation = v2->annotation;
        //if (v2->annotation) printf(
        //"surf2BoundaryLabels[%d].annotation = 0x%8.8X\n",
        //colorTabIndex2,v2->annotation);
      }
      else if (v2->annotation != 0)
      {
        // sanity check of our data structures
        if (surf2BoundaryLabels[colorTabIndex2].annotation != v2->annotation)
        {
          printf("ERROR:"
                 "surf2BoundaryLabels[%d].annotation=0x%8.8X != "
                 "v2->annotation=0x%8.8X\n",
                 colorTabIndex2,
                 surf2BoundaryLabels[colorTabIndex2].annotation, 
                 v2->annotation);
          exit(1);
        }
      }
      // flag this vertex number as occupied, if it is a boundary vertex
      if (v2->border) surf2BoundaryLabels[colorTabIndex2].vnos[n] = 1;
    }

    /*
     * create an overlap file, where vertices where labels are equal
     * are relabelled as unknown, thus leaving unequal labels labelled
     * as-is, allowing visualization of overlap failure.
     */
    if (colorTabIndex1 == colorTabIndex2)
    {
      VERTEX *v = &overlapSurf->vertices[n];
      v->annotation = unknown_annot;
    }
  }

  /*
   * print Dice results, and write-out the ?h.overlap.annot debug file
   */

  if (debug_overlap)
  {
    sprintf(tmpstr,"%s/%s/label/%s.overlap.annot",
            SUBJECTS_DIR,subject,hemi);
    printf("Writing %s...",tmpstr);
    MRISwriteAnnotation(overlapSurf,tmpstr);
    printf("done.\n");
  }

  printf("Found %d overlaps (%d mismatches) out of %d checked vertices\n",
         dice_overlap,mismatchCount,dice_union);
  printf("\n"
         "----------------------\n"
         "Overall Dice = %1.4f \n"
         "----------------------\n",
         dice_overlap*2.0/(float)(dice_surf1 + dice_surf2));

  if (logfile)
    fprintf(logfile, "%f\n", dice_overlap*2.0/(float)(dice_surf1 + dice_surf2)) ;
  /*
   * Calc and print mean-distance results
   */
  calcMeanMinLabelDistances();

  if (logfile)
    fclose(logfile) ;
  exit(0);

}  /*  end main()  */


static int *excludedLabelsList=NULL;
static int numExcludedLabels;
static void addToExcludedLabelsList(COLOR_TABLE *ct, const char *labelToExclude)
{
  // first-time setup
  if (excludedLabelsList == NULL)
  {
    excludedLabelsList=(int *)malloc(10000*sizeof(int));
    numExcludedLabels=0;
  }
  if (excludedLabelsList == NULL)
  {
    printf("ERROR: failed to malloc memory for excludedLabels list\n");
    exit(1);
  }

  int index=0;
  CTABfindName(ct, labelToExclude, &index);
  if (index == -1)
  {
    //printf("INFO: could not retrieve index for label '%s'\n",labelToExclude);
  }
  else
  {
    //printf("Excluding label '%s' from measurements\n",labelToExclude);
    excludedLabelsList[numExcludedLabels] = index; // add to list
    if (++numExcludedLabels >= 10000)
    {
      printf("ERROR: exceeded max num excluded labels\n");
      exit(1);
    }
  }
}
static int isExcludedLabel(int colortabIndex)
{
  if (excludedLabelsList)
  {
    int i;
    for (i=0; i < numExcludedLabels; i++)
    {
      if (excludedLabelsList[i] == colortabIndex) return 1; // excluded label
    }
  }

  return 0;
}


static int *includedLabelsList=NULL;
static int numIncludedLabels;
static void addToIncludedLabelsList(COLOR_TABLE *ct, const char *labelToInclude)
{
  // first-time setup
  if (includedLabelsList == NULL)
  {
    includedLabelsList=(int *)malloc(10000*sizeof(int));
    numIncludedLabels=0;
  }
  if (includedLabelsList == NULL)
  {
    printf("ERROR: failed to malloc memory for includedLabels list\n");
    exit(1);
  }

  int index=0;
  CTABfindName(ct, labelToInclude, &index);
  if (index == -1)
  {
    //printf("INFO: could not retrieve index for label '%s'\n",labelToInclude);
  }
  else
  {
    //printf("Including label '%s' in measurements\n",labelToInclude);
    includedLabelsList[numIncludedLabels] = index; // add to list
    if (++numIncludedLabels >= 10000)
    {
      printf("ERROR: exceeded max num included labels\n");
      exit(1);
    }
  }
}
static int isIncludedLabel(int colortabIndex)
{
  if (includedLabelsList)
  {
    int i;
    for (i=0; i < numIncludedLabels; i++)
    {
      if (includedLabelsList[i] == colortabIndex) return 1; // included label
    }
  }

  return 0;
}


/*
 * A superior accuracy measure to the Dice coeffient is the following:
 *
 *  for every point on the boundary of label 1
 *    compute min distance to boundary of label 2.
 *  end
 *
 *  average these to get avg_label1_to_label2
 *
 *  for every point on the boundary of label 2
 *    compute min distance to boundary of label 1.
 *  end
 *
 *  average these to get avg_label2_to_label1
 * 
 *  avg = (avg_label2_to_label1 + avg_label1_to_label2)/2
 *
 * This doesn't preference big structures over small, and also yields a 
 * measure in mm that is much more interpretable. Also, the vertex-based way 
 * of computing DICE doesn't account for variable size triangles (it should 
 * probably be computed over faces not vertices, and scale by area).
 *
 *
 * Static global inputs: 
 *   MRIS* surface1,surface2
 *   LABEL_INFO surf1BoundaryLabels, surf2BoundaryLabels
 *
 * Output to stdout:
 *   mean min distance of each label
 *   overal mean min distance
 */
static void calcMeanMinLabelDistances(void)
{
  if (debug_boundaries)
  {
    // write-out annotation files where just the boundaries are labelled, 
    // useful for debug (to check if indeed the 'border' vertices are borders)
    int n;
    for (n = 0; n < surface1->nvertices; n++)
    {
      VERTEX *v1 = &surface1->vertices[n];
      if ( ! v1->border)
      {
        // not a boundary vertex, so make it black (unknown), so that only
        // boundary vertices are displayed
        v1->annotation = unknown_annot;
      }
    }
    sprintf(tmpstr,"%s/%s/label/%s.boundary1.annot",
            SUBJECTS_DIR,subject,hemi);
    printf("Writing %s...",tmpstr);
    MRISwriteAnnotation(surface1,tmpstr);
    printf("done.\n");

    // surface 2
    for (n = 0; n < surface2->nvertices; n++)
    {
      VERTEX *v2 = &surface2->vertices[n];
      if ( ! v2->border)
      {
        // not a boundary vertex, so make it black (unknown), so that only
        // boundary vertices are displayed
        v2->annotation = unknown_annot;
      }
    }
    sprintf(tmpstr,"%s/%s/label/%s.boundary2.annot",
            SUBJECTS_DIR,subject,hemi);
    printf("Writing %s...",tmpstr);
    MRISwriteAnnotation(surface2,tmpstr);
    printf("done.\n");
  }


  /*
   * for every point on the boundary of label on surface 1,
   * compute min distance to boundary of same label on surface 2
   */
  int cti; // colortable indexl (ie, label index)
  for (cti=0; cti < MAX_LABELS; cti++)
  {
    // skip excluded labels
    if ( isExcludedLabel(cti) ) continue;

    // if a label-list file was specified on command-line, 
    // then skip labels not found in that list
    if ((NULL != labelsfile) && ( ! isIncludedLabel(cti) )) continue;

    int boundaryVertexCount=0;
    // for each label...(valid label if nonzero annot id)
    if (surf1BoundaryLabels[cti].annotation != 0 &&
        surf1BoundaryLabels[cti].annotation != 0xFFFFFFFF)
    {
      int vno;
      // for each boundary vertex on surface1...
      for (vno=0; vno < surface1->nvertices; vno++)
      {
        if (surf1BoundaryLabels[cti].vnos[vno]) // found a boundary vertex
        { 
          VERTEX *v1 = &surface1->vertices[vno];
          // compute min distance to boundary of same label on surface 2
          float minDist=100000000;
          int vno2;
          for (vno2=0; vno2 < surface2->nvertices; vno2++)
          {
            if (surf2BoundaryLabels[cti].vnos[vno2])
            {// found boundary vertex
              VERTEX *v2 = &surface2->vertices[vno2];
#if 1
              // note: this is the Euclidian distance! not geodesic
              float dx = v2->x - v1->x;
              float dy = v2->y - v1->y;
              float dz = v2->z - v1->z;
              float dist = sqrt((dx*dx)+(dy*dy)+(dz*dz));
#endif
              //printf("dist=%f\n",dist);
              if (dist < minDist) minDist=dist;
              //if (dist==0) printf("cti=%d,vno=%d,vno2=%d\n",cti,vno,vno2);

              // sanity check: make sure annot id is what is expected,
              // (checking that memory hasn't been trampled in some manner)
              if (v1->annotation != surf2BoundaryLabels[cti].annotation)
              {
                printf("ERROR: v1->annotation=0x%8.8X != "
                       "surf2BoundaryLabels[%d].annotation=0x%8.8X\n",
                       v1->annotation,
                       cti,
                       surf2BoundaryLabels[cti].annotation);
                exit(1);
              }
              if (v2->annotation != surf1BoundaryLabels[cti].annotation)
              {
                printf("ERROR: v2->annotation=0x%8.8X != "
                       "surf1BoundaryLabels[%d].annotation=0x%8.8X\n",
                       v2->annotation,
                       cti,
                       surf1BoundaryLabels[cti].annotation);
                exit(1);
              }
            }
            // keep looping on finding the lowest distance from this surface1
            // label vertex to surface2's same label border vertices 
          }
          // sanity check: we should have found a minimum distance
          if (minDist==100000000) 
          {
	    minDist = abs(surface1->yhi-surface1->ylo) ;
            printf("ERROR: minDist==100000000 (%f) "
                   "(cti=%d,vno=%d,annotation1=0x%8.8X,annotation2=0x%8.8X)\n",
                   minDist,cti,vno,
                   surf1BoundaryLabels[cti].annotation,
                   surf2BoundaryLabels[cti].annotation);
//            exit(1);  in CD parcellation some folds may not occur
          }
          if (minDist < 0) // stupid sanity check
          {
            printf("ERROR: minDist<0 (cti=%d,vno=%d)\n",cti,vno);
            exit(1);
          }
          // save away the min distance info for this vertex
          surf1BoundaryLabels[cti].meanMinDist += minDist;
          boundaryVertexCount++;  // we'll need this to compute the average
        }
        // goto next border vertex for this label on surface1
      }
      // end of for each boundary vertices on surface1...
      // compute the average min distance of this label to its same label
      // on the other surface
      if (boundaryVertexCount==0) // sanity check
      {
        printf("ERROR: boundaryVertexCount==0\n");
        exit(1);
      }
      surf1BoundaryLabels[cti].meanMinDist /= boundaryVertexCount;
    }
  }

  /*
   * now do it all over again, this time from surface 2 to surface 1 (since
   * the mean minimum distance is not necessarily the same
   */
  for (cti=0; cti < MAX_LABELS; cti++)
  {
    // skip excluded labels
    if ( isExcludedLabel(cti) ) continue;

    // if a label-list file was specified on command-line, 
    // then skip labels not found in that list
    if ((NULL != labelsfile) && ( ! isIncludedLabel(cti) )) continue;

    int boundaryVertexCount=0;
    // for each label...(valid label if nonzero annot id)
    if (surf2BoundaryLabels[cti].annotation != 0 &&
        surf2BoundaryLabels[cti].annotation != 0xFFFFFFFF)
    {
      int vno;
      // for each boundary vertex on surface2...
      for (vno=0; vno < surface2->nvertices; vno++)
      {
        if (surf2BoundaryLabels[cti].vnos[vno]) // found a boundary vertex
        { 
          VERTEX *v2 = &surface2->vertices[vno];
          // compute min distance to boundary of same label on surface 1
          float minDist=100000000;
          int vno1;
          for (vno1=0; vno1 < surface1->nvertices; vno1++)
          {
            if (surf1BoundaryLabels[cti].vnos[vno1])
            {// found boundary vertex
              VERTEX *v1 = &surface1->vertices[vno1];
              float dx = v1->x - v2->x;
              float dy = v1->y - v2->y;
              float dz = v1->z - v2->z;
              float dist = sqrt((dx*dx)+(dy*dy)+(dz*dz));
              //printf("dist=%f\n",dist);
              if (dist < minDist) minDist=dist;
              //if (dist==0) printf("cti=%d,vno=%d,vno1=%d\n",cti,vno,vno1);

              // sanity check: make sure annot id is what is expected,
              // (checking that memory hasn't been trampled in some manner)
              if (v2->annotation != surf1BoundaryLabels[cti].annotation)
              {
                printf("ERROR: v2->annotation=0x%8.8X != "
                       "surf1BoundaryLabels[%d].annotation=0x%8.8X\n",
                       v2->annotation,
                       cti,
                       surf1BoundaryLabels[cti].annotation);
                exit(1);
              }
              if (v1->annotation != surf2BoundaryLabels[cti].annotation)
              {
                printf("ERROR: v1->annotation=0x%8.8X != "
                       "surf2BoundaryLabels[%d].annotation=0x%8.8X\n",
                       v1->annotation,
                       cti,
                       surf2BoundaryLabels[cti].annotation);
                exit(1);
              }
            }
          }
          // sanity check: we should have found a minimum distance
          if (minDist==100000000) 
          {
	    minDist = abs(surface1->yhi-surface1->ylo) ;
            printf("ERROR: minDist==100000000 (%f) "
                   "(cti=%d,vno=%d,annotation1=0x%8.8X,annotation2=0x%8.8X)\n",
                   minDist, cti,vno,
                   surf1BoundaryLabels[cti].annotation,
                   surf2BoundaryLabels[cti].annotation);
//            exit(1);  not an error - in CD parcellation some folds may not occur
          }
          // save away the min distance info for this vertex
          surf2BoundaryLabels[cti].meanMinDist += minDist;
          boundaryVertexCount++;  // we'll need this to compute the average
        }
      }
      // end of for each boundary vertex on surface2...
      // compute the average min distance of this label to its same label
      // on the other surface
      if (boundaryVertexCount==0) // sanity check
      {
        printf("ERROR: boundaryVertexCount==0\n");
        exit(1);
      }
      surf2BoundaryLabels[cti].meanMinDist /= boundaryVertexCount;
    }
  }

  /*
   * now, last but not least, average the two sets of mean min distances
   * for each label
   */
  float overallMeanMinDist=0;
  int ctiCount=0;
  printf("\n"
         "Average minimium distance error table:\n"
         "--------------------------------------------\n"
         "Label:                              Avg(mm):\n");
  for (cti=0; cti < MAX_LABELS; cti++)
  {
    if (isExcludedLabel(cti)) continue;  // skip excluded labels from calc

    // if a label-list file was specified on command-line, 
    // then skip labels not found in that list
    if ((NULL != labelsfile) && ( ! isIncludedLabel(cti) )) continue;

    // for each label...(valid label if nonzero annot id)
    if (surf1BoundaryLabels[cti].annotation != 0 &&
        surf1BoundaryLabels[cti].annotation != 0xFFFFFFFF)
    {
      // sanity check
      if (surf1BoundaryLabels[cti].annotation != 
          surf2BoundaryLabels[cti].annotation)
      {
	const char *n1, *n2 ;
	n1  = CTABgetAnnotationName(surface1->ct, surf1BoundaryLabels[cti].annotation) ;
	n2  = CTABgetAnnotationName(surface2->ct, surf2BoundaryLabels[cti].annotation) ;

        printf("ERROR: surf1BoundaryLabels[%d].annotation=0x%8.8X (%s) != "
               "surf2BoundaryLabels[%d].annotation=0x%8.8X (%s) \n",
               cti, surf1BoundaryLabels[cti].annotation, n1,
               cti, surf2BoundaryLabels[cti].annotation, n2);
//        exit(1);   Not an error - in CD parcellation some folds may not occur
      }
      // calc the average
      float avg = surf1BoundaryLabels[cti].meanMinDist;
      avg += surf2BoundaryLabels[cti].meanMinDist;
      avg /= 2;
      CTABcopyName(surface1->ct, cti, tmpstr, sizeof(tmpstr));
      padWhite(tmpstr,35);
      printf("%s %f\n", tmpstr, avg);
      // for overall average:
      overallMeanMinDist+=avg;
      ctiCount++;
    }
  }

  if (ctiCount==0) // sanity check
  {
    printf("ERROR: ctiCount==0\n");
    exit(1);
  }

  overallMeanMinDist /= ctiCount;
  printf("\n"
         "-------------------------------------------\n"
         "Overall mean min distance (mm) = %f\n"
         "-------------------------------------------\n",
         overallMeanMinDist);

  if (logfile)
    fprintf(logfile, "%f\n", overallMeanMinDist) ;

  // if any labels were skipped, print a warning...
  int cnt=0;
  for (cti=0; cti < MAX_SKIPPED_LABELS; cti++) if (skippedLabels[cti]) cnt++;

  if (cnt)
  {
    printf("\nWARNING!  The following labels were not included in the\n"
           "mean distance error table due to lack of memory space:\n"
           "index:\tname:\n");
    for (cti=0; cti < MAX_SKIPPED_LABELS; cti++) 
    {
      if (skippedLabels[cti])
      {
        CTABcopyName(surface1->ct, cti, tmpstr, sizeof(tmpstr));
        printf("%d\t%s\n",cti,tmpstr);
      }
    }    
  }
}



/* --------------------------------------------- */
static void usage(int exit_val)
{
  FILE *fout;
  char progname[] = "mris_compute_parc_overlap";

  fout = (exit_val ? stderr : stdout);

  fprintf
  (fout,
   "Compares two parcellated (annotated or labeled) surfaces\n"
   "and computes an overall Dice coefficient\n"
   "and mean minimum distances (mm).\n\n") ;
  fprintf
  (fout,
   "Usage:\n"
   "  %s --s subject --hemi hemi \\ \n"
   "    --annot1 annotfile --annot2 annotfile\n\n",
   progname);
  fprintf(fout, "Required:\n");
  fprintf(fout, "  --s subject              subject to check\n");
  fprintf(fout, "  --hemi hemi              hemisphere: rh or lh\n");
  fprintf(fout, "  and:\n");
  fprintf(fout, "    --annot1 annotfile     first .annot file\n");
  fprintf(fout, "    --annot2 annotfile     second .annot file\n");
  fprintf(fout, "  or:\n");
  fprintf(fout, "    --label1 labelfile     first .label file\n");
  fprintf(fout, "    --label2 labelfile     second .label file\n");
  fprintf(fout, "\nOptional:\n");
  fprintf(fout, "  --sd subj_dir            set SUBJECTS_DIR\n");
  fprintf(fout, "  --log filename           output the overall DICE and min dist to filename\n");
  fprintf(fout, "  --label-list file        file containing labels to \n"
                "                           check, one per line");
  fprintf(fout, "  --nocheck-label1-xyz     when loading label1 file, don't\n"
                "                           check x,y,z coords to surface\n"
                "                           default: check x,y,x\n");
  fprintf(fout, "  --nocheck-label2-xyz     ditto for label2\n");
  fprintf(fout, "  --nocheck-label-xyz      do not check label1 and label2\n");
  fprintf(fout, "  --use-label1-xyz         replace surface x,y,z coords\n"
                "                           with those in label1 file\n");
  fprintf(fout, "  --use-label2-xyz         ditto for label2\n");
  fprintf(fout, "  --use-label-xyz          use label1 and label2 coords\n");
  fprintf(fout, "  --debug-overlap          generate ?h.overlap.annot\n");
  fprintf(fout, "  --version                version info\n");
  fprintf(fout, "  --help                   this usage info\n");
  fprintf(fout, "\nExample 1:\n");
  fprintf
  (fout,
   "  %s --s bert --hemi lh \\"
   "\n    --annot1 aparc --annot2 aparc.ernie\n",
   progname);
  fprintf
  (fout,
   "\nIn this example, the annotation file named lh.aparc.annot, which is\n"
   "created by the utility mris_ca_label (executed during the -autorecon3\n"
   "stage of recon-all), is compared against the annotation file named \n"
   "lh.aparc.ernie.annot.  This second annotation file is created by the\n"
   "utility mri_surf2surf, which resamples one surface onto another. This\n"
   "resampling is necessary so that a vertex-by-vertex comparison is\n"
   "meaningful.  An example command-line is: \n"
   "  mri_surf2surf --srcsubject ernie --trgsubject bert --hemi lh \\ \n"
   "    --sval-annot $SUBJECTS_DIR/ernie/label/lh.aparc.annot \\ \n"
   "    --tval       $SUBJECTS_DIR/bert/label/lh.ernie.aparc.annot\n\n"
   "Note that the resampling output file, lh.ernie.annot, is deposited\n"
   "in the label directory of the subject (bert) supplied as the input\n"
   "to the %s utility.  Supply --help to\n"
   "mri_surf2surf for its usage information.\n\n"
   "There are two measures output by %s: an overall\n"
   "Dice coefficient, and a table of mean minimum distances between \n"
   "corresponding labels. A Dice value of 1 indicates perfect overlap.\n"
   "A mean minimum distance of 0 (mm) indicates perfect overlap.\n\n"
   "If --debug-overlap is specified, a file called ?h.overlap.annot is \n"
   "created in the subjects label directory, and is a copy of the annot2\n"
   "input, except wherever the labels are identical with annot1, that label\n"
   "is replaced with the label 'unknown', thus leaving any mismatches\n"
   "labeled as they were in the annot2 file.  If the Dice coefficient is \n"
   "less than 1, then this ile is a way to visualize the mismatches, \n"
   "which are typically around the label borders."
   "\n",progname,progname);

  fprintf(fout, "\nExample 2:\n");
  fprintf
  (fout,
   "  %s --s bert --hemi lh \\"
   "\n    --label1 precentral --label2 precentral\n",
   progname);
  fprintf
  (fout,
   "\nIn this example, two label files are specified.  Just those two\n"
   "labels are compared.\n");
  fprintf
  (fout,
   "\nWhen comparing the labels at two vertices, the following labels are\n"
   "excluded from measurements if both vertices have the same label:\n"
    "  unknown\n"
    "  corpuscallosum\n"
    "  Unknown\n"
    "  Corpus_callosum\n"
    );

  exit(exit_val);
}  /*  end usage()  */


/* --------------------------------------------- */
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage(1);

  nargc   = argc;
  pargv = argv;
  while (nargc > 0)
  {

    option = pargv[0];
    //printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug-overlap")) debug_overlap = 1;
    else if (!strcasecmp(option, "--nodebug-overlap")) debug_overlap = 0;
    else if (!strcasecmp(option, "--debug-boundaries")) debug_boundaries = 1;
    else if (!strcasecmp(option, "--nodebug-boundaries")) debug_boundaries = 0;
    else if (!strcasecmp(option, "--debug-labels")) debug_labels = 1;
    else if (!strcasecmp(option, "--nodebug-labels")) debug_labels = 0;
    else if (!strcasecmp(option, "--check-label1-xyz")) check_label1_xyz = 1;
    else if (!strcasecmp(option, "--nocheck-label1-xyz")) check_label1_xyz = 0;
    else if (!strcasecmp(option, "--check-label2-xyz")) check_label2_xyz = 1;
    else if (!strcasecmp(option, "--nocheck-label2-xyz")) check_label2_xyz = 0;
    else if (!strcasecmp(option, "--nocheck-xyz"))
    {
      check_label1_xyz = 0;
      check_label2_xyz = 0;
    }
    else if (!strcasecmp(option, "--use-label1-xyz")) use_label1_xyz = 1;
    else if (!strcasecmp(option, "--use-label2-xyz")) use_label2_xyz = 1;
    else if (!strcasecmp(option, "--use-label-xyz")) 
    {
      use_label1_xyz = 1;
      use_label2_xyz = 1;
    }
    else if (!strcmp(option, "--sd"))
    {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      setenv("SUBJECTS_DIR",SUBJECTS_DIR,1);
      nargsused = 1;
    }
    else if (!strcmp(option, "--log"))
    {
      if (nargc < 1) argnerr(option,1);
      logfile = fopen(pargv[0], "w");
      if (logfile == NULL)
	ErrorExit(ERROR_NOFILE, "could not open log file %s\n", pargv[0]) ;
      nargsused = 1;
    }
    else if (!strcmp(option, "--s"))
    {
      if (nargc < 1) argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--hemi"))
    {
      if (nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--annot1"))
    {
      if (nargc < 1) argnerr(option,1);
      annot1 = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--annot2"))
    {
      if (nargc < 1) argnerr(option,1);
      annot2 = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--label1"))
    {
      if (nargc < 1) argnerr(option,1);
      label1name = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--label2"))
    {
      if (nargc < 1) argnerr(option,1);
      label2name = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--label-list"))
    {
      if (nargc < 1) argnerr(option,1);
      labelsfile = pargv[0];
      nargsused = 1;
    }
    else
    {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}


/* --------------------------------------------- */
static void print_help(void)
{
  usage(1);
}


/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}


/* --------------------------------------------- */
static void check_options(void)
{
  if (subject == NULL)
  {
    printf("ERROR: must specify a subject\n");
    exit(1);
  }
  if (hemi == NULL)
  {
    printf("ERROR: must specify a hemisphere (rh or lh)\n");
    exit(1);
  }
  if ((annot1 == NULL) && (annot2 == NULL) &&
      (label1name == NULL) && (label2name == NULL))
  {
    printf("ERROR: must specify a pair of annotation or label files\n");
    exit(1);
  }
  if ((annot1 || annot2) && (label1name || label2name))
  {
    printf("ERROR: cannot specify both annotation and label file pairs\n");
    exit(1);
  }
  return;
}


/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\nSUBJECTS_DIR: %s\n",SUBJECTS_DIR);
  fprintf(fp,"subject:      %s\n",subject);
  fprintf(fp,"hemi:         %s\n",hemi);
  if (annot1) fprintf(fp,"annot1:       %s\n",annot1);
  if (annot2) fprintf(fp,"annot2:       %s\n",annot2);
  if (label1name) fprintf(fp,"label1:       %s\n",label1name);
  if (label2name) fprintf(fp,"label2:       %s\n",label2name);
  return;
}


/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}

/*---------------------------------------------------------------*/
static void padWhite(char* str, int maxLen)
{
  // pad string with spaces, upto maxLen string length
  while (strlen(str) < maxLen)
  {
    strcat(str," ");
  }
}

/*  EOF  */
