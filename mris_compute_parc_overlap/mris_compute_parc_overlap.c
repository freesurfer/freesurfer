/**
 * @file  mris_compute_parc_overlap.c
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
 *   --version            version info
 *   --help               this usage info
 *
 * Example:
 *   mris_compute_parc_overlap --s bert --hemi lh \
 *     --annot1 aparc --annot2 aparc.ernie
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/03/28 20:06:09 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2007,
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

#include <stdio.h>
#include <stdlib.h>

#include "mrisurf.h"
#include "annotation.h"
#include "version.h"
#include "hipsu.h"
#include "error.h"

#define MAX_VNOS 150000
typedef struct _labelInfo
{
  int annotation;      /* its annotation identifier */
  char vnos[MAX_VNOS]; /* if array element is 1, then that array index, which
                          is a vertex numbers, is a border vertex of this 
                          label */
  float meanMinDist;   /* average of the minimum distance of every boundary
                          vertex to the other surfaces boundary vertices for
                          this same label */
}
LABEL_INFO;
#define MAX_LABELS 100
static LABEL_INFO surf1BoundaryLabels[MAX_LABELS];
static LABEL_INFO surf2BoundaryLabels[MAX_LABELS];

static void calcMeanMinLabelDistances(void);
static void usage(int exit_val);
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void dump_options(FILE *fp);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static int  singledash(char *flag);

char *Progname;
static char vcid[] =
  "$Id: mris_compute_parc_overlap.c,v 1.5 2007/03/28 20:06:09 nicks Exp $";
static char *SUBJECTS_DIR = NULL;
static char *subject = NULL;
static char *hemi = NULL;
static char *annot1 = NULL;
static char *annot2 = NULL;
static int unknown_annot=0;
static MRIS *surface1 = NULL;
static MRIS *surface2 = NULL;
static int debug = 0;
static char tmpstr[2000];


int main(int argc, char *argv[])
{
  int nargs,err,n;
  char annotfile1[2000];
  char annotfile2[2000];
  MRIS *overlapSurf = NULL;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc--;
  argv++;
  if (argc == 0) usage(1);
  if (argc < 8) usage(1);

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
  surface2 = MRISread(tmpstr);
  overlapSurf = MRISread(tmpstr); // a debugging surface

  /* Load subject's native annotation onto surface1 */
  sprintf(annotfile1,"%s/%s/label/%s.%s.annot",
          SUBJECTS_DIR,subject,hemi,annot1);
  printf("\nLoading %s annotations from %s\n",hemi,annotfile1);
  fflush(stdout);
  err = MRISreadAnnotation(surface1, annotfile1);
  if (err)
  {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile1);
    exit(1);
  }

  /* Load the resampled annotation (from mri_surf2surf) onto surface2 */
  sprintf(annotfile2,"%s/%s/label/%s.%s.annot",
          SUBJECTS_DIR,subject,hemi,annot2);
  printf("\nLoading rh annotations from %s\n",annotfile2);
  fflush(stdout);
  err = MRISreadAnnotation(surface2, annotfile2);
  if (err)
  {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile2);
    exit(1);
  }
  MRISreadAnnotation(overlapSurf, annotfile2);

  /* make sure files contain embedded color look-up table. */
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
    printf("ERROR: (surface1->nvertices=%d != surface2->nvertices=%d)",
           surface1->nvertices,surface2->nvertices);
    exit(1);
  }

  // don't want to overrun our global array of label info.
  if (surface1->nvertices >= MAX_VNOS)
  {
    printf("ERROR: (surface1->nvertices=%d >= MAX_VNOS=%d)",
           surface1->nvertices,MAX_VNOS);
    exit(1); // need to increase MAX_VNOS
  }

  /*
     get annot info on the label 'unknown', so that we can use it to
     create an .annot file containing vertices where perfect overlap occurs
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
  }//else printf("label 'unknown' has index %d\n",unknown_index);
  err=CTABannotationAtIndex(surface1->ct,unknown_index,&unknown_annot);
  if (err != NO_ERROR)
  {
    printf
    ("ERROR: could not retrieve annotation for label 'unknown'\n");
    exit(1);
  }//else printf("label 'unknown' has annotation 0x%8.8X\n",unknown_annot);

  /* ------- Do vertex-by-vertex comparison ------- */
  printf("checking %d vertices...\n",surface1->nvertices);
  fflush(stdout);
  fflush(stderr);
  err=0;
  int surfannot_overlap = 0;
  int surfannot1 = 0;
  int surfannot2 = 0;
  memset(surf1BoundaryLabels,0,sizeof(surf1BoundaryLabels));
  memset(surf2BoundaryLabels,0,sizeof(surf2BoundaryLabels));
  for (n = 0; n < surface1->nvertices; n++)
  {
    /*
     * surface 1: get the color and then the index.
     */
    int colorTabIndex1 = 0;
    VERTEX *v1 = &surface1->vertices[n];
    if ( 0 != v1->annotation )
    {
      //printf("v1->annotation=0x%8.8.8X\n",v1->annotation);
      CTABfindAnnotation (surface1->ct, v1->annotation, &colorTabIndex1);
    }
    /*
     * surface 2: get the color and then the index.
     */
    int colorTabIndex2 = 0;
    VERTEX *v2 = &surface2->vertices[n];
    if ( 0 != v2->annotation )
    {
      CTABfindAnnotation (surface2->ct, v2->annotation, &colorTabIndex2);
    }
    /*
     * compare
     */
    //printf("colorTabIndex1=%d, colorTabIndex2=%d\n",
    //colorTabIndex1, colorTabIndex2);
    if (colorTabIndex1 != colorTabIndex2)
    {
      if (debug)
        printf("colorTabIndex1=%d != colorTabIndex2=%d\n",
               colorTabIndex1, colorTabIndex2);
      err++; // hey!  a mismatch!  used later for Dice calc
    }
    if (colorTabIndex1 != 0)
    {
      if (colorTabIndex1 == colorTabIndex2) surfannot_overlap++; // union
      surfannot1++;  // used later for Dice calc
      surfannot2++;  // used later for Dice calc
    }

    /*
     * while we are looping in this Dice calc, we need to gather some
     * info for the other way of assessing accuracy (the mean-min-distances): 
     * determine if this vertex is a label boundary (outline).
     * if so, indicate by setting the 'border' field (we are stealing this
     * field for our own use).
     */
    v1->border=0;
    v2->border=0;
    int neighbor_vno;
    for (neighbor_vno = 0; neighbor_vno < v1->vnum; neighbor_vno++ )
    {
      VERTEX *v1Neighbor = &surface1->vertices[v1->v[neighbor_vno]];
      if (v1Neighbor->annotation != v1->annotation)
      {
        v1->border=1; // this neighbor is a foreigner! (a different label)
        v2->border=1; // surface2 is a copy of surface1
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
      printf("ERROR: Colortable1 index (%d) exceeds MAX_LABELS (%d)\n",
             colorTabIndex1,MAX_LABELS);
      exit(1);
    }
    if (surf1BoundaryLabels[colorTabIndex1].annotation == 0 )
    {
      // store-away the annotation identifier
      surf1BoundaryLabels[colorTabIndex1].annotation = v1->annotation;
    }
    else
    {
      // sanity check of our data structures
      if (surf1BoundaryLabels[colorTabIndex1].annotation != v1->annotation)
      {
        printf("ERROR: "
               "surf1BoundaryLabels[colorTabIndex1].annotation=0x%8.8X != "
               "v1->annotation=0x%8.8X\n",
               surf1BoundaryLabels[colorTabIndex1].annotation, v1->annotation);
        exit(1);
      }
    }
    // flag this vertex number as occupied, if it is a boundary vertex
    if (v1->border) surf1BoundaryLabels[colorTabIndex1].vnos[n] = 1;

    // repeat for the other surface
    if (colorTabIndex2 >= MAX_LABELS)
    {
      printf("ERROR: Colortable2 index (%d) exceeds MAX_LABELS (%d)\n",
             colorTabIndex2,MAX_LABELS);
      exit(1);
    }
    if (surf2BoundaryLabels[colorTabIndex2].annotation == 0 )
    {
      // store-away the annotation identifier
      surf2BoundaryLabels[colorTabIndex2].annotation = v2->annotation;
    }
    else
    {
      // sanity check of our data structures
      if (surf2BoundaryLabels[colorTabIndex2].annotation != v2->annotation)
      {
        printf("ERROR:"
               "surf2BoundaryLabels[colorTabIndex2].annotation=0x%8.8X != "
               "v2->annotation=0x%8.8X\n",
               surf2BoundaryLabels[colorTabIndex2].annotation, v2->annotation);
        exit(1);
      }
    }
    // flag this vertex number as occupied, if it is a boundary vertex
    if (v2->border) surf2BoundaryLabels[colorTabIndex2].vnos[n] = 1;

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
   * print Dice results
   */

  printf("Found %d mismatches out of %d vertices\n",err,surface1->nvertices);

  sprintf(tmpstr,"%s/%s/label/%s.overlap.annot",
          SUBJECTS_DIR,subject,hemi);
  printf("Writing %s\n",tmpstr);
  MRISwriteAnnotation(overlapSurf,tmpstr);

  printf("Overall Dice = %1.4f \n",
         surfannot_overlap*2.0/(float)(surfannot1 + surfannot2));

  /*
   * Calc and print mean-distance results
   */
  calcMeanMinLabelDistances();


  exit(0);

}  /*  end main()  */


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

#if 0
  // write-out a label file where just the boundaries are labelled, useful
  // for debug (to check if indeed the 'border' vertices are borders)
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
  sprintf(tmpstr,"%s/%s/label/%s.boundary.annot",
          SUBJECTS_DIR,subject,hemi);
  printf("Writing %s\n",tmpstr);
  MRISwriteAnnotation(surface1,tmpstr);
#endif

  /*
   * for every point on the boundary of label on surface 1,
   * compute min distance to boundary of same label on surface 2
   */
  int cti; // colortable indexl (ie, label index)
  for (cti=0; cti < MAX_LABELS; cti++)
  {
    int boundaryVertexNum=0;
    // for each label...(valid label if nonzero annot id, and skip unknown)
    if ((surf1BoundaryLabels[cti].annotation != 0) &&
        (surf1BoundaryLabels[cti].annotation != unknown_annot))
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
              float dx = v2->x - v1->x;
              float dy = v2->y - v1->y;
              float dz = v2->z - v1->z;
              float dist = sqrt((dx*dx)+(dy*dy)+(dz*dz));
              //printf("dist=%f\n",dist);
              if (dist < minDist) minDist=dist;
              //if (dist==0) printf("cti=%d,vno=%d,vno2=%d\n",cti,vno,vno2);

              // sanity check: make sure annot id is what is expected,
              // (checking that memory hasn't been trampled in some manner)
              if (v1->annotation != surf2BoundaryLabels[cti].annotation)
              {
                printf("ERROR: v1->annotation=0x%8.8X != "
                       "surf2BoundaryLabels[cti].annotation=0x%8.8X\n",
                       v1->annotation,
                       surf2BoundaryLabels[cti].annotation);
                exit(1);
              }
              if (v2->annotation != surf1BoundaryLabels[cti].annotation)
              {
                printf("ERROR: v2->annotation=0x%8.8X != "
                       "surf1BoundaryLabels[cti].annotation=0x%8.8X\n",
                       v2->annotation,
                       surf1BoundaryLabels[cti].annotation);
                exit(1);
              }
            }
          }
          // sanity check: we should have found a minimum distance
          if (minDist==100000000) 
          {
            printf("ERROR: minDist==100000000 (cti=%d,vno=%d)\n",cti,vno);
            exit(1);
          }
          // save away the min distance info for this vertex
          surf1BoundaryLabels[cti].meanMinDist += minDist;
          boundaryVertexNum++;  // we'll need this to compute the average
        }
      }
      // end of for each boundary vertex on surface1...
      // compute the average min distance of this label to its same label
      // on the other surface
      if (boundaryVertexNum==0) // sanity check
      {
        printf("ERROR: boundaryVertexNum==0\n");
        exit(1);
      }
      surf1BoundaryLabels[cti].meanMinDist /= boundaryVertexNum;
    }
  }


  /*
   * now do it all over again, this time from surface 2 to surface 1 (since
   * the mean minimum distance is not necessarily the same
   */
  for (cti=0; cti < MAX_LABELS; cti++)
  {
    int boundaryVertexNum=0;
    // for each label...(valid label if nonzero annot id, and skip unknown)
    if ((surf2BoundaryLabels[cti].annotation != 0) &&
        (surf2BoundaryLabels[cti].annotation != unknown_annot))
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
                       "surf1BoundaryLabels[cti].annotation=0x%8.8X\n",
                       v2->annotation,
                       surf1BoundaryLabels[cti].annotation);
                exit(1);
              }
              if (v1->annotation != surf2BoundaryLabels[cti].annotation)
              {
                printf("ERROR: v1->annotation=0x%8.8X != "
                       "surf2BoundaryLabels[cti].annotation=0x%8.8X\n",
                       v1->annotation,
                       surf2BoundaryLabels[cti].annotation);
                exit(1);
              }
            }
          }
          // sanity check: we should have found a minimum distance
          if (minDist==100000000) 
          {
            printf("ERROR: minDist==100000000 (cti=%d,vno=%d)\n",cti,vno);
            exit(1);
          }
          // save away the min distance info for this vertex
          surf2BoundaryLabels[cti].meanMinDist += minDist;
          boundaryVertexNum++;  // we'll need this to compute the average
        }
      }
      // end of for each boundary vertex on surface2...
      // compute the average min distance of this label to its same label
      // on the other surface
      if (boundaryVertexNum==0) // sanity check
      {
        printf("ERROR: boundaryVertexNum==0\n");
        exit(1);
      }
      surf2BoundaryLabels[cti].meanMinDist /= boundaryVertexNum;
    }
  }

  /*
   * now, last but not least, average the two sets of mean min distances
   * for each label
   */
  float overallMeanMinDist=0;
  int ctiCount=0;
  printf("\nLabel name:\tLabel-to-label mean minimum distance (mm):\n");
  for (cti=0; cti < MAX_LABELS; cti++)
  {
    // for each label...(valid label if nonzero annot id, and skip unknown)
    if ((surf1BoundaryLabels[cti].annotation != 0) &&
        (surf1BoundaryLabels[cti].annotation != unknown_annot))
    {
      // sanity check
      if (surf1BoundaryLabels[cti].annotation != 
          surf2BoundaryLabels[cti].annotation)
      {
        printf("ERROR: surf1BoundaryLabels[%d].annotation=0x%8.8X != "
               "surf2BoundaryLabels[%d].annotation=0x%8.8X \n",
               cti, surf1BoundaryLabels[cti].annotation,
               cti, surf2BoundaryLabels[cti].annotation);
        exit(1);
      }
      // calc the average
      float avg = surf1BoundaryLabels[cti].meanMinDist;
      avg += surf2BoundaryLabels[cti].meanMinDist;
      avg /= 2;
      CTABcopyName(surface1->ct, cti, tmpstr, sizeof(tmpstr));
      printf("%s:\t%f\n", tmpstr, avg);
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
  printf("Overall mean min distance (mm) = %f\n",overallMeanMinDist);

}



/* --------------------------------------------- */
static void usage(int exit_val)
{
  FILE *fout;
  char progname[] = "mris_compute_parc_overlap";

  fout = (exit_val ? stderr : stdout);

  fprintf
  (fout,
   "Compares two parcellated (annotated) surfaces\n"
   "and computes an overall Dice coefficient\n"
   "and mean minimum distances (mm).\n\n") ;
  fprintf
  (fout,
   "Usage:\n"
   "  %s --s subject --hemi hemi \\ \n"
   "    --annot1 annotfile --annot2 annotfile\n\n",
   progname);
  fprintf(fout, "Required:\n");
  fprintf(fout, "  --s subject          subject to check\n");
  fprintf(fout, "  --hemi hemi          hemisphere: rh or lh\n");
  fprintf(fout, "  --annot1 annotfile   first .annot file\n");
  fprintf(fout, "  --annot2 annotfile   second .annot file\n");
  fprintf(fout, "\nOptional:\n");
  fprintf(fout, "  --sd subj_dir        set SUBJECTS_DIR\n");
  fprintf(fout, "  --version            version info\n");
  fprintf(fout, "  --help               this usage info\n");
  fprintf(fout, "\nExample:\n");
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
   "The output of %s is an overall Dice coefficient and mean min distance.\n"
   "A Dice value of 1 indicates perfect overlap.\n\n"
   "A mean minimum distance of 0 (mm) indicates perfect overlap.\n\n"
   "A file called ?h.overlap.annot is created in the subjects label\n"
   "directory, and is a copy of the annot2 input, except wherever the\n"
   "labels are identical with annot1, that label is replaced with the\n"
   "label 'unknown', thus leaving any mismatches labeled as they were\n"
   "in the annot2 file.  If the Dice coefficient is less than 1, then this\n"
   "file is a handy way to visualize the mismatches, which are typically\n"
   "around the label borders."
   "\n",progname,progname);

  exit(exit_val);
}  /*  end usage()  */


/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
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
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcmp(option, "--sd"))
    {
      if (nargc < 1) argnerr(option,1);
      SUBJECTS_DIR = pargv[0];
      setenv("SUBJECTS_DIR",SUBJECTS_DIR,1);
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
  if (annot1 == NULL)
  {
    printf("ERROR: must specify an annotation file\n");
    exit(1);
  }
  if (annot2 == NULL)
  {
    printf("ERROR: must specify the second annotation file\n");
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
  fprintf(fp,"annot1:       %s\n",annot1);
  fprintf(fp,"annot2:       %s\n",annot2);
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

/*  EOF  */
