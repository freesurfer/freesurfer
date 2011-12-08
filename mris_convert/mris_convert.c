/**
 * @file  mris_convert.c
 * @brief Format conversions of surface files and scalar overlay files
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2011/12/08 20:42:03 $
 *    $Revision: 1.40 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "macros.h"
#include "fio.h"
#include "version.h"
#include "matrix.h"
#include "transform.h"
#include "gifti_local.h"


//------------------------------------------------------------------------
static char vcid[] =
"$Id: mris_convert.c,v 1.40 2011/12/08 20:42:03 greve Exp $";

/*-------------------------------- CONSTANTS -----------------------------*/
// this mini colortable is used when .label file gets converted to gifti
static const COLOR_TABLE_ENTRY unknown = 
{"unknown", 0,0,0,0, 0,0,0,0};
static COLOR_TABLE_ENTRY userLabel = 
{ "user label name gets copied here                   ", 
  220,20,20,255, 0.8,0.08,0.08,1};
static const CTE *entries[2] = {&unknown, &userLabel};
static const COLOR_TABLE miniColorTable = 
{(CTE**)entries, 2, "miniColorTable", 2};

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int convertToWFile(char *in_fname, char *out_fname) ;
static int convertFromWFile(char *in_fname, char *out_fname) ;
static int writeAsciiCurvFile(MRI_SURFACE *mris, char *out_fname) ;

/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static int talairach_flag = 0 ;
static char *talxfmsubject = NULL;
static int patch_flag = 0 ;
static int read_orig_positions = 0 ;
static int w_file_dst_flag = 0 ;
static int w_file_src_flag = 0 ;
static int curv_file_flag = 0 ;
static char *curv_fname ;
static int func_file_flag = 0 ;
static char *func_fname ;
static int annot_file_flag = 0 ;
static char *annot_fname ;
static int gifti_da_num = -1;
static int label_file_flag = 0 ;
static char *label_fname ;
static char *label_name ;
static int labelstats_file_flag = 0 ;
static char *labelstats_fname ;
static int parcstats_file_flag = 0 ;
static char *parcstats_fname ;
static char *orig_surf_name = NULL ;
static double scale=0;
static int rescale=0;  // for rescaling group average surfaces
static int output_normals=0;
static int PrintXYZOnly = 0;
static MATRIX *XFM=NULL;
static int write_vertex_neighbors = 0;
static int combinesurfs_flag = 0;
int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname);

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  MRI_SURFACE  *mris ;
  char **av, *in_fname, *out_fname, fname[STRLEN], hemi[10],
    *cp, path[STRLEN], *dot, ext[STRLEN] ;
  int ac, nargs,nthvtx ;
  FILE *fp=NULL;
  char *in2_fname=NULL;
  MRI_SURFACE *mris2=NULL;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_convert.c,v 1.40 2011/12/08 20:42:03 greve Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
    //printf("argc:%d, argv[1]:%s, argv[2]:%s, argv[3]:%s\n",
    //     argc,argv[1],argv[2],argv[3]);
  }

  // confirm that all options were eaten (this catches the case where options
  // were included at the end of the command string)
  if (combinesurfs_flag) {
    if (argc != 4) usage_exit() ;
  } 
  else {
    if (argc != 3) usage_exit() ;
  }

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  if (combinesurfs_flag) {
    in2_fname = argv[2];
    out_fname = argv[3];
  }

  if (talxfmsubject && curv_file_flag) {
    printf("ERROR: cannot specify -t and -c\n");
    exit(1);
  }

  if (labelstats_file_flag && ! label_file_flag) {
    printf("ERROR: cannot specify --labelstats without --label\n");
    exit(1);
  }

  if (parcstats_file_flag && ! annot_file_flag) {
    printf("ERROR: cannot specify --parcstats without --annot\n");
    exit(1);
  }

  // check whether output is a .w file
  dot = strrchr(out_fname, '.') ;
  if (dot) {
    strcpy(ext, dot+1) ;
    if (!stricmp(ext, "W")) w_file_dst_flag = 1 ;
  }

  if (w_file_dst_flag) {
    convertToWFile(in_fname, out_fname) ; //???
    exit(0) ;
  }

  // check whether input is a .w file
  dot = strrchr(in_fname, '.') ;
  if (dot) {
    strcpy(ext, dot+1) ;
    if (!stricmp(ext, "W")) w_file_src_flag = 1 ;
  }

  if (w_file_src_flag) {
    convertFromWFile(in_fname, out_fname) ; //???
    exit(0) ;
  }

  if (patch_flag) {   /* read in orig surface before reading in patch */
    char name[100] ;

    FileNamePath(in_fname, path) ;
    FileNameOnly(in_fname, name) ;
    cp = strchr(name, '.') ;
    if (cp) {
      strncpy(hemi, cp-2, 2) ;
      hemi[2] = 0 ;
    } else
      strcpy(hemi, "lh") ;

    sprintf(fname, "%s/%s.orig", path, hemi) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    if (MRISreadPatch(mris, in_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, in_fname) ;
    if (read_orig_positions) {
      if (MRISreadVertexPositions(mris, orig_surf_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, orig_surf_name) ;
    }
  } else {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
    if (combinesurfs_flag) {
      mris2 = MRISread(in2_fname) ;
      if (!mris2)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, in2_fname) ;
    }
  }

  if (talxfmsubject) {
    XFM = DevolveXFM(talxfmsubject, NULL, NULL);
    if (XFM == NULL) exit(1);
    printf("Applying talairach transform\n");
    MatrixPrint(stdout,XFM);
    MRISmatrixMultiply(mris,XFM);
  }

  if (rescale) {
    if (mris->group_avg_surface_area == 0) {
      printf("ERROR: cannot rescale a non-group surface\n");
      exit(1);
    }
    scale = sqrt((double)mris->group_avg_surface_area/mris->total_area);
  }

  if (scale > 0) {
    printf("scale = %lf\n",scale);
    MRISscale(mris,scale);
    MRIScomputeMetricProperties(mris);
  }

  if (curv_file_flag) {
    int type ;

    if (MRISreadCurvatureFile(mris, curv_fname) != NO_ERROR) exit(1);
    type = MRISfileNameType(out_fname) ;
    if (type == MRIS_ASCII_FILE)
      writeAsciiCurvFile(mris, out_fname) ;
    else if (type == MRIS_GIFTI_FILE)
      MRISwriteGIFTI(mris, NIFTI_INTENT_SHAPE, out_fname, curv_fname);
    else
      MRISwriteCurvature(mris, out_fname) ;
  } 
  else if (annot_file_flag) {
    // first read the annotation/gifti label data...
    int type = MRISfileNameType(annot_fname);
    if (type == MRIS_ANNOT_FILE) {
      if (MRISreadAnnotation(mris, annot_fname) != NO_ERROR) exit(1);
    } else if (type == MRIS_GIFTI_FILE) {
      if (NULL == mrisReadGIFTIdanum(annot_fname, mris, gifti_da_num)) exit(1);
    } else {
      printf("ERROR: unknown file annot file type specified for --annot: "
             "%s\n",annot_fname);
      exit(1);
    }
    // read parcstats text file (pairs of parc labels and stat values) and
    // save value associated with that parc label into the vertex with that
    // parc (annot) label
    if (parcstats_file_flag) {
      FILE* fp;
      if ((fp = fopen(parcstats_fname, "r")) == NULL)
      {
        errno = 0;
        ErrorExit(ERROR_BADFILE, "ERROR: can't open file %s", parcstats_fname);
      }
      char line[STRLEN];
      while (fgets(line, STRLEN, fp) != NULL) {
        char label[STRLEN];
        float val;
        sscanf(line,"%s %f",label,&val);
        // get the annotation value for this label from the colortable
        int annot = CTABentryNameToAnnotation(label, mris->ct);
        int vno;
        int doprint=1;
        for (vno=0; vno < mris->nvertices; vno++) {
          if (annot == mris->vertices[vno].annotation) {
            mris->vertices[vno].curv = val;
            if (doprint) {
              printf("label: %s, val: %9.9f\n",label,val);
              doprint = 0;
            }
          }
        }
      }
      // now write the 'curv' data (the parc stats we just assigned) to file
      type = MRISfileNameType(out_fname) ;
      if (type == MRIS_ASCII_FILE) {
        writeAsciiCurvFile(mris, out_fname) ;
      }
      else if (type == MRIS_GIFTI_FILE) {
        MRISwriteGIFTI(mris, NIFTI_INTENT_SHAPE, out_fname, curv_fname);
      }
      else {
        MRISwriteCurvature(mris, out_fname) ;
      }
      exit(0);
    }

    // if fall through, then write annot file
    type = MRISfileNameType(out_fname);
    if (type == MRIS_ANNOT_FILE) {
      if (MRISwriteAnnotation(mris, out_fname) != NO_ERROR) exit(1);
    } else if (type == MRIS_GIFTI_FILE) {
      if (MRISwriteGIFTI(mris,NIFTI_INTENT_LABEL,out_fname,NULL) != NO_ERROR) {
        exit(1);
      }
    } else {
      printf("ERROR: unknown file annot file type specified for output: "
             "%s\n",out_fname);
      exit(1);
    }
  }
  else if (label_file_flag) {
    // first read the freesurfer .label file...
    LABEL* label = LabelRead(NULL, label_fname);
    if (NULL == label) {
      printf("ERROR: reading .label file specified for --label: "
             "%s\n",label_fname);
      exit(1);
    }
    // give this label its own colortable:
    mris->ct = (CT*)&miniColorTable;
    strcpy(userLabel.name,label_name);
    // try to find this label in the FreeSurferColorLUT, so we have a unique
    // color (annotation) for it (otherwise, just use default miniColorTable)
    COLOR_TABLE *ct0;
    char ctabfile[2000];
    sprintf(ctabfile,"%s/FreeSurferColorLUT.txt",getenv("FREESURFER_HOME"));
    ct0 = CTABreadASCII(ctabfile);
    if (ct0) {
      int cno;
      for (cno=0; cno < ct0->nentries; cno++) {
        if (ct0->entries[cno]) {
          if (ct0->entries[cno]->name) {
            if (0==strcmp(label_name,ct0->entries[cno]->name)) {
              // we found this label! so update local colortable with info
              memcpy(miniColorTable.entries[1],
                     ct0->entries[cno],
                     sizeof(COLOR_TABLE_ENTRY));
              break;
            }
          }
        }
      }
    }
    // assign annotation to each label vertex (and while we're at it, stats)
    int annotation = CTABrgb2Annotation(miniColorTable.entries[1]->ri,
                                        miniColorTable.entries[1]->gi,
                                        miniColorTable.entries[1]->bi);
    int lno;
    for (lno=0; lno < label->n_points; lno++) {
      int vno = label->lv[lno].vno;
      mris->vertices[vno].annotation = annotation;
      mris->vertices[vno].stat = label->lv[lno].stat; // in case --labelstats
    }

    // now write the annot file (either in .annot format, or gifti LabelTable)
    int type = MRISfileNameType(out_fname);
    if (type == MRIS_ANNOT_FILE) {
      if (MRISwriteAnnotation(mris, out_fname) != NO_ERROR) exit(1);
    } else if (type == MRIS_GIFTI_FILE) {
      if (MRISwriteGIFTI(mris,NIFTI_INTENT_LABEL,out_fname,NULL) != NO_ERROR) {
        exit(1);
      }
    } else {
      printf("ERROR: unknown file annot file type specified for output: "
             "%s\n",out_fname);
      exit(1);
    }

    mris->ct = NULL; // to avoid calling CTABfree (our table is static memory)

    // if --labelstats was given, then we want to write-out the stats values
    // found in the .label file to a file
    if (labelstats_file_flag) {
      int type = MRISfileNameType(labelstats_fname);
      if (type == MRIS_GIFTI_FILE) {
        if (MRISwriteGIFTI(mris,NIFTI_INTENT_UNIFORM,labelstats_fname,NULL) 
            != NO_ERROR) {
          exit(1);
        }
      } else {
        printf("ERROR: unknown file file type specified for --labelstats: "
               "%s\n",labelstats_fname);
        exit(1);
      }
    }
  }
  else if (func_file_flag) {
    MRI* mri = MRIread( func_fname );
    if (NULL == mri) {
      printf("ERROR: unable to to read %s\n",func_fname);
      exit(1);
    }
    MRIwrite(mri, out_fname);
    exit(0);
  }
  else if (mris->patch) {
    int type = MRISfileNameType(out_fname) ;
    if (type == MRIS_GIFTI_FILE) {
      MRISwrite(mris, out_fname);
    }
    else {
      MRISwritePatch(mris, out_fname) ;
    }
  }
  else if (output_normals)
    MRISwriteNormalsAscii(mris, out_fname) ;
  else if (write_vertex_neighbors)
    MRISwriteVertexNeighborsAscii(mris, out_fname) ;
  else if(PrintXYZOnly){
    printf("Printing only XYZ to ascii file\n");
    fp = fopen(out_fname,"w");
    for(nthvtx = 0; nthvtx < mris->nvertices; nthvtx++){
      fprintf(fp,"%9.4f %9.4f %9.4f\n",mris->vertices[nthvtx].x,
              mris->vertices[nthvtx].y,mris->vertices[nthvtx].z);
    }
    fclose(fp);
  }
  else if (combinesurfs_flag) {
    MRI_SURFACE *mris3 = MRISalloc(mris->nvertices+mris2->nvertices,
                                   mris->nfaces+mris2->nfaces);
    int vno,vno2,vno3;
    for (vno=0,vno3=0; vno < mris->nvertices; vno++, vno3++) {
      mris3->vertices[vno3].x = mris->vertices[vno].x;
      mris3->vertices[vno3].y = mris->vertices[vno].y;
      mris3->vertices[vno3].z = mris->vertices[vno].z;
    }
    for (vno2=0; vno2 < mris2->nvertices; vno2++, vno3++) {
      mris3->vertices[vno3].x = mris2->vertices[vno2].x;
      mris3->vertices[vno3].y = mris2->vertices[vno2].y;
      mris3->vertices[vno3].z = mris2->vertices[vno2].z;
    }
    int fno,fno2,fno3;
    for (fno=0,fno3=0; fno < mris->nfaces; fno++, fno3++) {
      mris3->faces[fno3].v[0] = mris->faces[fno].v[0];
      mris3->faces[fno3].v[1] = mris->faces[fno].v[1];
      mris3->faces[fno3].v[2] = mris->faces[fno].v[2];
    }
    int offset = mris->nvertices;
    for (fno2=0; fno2 < mris2->nfaces; fno2++, fno3++) {
      mris3->faces[fno3].v[0] = mris2->faces[fno2].v[0] + offset;
      mris3->faces[fno3].v[1] = mris2->faces[fno2].v[1] + offset;
      mris3->faces[fno3].v[2] = mris2->faces[fno2].v[2] + offset;
    }
    MRISwrite(mris3, out_fname);
    MRISfree(&mris2);
    MRISfree(&mris3);
  }
  else {
    if(MRISfileNameType(out_fname) == MRIS_VOLUME_FILE){
      printf("Saving surface xyz %s as a volume format\n",out_fname);
      MRI *vol = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);
      MRIScopyMRI(mris,vol,0,"x");
      MRIScopyMRI(mris,vol,1,"y");
      MRIScopyMRI(mris,vol,2,"z");
      MRIwrite(vol,out_fname);
      MRIfree(&vol);
    } else {
      // default output: 
      printf("Saving %s as a surface\n",out_fname);
      MRISwrite(mris, out_fname) ;
    }
  }

  MRISfree(&mris) ;

  exit(0) ;

  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "-annot")) {
    annot_fname = argv[2] ;
    annot_file_flag = 1;
    nargs = 1 ;
  } else if (!stricmp(option, "-da_num")) {
    sscanf(argv[2],"%d",&gifti_da_num);
    nargs = 1 ;
  } else if (!stricmp(option, "-label")) {
    label_fname = argv[2] ;
    label_name = argv[3] ;
    label_file_flag = 1;
    nargs = 2 ;
  } else if (!stricmp(option, "-labelstats")) {
    labelstats_fname = argv[2] ;
    labelstats_file_flag = 1;
    nargs = 1 ;
  } else if (!stricmp(option, "-parcstats")) {
    parcstats_fname = argv[2] ;
    parcstats_file_flag = 1;
    nargs = 1 ;
  } else if (!stricmp(option, "-combinesurfs")) {
    combinesurfs_flag = 1;
  } else switch (toupper(*option)) {
  case 'A':
    PrintXYZOnly = 1;
    break ;
  case 'F':
    func_file_flag = 1 ;
    func_fname = argv[2] ;
    nargs = 1 ;
    break ;
  case 'C':
    curv_file_flag = 1 ;
    curv_fname = argv[2] ;
    nargs = 1 ;
    break ;
  case 'N':
    output_normals = 1;
    break ;
  case 'O':
    read_orig_positions = 1 ;
    orig_surf_name = argv[2] ;
    nargs = 1 ;
    break ;
  case 'S':
    sscanf(argv[2],"%lf",&scale);
    nargs = 1 ;
    break ;
  case 'R':
    rescale = 1;
    break ;
  case 'P':
    patch_flag = 1 ;
    nargs = 0 ;
    break ;
  case 'T':
    talairach_flag = 1 ;
    talxfmsubject = argv[2] ;
    nargs = 1 ;
    break ;
  case 'V':
    write_vertex_neighbors = 1;
    nargs = 0 ;
    break ;
  case '?':
  case 'U':
  case 'H':
    print_help() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "Usage: %s [options] <input file> <output file>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  printf(
    "\nThis program will convert MRI-surface data formats.\n") ;
  printf( "\nValid options are:\n") ;
  printf( "  -p                input is a patch, not a full surface\n") ;
  printf( "  -c <scalar file>  input is scalar curv overlay file (must still\n"
          "                    specify surface)\n") ;
  printf( "  -f <scalar file>  input is functional time-series or other\n"
          "                    multi-frame data (must specify surface)\n") ;
  printf( "  --annot <annot file> input is annotation or gifti label data\n") ;
  printf( "  --parcstats <infile>  infile is name of text file containing\n") ;
  printf( "                    label/val pairs, where label is an annot name\n") ;
  printf( "                    and val is a value associated with that label.\n") ;
  printf( "                    The output file will be a scalar file.\n") ;
  printf( "  --da_num <num>    if input is gifti, 'num' specifies which\n"
          "                    data array to use\n");
  printf( "  --label <infile> <label>  infile is .label file\n") ;
  printf( "                    label is name of this label\n") ;
  printf( "  --labelstats <outfile>  outfile is name of gifti file\n") ;
  printf( "                    to which label stats will be written\n") ;
  printf( "  -o origname       read orig positions\n") ;
  printf( "  -s scale          scale vertex xyz by scale\n") ;
  printf( "  -r                rescale vertex xyz so total area is\n"
          "                    same as group average\n") ;
  printf( "  -t subject        apply talairach xfm of subject to\n"
          "                    vertex xyz\n");
  printf( "  -n                output is an ascii file where vertex data\n") ;
  printf( "                    is the surface normal vector\n") ;
  printf( "  -v Writes out neighbors of a vertex in each row. The first\n");
  printf( "     column is the vertex number, the 2nd col is the number of neighbors,\n");
  printf( "     the remaining cols are the vertex numbers of the neighbors.  \n");
  printf( "     Note: there can be a different number of neighbors for each vertex.\n") ;
  printf( "  -a                Print only surface xyz to ascii file\n") ;
  printf( "  --combinesurfs <infile> <in2file> <outfile>\n") ;
  printf( "\n") ;
  printf( "These file formats are supported:\n") ;
  printf( "  ASCII:       .asc\n");
  printf( "  ICO:         .ico, .tri\n");
  printf( "  GEO:         .geo\n");
  printf( "  STL:         .stl\n");
  printf( "  VTK:         .vtk\n");
  printf( "  GIFTI:       .gii\n");
  printf( "  MGH surface-encoded 'volume': .mgh, .mgz\n");
  printf( "Freesurfer binary format assumed for all other extensions.\n") ;
  printf( "\n") ;
  printf( "EXAMPLES:\n") ;
  printf( "\n");
  printf( "Convert a surface file to ascii:\n");
  printf( "  mris_convert lh.white lh.white.asc\n") ;
  printf( "\n");
  printf( "Write vertex neighbors to ascii:\n");
  printf( "  mris_convert -v lh.white lh.white.neighbors.asc\n") ;
  printf( "\n");
  printf( "Convert a surface file to ascii (vertices are surface normals):\n");
  printf( "  mris_convert -n lh.white lh.white.normals.asc\n") ;
  printf( "\n");
  printf( "Apply talairach xfm to white surface, save as binary:\n");
  printf( "  mris_convert -t bert lh.white lh.white.tal\n") ;
  printf( "\n");
  printf( "Convert a scalar overlay file to ascii:\n");
  printf( "  mris_convert -c lh.thickness lh.white lh.thickness.asc\n") ;
  printf( "\n") ;
  printf( "Convert a .annot file to Gifti label file:\n");
  printf( "  mris_convert --annot lh.aparc.annot lh.white lh.aparc.gii\n") ;
  printf( "\n") ;
  printf( "Convert a Gifti label file to .annot:\n");
  printf( "  mris_convert --annot lh.aparc.gii lh.white.gii lh.aparc.annot\n");
  printf( "\n") ;
  printf( "Convert a Freesurfer .label file to Gifti label format:\n");
  printf( "  mris_convert --label lh.V1.label V1 lh.white lh.V1.label.gii\n") ;
  printf( "\n") ;
  printf( "Create a scalar overlay file where each parcellation region\n") ;
  printf( "contains a single value:\n") ;
  printf( "  mris_convert --annot lh.aparc.annot --parcstats lh.parcstats.txt\\ \n") ;
  printf( "               lh.white lh.parcstats\n") ;
  printf( "\n") ;
  printf( "See also mri_surf2surf\n") ;
  exit(1) ;
}

static void
print_version(void) {
  printf( "%s\n", vcid) ;
  exit(1) ;
}

static int
convertToWFile(char *in_fname, char *out_fname) {
  FILE   *infp, *outfp ;
  char   line[300], *cp ;
  int    vno, l = 0, num, ilat ;
  float  val ;

  fprintf(stderr, "writing w file %s...\n", out_fname) ;
  outfp = fopen(out_fname,"wb");
  if (outfp==NULL)
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,out_fname) ;

  infp = fopen(in_fname,"rb");
  if (infp==NULL)
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,in_fname) ;


  cp = fgetl(line, 299, infp) ;
  ilat = atoi(cp) ; /* not used at the moment */
  cp = fgetl(line, 299, infp) ;
  num = atoi(cp) ;
  fwrite2(0,outfp);
  fwrite3(num,outfp);

  while ((cp = fgetl(line, 299, infp)) != NULL) {
    l++ ;
    if (sscanf(cp, "%d %f", &vno, &val) != 2) {
      ErrorPrintf(ERROR_BADFILE,
                  "%s: could not scan parms from line %d: %s.\n",
                  Progname, l, cp) ;
      val = 0.0f ;   /* don't know what it is... */
    }
    fwrite3(vno,outfp);
    fwriteFloat(val, outfp) ;
  }
  fclose(infp) ;
  fclose(outfp) ;
  return(NO_ERROR) ;
}


static int
convertFromWFile(char *in_fname, char *out_fname) {
  FILE   *infp, *outfp ;
  int    vno, num, ilat, i ;
  float  val, lat ;

  fprintf(stderr, "writing ascii w file %s...\n", out_fname) ;
  outfp = fopen(out_fname,"wb");
  if (outfp==NULL)
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,out_fname) ;

  infp = fopen(in_fname,"rb");
  if (infp==NULL)
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,in_fname) ;

  fread2(&ilat,infp);
  lat = (float)ilat/10.0;
  fprintf(outfp, "%2.3f\n", lat) ;
  fread3(&num,infp);
  fprintf(outfp, "%d\n", num) ;
  for (i=0;i<num;i++) {
    fread3(&vno,infp);
    val = freadFloat(infp) ;
    fprintf(outfp, "%d  %f\n", vno, val) ;
  }
  fclose(outfp);
  fclose(infp) ;

  return(NO_ERROR) ;
}


static int
writeAsciiCurvFile(MRI_SURFACE *mris, char *out_fname) {
  FILE   *fp ;
  int    vno ;
  VERTEX *v ;

  fp = fopen(out_fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_BADFILE, "%s could not open output file %s.\n",
              Progname, out_fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    fprintf(fp, "%3.3d %2.5f %2.5f %2.5f %2.5f\n",
            vno, v->x, v->y, v->z, v->curv) ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

/*
  \fn int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname)
  \brief Writes out neighbors of a vertex in each row. The first
  column is the vertex number, the 2nd col is the number of neighbors,
  the remaining cols are the vertex numbers of the neighbors.
*/

int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname)
{
  int vno, nnbrs, nbrvno;
  FILE *fp;

  fp = fopen(out_fname,"w");
  if(fp == NULL){
    printf("ERROR: opening %s\n",out_fname);
    exit(1);
  }

  for(vno=0; vno < mris->nvertices; vno++){
    nnbrs = mris->vertices[vno].vnum;
    fprintf(fp,"%6d %2d   ",vno,nnbrs);
    for (nbrvno = 0; nbrvno < nnbrs; nbrvno++)
      fprintf(fp,"%6d ",mris->vertices[vno].v[nbrvno]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  return(0);
}

