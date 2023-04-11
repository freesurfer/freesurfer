#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/**
 * @brief Format conversions of surface files and scalar overlay files
 *
 */
/*
 * Original Author: Bruce Fischl
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
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

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
#include "gifti.h"
#include "mri_identify.h"
#include "fsenv.h"

#include "compilerdefs.h"

#include "MRISurfOverlay.h"

#define __COMBINESURFS_TAKE_INFILE2 0

//------------------------------------------------------------------------

/*-------------------------------- CONSTANTS -----------------------------*/
// this mini colortable is used when .label file gets converted to gifti
#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
#endif

static const COLOR_TABLE_ENTRY unknown = {"unknown", 0,0,0,0, 0,0,0,0};
static COLOR_TABLE_ENTRY userLabel =
    {
      "user label name gets copied here                   ",
      220,20,20,255, 0.8,0.08,0.08,1
    };
static const CTE *entries[2] = {&unknown, &userLabel};
static const COLOR_TABLE miniColorTable =
    {(CTE**)entries, 2, "miniColorTable", 2};

#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic pop
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic pop
#endif

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static MRIS *__splitGIFTI(const char *fgifti, const char *outdir, const char *fout=NULL);
static void __convertCurvatureFile(MRIS *mris, int noverlay, const char **foverlays, char *out_fname);
static void __convertLabelFile(MRIS *mris, const char *flabel, const char *label, const char *labelstats, char *out_fname);
static void __convertAnnotFile(MRIS *mris, const char *fannot, int giftiDaNum, const char *parcstats, char *out_fname);
static void __convertFuncFile(const char *ffunc, char *out_fname);
static void __convertNormals(MRIS *mris, char *out_fname);
static void __convertMRISPatch(MRIS *mris, char *out_fname);

static int  get_option(int argc, char *argv[]) ;
static void check_options(void);
static int isflag(char *flag);
static int nth_is_arg(int nargc, char **argv, int nth);
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int convertToWFile(char *in_fname, char *out_fname) ;
static int convertFromWFile(char *in_fname, char *out_fname) ;
//static int writeAsciiCurvFile(MRI_SURFACE *mris, char *out_fname) ;
static MRI* computeAngles(MRIS* surf) ;
static int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname);

/*-------------------------------- DATA ----------------------------*/

const char *Progname ;

static bool doTkrRASConvert = false;
static int center_surface=0 ;
static int talairach_flag = 0 ;
static char *talxfmsubject = NULL;
static int patch_flag = 0 ;
static int read_orig_positions = 0 ;
static int w_file_dst_flag = 0 ;
static int w_file_src_flag = 0 ;
static int curv_file_flag = 0 ;
//static char *curv_fname = NULL;
static const char **arr_fcurv = NULL;
static int  nfcurv = 0;
static int func_file_flag = 0 ;
static char *func_fname = NULL;
static int annot_file_flag = 0 ;
static char *annot_fname = NULL;
static int gifti_da_num = -1;
static int label_file_flag = 0 ;
static char *label_fname = NULL;
static char *label_name = NULL;
static int labelstats_file_flag = 0 ;
static char *labelstats_fname = NULL ;
static int parcstats_file_flag = 0 ;
static char *parcstats_fname = NULL;
static char *orig_surf_name = NULL ;
static double scale=0;
static int rescale=0;  // for rescaling group average surfaces
static int output_normals=0;
static int PrintXYZOnly = 0;
static MATRIX *XFM=NULL;
static int write_vertex_neighbors = 0;
static int combinesurfs_flag = 0;
static int mergegifti_flag = 0;
static int splitgifti_flag = 0;
static char *giftioutdir = NULL;
static int userealras_flag = 0;
static int usesurfras_flag = 0;
static MRI *VolGeomMRI=NULL;
static int RemoveVolGeom = 0;
static int cras_add = 0;
static int cras_subtract = 0;
static int ToScanner = 0;
static int ToTkr = 0;
static int ToSurfCoords = 0;
MRIS *SurfCoords = NULL;
int WriteArea = 0;
int nUpsample = 0,UpsampleSortType=0;
int DeleteCommands = 0;

static char *infile = NULL, *infile2 = NULL, *outfile = NULL; 


/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[])
{
  MRI_SURFACE  *mris ;
  char **av, *in_fname, *out_fname, fname[STRLEN], hemi[10],
       *cp, path[STRLEN], *dot, ext[STRLEN] ;
  int ac, nargs,nthvtx,n;
  FILE *fp=NULL;
  char *in2_fname=NULL;

  std::string cmdline = getAllInfo(argc, argv, "mris_convert");

  nargs = handleVersionOption(argc, argv, "mris_convert");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
    //printf("argc:%d, argv[1]:%s, argv[2]:%s, argv[3]:%s\n",
    //     argc,argv[1],argv[2],argv[3]);
  }

  if (argc == 1) // no input parameters provided
    usage_exit();

#if !__COMBINESURFS_TAKE_INFILE2
  // confirm that all options were eaten (this catches the case where options
  // were included at the end of the command string)
  if (combinesurfs_flag)
  {
    if (argc != 4)
    {
      printf("--combinesurfs option is specified.\n");
      printf("ERROR in-file in2-file out-file are required!!!\n");
      exit(1);
      //usage_exit() ;
    }
  }
  else
  {
    if (argc != 3)
    {
      printf("ERROR in-file out-file are required!!!\n");
      exit(1);
      //usage_exit() ;
    }
  }
#else
  if (argc != 3)
  {
    usage_exit() ;
  }
#endif

  in_fname = argv[1] ;
  out_fname = argv[2] ;

#if !__COMBINESURFS_TAKE_INFILE2
  if (combinesurfs_flag)
  {
    in2_fname = argv[2];
    out_fname = argv[3];
  }
#endif

  infile  = in_fname;
  infile2 = in2_fname;
  outfile = out_fname;

  check_options();

  // check whether output is a .w file
  dot = strrchr(out_fname, '.') ;
  if (dot)
  {
    strcpy(ext, dot+1) ;
    if (!stricmp(ext, "W"))
    {
      w_file_dst_flag = 1 ;
    }
  }

  if (w_file_dst_flag)
  {
    convertToWFile(in_fname, out_fname) ; //???
    exit(0) ;
  }

  // check whether input is a .w file
  dot = strrchr(in_fname, '.') ;
  if (dot)
  {
    strcpy(ext, dot+1) ;
    if (!stricmp(ext, "W"))
    {
      w_file_src_flag = 1 ;
    }
  }

  if (w_file_src_flag)
  {
    convertFromWFile(in_fname, out_fname) ; //???
    exit(0) ;
  }

  if (patch_flag)     /* read in orig surface before reading in patch */
  {
    char name[100] ;

    FileNamePath(in_fname, path) ;
    FileNameOnly(in_fname, name) ;
    cp = strchr(name, '.') ;
    if (cp)
    {
      strncpy(hemi, cp-2, 2) ;
      hemi[2] = 0 ;
    }
    else
    {
      strcpy(hemi, "lh") ;
    }

    int req = snprintf(fname, STRLEN, "%s/%s.orig", path, hemi) ;    
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    
    // why do we need to read .orig before reading the patch ???
    mris = MRISread(fname, doTkrRASConvert) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    if (MRISreadPatch(mris, in_fname, doTkrRASConvert) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, in_fname) ;
    if (read_orig_positions)
    {
      if (MRISreadVertexPositions(mris, orig_surf_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, orig_surf_name) ;
    }
  }
  else
  {
    if (splitgifti_flag)
    {
      mris = __splitGIFTI(in_fname, giftioutdir);
    }
    else
    {
      mris = MRISread(in_fname, doTkrRASConvert) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
      if (center_surface)
        MRIScenter(mris, mris) ;
    }
  }

  if(nUpsample > 0){
    printf("Upsampling %d times, SortType %d\n",nUpsample,UpsampleSortType);
    MRIS *tmpmris = MRISupsampleSplit(mris,nUpsample,UpsampleSortType);
    MRISfree(&mris);
    mris = tmpmris;
  }

  if(DeleteCommands)
  {
    fprintf(stderr,"Deleting %d commands from surface header\n",mris->ncmds);
    for(n=0; n<mris->ncmds; n++)
    {
      free(mris->cmdlines[n]);
      mris->cmdlines[n] = NULL;
    }
    mris->ncmds = 0;
  }

  MRISaddCommandLine(mris, cmdline);

  if(cras_add){
    printf("Adding scanner CRAS to surface xyz\n");
    MRISshiftCRAS(mris, 1);
  }

  if(cras_subtract){
    printf("Subtracting scanner CRAS from surface xyz\n");
    MRISshiftCRAS(mris, -1);
  }

  if(ToSurfCoords){
    printf("Converting to coords of --to-surf surface\n");
    int err = MRIScopyCoords(mris,SurfCoords);
    if(err) exit(1);
  }

  if(ToScanner || userealras_flag){
    if (mris->useRealRAS)
      printf("Surface XYZ coordinates are already in scanner space. No conversion needed.\n");
    else
    {
      if (mris->vg.valid)
      {
        printf("Converting surface XYZ coordinates from tkr to scanner space.\n");
        MRIStkr2Scanner(mris);
      }
      else
      {
        printf("Invalid volume geometry. Cannot convert surface XYZ coordinates from tkr to scanner space.\n"); 
        mris->useRealRAS = 0;
      }
    }
  }

  if(ToTkr || usesurfras_flag){
    if (!mris->useRealRAS)
      printf("Surface XYZ coordinates are already in tkr space. No conversion needed.\n");
    else
    {
      if (mris->vg.valid)
      {
        printf("Converting surface XYZ coordinates from scanner to tkr space.\n");
        MRISscanner2Tkr(mris);
      }
      else
      {
        printf("Invalid volume geometry. Cannot convert surface XYZ coordinates from scanner to tkr space.\n");
        mris->useRealRAS = 1;
      }
    }
  }

  if (talxfmsubject)
  {
    XFM = DevolveXFM(talxfmsubject, NULL, NULL);
    if (XFM == NULL)
    {
      exit(1);
    }
    printf("Applying talairach transform\n");
    MatrixPrint(stdout,XFM);
    MRISmatrixMultiply(mris,XFM);
  }

  if (rescale)
  {
    if (mris->group_avg_surface_area == 0)
    {
      printf("ERROR: cannot rescale a non-group surface\n");
      exit(1);
    }
    scale = sqrt((double)mris->group_avg_surface_area/mris->total_area);
  }

  if (scale > 0)
  {
    printf("scale = %lf\n",scale);
    MRISscale(mris,scale);
    MRIScomputeMetricProperties(mris);
  }

  if (curv_file_flag)
    __convertCurvatureFile(mris, nfcurv, arr_fcurv, out_fname);
  else if (annot_file_flag)
    __convertAnnotFile(mris, annot_fname, gifti_da_num, parcstats_fname, out_fname);
  else if (label_file_flag)
    __convertLabelFile(mris, label_fname, label_name, labelstats_fname, out_fname);
  else if (func_file_flag)
  {
    __convertFuncFile(func_fname, out_fname);
    exit(0);
  }
  else if (mris->patch)
    __convertMRISPatch(mris, out_fname);
  else if (output_normals)
    __convertNormals(mris, out_fname);
  else if (write_vertex_neighbors)
    MRISwriteVertexNeighborsAscii(mris, out_fname) ;
  else if(PrintXYZOnly)
  {
    fp = fopen(out_fname, "w");
    printf("Printing only XYZ to ascii file %s\n", out_fname);

    for(nthvtx = 0; nthvtx < mris->nvertices; nthvtx++)
    {
      fprintf(fp, "%9.4f %9.4f %9.4f\n",
                  mris->vertices[nthvtx].x,
                  mris->vertices[nthvtx].y,
                  mris->vertices[nthvtx].z);
    }

    fclose(fp);
  }
  else if (combinesurfs_flag)
  {
    MRI_SURFACE *mris2 = MRISread(in2_fname, doTkrRASConvert) ;
    if (!mris2)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, in2_fname) ;

    MRIS* mris3 = MRISunion(mris, mris2);
    
    if (mris3 != NULL)
    {
      MRISwrite(mris3, out_fname);
      MRISfree(&mris3);
    }
    MRISfree(&mris2);
  }
  else
  {
    if(RemoveVolGeom){
      printf("Removing Vol Geom\n");
      mris->vg.valid = 0;
    }

    if(MRISfileNameType(out_fname) == MRIS_VOLUME_FILE) {
      printf("Saving surface xyz %s as a volume format\n",out_fname);
      MRI *vol = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);
      MRIcopyMRIS(vol,mris,0,"x");
      MRIcopyMRIS(vol,mris,1,"y");
      MRIcopyMRIS(vol,mris,2,"z");
      MRIwrite(vol,out_fname);
      MRIfree(&vol);
    }
    else
    {
      // default output:
      printf("Saving %s as a surface in %s space\n", out_fname, (mris->useRealRAS) ? "SCANNER" : "TKREGISTER");
      if(VolGeomMRI) {
	printf("Adding Vol Geom\n");
        getVolGeom(VolGeomMRI,&mris->vg);
      }
      MRISwrite(mris, out_fname) ;
    }
  }

  MRISfree(&mris) ;

  exit(0) ;

  return(0) ;  /* for ansi */
}


static MRIS *__splitGIFTI(const char *fgifti, const char *outdir, const char *fout)
{
  MRIS *mris = NULL;

  MRISurfOverlay *fsOverlay = new MRISurfOverlay(fgifti);
  mris = fsOverlay->readGIFTICombined(fgifti, mris);
  if (mris == NULL)
    return NULL;

  fsOverlay->separateGIFTIDataArray(mris, outdir);

  return mris;
}


static void __convertCurvatureFile(MRIS *mris, int noverlay, const char **foverlays, char *out_fname)
{
  MRISurfOverlay *fsOverlay = new MRISurfOverlay(mris, noverlay, foverlays);
  int error = fsOverlay->read(TRUE, mris);
  if (error != NO_ERROR)
    return;

  fsOverlay->write(out_fname, mris, (mergegifti_flag) ? true : false);
}


static void __convertLabelFile(MRIS *mris, const char *flabel, const char *labelname, const char *labelstats, char *out_fname)
{
#if 0  // for some reasons, defining the variables here doesn't work
  const COLOR_TABLE_ENTRY unknown = {"unknown", 0,0,0,0, 0,0,0,0};
  COLOR_TABLE_ENTRY userLabel =
    {
      "user label name gets copied here                   ",
      220,20,20,255, 0.8,0.08,0.08,1
    };
  const CTE *entries[2] = {&unknown, &userLabel};
  const COLOR_TABLE miniColorTable =
    {(CTE**)entries, 2, "miniColorTable", 2};
#endif

  // first read the freesurfer .label file...
  LABEL* label = LabelRead(NULL, flabel);
  if (NULL == labelname)
  {
    printf("ERROR: reading .label file specified for --label: %s\n", flabel);
    exit(1);
  }

  // give this label its own colortable:
  mris->ct = (CT*)&miniColorTable;
  strcpy(userLabel.name, labelname);
  // try to find this label in the FreeSurferColorLUT, so we have a unique
  // color (annotation) for it (otherwise, just use default miniColorTable)
  COLOR_TABLE *ct0;
  char ctabfile[2000];
  sprintf(ctabfile,"%s/FreeSurferColorLUT.txt", getenv("FREESURFER_HOME"));
  ct0 = CTABreadASCII(ctabfile);
  if (ct0)
  {
    int cno;
    for (cno=0; cno < ct0->nentries; cno++) {
      if ((ct0->entries[cno]) && (strcmp(labelname, ct0->entries[cno]->name) == 0)) {
        // we found this label! so update local colortable with info
        memcpy(miniColorTable.entries[1], ct0->entries[cno], sizeof(COLOR_TABLE_ENTRY));
        break;
      }
    }
  }
  // assign annotation to each label vertex (and while we're at it, stats)
  int annotation = CTABrgb2Annotation(miniColorTable.entries[1]->ri,
                                      miniColorTable.entries[1]->gi,
                                      miniColorTable.entries[1]->bi);
  int lno;
  for (lno=0; lno < label->n_points; lno++)
  {
    int vno = label->lv[lno].vno;
    mris->vertices[vno].annotation = annotation;
    mris->vertices[vno].stat = label->lv[lno].stat; // in case --labelstats
  }

  // now write the annot file (either in .annot format, or gifti LabelTable)
  int type = MRISfileNameType(out_fname);
  if (type == MRIS_ANNOT_FILE)
  {
    if (MRISwriteAnnotation(mris, out_fname) != NO_ERROR)
    {
      exit(1);
    }
  }
  else if (type == MRIS_GIFTI_FILE)
  {
    if (MRISwriteGIFTI(mris,NIFTI_INTENT_LABEL, out_fname, NULL) != NO_ERROR)
    {
      exit(1);
    }
  }
  else
  {
    printf("ERROR: unknown file annot file type specified for output: "
           "%s\n", out_fname);
    exit(1);
  }

  mris->ct = NULL; // to avoid calling CTABfree (our table is static memory)

  // if --labelstats was given, then we want to write-out the stats values
  // found in the .label file to a file
  if (labelstats != NULL)
  {
    int type = MRISfileNameType(labelstats);
    if (type == MRIS_GIFTI_FILE)
    {
      if (MRISwriteGIFTI(mris, NIFTI_INTENT_UNIFORM, labelstats, NULL) != NO_ERROR)
      {
        exit(1);
      }
    }
    else
    {
      printf("ERROR: unknown file file type specified for --labelstats: %s\n", labelstats);
      exit(1);
    }
  }
}


static void __convertAnnotFile(MRIS *mris, const char *fannot, int giftiDaNum, const char *parcstats, char *out_fname)
{
    if (MRISreadAnnotation(mris, fannot) != NO_ERROR)
    {
      printf("ERROR: failed to read annotation %s\n", fannot);
      exit(1);
    }

    // read parcstats text file (pairs of parc labels and stat values) and
    // save value associated with that parc label into the vertex with that
    // parc (annot) label
    if (parcstats != NULL)
    {
      FILE* fp;
      if ((fp = fopen(parcstats, "r")) == NULL)
      {
        errno = 0;
        ErrorExit(ERROR_BADFILE, "ERROR: can't open file %s", parcstats);
      }
      char line[STRLEN];
      while (fgets(line, STRLEN, fp) != NULL)
      {
        char label[STRLEN];
        float val;
        sscanf(line,"%s %f", label, &val);
        // get the annotation value for this label from the colortable
        int annot = CTABentryNameToAnnotation(label, mris->ct);
        int vno;
        int doprint=1;
        for (vno=0; vno < mris->nvertices; vno++)
        {
          if (annot == mris->vertices[vno].annotation)
          {
            mris->vertices[vno].curv = val;
            if (doprint)
            {
              printf("label: %s, val: %9.9f\n", label, val);
              doprint = 0;
            }
          }
        }
      }
      // now write the 'curv' data (the parc stats we just assigned) to file
      int type = MRISfileNameType(out_fname) ;
      if (type == MRIS_ASCII_FILE)
      {
        mrisWriteAsciiCurvatureFile(mris, out_fname);
        //writeAsciiCurvFile(mris, out_fname) ;
      }
      else if (type == MRIS_GIFTI_FILE)
      {
        MRISwriteGIFTI(mris, NIFTI_INTENT_SHAPE, out_fname, NULL);
      }
      else
      {
        MRISwriteCurvature(mris, out_fname) ;
      }
      exit(0);
    } // if (parcstats != NULL)

    // if fall through, then write annot file
    if (MRISwriteAnnotation(mris, out_fname) != NO_ERROR)
    {
      printf("ERROR: failed to write annotation %s\n", out_fname);
      exit(1);
    }
}


static void __convertFuncFile(const char *ffunc, char *out_fname)
{
    MRI* mri = MRIread( ffunc );
    if (NULL == mri)
    {
      printf("ERROR: unable to to read %s\n", ffunc);
      exit(1);
    }
    MRIwrite(mri, out_fname);
}


static void __convertNormals(MRIS *mris, char *out_fname)
{
    if (MRISfileNameType(out_fname) == MRIS_ASCII_TRIANGLE_FILE)
      MRISwriteNormalsAscii(mris, out_fname) ;
    else
      MRISwriteNormals(mris, out_fname) ;
}


static void __convertMRISPatch(MRIS *mris, char *out_fname)
{
    int type = MRISfileNameType(out_fname) ;
    if (type == MRIS_GIFTI_FILE)
    {
      MRISwrite(mris, out_fname);
    }
    else
    {
      MRISwritePatch(mris, out_fname) ;
    }
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "center") || !stricmp(option, "-center"))
  {
    center_surface = 1 ;
    printf("centering surface\n") ;
  }
  else if (!stricmp(option, "-annot"))
  {
    annot_fname = argv[2] ;
    annot_file_flag = 1;
    nargs = 1 ;
  }
  else if (!stricmp(option, "-label2mask")) {
    MRIS *surf = MRISread(argv[2]);
    if(surf==NULL) exit(1);
    LABEL *srclabel = LabelRead(NULL, argv[3]);
    if(srclabel==NULL) exit(1);
    MRI *outmask = MRISlabel2Mask(surf,srclabel,NULL);
    if(outmask==NULL) exit(1);
    int err = MRIwrite(outmask,argv[4]);
    MRIfree(&outmask); MRISfree(&surf); LabelFree(&srclabel) ;
    exit(err);
  }
  else if (!stricmp(option, "-area")) {
    // This little bit of code is self-contained, run like
    // mris_convert --area surface area.mgz
    MRIS *surf = MRISread(argv[2]);
    if(surf==NULL) exit(1);
    MRI *SrcVals = MRIcopyMRIS(NULL, surf, 0, "area");
    if(surf->group_avg_surface_area > 0) {
      double val = surf->group_avg_surface_area / surf->total_area;
      printf("group surface, scaling area by %g\n",val);
      MRIscalarMul(SrcVals,SrcVals,val);
    }
    printf("Writing vertex area to %s\n",argv[3]);
    int err = MRIwrite(SrcVals,argv[3]);
    MRIfree(&SrcVals); MRISfree(&surf);
    exit(err);
  }
  else if (!stricmp(option, "-angle")) {
    // This little bit of code is self-contained. Run like:
    // mris_convert --angle <surface> <angles.mgz>
    MRIS *surf = MRISread(argv[2]);
    if (!surf) exit(1);
    MRI *angles = computeAngles(surf);
    int error = MRIwrite(angles, argv[3]);
    MRIfree(&angles);
    MRISfree(&surf);
    exit(error);
  }
  else if (!stricmp(option, "-volume")){
    // This little bit of code is self-contained, run like
    // mris_convert --volume subject hemi outcurv
    // to compute vertex-wise volume. Better than using
    // thickness * midarea
    int err;
    err = ComputeMRISvolumeTH3(argv[2],argv[3],1,argv[4]);
    exit(err);
  }
  else if (!stricmp(option, "-da_num"))
  {
    sscanf(argv[2],"%d",&gifti_da_num);
    nargs = 1 ;
  }
  else if (!stricmp(option, "-label"))
  {
    label_fname = argv[2] ;
    label_name = argv[3] ;
    label_file_flag = 1;
    nargs = 2 ;
  }
  else if (!stricmp(option, "-labelstats"))
  {
    labelstats_fname = argv[2] ;
    labelstats_file_flag = 1;
    nargs = 1 ;
  }
  else if (!stricmp(option, "-parcstats"))
  {
    parcstats_fname = argv[2] ;
    parcstats_file_flag = 1;
    nargs = 1 ;
  }
  else if (!stricmp(option, "-combinesurfs"))
  {
    combinesurfs_flag = 1;
#if __COMBINESURFS_TAKE_INFILE2
    in2_fname = argv[2];
    nargs = 1 ;
#endif
  }
  else if (!stricmp(option, "-mergegifti"))
  {
    mergegifti_flag = 1;
  }
  else if (!stricmp(option, "-splitgifti"))
  {
    splitgifti_flag = 1;
  }
  else if (!stricmp(option, "-giftioutdir"))
  {
    giftioutdir = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "-delete-cmds"))
  {
    DeleteCommands = 1;
  }
  else if (!stricmp(option, "-userealras"))
  {
    userealras_flag = 1;
    //usesurfras_flag = 0;
  }
  else if (!stricmp(option, "-usesurfras"))
  {
    //userealras_flag = 0;
    usesurfras_flag = 1;
  }
  else if (!stricmp(option, "-cras_correction") || !stricmp(option, "-cras_add")) {
    cras_add = 1;
    cras_subtract = 0;
    ToScanner = 0;
    ToTkr = 0;
  }
  else if (!stricmp(option, "-cras_remove") || !stricmp(option, "-cras_subtract")) {
    cras_add = 0;
    cras_subtract = 1;
    ToScanner = 0;
    ToTkr = 0;
  }
  else if (!stricmp(option, "-to-scanner")) {
    ToScanner = 1;
    //ToTkr = 0;
    cras_add = 0;
    cras_subtract = 0;
  }
  else if (!stricmp(option, "-to-tkr")) {
    //ToScanner = 0;
    ToTkr = 1;
    cras_add = 0;
    cras_subtract = 0;
  }
  else if (!stricmp(option, "-upsample")) {
    sscanf(argv[2],"%d",&nUpsample);
    sscanf(argv[3],"%d",&UpsampleSortType);
    nargs = 2 ;
  }
  else if (!stricmp(option, "-to-surf")) {
    ToSurfCoords = 1;
    SurfCoords = MRISread(argv[2], doTkrRASConvert);
    if(SurfCoords == NULL) exit (1);
    nargs = 1 ;
  }
  else if (!stricmp(option, "-vol-geom"))
  {
    VolGeomMRI = MRIreadHeader(argv[2],MRI_VOLUME_TYPE_UNKNOWN);
    if(VolGeomMRI == NULL)
    {
      exit(1);
    }
    nargs = 1 ;
  }
  else if (!stricmp(option, "-remove-vol-geom"))
  {
    RemoveVolGeom = 1;
  }
  else switch (toupper(*option))
    {
    case 'A':
      PrintXYZOnly = 1;
      break ;
    case 'F':
      func_file_flag = 1 ;
      func_fname = argv[2] ;
      nargs = 1 ;
      break ;
    case 'C':
    {
      curv_file_flag = 1 ;
      //curv_fname = argv[2] ;
      nargs = 1 ;

      int argc0 = 2;
      int nth = 3;
      // get additional input curvature files
      // minus 2 positional arguments
      while (nth_is_arg(argc-nth-2, argv, nth))
      {
        nargs++; nth++;
      }

      arr_fcurv = new const char*[nargs];
      for (int n = 0; n < nargs; n++)
      {
        arr_fcurv[n] = argv[argc0];
        argc0++;
      }

      nfcurv = nargs;
      break ;
    }
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

static void check_options(void)
{
  if (nfcurv > 1 && !mergegifti_flag)
  {
    printf("ERROR: more than one overlay are specified. Use --mergegifti to combine them into single GIFTI output.\n");
    exit(1);
  }

  if (mergegifti_flag && MRISurfOverlay::getFileFormat(outfile) != GIFTI_FILE)
  {
    printf("ERROR: --mergegifti only works with GIFTI output\n");
    exit(1);
  }

  if (splitgifti_flag && MRISurfOverlay::getFileFormat(infile) != GIFTI_FILE)
  {
    printf("ERROR: --splitgifti only works with GIFTI input\n");
    exit(1);
  }

  if (splitgifti_flag && giftioutdir == NULL)
  {
    printf("ERROR: use --giftioutdir to specify output directory for generated GIFTI files\n");
    exit(1);
  }

  if (combinesurfs_flag && curv_file_flag)
  {
    printf("ERROR: --combinesurf and -c are mutually exclusive. \n");
    exit(1);
  }

  // --to-scanner/--userealras and --to-tkr/--usesurfras are mutually exclusive
  if ( (ToScanner || userealras_flag) && 
       (ToTkr     || usesurfras_flag) )
  {
    printf("ERROR: --to-scanner/--userealras and --to-tkr/--usesurfras are mutually exclusive. \n");
    exit(1);
  }

  if (talxfmsubject && curv_file_flag)
  {
    printf("ERROR: cannot specify -t and -c\n");
    exit(1);
  }

  if (labelstats_file_flag && ! label_file_flag)
  {
    printf("ERROR: cannot specify --labelstats without --label\n");
    exit(1);
  }

  if (parcstats_file_flag && ! annot_file_flag)
  {
    printf("ERROR: cannot specify --parcstats without --annot\n");
    exit(1);
  }
}

/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}

/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth) {
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if (nargc <= 0) return(0);

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);

  return(1);
}

static void 
usage_exit(void)
{
  print_help() ;
  exit(1) ;
}


#include "mris_convert.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_convert_help_xml, mris_convert_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  printf( "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static int
convertToWFile(char *in_fname, char *out_fname)
{
  FILE   *infp, *outfp ;
  char   line[300], *cp ;
  int    vno, l = 0, num, ilat ;
  float  val ;

  fprintf(stderr, "writing w file %s...\n", out_fname) ;
  outfp = fopen(out_fname,"wb");
  if (outfp==NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,out_fname) ;
  }

  infp = fopen(in_fname,"rb");
  if (infp==NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,in_fname) ;
  }


  cp = fgetl(line, 299, infp) ;
  ilat = atoi(cp) ; /* not used at the moment */
  cp = fgetl(line, 299, infp) ;
  num = atoi(cp) ;
  fwrite2(0,outfp);
  fwrite3(num,outfp);

  while ((cp = fgetl(line, 299, infp)) != NULL)
  {
    l++ ;
    if (sscanf(cp, "%d %f", &vno, &val) != 2)
    {
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
convertFromWFile(char *in_fname, char *out_fname)
{
  FILE   *infp, *outfp ;
  int    vno, num, ilat, i ;
  float  val, lat ;

  fprintf(stderr, "writing ascii w file %s...\n", out_fname) ;
  outfp = fopen(out_fname,"wb");
  if (outfp==NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,out_fname) ;
  }

  infp = fopen(in_fname,"rb");
  if (infp==NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,in_fname) ;
  }

  fread2(&ilat,infp);
  lat = (float)ilat/10.0;
  fprintf(outfp, "%2.3f\n", lat) ;
  fread3(&num,infp);
  fprintf(outfp, "%d\n", num) ;
  for (i=0; i<num; i++)
  {
    fread3(&vno,infp);
    val = freadFloat(infp) ;
    fprintf(outfp, "%d  %f\n", vno, val) ;
  }
  fclose(outfp);
  fclose(infp) ;

  return(NO_ERROR) ;
}


// mrisurf_io.cpp now has a copy of this - mrisWriteAsciiCurvatureFile()
//static int
//writeAsciiCurvFile(MRI_SURFACE *mris, char *out_fname)
//{
//  FILE   *fp ;
//  int    vno ;
//  VERTEX *v ;
//
//  fp = fopen(out_fname, "w") ;
//  if (!fp)
//    ErrorExit(ERROR_BADFILE, "%s could not open output file %s.\n",
//              Progname, out_fname) ;
//
//  for (vno = 0 ; vno < mris->nvertices ; vno++)
//  {
//    v = &mris->vertices[vno] ;
//    if (v->ripflag)
//    {
//      continue ;
//    }
//    fprintf(fp, "%3.3d %2.5f %2.5f %2.5f %2.5f\n",
//            vno, v->x, v->y, v->z, v->curv) ;
//  }
//
//  fclose(fp) ;
//  return(NO_ERROR) ;
//}

/*
  \fn int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname)
  \brief Writes out neighbors of a vertex in each row. The first
  column is the vertex number, the 2nd col is the number of neighbors,
  the remaining cols are the vertex numbers of the neighbors.
*/

static int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname)
{
  int vno, nnbrs, nbrvno;
  FILE *fp;

  fp = fopen(out_fname,"w");
  if(fp == NULL)
  {
    printf("ERROR: opening %s\n",out_fname);
    exit(1);
  }

  for(vno=0; vno < mris->nvertices; vno++)
  {
    nnbrs = mris->vertices_topology[vno].vnum;
    fprintf(fp,"%6d %2d   ",vno,nnbrs);
    for (nbrvno = 0; nbrvno < nnbrs; nbrvno++)
    {
      fprintf(fp,"%6d ",mris->vertices_topology[vno].v[nbrvno]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  return(0);
}


MRI* computeAngles(MRIS* surf)
{
  // added by OVOX (see pull request #707)

  printf("************************************* CAUTION **************************************\n"); 
  printf("If you are running \"mris_convert --angle\" on an original surface (i.e. has not\n"); 
  printf("been transformed by mris_surf2surf), then you can ignore this message.\n");
  printf("However, if you need to calculate orientation angles of a TRANSFOMRED surface you \n");
  printf("will need to convert the .lta/.dat registration file first and force the original \n");
  printf("volume as the source. Otherwise important header information will be lost and the \n");
  printf("angles will be incorrect (do not be fooled, they might appear close, but they won't\n");
  printf("be exact). A detailed description and example can be found on the freesurfer wiki.\n");
  printf("************************************************************************************\n");

  // this will enable angle-weighted surface normals
  UnitizeNormalFace = 0;
  MRIScomputeNormals(surf);
  UnitizeNormalFace = 1;

  MATRIX* vox2ras = vg_i_to_r(&surf->vg);  // scanner space vox2ras from volume geometry
  MATRIX* ras2vox = MatrixInverse(vox2ras, nullptr);  // inverse from imaging reference frame to scanner coordinate system

  // imaging volume axes
  double vi_x[3] = { ras2vox->rptr[1][1], ras2vox->rptr[2][1], ras2vox->rptr[3][1] };
  double vi_y[3] = { ras2vox->rptr[1][2], ras2vox->rptr[2][2], ras2vox->rptr[3][2] };
  double vi_z[3] = { ras2vox->rptr[1][3], ras2vox->rptr[2][3], ras2vox->rptr[3][3] };

  // normalize the axes
  auto normalize = [](double* vec) {
    double len = sqrt(pow(vec[0], 2) +  pow(vec[1], 2) + pow(vec[2], 2));
    vec[0] /= len;
    vec[1] /= len;
    vec[2] /= len;
  };

  normalize(vi_x);
  normalize(vi_y);
  normalize(vi_z);

  // the angle between normal and volume axes is simply the acos of the normal's coordinate, 
  // the angle to the scanner axes is the acos of the dot product between the normal and 
  // the scanner directions. Store all angles in one overlay file. The first three will be w.r.t.
  // B0, the other w.r.t. the scanner volume
  MRI *angles = MRIallocSequence(surf->nvertices, 1, 1, MRI_FLOAT, 6);

  for (int vno = 0; vno < surf->nvertices; vno++) {
    VERTEX* v = &surf->vertices[vno];
    MRIsetVoxVal(angles, vno, 0, 0, 0, acos((vi_x[0] * v->nx + vi_x[1] * v->ny + vi_x[2] * v->nz)) * 180 / PI); 
    MRIsetVoxVal(angles, vno, 0, 0, 1, acos((vi_y[0] * v->nx + vi_y[1] * v->ny + vi_y[2] * v->nz)) * 180 / PI); 
    MRIsetVoxVal(angles, vno, 0, 0, 2, acos((vi_z[0] * v->nx + vi_z[1] * v->ny + vi_z[2] * v->nz)) * 180 / PI); 
    MRIsetVoxVal(angles, vno, 0, 0, 3, acos(v->nx) * 180 / PI) ;
    MRIsetVoxVal(angles, vno, 0, 0, 4, acos(v->ny) * 180 / PI); 
    MRIsetVoxVal(angles, vno, 0, 0, 5, acos(v->nz) * 180 / PI);
  }

  return angles;
}
