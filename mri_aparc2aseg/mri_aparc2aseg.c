/**
 * @file  mri_aparc2aseg.c
 * @brief Maps aparc labels to aseg
 *
 * Maps the cortical labels from the automatic cortical parcellation (aparc)
 * to the automatic segmentation volume (aseg). The result can be used as
 * the aseg would. The algorithm is to find each aseg voxel labeled as
 * cortex (3 and 42) and assign it the label of the closest cortical vertex.
 * If the voxel is not in the ribbon (as defined by mri/lh.ribbon and
 * rh.ribbon), then the voxel is marked as unknown (0). This can be turned
 * off with --noribbon. The cortical parcellation is obtained from
 * subject/label/hemi.aparc.annot which should be based on the
 * curvature.buckner40.filled.desikan_killiany.gcs atlas
 * The aseg is obtained from subject/mri/aseg.mgz and should be based on
 * the RB40_talairach_2005-07-20.gca atlas. If these atlases are used, then
 * the segmentations can be viewed with tkmedit and the FreeSurferColorLUT.txt
 * color table found in $FREESURFER_HOME. These are the default atlases
 * used by recon-all.
 */
/*
 * Original Author: Doug Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/11/15 16:34:15 $
 *    $Revision: 1.43 $
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
#include <math.h>
#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "annotation.h"
#include "version.h"
#include "mrisegment.h"
#include "cma.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
int FindClosestLRWPVertexNo(int c, int r, int s,
                            int *lhwvtx, int *lhpvtx,
                            int *rhwvtx, int *rhpvtx,
                            MATRIX *Vox2RAS,
                            MRIS *lhwite,  MRIS *lhpial,
                            MRIS *rhwhite, MRIS *rhpial,
                            MHT *lhwhite_hash, MHT *lhpial_hash,
                            MHT *rhwhite_hash, MHT *rhpial_hash);
int CCSegment(MRI *seg, int segid, int segidunknown);

int main(int argc, char *argv[]) ;

static char vcid[] =
  "$Id: mri_aparc2aseg.c,v 1.43 2012/11/15 16:34:15 greve Exp $";
char *Progname = NULL;
static char *SUBJECTS_DIR = NULL;
static char *subject = NULL;
static char *OutASegFile = NULL;
static char *OutAParcFile = NULL;
static char *OutDistFile = NULL;
static int debug = 0;
static int UseRibbon = 0;
static int UseNewRibbon = 1;
static MRI *ASeg, *filled, *mritmp;
static MRI *AParc;
static MRI *Dist;
static MRI *lhRibbon,*rhRibbon,*RibbonSeg;
static MRIS *lhwhite, *rhwhite;
static MRIS *lhpial, *rhpial;
static MHT *lhwhite_hash, *rhwhite_hash;
static MHT *lhpial_hash, *rhpial_hash;
static VERTEX vtx;
static int  lhwvtx, lhpvtx, rhwvtx, rhpvtx;
static MATRIX *Vox2RAS, *CRS, *RAS;
static float dlhw, dlhp, drhw, drhp;
static float dmaxctx = 5.0;
static int LabelWM=0;
static int LabelHypoAsWM=0;
static int RipUnknown = 0;

static char tmpstr[2000];
static char annotfile[1000];
static char *annotname = "aparc";
static char *asegname = "aseg";
static int baseoffset = 0;
static float hashres = 16;

int crsTest = 0, ctest=0, rtest=0, stest=0;
int UseHash = 1;

char *CtxSegFile = NULL;
MRI *CtxSeg = NULL;

int FixParaHipWM = 1;

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, err, asegid, c, r, s, nctx, annot,vtxno,nripped;
  int annotid, IsCortex=0, IsWM=0, IsHypo=0, hemi=0, segval=0;
  int RibbonVal=0,nbrute=0;
  float dmin=0.0, lhRibbonVal=0, rhRibbonVal=0;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0)
  {
    usage_exit();
  }

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL)
  {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  /* ------ Load subject's lh white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
  printf("\nReading lh white surface \n %s\n",tmpstr);
  lhwhite = MRISread(tmpstr);
  if (lhwhite == NULL)
  {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  /* ------ Load subject's lh pial surface ------ */
  sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,subject);
  printf("\nReading lh pial surface \n %s\n",tmpstr);
  lhpial = MRISread(tmpstr);
  if (lhpial == NULL)
  {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  if (lhwhite->nvertices != lhpial->nvertices)
  {
    printf("ERROR: lh white and pial have a different number of "
           "vertices (%d,%d)\n",
           lhwhite->nvertices,lhpial->nvertices);
    exit(1);
  }

  /* ------ Load lh annotation ------ */
  sprintf(annotfile,"%s/%s/label/lh.%s.annot",SUBJECTS_DIR,subject,annotname);
  printf("\nLoading lh annotations from %s\n",annotfile);
  err = MRISreadAnnotation(lhwhite, annotfile);
  if (err)
  {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
    exit(1);
  }

  /* ------ Load subject's rh white surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
  printf("\nReading rh white surface \n %s\n",tmpstr);
  rhwhite = MRISread(tmpstr);
  if (rhwhite == NULL)
  {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  /* ------ Load subject's rh pial surface ------ */
  sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,subject);
  printf("\nReading rh pial surface \n %s\n",tmpstr);
  rhpial = MRISread(tmpstr);
  if (rhpial == NULL)
  {
    fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
    exit(1);
  }
  if (rhwhite->nvertices != rhpial->nvertices)
  {
    printf("ERROR: rh white and pial have a different "
           "number of vertices (%d,%d)\n",
           rhwhite->nvertices,rhpial->nvertices);
    exit(1);
  }

  /* ------ Load rh annotation ------ */
  sprintf(annotfile,"%s/%s/label/rh.%s.annot",SUBJECTS_DIR,subject,annotname);
  printf("\nLoading rh annotations from %s\n",annotfile);
  err = MRISreadAnnotation(rhwhite, annotfile);
  if (err)
  {
    printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
    exit(1);
  }

  if (lhwhite->ct)
  {
    printf("Have color table for lh white annotation\n");
  }
  if (rhwhite->ct)
  {
    printf("Have color table for rh white annotation\n");
  }
  //print_annotation_table(stdout);

  if (UseRibbon)
  {
    sprintf(tmpstr,"%s/%s/mri/lh.ribbon.mgz",SUBJECTS_DIR,subject);
    printf("Loading lh ribbon mask from %s\n",tmpstr);
    lhRibbon = MRIread(tmpstr);
    if (lhRibbon == NULL)
    {
      printf("ERROR: loading %s\n",tmpstr);
      exit(1);
    }
    sprintf(tmpstr,"%s/%s/mri/rh.ribbon.mgz",SUBJECTS_DIR,subject);
    printf("Loading rh ribbon mask from %s\n",tmpstr);
    rhRibbon = MRIread(tmpstr);
    if (rhRibbon == NULL)
    {
      printf("ERROR: loading  %s\n",tmpstr);
      exit(1);
    }
  }

  if (UseNewRibbon)
  {
    sprintf(tmpstr,"%s/%s/mri/ribbon.mgz",SUBJECTS_DIR,subject);
    printf("Loading ribbon segmentation from %s\n",tmpstr);
    RibbonSeg = MRIread(tmpstr);
    if (RibbonSeg == NULL)
    {
      printf("ERROR: loading %s\n",tmpstr);
      exit(1);
    }
  }

  if (LabelHypoAsWM)
  {
    sprintf(tmpstr,"%s/%s/mri/filled.mgz",SUBJECTS_DIR,subject);
    printf("Loading filled from %s\n",tmpstr);
    filled = MRIread(tmpstr);
    if (filled == NULL)
    {
      printf("ERROR: loading filled %s\n",tmpstr);
      exit(1);
    }
  }

  // ------------ Rip -----------------------
  if (RipUnknown)
  {
    printf("Ripping vertices labeled as unkown\n");
    nripped = 0;
    for (vtxno = 0; vtxno < lhwhite->nvertices; vtxno++)
    {
      annot = lhwhite->vertices[vtxno].annotation;
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
      // Sometimes the annotation will be "none" indicated by
      // annotid = -1. We interpret this as "unknown".
      if (annotid == 0 || annotid == -1)
      {
        lhwhite->vertices[vtxno].ripflag = 1;
        lhpial->vertices[vtxno].ripflag = 1;
        nripped++;
      }
    }
    printf("Ripped %d vertices from left hemi\n",nripped);
    nripped = 0;
    for (vtxno = 0; vtxno < rhwhite->nvertices; vtxno++)
    {
      annot = rhwhite->vertices[vtxno].annotation;
      CTABfindAnnotation(rhwhite->ct, annot, &annotid);
      if (annotid == 0 || annotid == -1)
      {
        rhwhite->vertices[vtxno].ripflag = 1;
        rhpial->vertices[vtxno].ripflag = 1;
        nripped++;
      }
    }
    printf("Ripped %d vertices from right hemi\n",nripped);
  }


  printf("\n");
  printf("Building hash of lh white\n");
  lhwhite_hash = MHTfillVertexTableRes(lhwhite, NULL,CURRENT_VERTICES,hashres);
  printf("\n");
  printf("Building hash of lh pial\n");
  lhpial_hash = MHTfillVertexTableRes(lhpial, NULL,CURRENT_VERTICES,hashres);
  printf("\n");
  printf("Building hash of rh white\n");
  rhwhite_hash = MHTfillVertexTableRes(rhwhite, NULL,CURRENT_VERTICES,hashres);
  printf("\n");
  printf("Building hash of rh pial\n");
  rhpial_hash = MHTfillVertexTableRes(rhpial, NULL,CURRENT_VERTICES,hashres);

  /* ------ Load ASeg ------ */
  sprintf(tmpstr,"%s/%s/mri/%s.mgz",SUBJECTS_DIR,subject,asegname);
  if (!fio_FileExistsReadable(tmpstr))
  {
    sprintf(tmpstr,"%s/%s/mri/%s.mgh",SUBJECTS_DIR,subject,asegname);
    if (!fio_FileExistsReadable(tmpstr))
    {
      sprintf(tmpstr,"%s/%s/mri/aseg/COR-.info",SUBJECTS_DIR,subject);
      if (!fio_FileExistsReadable(tmpstr))
      {
        printf("ERROR: cannot find aseg\n");
        exit(1);
      }
      else
      {
        sprintf(tmpstr,"%s/%s/mri/aseg/",SUBJECTS_DIR,subject);
      }
    }
  }

  printf("\nLoading aseg from %s\n",tmpstr);
  ASeg = MRIread(tmpstr);
  if (ASeg == NULL)
  {
    printf("ERROR: loading aseg %s\n",tmpstr);
    exit(1);
  }
  mritmp = MRIchangeType(ASeg,MRI_INT,0,0,1);
  MRIfree(&ASeg);
  ASeg = mritmp;

  if (CtxSegFile)
  {
    printf("Loading Ctx Seg File %s\n",CtxSegFile);
    CtxSeg = MRIread(CtxSegFile);
    if (CtxSeg == NULL)
    {
      exit(1);
    }
  }

  AParc = MRIclone(ASeg,NULL);
  if (OutDistFile != NULL)
  {
    Dist = MRIclone(ASeg,NULL);
    mritmp = MRIchangeType(Dist,MRI_FLOAT,0,0,0);
    if (mritmp == NULL)
    {
      printf("ERROR: could change type\n");
      exit(1);
    }
    MRIfree(&Dist);
    Dist = mritmp;
  }

  Vox2RAS = MRIxfmCRS2XYZtkreg(ASeg);
  printf("ASeg Vox2RAS: -----------\n");
  MatrixPrint(stdout,Vox2RAS);
  printf("-------------------------\n");
  CRS = MatrixAlloc(4,1,MATRIX_REAL);
  CRS->rptr[4][1] = 1;
  RAS = MatrixAlloc(4,1,MATRIX_REAL);
  RAS->rptr[4][1] = 1;

  if (crsTest)
  {
    printf("Testing point %d %d %d\n",ctest,rtest,stest);
    err = FindClosestLRWPVertexNo(ctest,rtest,stest,
                                  &lhwvtx, &lhpvtx,
                                  &rhwvtx, &rhpvtx, Vox2RAS,
                                  lhwhite,  lhpial,
                                  rhwhite, rhpial,
                                  lhwhite_hash, lhpial_hash,
                                  rhwhite_hash, rhpial_hash);

    printf("Result: err = %d\n",err);
    exit(err);
  }

  printf("\nLabeling Slice\n");
  nctx = 0;
  annot = 0;
  annotid = 0;
  nbrute = 0;

  // Go through each voxel in the aseg
  for (c=0; c < ASeg->width; c++)
  {
    printf("%3d ",c);
    if (c%20 ==19)
    {
      printf("\n");
    }
    fflush(stdout);
    for (r=0; r < ASeg->height; r++)
    {
      for (s=0; s < ASeg->depth; s++)
      {

        asegid = MRIgetVoxVal(ASeg,c,r,s,0);
        if (asegid == 3 || asegid == 42)
        {
          IsCortex = 1;
        }
        else
        {
          IsCortex = 0;
        }
        if (asegid >= 77 && asegid <= 82)
        {
          IsHypo = 1;
        }
        else
        {
          IsHypo = 0;
        }
        if (asegid == 2 || asegid == 41)
        {
          IsWM = 1;
        }
        else
        {
          IsWM = 0;
        }
        if (IsHypo && LabelHypoAsWM && MRIgetVoxVal(filled,c,r,s,0))
        {
          IsWM = 1;
        }

        // integrate surface information
        //
        // Only Do This for GM,WM or Unknown labels in the ASEG !!!
        //
        // priority is given to the ribbon computed from the surface
        // namely
        //  ribbon=GM => GM
        //  aseg=GM AND ribbon=WM => WM
        //  ribbon=UNKNOWN => UNKNOWN
        if (UseNewRibbon && ( IsCortex || IsWM || asegid==0 ) )
        {
          RibbonVal = MRIgetVoxVal(RibbonSeg,c,r,s,0);
          MRIsetVoxVal(ASeg,c,r,s,0, RibbonVal);
          if (RibbonVal==2 || RibbonVal==41)
          {
            IsWM = 1;
            IsCortex = 0;
          }
          else if (RibbonVal==3 || RibbonVal==42)
          {
            IsWM = 0;
            IsCortex = 1;
          }
          if (RibbonVal==0)
          {
            IsWM = 0;
            IsCortex = 0;
          }
        }

        // If it's not labeled as cortex or wm in the aseg, skip
        if (!IsCortex && !IsWM)
        {
          continue;
        }

        // If it's wm but not labeling wm, skip
        if (IsWM && !LabelWM)
        {
          continue;
        }

        // Check whether this point is in the ribbon
        if (UseRibbon)
        {
          lhRibbonVal = MRIgetVoxVal(lhRibbon,c,r,s,0);
          rhRibbonVal = MRIgetVoxVal(rhRibbon,c,r,s,0);
          if (IsCortex)
          {
            // ASeg says it's in cortex
            if (lhRibbonVal < 0.5 && rhRibbonVal < 0.5)
            {
              // but it is not part of the ribbon,
              // so set it to unknown (0) and go to the next voxel.
              MRIsetVoxVal(ASeg,c,r,s,0,0);
              continue;
            }
          }
        }

        // Convert the CRS to RAS
        CRS->rptr[1][1] = c;
        CRS->rptr[2][1] = r;
        CRS->rptr[3][1] = s;
        RAS = MatrixMultiply(Vox2RAS,CRS,RAS);
        vtx.x = RAS->rptr[1][1];
        vtx.y = RAS->rptr[2][1];
        vtx.z = RAS->rptr[3][1];

        // Get the index of the closest vertex in the
        // lh.white, lh.pial, rh.white, rh.pial
        if (UseHash)
        {
          lhwvtx = MHTfindClosestVertexNo(lhwhite_hash,lhwhite,&vtx,&dlhw);
          lhpvtx = MHTfindClosestVertexNo(lhpial_hash, lhpial, &vtx,&dlhp);
          rhwvtx = MHTfindClosestVertexNo(rhwhite_hash,rhwhite,&vtx,&drhw);
          rhpvtx = MHTfindClosestVertexNo(rhpial_hash, rhpial, &vtx,&drhp);
          if (lhwvtx < 0 && lhpvtx < 0 && rhwvtx < 0 && rhpvtx < 0)
          {
            /*
            printf("  Could not map to any surface with hash table:\n");
            printf("  crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
            c,r,s,vtx.x,vtx.y,vtx.z);
            printf("  Using brute force search %d ... \n",nbrute);
            fflush(stdout);
            */
            lhwvtx = MRISfindClosestVertex(lhwhite,vtx.x,vtx.y,vtx.z,&dlhw);
            lhpvtx = MRISfindClosestVertex(lhpial,vtx.x,vtx.y,vtx.z,&dlhp);
            rhwvtx = MRISfindClosestVertex(rhwhite,vtx.x,vtx.y,vtx.z,&drhw);
            rhpvtx = MRISfindClosestVertex(rhpial,vtx.x,vtx.y,vtx.z,&drhp);
            nbrute ++;
            //exit(1);
          }
        }
        else
        {
          lhwvtx = MRISfindClosestVertex(lhwhite,vtx.x,vtx.y,vtx.z,&dlhw);
          lhpvtx = MRISfindClosestVertex(lhpial,vtx.x,vtx.y,vtx.z,&dlhp);
          rhwvtx = MRISfindClosestVertex(rhwhite,vtx.x,vtx.y,vtx.z,&drhw);
          rhpvtx = MRISfindClosestVertex(rhpial,vtx.x,vtx.y,vtx.z,&drhp);
        }

        if (lhwvtx < 0)
        {
          dlhw = 1000000000000000.0;
        }
        if (lhpvtx < 0)
        {
          dlhp = 1000000000000000.0;
        }
        if (rhwvtx < 0)
        {
          drhw = 1000000000000000.0;
        }
        if (rhpvtx < 0)
        {
          drhp = 1000000000000000.0;
        }

        if (dlhw <= dlhp && dlhw < drhw && dlhw < drhp && lhwvtx >= 0)
        {
          annot = lhwhite->vertices[lhwvtx].annotation;
          hemi = 1;
          if (lhwhite->ct)
          {
            CTABfindAnnotation(lhwhite->ct, annot, &annotid);
          }
          else
          {
            annotid = annotation_to_index(annot);
          }
          dmin = dlhw;
        }
        if (dlhp < dlhw && dlhp < drhw && dlhp < drhp && lhpvtx >= 0)
        {
          annot = lhwhite->vertices[lhpvtx].annotation;
          hemi = 1;
          if (lhwhite->ct)
          {
            CTABfindAnnotation(lhwhite->ct, annot, &annotid);
          }
          else
          {
            annotid = annotation_to_index(annot);
          }
          dmin = dlhp;
        }

        if (drhw < dlhp && drhw < dlhw && drhw <= drhp && rhwvtx >= 0)
        {
          annot = rhwhite->vertices[rhwvtx].annotation;
          hemi = 2;
          if (rhwhite->ct)
          {
            CTABfindAnnotation(rhwhite->ct, annot, &annotid);
          }
          else
          {
            annotid = annotation_to_index(annot);
          }
          dmin = drhw;
        }
        if (drhp < dlhp && drhp < drhw && drhp < dlhw && rhpvtx >= 0)
        {
          annot = rhwhite->vertices[rhpvtx].annotation;
          hemi = 2;
          if (rhwhite->ct)
          {
            CTABfindAnnotation(rhwhite->ct, annot, &annotid);
          }
          else
          {
            annotid = annotation_to_index(annot);
          }
          dmin = drhp;
        }

        // Sometimes the annotation will be "none" indicated by
        // annotid = -1. We interpret this as "unknown".
        if (annotid == -1)
        {
          annotid = 0;
        }

	/* If the cortical label is "unkown", it is difficult to
	   determine what to put here. If the aseg says it is WM, then
	   that is kept. If the aseg says it is GM, then it is given
	   "ctx-?h-unknown". These voxels can show up in funny places
	   (eg, between hippo and amyg), so this is really just a
	   hack. The real fix should be the surface creation or the
	   aseg. */
	if(annotid == 0 && !LabelWM){
	  if(asegid == Left_Cerebral_Cortex)  MRIsetVoxVal(ASeg,c,r,s,0,1000);
	  else if(asegid == Right_Cerebral_Cortex) MRIsetVoxVal(ASeg,c,r,s,0,2000);
	  else MRIsetVoxVal(ASeg,c,r,s,0,asegid);
	  continue;
	  //{if(hemi == 1) MRIsetVoxVal(ASeg,c,r,s,0,Left_Cerebral_White_Matter);
	  //if(hemi == 2) MRIsetVoxVal(ASeg,c,r,s,0,Right_Cerebral_White_Matter);
	  //continue;
	}

        // why was this here in the first place?
        /*
               if (annotid == 0 &&
                   lhwvtx >= 0 &&
                   lhpvtx >= 0 &&
                   rhwvtx >= 0 &&
                   rhpvtx >= 0) {
                 printf("%d %d %d %d\n",
                        lhwhite->vertices[lhwvtx].ripflag,
                        lhpial->vertices[lhpvtx].ripflag,
                        rhwhite->vertices[rhwvtx].ripflag,
                        rhpial->vertices[rhpvtx].ripflag);
          } */

        if ( IsCortex && hemi == 1)
        {
          segval = annotid+1000 + baseoffset;  //ctx-lh
        }
        if ( IsCortex && hemi == 2)
        {
          segval = annotid+2000 + baseoffset;  //ctx-rh
        }
        if (!IsCortex && hemi == 1)
        {
          segval = annotid+3000 + baseoffset;  // wm-lh
        }
        if (!IsCortex && hemi == 2)
        {
          segval = annotid+4000 + baseoffset;  // wm-rh
        }

        if (!IsCortex && dmin > dmaxctx && hemi == 1)
        {
          segval = 5001;
        }
        if (!IsCortex && dmin > dmaxctx && hemi == 2)
        {
          segval = 5002;
        }

        // This is a hack for getting the right cortical seg with --rip-unknown
        // The aparc+aseg should be passed as CtxSeg. Used with WMParc
        if (IsCortex && CtxSeg)
        {
          segval = MRIgetVoxVal(CtxSeg,c,r,s,0);
        }

        MRIsetVoxVal(ASeg,c,r,s,0,segval);
        MRIsetVoxVal(AParc,c,r,s,0,annot);
        if (OutDistFile != NULL)
        {
          MRIsetVoxVal(Dist,c,r,s,0,dmin);
        }

        if (debug || annotid == -1)
        {
          // Gets here when there is no label at the found vertex.
          // This is different than having a vertex labeled as "unknown"
          if (!debug)
          {
            continue;
          }
          printf("\n");
          printf("Found closest vertex, but it has no label.\n");
          printf("aseg id = %d\n",asegid);
          printf("crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
                 c,r,s,vtx.x,vtx.y,vtx.z);
          if (lhwvtx > 0) printf("lhw  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   lhwvtx, dlhw,
                                   lhwhite->vertices[lhwvtx].x,
                                   lhwhite->vertices[lhwvtx].y,
                                   lhwhite->vertices[lhwvtx].z);
          if (lhpvtx > 0) printf("lhp  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   lhpvtx, dlhp,
                                   lhpial->vertices[lhpvtx].x,
                                   lhpial->vertices[lhpvtx].y,
                                   lhpial->vertices[lhpvtx].z);
          if (rhwvtx > 0) printf("rhw  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   rhwvtx, drhw,
                                   rhwhite->vertices[rhwvtx].x,
                                   rhwhite->vertices[rhwvtx].y,
                                   rhwhite->vertices[rhwvtx].z);
          if (rhpvtx > 0) printf("rhp  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   rhpvtx, drhp,
                                   rhpial->vertices[rhpvtx].x,
                                   rhpial->vertices[rhpvtx].y,
                                   rhpial->vertices[rhpvtx].z);
          printf("annot = %d, annotid = %d\n",annot,annotid);
          CTABprintASCII(lhwhite->ct,stdout);
          continue;
        }

        nctx++;
      }
    }
  }
  printf("nctx = %d\n",nctx);
  printf("Used brute-force search on %d voxels\n",nbrute);

  if (FixParaHipWM)
  {
    /* This is a bit of a hack. There are some vertices that have been
       ripped because they are "unkown". When the above alorithm finds
       these, it searches for the closest known vertex. If this is
       less than dmax away, then the wm voxel gets labeled
       accordingly.  However, there are often some voxels near
       ventralDC that are just close enough in 3d space to parahip to
       get labeled even though they are very far away along the
       surface. These voxels end up forming an island. CCSegment()
       will eliminate any islands. Unforunately, CCSegment() uses
       6-neighbor (face) definition of connectedness, so some voxels
       may be eliminated.
     */
    printf("Fixing Parahip LH WM\n");
    CCSegment(ASeg, 3016, 5001); //3016 = lhphwm, 5001 = unsegmented WM left
    printf("Fixing Parahip RH WM\n");
    CCSegment(ASeg, 4016, 5002); //4016 = rhphwm, 5002 = unsegmented WM right
  }

  printf("Writing output aseg to %s\n",OutASegFile);
  MRIwrite(ASeg,OutASegFile);

  if (OutAParcFile != NULL)
  {
    printf("Writing output aparc to %s\n",OutAParcFile);
    MRIwrite(AParc,OutAParcFile);
  }
  if (OutDistFile != NULL)
  {
    printf("Writing output dist file to %s\n",OutDistFile);
    MRIwrite(Dist,OutDistFile);
  }


  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1)
  {
    usage_exit();
  }

  nargc   = argc;
  pargv = argv;
  while (nargc > 0)
  {
    option = pargv[0];
    if (debug)
    {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help")||
        !strcasecmp(option, "-h")||
        !strcasecmp(option, "--usage")||
        !strcasecmp(option, "-u"))
    {
      print_help() ;
    }
    else if (!strcasecmp(option, "--version"))
    {
      print_version() ;
    }
    else if (!strcasecmp(option, "--debug"))
    {
      debug = 1;
    }
    // This was --ribbon, but changed to --old-ribbon 4/17/08 DNG
    else if (!strcasecmp(option, "--old-ribbon"))
    {
      UseRibbon = 1;
      UseNewRibbon = 0;
    }
    else if (!strcasecmp(option, "--volmask") ||
             !strcasecmp(option, "--new-ribbon"))
    {
      UseNewRibbon = 1;
    }
    else if (!strcasecmp(option, "--noribbon"))
    {
      UseRibbon = 0;
      UseNewRibbon = 0;
    }
    else if (!strcasecmp(option, "--labelwm"))
    {
      LabelWM = 1;
    }
    else if (!strcasecmp(option, "--fix-parahipwm"))
    {
      FixParaHipWM = 1;
    }
    else if (!strcasecmp(option, "--no-fix-parahipwm"))
    {
      FixParaHipWM = 0;
    }
    else if (!strcasecmp(option, "--hypo-as-wm"))
    {
      LabelHypoAsWM = 1;
    }
    else if (!strcasecmp(option, "--rip-unknown"))
    {
      RipUnknown = 1;
    }
    else if (!strcasecmp(option, "--no-hash"))
    {
      UseHash = 0;
    }
    else if (!strcmp(option, "--sd"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      SUBJECTS_DIR = pargv[0];
      setenv("SUBJECTS_DIR",SUBJECTS_DIR,1);
      nargsused = 1;
    }
    else if (!strcmp(option, "--s"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--oaseg") || !strcmp(option, "--o"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      OutASegFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--a2005s"))
    {
      annotname = "aparc.a2005s";
      baseoffset = 100;
    }
    else if (!strcmp(option, "--a2009s"))
    {
      annotname = "aparc.a2009s";
      baseoffset = 10100;
    }
   else if (!strcmp(option, "--aseg"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      asegname = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--annot"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      annotname = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--annot-table"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      // annotation_table_file is declared in annotation.h
      // default is $FREESURFER_HOME/Simple_surface_labels2009.txt
      annotation_table_file = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--oaparc"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      OutAParcFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--ctxseg"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      CtxSegFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--dist"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      OutDistFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--hashres"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%f",&hashres);
      nargsused = 1;
    }
    else if (!strcmp(option, "--wmparc-dmax"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%f",&dmaxctx);
      nargsused = 1;
    }
    else if (!strcmp(option, "--crs-test"))
    {
      if (nargc < 3)
      {
        argnerr(option,3);
      }
      sscanf(pargv[0],"%d",&ctest);
      sscanf(pargv[1],"%d",&rtest);
      sscanf(pargv[2],"%d",&stest);
      crsTest = 1;
      nargsused = 3;
    }
    else
    {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
      {
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      }
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
#include "mri_aparc2aseg.help.xml.h"
static void print_usage(void)
{
  outputHelpXml(mri_aparc2aseg_help_xml,mri_aparc2aseg_help_xml_len);
}
/* --------------------------------------------- */
static void print_help(void)
{
  outputHelpXml(mri_aparc2aseg_help_xml,mri_aparc2aseg_help_xml_len);
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1)
  {
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  }
  else
  {
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  }
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
  if (OutASegFile == NULL)
  {
    sprintf(tmpstr,"%s/%s/mri/%s+aseg.mgz",SUBJECTS_DIR,subject,annotname);
    OutASegFile = strcpyalloc(tmpstr);
  }
  if (UseRibbon && UseNewRibbon)
  {
    printf("ERROR: cannot --old-ribbon and --new-ribbon\n");
    exit(1);
  }
  if (CtxSegFile && ! RipUnknown)
  {
    printf("ERROR: can only use --ctxseg with --rip-unknown\n");
    exit(1);
  }
  if (CtxSegFile)
  {
    if (!fio_FileExistsReadable(CtxSegFile))
    {
      sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,CtxSegFile);
      if (! fio_FileExistsReadable(tmpstr))
      {
        printf("ERROR: cannot find %s or %s\n",CtxSegFile,tmpstr);
        exit(1);
      }
      CtxSegFile = strcpyalloc(tmpstr);
    }
  }
  if (FixParaHipWM && ! LabelWM)
  {
    FixParaHipWM  = 0;
  }

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"subject %s\n",subject);
  fprintf(fp,"outvol %s\n",OutASegFile);
  fprintf(fp,"useribbon %d\n",UseRibbon);
  fprintf(fp,"baseoffset %d\n",baseoffset);
  if (LabelWM)
  {
    printf("labeling wm\n");
    if (LabelHypoAsWM)
    {
      printf("labeling hypo-intensities as wm\n");
    }
    printf("dmaxctx %f\n",dmaxctx);
  }
  fprintf(fp,"RipUnknown %d\n",RipUnknown);
  if (CtxSegFile)
  {
    fprintf(fp,"CtxSeg %s\n",CtxSegFile);
  }
  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if (len < 2)
  {
    return(0);
  }

  if (flag[0] == '-' && flag[1] != '-')
  {
    return(1);
  }
  return(0);
}

/*---------------------------------------------------------------*/
int FindClosestLRWPVertexNo(int c, int r, int s,
                            int *lhwvtx, int *lhpvtx,
                            int *rhwvtx, int *rhpvtx,
                            MATRIX *Vox2RAS,
                            MRIS *lhwite,  MRIS *lhpial,
                            MRIS *rhwhite, MRIS *rhpial,
                            MHT *lhwhite_hash, MHT *lhpial_hash,
                            MHT *rhwhite_hash, MHT *rhpial_hash)
{
  static MATRIX *CRS = NULL;
  static MATRIX *RAS = NULL;
  static VERTEX vtx;
  static float dlhw, dlhp, drhw, drhp,dmin;
  int annot, hemi, annotid;

  if (CRS == NULL)
  {
    CRS = MatrixAlloc(4,1,MATRIX_REAL);
    CRS->rptr[4][1] = 1;
    RAS = MatrixAlloc(4,1,MATRIX_REAL);
    RAS->rptr[4][1] = 1;
  }

  CRS->rptr[1][1] = c;
  CRS->rptr[2][1] = r;
  CRS->rptr[3][1] = s;
  RAS = MatrixMultiply(Vox2RAS,CRS,RAS);
  vtx.x = RAS->rptr[1][1];
  vtx.y = RAS->rptr[2][1];
  vtx.z = RAS->rptr[3][1];

  *lhwvtx = MHTfindClosestVertexNo(lhwhite_hash,lhwhite,&vtx,&dlhw);
  *lhpvtx = MHTfindClosestVertexNo(lhpial_hash, lhpial, &vtx,&dlhp);
  *rhwvtx = MHTfindClosestVertexNo(rhwhite_hash,rhwhite,&vtx,&drhw);
  *rhpvtx = MHTfindClosestVertexNo(rhpial_hash, rhpial, &vtx,&drhp);

  printf("lh white: %d %g\n",*lhwvtx,dlhw);
  printf("lh pial:  %d %g\n",*lhpvtx,dlhp);
  printf("rh white: %d %g\n",*rhwvtx,drhw);
  printf("rh pial:  %d %g\n",*rhpvtx,drhp);

  hemi = 0;
  dmin = -1;
  if (*lhwvtx < 0 && *lhpvtx < 0 && *rhwvtx < 0 && *rhpvtx < 0)
  {
    printf("ERROR2: could not map to any surface.\n");
    printf("crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
           c,r,s,vtx.x,vtx.y,vtx.z);
    printf("Using Bruce Force\n");
    *lhwvtx = MRISfindClosestVertex(lhwhite,vtx.x,vtx.y,vtx.z,&dlhw);
    *lhpvtx = MRISfindClosestVertex(lhpial,vtx.x,vtx.y,vtx.z,&dlhp);
    *rhwvtx = MRISfindClosestVertex(rhwhite,vtx.x,vtx.y,vtx.z,&drhw);
    *rhpvtx = MRISfindClosestVertex(rhpial,vtx.x,vtx.y,vtx.z,&drhp);
    printf("lh white: %d %g\n",*lhwvtx,dlhw);
    printf("lh pial:  %d %g\n",*lhpvtx,dlhp);
    printf("rh white: %d %g\n",*rhwvtx,drhw);
    printf("rh pial:  %d %g\n",*rhpvtx,drhp);
    return(1);
  }
  if (dlhw <= dlhp && dlhw < drhw && dlhw < drhp && lhwvtx >= 0)
  {
    annot = lhwhite->vertices[*lhwvtx].annotation;
    hemi = 1;
    if (lhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = dlhw;
  }
  if (dlhp < dlhw && dlhp < drhw && dlhp < drhp && lhpvtx >= 0)
  {
    annot = lhwhite->vertices[*lhpvtx].annotation;
    hemi = 1;
    if (lhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = dlhp;
  }

  if (drhw < dlhp && drhw < dlhw && drhw <= drhp && rhwvtx >= 0)
  {
    annot = rhwhite->vertices[*rhwvtx].annotation;
    hemi = 2;
    if (rhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = drhw;
  }
  if (drhp < dlhp && drhp < drhw && drhp < dlhw && rhpvtx >= 0)
  {
    annot = rhwhite->vertices[*rhpvtx].annotation;
    hemi = 2;
    if (rhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = drhp;
  }

  printf("hemi = %d, annotid = %d, dist = %g\n",hemi,annotid,dmin);

  return(0);
}

/*!
  \fn int CCSegment(MRI *seg, int segid, int segidunknown)
  Constraints a sementation ID to consist of voxels that
  are spatially contiguous (6 face neighbors, not edge
  or corner). The voxels in the largest cluster are not
  changed. The voxels in the other clusters are set to
  segidunknown.
*/
int CCSegment(MRI *seg, int segid, int segidunknown)
{
  MRI_SEGMENTATION *sgmnt;
  int k,kmax,index,c,r,s;

  sgmnt = MRIsegment(seg,segid-.5,segid+.5);
  printf("  Found %d clusters\n",sgmnt->nsegments);

  kmax = 0;
  for (k=0; k < sgmnt->nsegments; k++)
    if (sgmnt->segments[k].nvoxels > sgmnt->segments[kmax].nvoxels)
    {
      kmax = k;
    }

  for (k=0; k < sgmnt->nsegments; k++)
  {
    printf("     %d k %f\n",k,sgmnt->segments[k].area);
    if (k==kmax)
    {
      continue;
    }
    for (index = 0; index < sgmnt->segments[k].nvoxels; index++)
    {
      c = sgmnt->segments[k].voxels[index].x;
      r = sgmnt->segments[k].voxels[index].y;
      s = sgmnt->segments[k].voxels[index].z;
      MRIsetVoxVal(seg,c,r,s,0,segidunknown);
    }
  }
  MRIsegmentFree(&sgmnt);
  return(0);
}
