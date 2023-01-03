#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED

/*
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_io.h"
#include "mrisp.h"

#include "mrisurf_base.h"

#include "mrisutils.h"


#define QUAD_FILE_MAGIC_NUMBER (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER (-3 & 0x00ffffff)

#define START_Y (-128)
#define SLICE_THICKNESS 1

static int mrisReadGeoFilePositions     (MRI_SURFACE *mris, const char *fname);
static int mrisReadTriangleFilePositions(MRI_SURFACE *mris, const char *fname);
static MRI_SURFACE *mrisReadTriangleFile(const char *fname, double nVFMultiplier);
static SMALL_SURFACE *mrisReadTriangleFileVertexPositionsOnly(const char *fname);


static int mris_readval_frame = -1;

/*--------------------------------------------------------
  MRISsetReadFrame() - sets frame to read when loading a "volume"
  with MRISreadValues().
  --------------------------------------------------------*/
void MRISsetReadFrame(int frame) { mris_readval_frame = frame; }

/*--------------------------------------------------------
  MRISgetReadFrame() - gets frame to read when loading a "volume"
  with MRISreadValues().
  --------------------------------------------------------*/
int MRISgetReadFrame(void)
{
  int frame = 0;
  char *envframe;
  if (mris_readval_frame >= 0) {
    return (mris_readval_frame);
  }
  envframe = getenv("MRIS_READVAL_FRAME");
  if (envframe == NULL) {
    return (0);
  }
  sscanf(envframe, "%d", &frame);
  return (frame);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadTriangleProperties(MRI_SURFACE *mris, const char *mris_fname)
{
  int ano, vnum, fnum, fno, vno;
  FACE *face;

  float f;
  FILE *fp;
  char fname[STRLEN], fpref[STRLEN], hemi[20];
  const char *cp;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "reading triangle files...");
  }

  cp = strrchr(mris_fname, '.');
  if (!cp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "MRISreadTriangleProperties(%s): could not scan"
                 " hemisphere from fname",
                 mris_fname));
  strncpy(hemi, cp - 2, 2);
  hemi[2] = 0;
  FileNamePath(mris_fname, fpref);

  int req = snprintf(fname, STRLEN, "%s/%s.triangle_area", fpref, hemi);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }


  fp = fopen(fname, "r");
  if (fp == NULL) {
    fprintf(stdout, "\nno precomputed triangle areas and angles - computing...\n");
    return (1); /* doesn't exist */
  }

  fread4(&f, fp);
  vnum = (int)f;
  fread4(&f, fp);
  fnum = (int)f;
  if (vnum != mris->nvertices) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadTriangleProperties: incompatible vertex "
                 "number in file %s",
                 fname));
  }

  mris->orig_area = 0.0f;
  for (fno = 0; fno < fnum; fno++) {
    face = &mris->faces[fno];
    f = freadFloat(fp);
    setFaceOrigArea(mris, fno, f);
    mris->orig_area += f;
  }

  /* compute original vertex areas from faces */
  for (vno = 0; vno < vnum; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    v->origarea = 0.0f;
    for (fno = 0; fno < vt->num; fno++) {
      v->origarea += getFaceOrigArea(mris, vt->f[fno]);
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "total area = %2.0f.\n", mris->orig_area);
  }

  /* now open and read the angle file */
  req = snprintf(fname, STRLEN, "%s/%s.triangle_angle", fpref, hemi);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  fp = fopen(fname, "r");
  if (fp == NULL) {
    return (1); /* doesn't exist */
  }

  fread4(&f, fp);
  vnum = (int)f;
  fread4(&f, fp);
  fnum = (int)f;
  if (vnum != mris->nvertices) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadTriangleProperties: incompatible vertex "
                 "number in file %s",
                 fname));
  }

  for (fno = 0; fno < fnum; fno++) {
    face = &mris->faces[fno];
    for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
      f = freadFloat(fp);
      face->orig_angle[ano] = f;
    }
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteTriangleProperties(MRI_SURFACE *mris, const char *mris_fname)
{
  int fno, ano, vno;
  FACE *face;
  FILE *fp;
  char fname[STRLEN], fpref[STRLEN], hemi[20];
  const char *cp;

  MRIScomputeTriangleProperties(mris);

  cp = strrchr(mris_fname, '.');
  if (!cp)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "MRISwriteTriangleProperties(%s): could not scan"
                 "hemisphere from fname",
                 mris_fname));
  strncpy(hemi, cp - 2, 2);
  hemi[2] = 0;
  FileNamePath(mris_fname, fpref);

  int req = snprintf(fname, STRLEN, "%s/%s.triangle_area", fpref, hemi); 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  fp = fopen(fname, "wb");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRISwriteTriangleProperties: could not open %s", fname));

  /* write out the distances to all neighboring vertices */
  fwrite4(mris->nvertices, fp);
  fwrite4(mris->nfaces, fp);
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    fwriteFloat(face->area, fp);
  }

  fclose(fp);

  req = snprintf(fname, STRLEN, "%s/%s.triangle_angle", fpref, hemi);    
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  fp = fopen(fname, "wb");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRISwriteTriangleProperties: could not open %s", fname));

  /* write out the area of all the triangles */
  fwrite4(mris->nvertices, fp);
  fwrite4(mris->nfaces, fp);
  for (vno = 0; vno < mris->nfaces; vno++) {
    face = &mris->faces[fno];
    for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
      fwriteFloat(face->angle[ano], fp);
    }
  }

  fclose(fp);

  return (NO_ERROR);
}

/*!
\fn int MRISwriteCurvature(MRI_SURFACE *mris, const char *sname)
\brief Writes the curvature field to the give file. If the file
is a volume format (eg, mgz), then it will use MRIwrite(). Can
save in VTK as well. 
*/
int MRISwriteCurvature(MRI_SURFACE *mris, const char *sname)
{
  int k, mritype;
  float curv;
  char fname[STRLEN], path[STRLEN], name[STRLEN];
  const char *hemi;
  const char *cp;
  FILE *fp;

  switch (mris->hemisphere) {
    case LEFT_HEMISPHERE:
      hemi = "lh";
      break;
    case BOTH_HEMISPHERES:
      hemi = "both";
      break;
    case RIGHT_HEMISPHERE:
      hemi = "rh";
      break;
    default:
      hemi = "unknown";
      break;
  }
  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path);
    cp = strchr(sname, '.');
    if (!cp || ((cp - sname) != 2) || *(cp - 1) != 'h' || ((*(cp - 2) != 'l' && *(cp - 2) != 'r'))) {
      if (getenv("FS_POSIX")) {
        // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
        int req = snprintf(fname, STRLEN, "./%s.%s", hemi, sname);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
        int req = snprintf(fname, STRLEN, "%s/%s.%s", path, hemi, sname);  
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
    }
    else {
      if (getenv("FS_POSIX")) {
        // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
        int req = snprintf(fname, STRLEN, "./%s", sname); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
        int req = snprintf(fname, STRLEN, "%s/%s", path, sname);   
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
    }
  }
  else {
    FileNamePath(sname, path);
    FileNameOnly(sname, name);
    cp = strchr(name, '.');
    if (!cp || ((cp - name) != 2) || *(cp - 1) != 'h' || ((*(cp - 2) != 'l' && *(cp - 2) != 'r'))) {
      int req = snprintf(fname, STRLEN, "%s/%s.%s", path, hemi, name); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else {
      int req = snprintf(fname, STRLEN, "%s/%s", path, name); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "writing curvature file %s\n", fname);
  }

  mritype = mri_identify(sname);

  if (mritype == MRI_MGH_FILE) {
    int vno;
    MRI *TempMRI;
    VERTEX *v;

    TempMRI = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT);
    if (TempMRI == NULL) {
      return (ERROR_NOMEMORY);
    }
    vno = 0;
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      MRIsetVoxVal(TempMRI, vno, 0, 0, 0, v->curv);
    }

    MRIwrite(TempMRI, fname);
    MRIfree(&TempMRI);
    return (NO_ERROR);
  }
  else if (mritype == VTK_FILE) {
    MRISwriteVTK(mris, fname);
    MRISwriteCurvVTK(mris, fname);
    return (NO_ERROR);
  }

  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteCurvature: could not open %s", fname));

  fwrite3(-1, fp); /* same old trick - mark it as new format */
  fwriteInt(mris->nvertices, fp);
  fwriteInt(mris->nfaces, fp);
  fwriteInt(1, fp); /* 1 value per vertex */
  for (k = 0; k < mris->nvertices; k++) {
    curv = mris->vertices[k].curv;
    fwriteFloat(curv, fp);
  }
  fclose(fp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteDists(MRI_SURFACE *mris, const char *sname)
{
  int k, i;
  float dist;
  char fname[STRLEN], path[STRLEN];
  const char *cp;
  FILE *fp;

  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path);
    cp = strchr(sname, '.');
    if (!cp)
      if (getenv("FS_POSIX")) {
        // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
        int req = snprintf(fname, STRLEN, "./%s.%s", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname);  
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
        // PW 2017/05/15: Legacy behaviour
        int req = snprintf(fname, STRLEN, "%s/%s.%s",
			   path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
    else {
      if (getenv("FS_POSIX")) {
        // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
        int req = snprintf(fname, STRLEN, "./%s", sname); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
        // PW 2017/05/15: Legacy behaviour
        int req = snprintf(fname, STRLEN, "%s/%s", path, sname);  
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
    }
  }
  else {
    strcpy(fname, sname); /* path specified explcitly */
  }
  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteDists: could not open %s", fname));

  fwrite3(mris->nvertices, fp);
  fwrite3(mris->nfaces, fp);
  for (k = 0; k < mris->nvertices; k++) {
    dist = mris->vertices[k].d;
    i = nint(dist * 100.0);
    fwrite2((int)i, fp);
  }
  fclose(fp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *MRISwriteIntoVolume(MRI_SURFACE *mris, MRI *mri, int which)
{
  int vno;
  VERTEX *v;
  float val = 0;

  if (mri == NULL) {
    mri = MRIalloc(1, 1, mris->nvertices, MRI_FLOAT);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    switch (which) {
      case VERTEX_DX:
        val = v->dx;
        break;
      case VERTEX_DY:
        val = v->dy;
        break;
      case VERTEX_DZ:
        val = v->dz;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VAL:
        val = v->val;
        break;
      case VERTEX_CURV:
        val = v->curv;
        break;
      case VERTEX_ANNOTATION:
        val = v->annotation;
        break;
      case VERTEX_AREA:
        val = v->area;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRISwriteIntoVolume: unsupported type %d", which);
        break;
    }
    if (mri->width == mris->nvertices)
      MRIsetVoxVal(mri, vno, 0, 0, 0, val);
    else if (mri->nframes == mris->nvertices)
      MRIsetVoxVal(mri, 0, 0, 0, vno, val);
    else {
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "MRISwriteIntoVolume: nvertices %d doesn't match with in MRI %d x %d x %d x %d\n",
                   mris->nvertices,
                   mri->width,
                   mri->height,
                   mri->depth,
                   mri->nframes));
    }
  }
  return (mri);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI_SURFACE *MRISreadFromVolume(MRI *mri, MRI_SURFACE *mris, int which)
{
  int vno;
  VERTEX *v;
  float val;

  if (mris->nvertices != mri->width)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MRISreadFromVolume: nvertices %d doesn't match with in MRI %d x %d x %d x %d\n",
                 mris->nvertices,
                 mri->width,
                 mri->height,
                 mri->depth,
                 mri->nframes));
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;
    val = MRIgetVoxVal(mri, vno, 0, 0, 0);
    switch (which) {
      case VERTEX_DX:
        v->dx = val;
        break;
      case VERTEX_DY:
        v->dy = val;
        break;
      case VERTEX_DZ:
        v->dz = val;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VAL:
        v->val = val;
        break;
      case VERTEX_CURV:
        v->curv = val;
        break;
      case VERTEX_ANNOTATION:
        v->annotation = val;
        break;
      case VERTEX_AREA:
        v->area = val;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRISreadFromVolume: unsupported type %d", which);
        break;
    }
  }
  return (mris);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteAreaError(MRI_SURFACE *mris, const char *name)
{
  int vno, fi, i;
  float area, orig_area;
  FACE *face;
  FILE *fp;
  char fname[STRLEN];

  strcpy(fname, name);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "writing area error file %s...", fname);
  }

  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteAreaError: could not open %s", fname));

  fwrite3(mris->nvertices, fp);
  fwrite3(mris->nfaces, fp);
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertex = &mris->vertices_topology[vno];
    area = orig_area = 0.0f;

    /* use average area of all faces this vertex is part of -
       this is not really correct, but should be good enough for
       visualization purposes.
    */
    for (fi = 0; fi < vertex->num; fi++) {
      face = &mris->faces[vertex->f[fi]];
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, vertex->f[fi]);
      area += face->area;
      orig_area += fNorm->orig_area;
    }
    i = nint((area - orig_area) * 100.0f / (float)(vertex->num));
    fwrite2((int)i, fp);
  }
  fclose(fp);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "done.\n");
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteAreaErrorToValFile(MRI_SURFACE *mris, const char *name)
{
  int vno, fno;
  float area, orig_area;
  FACE *face;

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    area = orig_area = 0.0f;

    /* use average area of all faces this vertex is part of -
       this is not really correct, but should be good enough for
       visualization purposes.
    */
    for (fno = 0; fno < vt->num; fno++) {
      face = &mris->faces[vt->f[fno]];
      area += face->area;
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
      orig_area += fNorm->orig_area;
    }
    v->val = (area - orig_area) / (float)(vt->num);
  }

  MRISwriteValues(mris, name);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteAngleError(MRI_SURFACE *mris, const char *fname)
{
  int vno, fno, ano, i;
  float error;
  FILE *fp;
  FACE *face;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "writing angular error file %s...", fname);
  }

  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteAngleError: could not open %s", fname));

  fwrite3(mris->nvertices, fp);
  fwrite3(mris->nfaces, fp);
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];
    error = 0.0f;
    for (fno = 0; fno < v->num; fno++) {
      face = &mris->faces[v->f[fno]];
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        error += fabs(deltaAngle(face->angle[ano], face->orig_angle[ano]));
      }
      error /= (float)(v->num * ANGLES_PER_TRIANGLE);
    }
    i = DEGREES(error) * 100.0;
    fwrite2((int)i, fp);
  }
  fclose(fp);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "done.\n");
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteCurvatureToWFile(MRI_SURFACE *mris, const char *fname)
{
  int k, num; /* loop counters */
  float f;
  FILE *fp;
  double sum = 0, sum2 = 0, max = -1000, min = 1000;

  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "writing out surface values to %s.\n", fname);
  }

  fp = fopen(fname, "wb");
  if (fp == NULL) {
    ErrorExit(ERROR_NOFILE, "Can't create file %s\n", fname);
  }

  num = mris->nvertices;
  fwrite2(0, fp);
  fwrite3(num, fp);
  for (k = 0; k < mris->nvertices; k++) {
    fwrite3(k, fp);
    f = mris->vertices[k].curv;
    if (!std::isfinite(f))
      ErrorPrintf(ERROR_BADPARM,
                  "MRISwriteCurvatureToWFile(%s): val at vertex %d is not"
                  "finite",
                  fname,
                  k);

    fwriteFloat(f, fp);
    sum += f;
    sum2 += f * f;
    if (f > max) {
      max = f;
    }
    if (f < min) {
      min = f;
    }
  }
  fclose(fp);
  sum /= num;
  sum2 = sqrt(sum2 / num - sum * sum);
  if (Gdiag & DIAG_SHOW) fprintf(stdout, "avg = %2.3f, stdev = %2.3f, min = %2.3f, max = %2.3f\n", sum, sum2, min, max);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteValues(MRI_SURFACE *mris, const char *sname)
{
  int k, num; /* loop counters */
  float f;
  char fname[STRLEN], *cp;
  FILE *fp;
  double sum = 0, sum2 = 0, max = -1000, min = 1000;
  int ftype, err;
  MRI *TempMRI;

  // Try saving it in a "volume" format -- but not img or nifti
  // as they use shorts for the number of vertices. Should add
  // a reshape.
  ftype = mri_identify(sname);
  if (ftype != MRI_VOLUME_TYPE_UNKNOWN) {
    TempMRI = MRIcopyMRIS(NULL, mris, 0, "val");
    if (TempMRI == NULL) {
      printf("ERROR: MRISwriteValues: could not alloc MRI\n");
      return (1);
    }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      printf("Saving surf vals in 'volume' format to %s\n", sname);
    }
    err = MRIwrite(TempMRI, sname);
    return (err);
  }

  strcpy(fname, sname);

  cp = strrchr(fname, '.');
  if (!cp || *(cp + 1) != 'w') {
    strcat(fname, ".w");
  }
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "writing surface values to %s.\n", fname);
  }

  fp = fopen(fname, "wb");
  if (fp == NULL) {
    ErrorExit(ERROR_NOFILE, "Can't create file %s\n", fname);
  }

  for (k = 0, num = 0; k < mris->nvertices; k++)
    if (mris->vertices[k].val != 0) {
      num++;
    }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("num = %d\n", num);
  }
  fwrite2(0, fp);
  fwrite3(num, fp);

  for (k = 0; k < mris->nvertices; k++) {
    if (mris->vertices[k].val != 0) {
      fwrite3(k, fp);
      f = mris->vertices[k].val;
      if (!std::isfinite(f)) ErrorPrintf(ERROR_BADPARM, "MRISwriteValues(%s): val at vertex %d is not finite", fname, k);

      fwriteFloat(f, fp);
      sum += f;
      sum2 += f * f;
      if (f > max) {
        max = f;
      }
      if (f < min) {
        min = f;
      }
    }
  }
  fclose(fp);

  if (num > 0) {
    sum /= num;
    sum2 = (sum2 / num - sum * sum);
    if (sum2 > 0) {
      sum2 = sqrt(sum2);
    }
    else {
      sum2 = 0;
    }
    printf("avg = %2.3f, stdev = %2.3f, min = %2.3f, max = %2.3f\n", sum, sum2, min, max);
  }
  else {
    printf("Warning: all vertex values are zero\n");
  }

  return (NO_ERROR);
}

int MRISwriteD(MRI_SURFACE *mris, const char *sname)
{
  float *curv_save;

  curv_save = (float *)calloc(mris->nvertices, sizeof(float));
  if (!curv_save) ErrorExit(ERROR_NOMEMORY, "MRISwriteD: could not alloc %d vertex curv storage", mris->nvertices);

  MRISextractCurvatureVector(mris, curv_save);
  MRISdToCurv(mris);
  MRISwriteCurvature(mris, sname);
  MRISimportCurvatureVector(mris, curv_save);
  free(curv_save);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadInflatedCoordinates(MRI_SURFACE *mris, const char *sname)
{
  if (!sname) {
    sname = "inflated";
  }
  MRISsaveVertexPositions(mris, TMP_VERTICES);
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR) {
    return (Gerror);
  }
  MRISsaveVertexPositions(mris, INFLATED_VERTICES);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadFlattenedCoordinates(MRI_SURFACE *mris, const char *sname)
{
  int vno, fno;
  VERTEX *v;
  FACE *f;

  if (!sname) {
    sname = "patch";
  }
  MRISsaveVertexPositions(mris, TMP_VERTICES);
  if (MRISreadPatchNoRemove(mris, sname) != NO_ERROR) {
    return (Gerror);
  }
  MRISsaveVertexPositions(mris, FLATTENED_VERTICES);
  
  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      MRISsetXYZ(mris,vno, v->x, v->y, -10000);
      v->ripflag = 0;
    }
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->ripflag) {
      f->ripflag = 0;
    }
  }
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  MRIScomputeMetricProperties(mris);

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadWhiteCoordinates(MRI_SURFACE *mris, const char *sname)
{
  if (!sname) {
    sname = "white";
  }
  MRISsaveVertexPositions(mris, TMP_VERTICES);
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR) {
    return (Gerror);
  }
  MRISsaveVertexPositions(mris, WHITE_VERTICES);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadPialCoordinates(MRI_SURFACE *mris, const char *sname)
{
  if (!sname) {
    sname = "pial";
  }
  MRISsaveVertexPositions(mris, TMP_VERTICES);
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR) {
    return (Gerror);
  }
  MRISsaveVertexPositions(mris, PIAL_VERTICES);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadCanonicalCoordinates(MRI_SURFACE *mris, const char *sname)
{
  MRISsaveVertexPositions(mris, TMP_VERTICES);
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR) {
    return (Gerror);
  }
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  MRIScomputeCanonicalCoordinates(mris);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadPatchNoRemove(MRI_SURFACE *mris, const char *pname)
{
  char fname[STRLEN];

  MRISbuildFileName_read(mris, pname, fname);

  int const type = MRISfileNameType(fname); /* using extension to get type */

  int ix, iy, iz, k, i, j, npts;
  double rx, ry, rz;
  FILE *fp = NULL;
  char line[256];
  char *cp;

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  // check whether the patch file is ascii or binary
  if (type == MRIS_GIFTI_FILE)    /* .gii */
  {
    mris = MRISread(fname);
  }
  else if (type == MRIS_ASCII_TRIANGLE_FILE) /* .ASC */
  {
    fp = fopen(fname, "r");
    if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadPatchNoRemove(%s): could not open file", fname));
    cp = fgetl(line, 256, fp);    // this would skip # lines
    sscanf(cp, "%d %*s", &npts);  // get points
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout,
              "reading .ASC patch %s with %d vertices (%2.1f%% of total)\n",
              pname,
              npts,
              100.0f * (float)npts / (float)mris->nvertices);

    // set all vertices ripflag to be true
    for (k = 0; k < mris->nvertices; k++) {
      mris->vertices[k].ripflag = TRUE;
    }

    // go through points
    for (j = 0; j < npts; j++) {
      // read int
      if ((cp = fgetl(line, 256, fp))) {
        sscanf(cp, "%d %*s", &i);
      }
      else
        ErrorReturn(ERROR_BADPARM, (ERROR_BAD_PARM, "MRISreadPatchNoRemove(%s): could not read line for point %d\n", fname, j));

      // if negative, flip it
      if (i < 0) {
        k = -i - 1;  // convert it to zero based number
        if (k < 0 || k >= mris->nvertices)
          ErrorExit(ERROR_BADFILE, "MRISreadPatch: bad vertex # (%d) found in patch file", k);
        // negative -> set the border to be true
        mris->vertices[k].border = TRUE;
      }
      // if positive
      else {
        k = i - 1;  // vertex number is zero based
        if (k < 0 || k >= mris->nvertices)
          ErrorExit(ERROR_BADFILE, "MRISreadPatch: bad vertex # (%d) found in patch file", k);
        // positive -> set the border to be false
        mris->vertices[k].border = FALSE;
      }
      // rip flag for this vertex to be false
      mris->vertices[k].ripflag = FALSE;
      // read 3 positions
      if ((cp = fgetl(line, 256, fp))) {
        sscanf(cp, "%lf %lf %lf", &rx, &ry, &rz);
      }
      else
        ErrorReturn(ERROR_BADPARM,
                    (ERROR_BADPARM,
                     "MRISreadPatch(%s): could not read 3 floating "
                     "values line for point %d\n",
                     fname,
                     j));

      // convert it to mm, i.e. change the vertex position
      MRISsetXYZ(mris,k,
        rx,
        ry,
        rz);
      if (k == Gdiag_no && Gdiag & DIAG_SHOW)
        fprintf(stdout,
                "vertex %d read @ (%2.2f, %2.2f, %2.2f)\n",
                k,
                mris->vertices[k].x,
                mris->vertices[k].y,
                mris->vertices[k].z);
    }
  }
  /////////////////////////////////////////////////////////////////////////
  // here file was binary
  /////////////////////////////////////////////////////////////////////////
  else {
    fp = fopen(fname, "rb");
    if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadPatchNoRemove(%s): could not open file", fname));

    // read number of vertices
    npts = freadInt(fp);
    if (npts >= 0) {
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout,
                "reading binary patch %s with %d vertices (%2.1f%% of total)\n",
                pname,
                npts,
                100.0f * (float)npts / (float)mris->nvertices);
      // set all vertices ripflag to be true
      for (k = 0; k < mris->nvertices; k++) {
        mris->vertices[k].ripflag = TRUE;
      }

      // go through points
      for (j = 0; j < npts; j++) {
        // read int
        i = freadInt(fp);
        // if negative, flip it
        if (i < 0) {
          k = -i - 1;  // convert it to zero based number
          if (k < 0 || k >= mris->nvertices)
            ErrorExit(ERROR_BADFILE, "MRISreadPatch: bad vertex # (%d) found in patch file", k);
          // negative -> set the border to be true
          mris->vertices[k].border = TRUE;
        }
        // if positive
        else {
          k = i - 1;  // vertex number is zero based
          if (k < 0 || k >= mris->nvertices)
            ErrorExit(ERROR_BADFILE, "MRISreadPatch: bad vertex # (%d) found in patch file", k);
          // positive -> set the border to be false
          mris->vertices[k].border = FALSE;
        }
        // rip flag for this vertex to be false
        mris->vertices[k].ripflag = FALSE;
        // read 3 positions
        fread2(&ix, fp);
        fread2(&iy, fp);
        fread2(&iz, fp);
        // convert it to mm, i.e. change the vertex position
        MRISsetXYZ(mris,k,
          ix / 100.0,
          iy / 100.0,
          iz / 100.0);
        if (k == Gdiag_no && Gdiag & DIAG_SHOW)
          fprintf(stdout,
                  "vertex %d read @ (%2.2f, %2.2f, %2.2f)\n",
                  k,
                  mris->vertices[k].x,
                  mris->vertices[k].y,
                  mris->vertices[k].z);
      }
    }
    else  // new surface format
    {
      // read number of vertices
      npts = freadInt(fp);
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout,
                "reading new surface format patch %s with %d vertices (%2.1f%% of total)\n",
                pname,
                npts,
                100.0f * (float)npts / (float)mris->nvertices);
      // set all vertices ripflag to be true
      for (k = 0; k < mris->nvertices; k++) {
        mris->vertices[k].ripflag = TRUE;
      }

      // go through points
      for (j = 0; j < npts; j++) {
        // read int
        i = freadInt(fp);
        // if negative, flip it
        if (i < 0) {
          k = -i - 1;  // convert it to zero based number
          if (k < 0 || k >= mris->nvertices)
            ErrorExit(ERROR_BADFILE, "MRISreadPatch: bad vertex # (%d) found in patch file", k);
          // negative -> set the border to be true
          mris->vertices[k].border = TRUE;
        }
        // if positive
        else {
          k = i - 1;  // vertex number is zero based
          if (k < 0 || k >= mris->nvertices)
            ErrorExit(ERROR_BADFILE, "MRISreadPatch: bad vertex # (%d) found in patch file", k);
          // positive -> set the border to be false
          mris->vertices[k].border = FALSE;
        }
        // rip flag for this vertex to be false
        mris->vertices[k].ripflag = FALSE;
        // read 3 positions
        // convert it to mm, i.e. change the vertex position
        float x = freadFloat(fp);
        float y = freadFloat(fp);
        float z = freadFloat(fp);
        MRISsetXYZ(mris,k,x,y,z);
        if (k == Gdiag_no && Gdiag & DIAG_SHOW)
          fprintf(stdout,
                  "vertex %d read @ (%2.2f, %2.2f, %2.2f)\n",
                  k,
                  mris->vertices[k].x,
                  mris->vertices[k].y,
                  mris->vertices[k].z);
      }
    }
  }
  if (fp) {
    fclose(fp);
  }
  for (k = 0; k < mris->nvertices; k++)
    if (mris->vertices_topology[k].num == 0 || mris->vertices_topology[k].vnum == 0) {
      mris->vertices[k].ripflag = 1;
    }

  // remove ripflag set vertices
  MRISremoveRipped(mris);
  mrisComputeSurfaceDimensions(mris);
  // set the patch flag
  mris->patch = 1;
  mris->status = MRIS_CUT;
  // recalculate properties
  mrisComputeBoundaryNormals(mris);
  mrisSmoothBoundaryNormals(mris, 10);
  MRISupdateSurface(mris);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadPatch(MRI_SURFACE *mris, const char *pname)
{
  int ret;

  // update the vertices in patch file
  ret = MRISreadPatchNoRemove(mris, pname);
  if (ret != NO_ERROR) {
    return (ret);
  }
  // remove ripflag set vertices
  MRISremoveRipped(mris);

  {
    int k;
    for (k = 0; k < mris->nvertices; k++)
      if (mris->vertices_topology[k].num == 0 || mris->vertices_topology[k].vnum == 0) {
        mris->vertices[k].ripflag = 1;
      }
  }

  // recalculate properties
  MRISupdateSurface(mris);

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwritePatch(MRI_SURFACE *mris, const char *fname)
{
  int k, i, npts, type;
  float x, y, z;
  FILE *fp;

  type = MRISfileNameType(fname);
  if (type == MRIS_ASCII_TRIANGLE_FILE)  // extension is ASC
  {
    return (MRISwritePatchAscii(mris, fname));
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE)  // extension is GEO
  {
    return (MRISwriteGeo(mris, fname));
  }
  else if (type == MRIS_STL_FILE) {
    return (mrisWriteSTL(mris, fname));
  }

  // binary file write
  // count number of points
  npts = 0;
  for (k = 0; k < mris->nvertices; k++)
    if (!mris->vertices[k].ripflag) {
      npts++;
    }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "writing surface patch with npts=%d to %s\n", npts, fname);
  }
  // binary write
  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwritePatch: can't create file %s\n", fname));
  // write num points
  fwriteInt(-1, fp);  // "version" #
  fwriteInt(npts, fp);
  // go through all points
  for (k = 0; k < mris->nvertices; k++)
    if (!mris->vertices[k].ripflag) {
      i = (mris->vertices[k].border) ? (-(k + 1)) : (k + 1);
      fwriteInt(i, fp);
      x = mris->vertices[k].x;
      y = mris->vertices[k].y;
      z = mris->vertices[k].z;
      fwriteFloat(x, fp);
      fwriteFloat(y, fp);
      fwriteFloat(z, fp);
      /*
        printf("k=%d, i=%d\n",k,i);
      */
    }
  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
static int mrisLabelVertices(MRI_SURFACE *mris, float cx, float cy, float cz, int label, float radius)
{
  VERTEX *v;
  float xd, yd, zd, d;
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    xd = cx - v->x;
    yd = cy - v->y;
    zd = cz - v->z;
    d = sqrt(xd * xd + yd * yd + zd * zd);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
int MRISreadTetherFile(MRI_SURFACE *mris, const char *fname, float radius)
{
  int l;
  float cx, cy, cz;
  FILE *fp;
  char line[STRLEN], *cp;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadTetherFile(%s): could not open file", fname));

  /* first count # of labels in file */
  mris->nlabels = 0;
  while ((cp = fgetl(line, 199, fp)) != NULL) {
    mris->nlabels++;
  }
  mris->labels = (MRIS_AREA_LABEL *)calloc(mris->nlabels, sizeof(MRIS_AREA_LABEL));
  if (!mris->labels) ErrorExit(ERROR_NOMEMORY, "MRISreadTetherFile: could not allocate %d labels", mris->nlabels);

  for (l = 0; l < mris->nlabels; l++) {
    cp = fgetl(line, 199, fp);
    if (sscanf(cp, "%s %f %f %f", mris->labels[l].name, &cx, &cy, &cz) != 4) {
      fclose(fp);
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadTetherFile(%s): could not scan parms from %s", fname, line));
    }
    mris->labels[l].cx = cx;
    mris->labels[l].cy = cy;
    mris->labels[l].cz = cz;
    mrisLabelVertices(mris, cx, cy, cz, l, radius); // weirdly this is currently a no-op!
  }

  fclose(fp);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadAnnotation(MRI_SURFACE *mris, const char *sname)
{
  int vno, need_hemi;
  int return_code;
  char fname[STRLEN], path[STRLEN], fname_no_path[STRLEN];
  const char *cp;
  int *array;

  // first attempt to read as gifti file
  int mritype = mri_identify(sname);
  if (mritype == GIFTI_FILE) {
    mris = mrisReadGIFTIfile(sname, mris);
    if (mris) {
      return (NO_ERROR);
    }
    else {
      return (ERROR_BADFILE);
    }
  }
  // else fall-thru with default .annot processing...

  cp = strchr(sname, '/');
  if (!cp) {
    /* no path - use same one as mris was read from */
    FileNameOnly(sname, fname_no_path);
    cp = strstr(fname_no_path, ".annot");
    if (!cp) {
      strcat(fname_no_path, ".annot");
    }
      
    need_hemi = stricmp(fname_no_path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
      
    FileNamePath(mris->fname, path);
    if (!need_hemi) {
      int req = snprintf(fname, STRLEN, "%s/../label/%s", path, fname_no_path); 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else {
      /* no hemisphere specified */
      int req = snprintf(fname, STRLEN, "%s/../label/%s.%s", 
			 path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", fname_no_path); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }
  else {
    strcpy(fname, sname); /* full path specified */
    cp = strstr(fname, ".annot");
    if (!cp) {
      strcat(fname, ".annot");
    }
  }
  
  // As a last resort, just assume the sname is the path
  if (!fio_FileExistsReadable(fname) && fio_FileExistsReadable(sname)) {
    int req = snprintf(fname, STRLEN, "%s", sname);  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }

  /* Try to read it into an array. */
  return_code = MRISreadAnnotationIntoArray(fname, mris->nvertices, &array);
  if (NO_ERROR != return_code) {
    return return_code;
  }

  /* If we got an array, fill in our annotation values. */
  MRISclearAnnotations(mris);
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    mris->vertices[vno].annotation = array[vno];
  }

  /* Try to read in a color table. If we read one, it will be
     allocated, otherwise it will stay NULL. */
  return_code = MRISreadCTABFromAnnotationIfPresent(fname, &mris->ct);
  if (NO_ERROR != return_code) {
    return return_code;
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------

  Parameters: fname is the file name of the annotation. in_array_size
    is the size of the array to allocate, but also a 'hint' as to the
    number of vertices expected in the file. This will return the
    allocated array in *out_array.

  Returns value: Standard error code.

  Description: Opens an annotation file and attempts to read label
    indices from it. If successful, will allocate an array of ints,
    populate it, and return a pointer to it. Will return an error code
    if not succesful. This could be due to an invalid file, or an
    incorrect number of vertices in the annotation file.

 ------------------------------------------------------*/
int MRISreadAnnotationIntoArray(const char *fname, int in_array_size, int **out_array)
{
  int i, j, vno, num;
  FILE *fp;
  int *array = NULL;

  if (fname == NULL || out_array == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "Parameter was NULL."));
  if (in_array_size < 0) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "in_array_size was negative."));

  /* Initialize our array. Note that we do it to the size that we got
     passed in, which might be different than the number of values in
     the file. */
  array = (int *)calloc(in_array_size, sizeof(int));
  if (NULL == array)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "ERROR: Reading annot value file %s\n"
                 "  Couldn't allocate array of size %s\n",
                 fname,
                 in_array_size));

  /* Open the file. */
  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "could not read annot file %s", fname));

  /* First int is the number of elements. */
  num = freadInt(fp);

  /* For each one, read in a vno and an int for the annotation
     value. Check the vno. */
  for (j = 0; j < num; j++) {
    vno = freadInt(fp);
    i = freadInt(fp);
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    /* Check the index we read to make sure it's less than the size
       we're expecting. */
    if (vno >= in_array_size || vno < 0) {
      fprintf(stderr,
              "MRISreadAnnotationIntoArray: vertex index out of range: "
              "%d i=%8.8X, in_array_size=%d\n",
              vno,
              i,
              in_array_size);
      fprintf(stderr, "    annot file: %s\n", fname);
      
      static int vertexIndexOutOfRangeCounter = 0;
      if (++vertexIndexOutOfRangeCounter > 200000) {
        // this check prevents creating 100GB error files
        ErrorExit(ERROR_BADFILE,
                  "ERROR: Too many out-of-range vertex indices in "
                  "MRISreadAnnotationIntoArray!\n");
      }
    }
    else {
      array[vno] = i;
    }
  }

  fclose(fp);

  /* Set the outgoing array. */
  *out_array = array;

  return (NO_ERROR);
}

/*-----------------------------------------------------

  Parameters: fname is the file name of the annotation. This will
    return the allocated table in *out_table.

  Returns value: Standard error code. If the file is fine but no color
    table is present, it will return NO_ERROR, but out_table will not
    be set.

  Description: Opens an annotation file and skips to the tag section
    of the file, after the values. Looks for the TAG_OLD_COLORTABLE
    value, and if present, reads in a color table.

*/
int MRISreadCTABFromAnnotationIfPresent(const char *fname, COLOR_TABLE **out_table)
{
  int skip, j, num;
  FILE *fp;
  COLOR_TABLE *ctab = NULL;
  int tag;

  if (fname == NULL || out_table == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "Parameter was NULL."));

  /* Open the file. */
  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "could not read annot file %s", fname));

  /* First int is the number of elements. */
  num = freadInt(fp);

  /* Skip two ints per num, to the end of the values section, where
     the tags are. */
  for (j = 0; j < num; j++) {
    skip = freadInt(fp);
    skip = freadInt(fp);
  }

  /* Check for tags. Right now we only have one possibility for tags,
     but if we add in more tags, we'll have to skip past other tags
     and their data sections here. Of course, we don't have an easy
     way to determine the length of the data section for the tag. If
     we hit EOF here, there is no tag. */
  tag = freadInt(fp);
  if (feof(fp)) {
    fclose(fp);
    return ERROR_NONE;
  }

  if (TAG_OLD_COLORTABLE == tag) {
    /* We have a color table, read it with CTABreadFromBinary. If it
    fails, it will print its own error message. */
    if (DIAG_VERBOSE_ON)
      fprintf(stdout, "reading colortable from annotation file...\n");
    ctab = CTABreadFromBinary(fp);
    if ((NULL != ctab) && DIAG_VERBOSE_ON) fprintf(stdout, "colortable with %d entries read (originally %s)\n", ctab->nentries, ctab->fname);
  }

  fclose(fp);

  /* Return the table if we got one. */
  if (NULL != ctab) {
    *out_table = ctab;
  }

  return ERROR_NONE;
}

/*-----------------------------------------------------

  Parameters: fname is the file name of the annotation. On successful
    return, *present will be set to 1 if the CTAB has an annotation,
    and 0 if not.

  Returns value: Standard error code.

  Description: Opens an annotation file and skips to the tag section
    of the file, after the values. Looks for the TAG_OLD_COLORTABLE
    value, and if present, sets *present to 1, or 0 if not present.

*/
int MRISisCTABPresentInAnnotation(const char *fname, int *present)
{
  int skip, j, num;
  FILE *fp;
  int tag;

  if (fname == NULL || present == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "Parameter was NULL."));

  /* Open the file. */
  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "could not read annot file %s", fname));

  /* First int is the number of elements. */
  num = freadInt(fp);

  /* Skip two ints per num, to the end of the values section, where
     the tags are. */
  for (j = 0; j < num; j++) {
    skip = freadInt(fp);
    skip = freadInt(fp);
  }

  /* No tag found yet. */
  *present = 0;

  /* Check for tags. Right now we only have one possibility for tags,
     but if we add in more tags, we'll have to skip past other tags
     and their data sections here. Of course, we don't have an easy
     way to determine the length of the data section for the tag. If
     we hit EOF here, there is no tag. */
  tag = freadInt(fp);
  if (feof(fp)) {
    fclose(fp);
    return ERROR_NONE;
  }

  if (TAG_OLD_COLORTABLE == tag) {
    *present = 1;
  }

  fclose(fp);

  return ERROR_NONE;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteAnnotation(MRI_SURFACE *mris, const char *sname)
{
  int i, vno, need_hemi;
  FILE *fp;
  const char *cp;
  char fname[STRLEN], path[STRLEN], fname_no_path[STRLEN];

  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    FileNameOnly(sname, fname_no_path);
    cp = strstr(fname_no_path, ".annot");
    if (!cp) {
      strcat(fname_no_path, ".annot");
    }

    need_hemi = strncmp(fname_no_path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", 2);

    FileNamePath(mris->fname, path);

    if (!need_hemi) {
      if (getenv("FS_POSIX")) {
        // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
        int req = snprintf(fname, STRLEN, "./%s", fname_no_path);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
        // PW 2017/05/15: Legacy behaviour
        int req = snprintf(fname, STRLEN, "%s/../label/%s", path, fname_no_path); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
    }
    else { /* no hemisphere specified */
      if (getenv("FS_POSIX")) {
        // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
        int req = snprintf(fname, STRLEN, "./%s.%s",
			   mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", fname_no_path);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
        // PW 2017/05/15: Legacy behaviour
        int req = snprintf(fname, STRLEN, "%s/../label/%s.%s",
			   path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", fname_no_path); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
    }
  }
  else {
    strcpy(fname, sname); /* full path specified */
    cp = strstr(fname, ".annot");
    if (!cp) {
      strcat(fname, ".annot");
    }
  }

  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "could not write annot file %s", fname));
  fwriteInt(mris->nvertices, fp);
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    i = mris->vertices[vno].annotation;
    fwriteInt(vno, fp);
    i = fwriteInt(i, fp);
  }

  if (mris->ct) /* also write annotation in */
  {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      printf("writing colortable into annotation file...\n");
    }
    fwriteInt(TAG_OLD_COLORTABLE, fp);
    CTABwriteIntoBinary(mris->ct, fp);
  }

  fclose(fp);

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadValues(MRI_SURFACE *mris, const char *sname)
{
  float *array = NULL;
  int return_code = 0;
  int vno;

  /* Try to read an array of values from this file. If we get an
     error, pass the code back up. */
  return_code = MRISreadValuesIntoArray(sname, mris->nvertices, &array);
  if (NO_ERROR != return_code) {
    return return_code;
  }

  /* Read all the values from the array into our val field. */
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    mris->vertices[vno].val = array[vno];
  }

  /* MRISreadValuesIntoArray allocates its own array, so free it
     here. */
  free(array);

  return (NO_ERROR);
}

int MRISreadValuesIntoArray(const char *sname, int in_array_size, float **out_array)
{
  int i, k, num, ilat, vno;
  float f = 0;
  float lat;
  FILE *fp;
  char *cp, fname[STRLEN];
  int type, frame, nv, c, r, s;
  MRI *TempMRI;
  float *array;
  int return_code;

  if (sname == NULL || out_array == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "Parameter was NULL."));
  if (in_array_size < 0) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "in_array_size was negative."));

  /* First try to load it as a volume */
  strncpy(fname, sname, sizeof(fname)-1);
  type = mri_identify(fname);
  if (type != MRI_VOLUME_TYPE_UNKNOWN) {
    frame = MRISgetReadFrame();
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      printf("MRISreadValues: frame=%d\n", frame);
    }
    TempMRI = MRIreadHeader(fname, type);
    if (TempMRI == NULL) {
      return (ERROR_BADFILE);
    }
    if (TempMRI->nframes <= frame) {
      MRIfree(&TempMRI);
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "ERROR: attempted to read frame %d from %s\n"
                   "but this file only has %d frames.\n",
                   frame,
                   fname,
                   TempMRI->nframes));
    }
    nv = TempMRI->width * TempMRI->height * TempMRI->depth;
    MRIfree(&TempMRI);

    /* Make sure the size is what we're expecting. */
    if (nv != in_array_size) {
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "ERROR: Reading %s as volume encoded scalar file,"
                   " but size doesn't match expected size of %d\n",
                   sname,
                   in_array_size));
    }

    /* Read the MRI of values. */
    TempMRI = MRIread(fname);
    if (TempMRI == NULL) {
      return (ERROR_BADFILE);
    }

    /* Initialize our array. */
    array = (float *)calloc(in_array_size, sizeof(float));
    if (NULL == array)
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "ERROR: Reading file value file %s\n"
                   "  Couldn't allocate array of size %s\n",
                   fname,
                   nv));

    /* Read all the values into the array. */
    vno = 0;
    for (s = 0; s < TempMRI->depth; s++) {
      for (r = 0; r < TempMRI->height; r++) {
        for (c = 0; c < TempMRI->width; c++) {
          array[vno] = MRIgetVoxVal(TempMRI, c, r, s, frame);
          vno++;
        }
      }
    }
    MRIfree(&TempMRI);

    /* Set the outgoing array and size. */
    *out_array = array;

    return (NO_ERROR);
  }

  /* Next, try reading a curvature file. If we don't get an error, the
     array is filled out, so we can return now.*/
  return_code = MRISreadCurvatureIntoArray(sname, in_array_size, out_array);
  if (ERROR_NONE == return_code) {
    printf("reading values from curvature-format file...\n");
    return ERROR_NONE;
  }

  /* Now we try an old .w file. If the file name doesn't have .w on
     it, append it now. */
  strcpy(fname, sname);
  cp = strrchr(fname, '.');
  if (!cp || *(cp + 1) != 'w') {
    strcat(fname, ".w");
  }

  /* Try opening the file. */
  fp = fopen(fname, "rb");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadValuesIntoArray: File %s not found\n", fname));

  /* Read something... Seems to be ignored now. */
  fread2(&ilat, fp);
  lat = ilat / 10.0;

  /* Read the number of values. */
  if (fread3(&num, fp) < 1) {
    fclose(fp);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadValues(%s): could not read # of vertices", fname));
  }

  /* The .w format, being sparse, doesn't require the number of values
     in the file to match the expected amount (usually
     mris->vertices), so we don't test that here. We'll test each
     index invidiually later. */

  /* Initialize our array. Note that we do it to the size that we got
     passed in, which might be different than num. */
  array = (float *)calloc(in_array_size, sizeof(float));
  if (NULL == array) {
    fclose(fp);
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "ERROR: Reading file value file %s\n"
                 "  Couldn't allocate array of size %s\n",
                 sname,
                 num));
  }

  /* Read the values into the array. */
  for (i = 0; i < num; i++) {
    /* Read the value index. */
    if (fread3(&k, fp) < 1)
      ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadValues(%s): could not read %dth vno", fname, i));

    /* Check the index we read to make sure it's less than the size
    we're expecting. */
    if (k >= in_array_size || k < 0) {
      printf("MRISreadValues: vertex index out of range: %d f=%f\n", k, f);
      continue;
    }

    if (k == Gdiag_no) {
      DiagBreak();
    }

    /* Read the value. */
    f = freadFloat(fp);

    /* Set it in the array. */
    array[k] = f;
  }

  fclose(fp);

  /* Set the outgoing array. */
  *out_array = array;

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadValuesBak(MRI_SURFACE *mris, const char *fname)
{
  int i, k, num, ilat;
  float f;
  float lat;
  FILE *fp;

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadValuesBak: File %s not found\n", fname));
  fread2(&ilat, fp);
  lat = ilat / 10.0;

  for (k = 0; k < mris->nvertices; k++) {
    mris->vertices[k].valbak = 0;
  }
  fread3(&num, fp);
  for (i = 0; i < num; i++) {
    fread3(&k, fp);
    f = freadFloat(fp);
    if (k >= mris->nvertices || k < 0) {
      printf("MRISreadValuesBak: vertex index out of range: %d f=%f\n", k, f);
    }
    else {
      mris->vertices[k].valbak = f;
    }
  }
  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadImagValues(MRI_SURFACE *mris, const char *fname)
{
  int i, k, num, ilat;
  float f;
  float lat;
  FILE *fp;

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadImagValues: File %s not found\n", fname));
  fread2(&ilat, fp);
  lat = ilat / 10.0;

  for (k = 0; k < mris->nvertices; k++) {
    mris->vertices[k].imag_val = 0;
  }
  fread3(&num, fp);
  for (i = 0; i < num; i++) {
    fread3(&k, fp);
    f = freadFloat(fp);
    if (k >= mris->nvertices || k < 0) {
      printf("MRISreadImagValues: vertex index out of range: %d f=%f\n", k, f);
    }
    /*
      else if (mris->vertices[k].dist!=0)
      printf("MRISreadImagValues: subsample and data file mismatch\n");
    */
    else {
      mris->vertices[k].imag_val = f;
      /*      mris->vertices[k].dist=0;*/
    }
  }
  fclose(fp);
  return (NO_ERROR);
}
int MRIScopyMarksToAnnotation(MRI_SURFACE *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->annotation = v->marked;
  }
  return (NO_ERROR);
}
int MRISmaxMarked(MRI_SURFACE *mris)
{
  int vno, max_marked;
  VERTEX *v;

  max_marked = mris->vertices[0].marked;
  for (vno = 1; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    if (v->marked > max_marked) {
      max_marked = v->marked;
    };
  }
  return (max_marked);
}
int MRIScopyAnnotationsToMarkedIndex(MRI_SURFACE *mris)
{
  int vno, index;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    CTABfindIndexFromAnnotation(mris->ct, v->annotation, &index);
    if (index < 0)
      ErrorPrintf(ERROR_BADPARM,
                  "%s: could not find index for vno %d annotation %x\n",
                  "MRIScopyAnnotationsToMarkedIndex",
                  vno,
                  v->annotation);

    v->marked = index;
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisWriteSnapshots(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t)
{
  char base_name[STRLEN];

  MRISsaveVertexPositions(mris, TMP_VERTICES);

  strcpy(base_name, parms->base_name);

  MRISrestoreVertexPositions(mris, PIAL_VERTICES);
  int req = snprintf(parms->base_name, STRLEN, "%s_pial", base_name); 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  mrisWriteSnapshot(mris, parms, t);

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  req = snprintf(parms->base_name, STRLEN, "%s_white", base_name); 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  mrisWriteSnapshot(mris, parms, t);

  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  MRIScomputeMetricProperties(mris);

  strcpy(parms->base_name, base_name);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisWriteSnapshot(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int t)
{
  char fname[STRLEN], path[STRLEN], base_name[STRLEN], *cp;
  const char *hemi;

  switch (mris->hemisphere) {
    case LEFT_HEMISPHERE:
      hemi = "lh";
      break;
    case BOTH_HEMISPHERES:
      hemi = "both";
      break;
    case RIGHT_HEMISPHERE:
      hemi = "rh";
      break;
    default:
      hemi = "unknown";
      break;
  }
  FileNamePath(mris->fname, path);
  int req = snprintf(base_name, STRLEN, "%s/%s.%s", path, hemi, parms->base_name);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  if ((cp = strstr(base_name, ".geo")) != NULL) {
    *cp = 0;
    int req = snprintf(fname, STRLEN, "%s%4.4d.geo", base_name, t); 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    *cp = '.';
  }
  else {
    int req = snprintf(fname, STRLEN, "%s%4.4d", base_name, t);   
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }
#if 1
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "writing %s", fname);
  }
  if (mris->patch) {
    char fname2[STRLEN] ;
    MRISwrite(mris, fname);
    strcpy(fname2, fname) ;
    strcat(fname2, ".patch") ;
    printf("\nwriting patch to %s\n", fname2) ;
    MRISwritePatch(mris, fname2);
    strcpy(fname2, fname) ;
    strcat(fname2, ".curv") ;
    MRISdToCurv(mris) ;
    printf("\nwriting curvature to %s\n", fname2) ;
    MRISwriteCurvature(mris, fname2);
  }
  else {
    if (parms->l_thick_normal > 0 || parms->l_thick_parallel > 0 || parms->l_thick_min > 0 || parms->l_thick_spring) {
      int vno;
      MRI *mri_vector;
      float dx, dy, dz, xp, yp, zp, xw, yw, zw;

      mri_vector = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);

      for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX *v;

        v = &mris->vertices[vno];
        if (vno == Gdiag_no) DiagBreak();

        if (v->ripflag) continue;

        MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
        MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);

        dx = xp - xw;
        dy = yp - yw;
        dz = zp - zw;
        MRIsetVoxVal(mri_vector, vno, 0, 0, 0, dx);
        MRIsetVoxVal(mri_vector, vno, 0, 0, 1, dy);
        MRIsetVoxVal(mri_vector, vno, 0, 0, 2, dz);
      }
      int req = snprintf(fname, STRLEN, "%s%4.4d.mgz", base_name, t);
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing vector field to %s\n", fname);
      MRIwrite(mri_vector, fname);
      MRIfree(&mri_vector);
    }
    else
      MRISwrite(mris, fname);
  }

  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "done.\n");
    fflush(stderr);
  }
#endif

  if (mris->status == MRIS_PARAMETERIZED_SPHERE && DIAG_VERBOSE_ON) {
    MRI_SP *mrisp = (MRI_SP *)mris->vp;

    int req = snprintf(fname, STRLEN, "%s%4.4d.hipl", parms->base_name, t);  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "writing %s\n", fname);
    }
    MRIStoParameterization(mris, mrisp, 1, 0);
    MRISPwrite(mrisp, fname);
  }
  if (!FZERO(parms->l_area) && DIAG_VERBOSE_ON) {
    int req = snprintf(fname, STRLEN, "%s%4.4d.area_error", base_name, t); 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, " %s...", fname);
    }
    MRISwriteAreaError(mris, fname);
  }

  if (!FZERO(parms->l_corr) && DIAG_VERBOSE_ON) {
    int req = snprintf(fname, STRLEN, "%s%4.4d.angle_error", base_name, t); 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, " %s...", fname);
    }
    MRISwriteAngleError(mris, fname);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadVertexPositions(MRI_SURFACE *mris, const char *name)
{
  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  char fname[STRLEN];
  int vno, nvertices, nfaces, magic, version, tmp, ix, iy, iz, n, type;
  FILE *fp;

  MRISbuildFileName_read(mris, name, fname);

  type = MRISfileNameType(name);
  if (type == MRIS_GEO_TRIANGLE_FILE) {
    return (mrisReadGeoFilePositions(mris, fname));
  }
  else if (type == MRIS_ICO_FILE) {
    return (ICOreadVertexPositions(mris, fname, CURRENT_VERTICES));
  }
  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadVertexPosition(%s): could not open file %s", name, fname));

  fread3(&magic, fp);
  if (magic == QUAD_FILE_MAGIC_NUMBER) {
    version = -1;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      fprintf(stdout, "new surface file format\n");
    }
  }
  else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) {
    version = -2;
  }
  else if (magic == TRIANGLE_FILE_MAGIC_NUMBER) {
    fclose(fp);
    if (mrisReadTriangleFilePositions(mris, fname) != NO_ERROR) {
      ErrorReturn(Gerror, (Gerror, "mrisReadTriangleFile failed.\n"));
    }
    version = -3;
  }
  else {
    rewind(fp);
    version = 0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      printf("surfer: old surface file format\n");
    }
  }

  if (version >= -2) {
    fread3(&nvertices, fp);
    fread3(&nfaces, fp);
    nfaces *= 2;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stdout, "reading %d vertices and %d faces.\n", nvertices, nfaces);

    if (nvertices != mris->nvertices || nfaces != mris->nfaces) {
      fclose(fp);
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "\nERROR: \n      MRISreadVertexPositions(%s): surfaces differ. "
                   "Main: %d verts %d faces, %s: %d verts %d faces\n"
                   "Surfaces may be out-of-date\n\n",
                   fname,
                   mris->nvertices,
                   mris->nfaces,
                   name,
                   nvertices,
                   nfaces));
    }

    for (vno = 0; vno < nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
      if (version == -1) {
        fread2(&ix, fp);
        fread2(&iy, fp);
        fread2(&iz, fp);
        MRISsetXYZ(mris,vno,
          ix / 100.0,
          iy / 100.0,
          iz / 100.0);
      }
      else /* version == -2 */
      {
        float x = freadFloat(fp);
        float y = freadFloat(fp);
        float z = freadFloat(fp);
        MRISsetXYZ(mris,vno, x,y,z);
      }
      if (version == 0) /* old surface format */
      {
        fread1(&tmp, fp);
        for (n = 0; n < vertext->num; n++) {
          fread3(&tmp, fp);
        }
      }
    }

    fclose(fp);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  read in vertex positions from the original file,
  and compute the metric properties of that surface.
  Note we must change the status of the mris as the
  original surface has no externally defined orientation.
  ------------------------------------------------------*/
int MRISreadOriginalProperties(MRI_SURFACE *mris, const char *sname)
{
  if (!sname) {
    sname = "smoothwm";
  }

  MRISsaveVertexPositions(mris, TMP_VERTICES);

  auto const old_status = mris->status;
  mris->status = MRIS_PATCH; /* so no orientating will be done */
  if (MRISreadVertexPositions(mris, sname) != NO_ERROR)
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadOriginalProperties: could not read surface file %s", sname));

  MRISsetOriginalXYZfromXYZ(mris);
  MRIScomputeMetricProperties(mris);
  MRIScomputeTriangleProperties(mris);
  MRISstoreMetricProperties(mris);
  mris->status = old_status;
  
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  
  MRIScomputeMetricProperties(mris);
  MRIScomputeTriangleProperties(mris);
  mrisOrientSurface(mris);
  mris->orig_area = mris->total_area;

  return (NO_ERROR);
}

/*-----------------------------------------------------
  int MRISfileNameType()

  Parameters: char *fname

  Returns value: int

  Description: return filetype using the extension
  default is MRIS_BINARY_QUADRANGLE_FILE
  ------------------------------------------------------*/
int MRISfileNameType(const char *fname)
{
  int type;
  char *dot, ext[STRLEN], str[STRLEN], *idext;

  // First check whether it is a volume format
  idext = IDextensionFromName(fname);
  if (idext != NULL) {
    free(idext);
    return (MRIS_VOLUME_FILE);
  }

  FileNameOnly(fname, str); /* remove path */

  ext[0] = 0;
  dot = strrchr(str, '@'); /* forces the type of the file */
  if (dot) {
    *dot = 0; /* remove it from file name so that fopen will work */
    strcpy(ext, dot + 1);
  }
  else /* no explicit type - check extension */
  {
    dot = strrchr(str, '.');
    if (dot) {
      strcpy(ext, dot + 1);
    }
  }
  StrUpper(ext);
  if (!strcmp(ext, "ASC")) {
    type = MRIS_ASCII_TRIANGLE_FILE;
  }
  else if (!strcmp(ext, "GEO")) {
    type = MRIS_GEO_TRIANGLE_FILE;
  }
  else if (!strcmp(ext, "TRI") || !strcmp(ext, "ICO")) {
    type = MRIS_ICO_FILE;
  }
  else if (!strcmp(ext, "VTK")) {
    type = MRIS_VTK_FILE;
  }
  else if (!strcmp(ext, "STL")) {
    type = MRIS_STL_FILE;
  }
  else if (!strcmp(ext, "GII")) {
    type = MRIS_GIFTI_FILE;
  }
  else if (!strcmp(ext, "MGH") || !strcmp(ext, "MGZ")) {
    type = MRI_MGH_FILE;  // surface-encoded volume
  }
  else if (!strcmp(ext, "ANNOT")) {
    type = MRIS_ANNOT_FILE;
  }
  else {
    type = MRIS_BINARY_QUADRANGLE_FILE;
  }

  return (type);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteAscii(MRI_SURFACE *mris, const char *fname)
{
  int vno, fno, n;
  VERTEX *v;
  FACE *face;
  FILE *fp;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteAscii: could not open file %s", fname));

  fprintf(fp, "#!ascii version of %s\n", mris->fname.data());
  fprintf(fp, "%d %d\n", mris->nvertices, mris->nfaces);

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    fprintf(fp, "%f  %f  %f  %d\n", v->x, v->y, v->z, v->ripflag);
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      fprintf(fp, "%d ", face->v[n]);
    }
    fprintf(fp, "%d\n", face->ripflag);
  }

  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Description: Same as MRISwriteAscii, but write surface
  normals instead of
  ------------------------------------------------------*/
int MRISwriteNormalsAscii(MRI_SURFACE *mris, const char *fname)
{
  int vno, fno, n;
  VERTEX *v;
  FACE *face;
  FILE *fp;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteNormalsAscii: could not open file %s", fname));

  fprintf(fp, "#!ascii version of %s (vertices are surface normals)\n", mris->fname.data());
  fprintf(fp, "%d %d\n", mris->nvertices, mris->nfaces);

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    fprintf(fp, "%f  %f  %f  %d\n", v->nx, v->ny, v->nz, v->ripflag);
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      fprintf(fp, "%d ", face->v[n]);
    }
    fprintf(fp, "%d\n", face->ripflag);
  }

  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Description: Same write the surface normals into a .mgz volume
  ------------------------------------------------------*/
int MRISwriteNormals(MRI_SURFACE *mris, const char *fname)
{
  int vno;
  VERTEX *v;
  MRI *mri;

  mri = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);
  if (mri == NULL) {
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteNormals(%s): could not allocate volume", fname));
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    MRIsetVoxVal(mri, vno, 0, 0, 0, v->nx);
    MRIsetVoxVal(mri, vno, 0, 0, 1, v->ny);
    MRIsetVoxVal(mri, vno, 0, 0, 2, v->nz);
  }
  MRIwrite(mri, fname);
  MRIfree(&mri);
  return (NO_ERROR);
}
int MRISreadNormals(MRI_SURFACE *mris, const char *fname)
{
  int vno;
  VERTEX *v;
  MRI *mri;

  mri = MRIread(fname);
  if (mri == NULL) {
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadNormals(%s): could not open volume", fname));
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->nx = MRIgetVoxVal(mri, vno, 0, 0, 0);
    v->ny = MRIgetVoxVal(mri, vno, 0, 0, 1);
    v->nz = MRIgetVoxVal(mri, vno, 0, 0, 2);
  }
  MRIfree(&mri);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Description: Same write the white matter surface normals into a .mgz volume
  ------------------------------------------------------*/
int MRISwriteWhiteNormals(MRI_SURFACE *mris, const char *fname)
{
  int vno;
  VERTEX *v;
  MRI *mri;

  mri = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);
  if (mri == NULL) {
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteWhiteNormals(%s): could not allocate volume", fname));
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    MRIsetVoxVal(mri, vno, 0, 0, 0, v->wnx);
    MRIsetVoxVal(mri, vno, 0, 0, 1, v->wny);
    MRIsetVoxVal(mri, vno, 0, 0, 2, v->wnz);
  }
  MRIwrite(mri, fname);
  MRIfree(&mri);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Description: write the principal curvature directions
  ------------------------------------------------------*/
int MRISwritePrincipalDirection(MRI_SURFACE *mris, int dir_index, const char *fname)
{
  int vno;
  VERTEX *v;
  MRI *mri;

  mri = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);
  if (mri == NULL) {
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwritePrincipalDirections(%s): could not allocate volume", fname));
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];

    if (dir_index == 1) {
      MRIsetVoxVal(mri, vno, 0, 0, 0, v->e1x);
      MRIsetVoxVal(mri, vno, 0, 0, 1, v->e1y);
      MRIsetVoxVal(mri, vno, 0, 0, 2, v->e1z);
    }
    else {
      MRIsetVoxVal(mri, vno, 0, 0, 0, v->e2x);
      MRIsetVoxVal(mri, vno, 0, 0, 1, v->e2y);
      MRIsetVoxVal(mri, vno, 0, 0, 2, v->e2z);
    }
  }
  MRIwrite(mri, fname);
  MRIfree(&mri);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteVTK(MRI_SURFACE *mris, const char *fname)
{
  int vno, fno, n;
  VERTEX *v;
  FACE *face;
  FILE *fp;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteVTK: could not open file %s", fname));

  fprintf(fp, "# vtk DataFile Version 1.0\n");
  fprintf(fp, "vtk output\nASCII\nDATASET POLYDATA\nPOINTS %d float\n", mris->nvertices);
  /*  fprintf(fp, "%d %d\n", mris->nvertices, mris->nfaces) ;*/

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    fprintf(fp, "%.9f  %.9f  %.9f\n", v->x, v->y, v->z);
  }
  fprintf(fp, "POLYGONS %d %d\n", mris->nfaces, mris->nfaces * 4);
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    fprintf(fp, "%d ", VERTICES_PER_FACE);
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      fprintf(fp, "%d ", face->v[n]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description: Append the curv data to 'fname', which
  is assumed to have been written with point and face
  data via the MRISwriteVTK routine.
  ------------------------------------------------------*/
int MRISwriteCurvVTK(MRI_SURFACE *mris, const char *fname)
{
  FILE *fp = fopen(fname, "a");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteScalarVTK: could not open file %s", fname));

  fprintf(fp, "POINT_DATA %d\n", mris->nvertices);
  fprintf(fp, "FIELD FieldData 1\n");
  fprintf(fp, "curv 1 %d double\n", mris->nvertices);

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    fprintf(fp, "%.9f\n", v->curv);
  }

  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI_SURFACE *MRISreadVTK(MRI_SURFACE *mris, const char *fname)
{
  char line[STRLEN], *cp = NULL;

  FILE *fp = fopen(fname, "r");
  if (!fp) {
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: could not open file %s\n", fname));
  }

  /* a valid file must have these three lines (see MRISwriteVTK):
  vtk output
  ASCII
  DATASET POLYDATA
  */
  int checks = 0;
  while ((cp = fgetl(line, STRLEN, fp))) {
    // handle both upper and lower case
    if (strncmp(cp, "VTK OUTPUT", 10) == 0) {
      checks++;
    }
    if (strncmp(cp, "vtk output", 10) == 0) {
      checks++;
    }
    if (strncmp(cp, "ASCII", 5) == 0) {
      checks++;
    }
    if (strncmp(cp, "ascii", 5) == 0) {
      checks++;
    }
    if (strncmp(cp, "DATASET POLYDATA", 16) == 0) {
      checks++;
    }
    if (strncmp(cp, "dataset polydata", 16) == 0) {
      checks++;
    }
    // printf("%s\n",cp);
    if (checks == 3) {
      break;
    }
  }
  if (!cp || (checks != 3)) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error reading file %s\n", fname));
  }

  /* now the next line should look like this (get the number of vertices):
  POINTS 142921 float
  */
  int nvertices = 0;
  cp = fgetl(line, STRLEN, fp);
  if (!cp) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error parsing POINTS from file %s\n", fname));
  }
  if ((sscanf(cp, "POINTS %d float\n", &nvertices)) || (sscanf(cp, "points %d float\n", &nvertices))) {
    if (nvertices <= 0) {
      fclose(fp);
      ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error reading file %s, invalid nvertices=%d\n", fname, nvertices));
    }
  }
  else {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error reading file %s, no POINTS field\n", fname));
  }

  /* if needed, alloc MRIS structure (we'll handle faces later in the file) */
  int newMris = 0;  // flag to indicate vertices and face data are empty
  if (NULL == mris) {
    mris = MRISalloc(nvertices, 0);  // note: nfaces=0 for now...
    mris->type = MRIS_TRIANGULAR_SURFACE;
    newMris = 1;  // flag indicating need to fill-in vertex and face data
  }
  else {
    // confirm that the mris structure we've been passed has the same
    // number of vertices
    if (mris->nvertices != nvertices) {
      fclose(fp);
      ErrorReturn(NULL,
                  (ERROR_NOFILE,
                   "MRISreadVTK: error reading file %s, mismatched nvertices=%d/%d\n",
                   fname,
                   mris->nvertices,
                   nvertices));
    }
  }

  /* read vertices... */
  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    float x,y,z;
    fscanf(fp, "%f %f %f", &x, &y, &z);
    if (newMris)  // if passed an mris struct, dont overwrite x,y,z
    {
      MRISsetXYZ(mris,vno,x,y,z);
    }
  }

  /* now, if we're lucky, we've read all the vertices, and the next line
     should be something like this:
  POLYGONS 285838 1143352
  */
  cp = fgetl(line, STRLEN, fp);
  if (!cp) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error parsing POLYGONS in file %s", fname));
  }
  int nfaces = 0;
  int somenum = 0;  // should be nfaces * 4
  if ((sscanf(cp, "POLYGONS %d %d\n", &nfaces, &somenum)) || (sscanf(cp, "polygons %d %d\n", &nfaces, &somenum))) {
    if (nfaces <= 0) {
      fclose(fp);
      ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error reading file %s, invalid nfaces=%d\n", fname, nfaces));
    }
    if (somenum != (nfaces * 4))  // check for 3-vertex face data
    {
      fclose(fp);
      ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error reading file %s, invalid POLYGON data\n", fname));
    }
  }
  else {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error parsing POLYGONS in file %s", fname));
  }
  // POLYGONS field looks valid....

  /* now we can book some face space (if we have a new mris structure) */
  if (newMris) {
    if (!MRISallocateFaces(mris, nfaces)) {
      fclose(fp);
      ErrorExit(ERROR_NO_MEMORY, "MRISreadVTK(%d, %d): could not allocate faces", nfaces, sizeof(FACE));
    }
  }
  else {
    // confirm that the mris structure passed to us has the expected num faces
    if (mris->nfaces != nfaces) {
      fclose(fp);
      ErrorReturn(
          NULL,
          (ERROR_NOFILE, "MRISreadVTK: error reading file %s, mismatched nfaces=%d/%d\n", fname, mris->nfaces, nfaces));
    }
  }

  /* finally, read the face data...*/
  int fno, facepoints = 0;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE *face = &mris->faces[fno];
    fscanf(fp, "%d ", &facepoints);
    if (facepoints != 3) {
      fclose(fp);
      ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error reading file %s, facepoints != 3\n", fname));
    }
    int n;
    for (n = 0; n < 3; n++) {
      vno = -1;
      fscanf(fp, "%d", &vno);
      if ((vno < 0) || (vno >= nvertices)) {
        fclose(fp);
        ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: invalid vertex num in %s, vno = %d\n", fname, vno));
      }
      if (newMris)  // dont overwrite face data if passed an mris struct
      {
        // fill-in the face data to our new mris struct
        face->v[n] = vno;
        mris->vertices_topology[vno].num++;
      }
      else {
        // confirm that the mris structure passed to us has the same face num
        if (face->v[n] != vno) {
          fclose(fp);
          ErrorReturn(
              NULL,
              (ERROR_NOFILE, "MRISreadVTK: error reading file %s, mismatched faces=%d/%d\n", fname, face->v[n], vno));
        }
      }
    }
  }
  if (fno != mris->nfaces) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: failure reading %s, fno=%d != nfaces\n", fname, fno, mris->nfaces));
  }

  /* at this point in the file, we're at the end, or there is possibly
     curv or scalar data fields to read */
  cp = fgetl(line, STRLEN, fp);
  // printf("%s\n",cp);
  int nvertices_data = 0;
  if (cp)  // it appears we have another line...
  {
    if ((sscanf(cp, "POINT_DATA %d\n", &nvertices_data)) || (sscanf(cp, "point_data %d\n", &nvertices_data))) {
      if (nvertices != nvertices_data) {
        fclose(fp);
        ErrorReturn(NULL,
                    (ERROR_NOFILE,
                     "MRISreadVTK: error reading file %s, invalid nvertices=%d "
                     "in POINT_DATA field, expected %d vertices\n",
                     fname,
                     nvertices_data,
                     nvertices));
      }
      int nfields = 0;
      cp = fgetl(line, STRLEN, fp);
      // printf("%s\n",cp);
      if ((sscanf(cp, "FIELD %*s %d\n", &nfields)) || (sscanf(cp, "field %*s %d\n", &nfields))) {
        int num_fields;
        for (num_fields = 0; num_fields < nfields; num_fields++) {
          // parse each data field
          char fieldName[STRLEN];
          int fieldNum = 0;
          cp = fgetl(line, STRLEN, fp);
          // printf("%s\n",cp);
          if (sscanf(cp, "%s %d %d %*s\n", fieldName, &fieldNum, &nvertices_data)) {
            if (nvertices != nvertices_data) {
              fclose(fp);
              ErrorReturn(NULL,
                          (ERROR_NOFILE,
                           "MRISreadVTK: error reading file %s, invalid nvertices=%d "
                           "in %s field, expected %d vertices\n",
                           fname,
                           nvertices_data,
                           fieldName,
                           nvertices));
            }
            // read either scalar or curv data
            int isCurvData = 0;
            if (strncmp(fieldName, "curv", 4) == 0) {
              isCurvData = 1;
            }
            if (strncmp(fieldName, "CURV", 4) == 0) {
              isCurvData = 1;
            }
            float curvmin = 10000.0f;
            float curvmax = -10000.0f; /* for compiler warnings */
            for (vno = 0; vno < mris->nvertices; vno++) {
              VERTEX *v = &mris->vertices[vno];
              float f = 0.0;
              if (fscanf(fp, "%f", &f)) {
                v->curv = f;  // fill-in both curvature and scalar data fields
                v->val = f;
                if (isCurvData) {
                  // printf("%f\n",v->curv);
                  if (vno == 0) {
                    curvmin = curvmax = v->curv;
                  }
                  if (v->curv > curvmax) {
                    curvmax = v->curv;
                  }
                  if (v->curv < curvmin) {
                    curvmin = v->curv;
                  }
                }
              }
              else {
                fclose(fp);
                ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error parsing %s FIELD in file %s", fieldName, fname));
              }
            }
            if (isCurvData) {
              mris->max_curv = curvmax;
              mris->min_curv = curvmin;
              // printf("maxcurv=%f, mincurv=%f\n",curvmax, curvmin);
            }
          }
          else {
            fclose(fp);
            ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadVTK: error parsing FIELDs in file %s", fname));
          }
        }
      }
    }
  }

  fclose(fp);
  
  return (mris);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteGeo(MRI_SURFACE *mris, const char *fname)
{
  int vno, fno, n, actual_vno, toggle, nfaces, nvertices, vnos[300000];
  VERTEX *v;
  FACE *face;
  FILE *fp;

  nfaces = mrisValidFaces(mris);
  nvertices = MRISvalidVertices(mris);
  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteGeo: could not open file %s", fname));

  fprintf(fp, "    1 %d %d %d\n", nvertices, nfaces, VERTICES_PER_FACE * nfaces);
  fprintf(fp, "    1 %d\n", nfaces);
  for (actual_vno = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    vnos[vno] = actual_vno++;
  }

  for (toggle = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    fprintf(fp, "%12.5e %12.5e %12.5e", v->x, v->y, v->z);
    if (toggle) {
      toggle = 0;
      fprintf(fp, "\n");
    }
    else {
      toggle = 1;
      fprintf(fp, " ");
    }
  }
  if (toggle) {
    fprintf(fp, "\n");
  }

  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    for (n = 0; n < VERTICES_PER_FACE - 2; n++) {
      fprintf(fp, "%d ", vnos[face->v[n]] + 1); /* 1-based */
    }

    /* swap order on output to conform to movie.byu convention */
    fprintf(fp, "%d ", vnos[face->v[VERTICES_PER_FACE - 1]] + 1);
    fprintf(fp, "-%d\n", vnos[face->v[VERTICES_PER_FACE - 2]] + 1);
  }

  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters: MRI_SURFACE *mris, char *fname

  Returns value: int

  Description: write ascii icosahedron data (vertices and face vertices info)
  ------------------------------------------------------*/
/* note that .tri or .ico file.  numbering is 1-based output.*/
int MRISwriteICO(MRI_SURFACE *mris, const char *fname)
{
  int vno, fno, nfaces, nvertices;
  int actual_fno, actual_vno;
  VERTEX *v;
  FACE *face;
  FILE *fp;
  // get the valid faces and vertices numbers
  nfaces = mrisValidFaces(mris);
  nvertices = MRISvalidVertices(mris);

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwriteICO: could not open file %s", fname));

  // write number of vertices
  fprintf(fp, "%8d\n", nvertices);
  // go over all vertices and ignoring bad ones
  actual_vno = 1;  // count from 1 (1-based output)
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    fprintf(fp, "%8d %8.4f %8.4f %8.4f\n", actual_vno, v->x, v->y, v->z);
    actual_vno++;
  }
  // write number of faces
  fprintf(fp, "%8d\n", nfaces);
  // go over all faces and ignoring bad ones
  actual_fno = 1;  // count from 1 (1-based)
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    // make the vertex number 1 based
    // the vertex ordering flipped to clockwise (see icosahedron.c)
    fprintf(fp, "%8d %8d %8d %8d\n", actual_fno, face->v[0] + 1, face->v[2] + 1, face->v[1] + 1);
    actual_fno++;
  }
  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwritePatchAscii(MRI_SURFACE *mris, const char *fname)
{
  FILE *fp;
  int vno, fno, n, nvertices, nfaces, type;
  VERTEX *v;
  FACE *face;
  int i;

  type = MRISfileNameType(fname);

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISwritePatchAscii: could not open file %s", fname));

  for (nvertices = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    nvertices++;
  }
  for (nfaces = fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    nfaces++;
  }
  fprintf(fp,
          "#!ascii version of patch %s. "
          "The 1st index is not a vertex number\n",
          mris->fname.data());
  fprintf(fp, "%d %d\n", nvertices, nfaces);
  fprintf(stdout, "nvertices=%d (valid=%d) nfaces=%d\n", nvertices, MRISvalidVertices(mris), nfaces);

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (!v->ripflag) {
      // patch file uses border to change vertex index written to a file
      i = (v->border) ? (-(vno + 1)) : (vno + 1);
      fprintf(fp, "%d vno=%d\n", i, vno);
      fprintf(fp, "%f  %f  %f\n", v->x, v->y, v->z);
    }
  }
  // face vertex info
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (!face->ripflag) {
      fprintf(fp, "%d\n", fno);
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        fprintf(fp, "%d ", face->v[n]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static MRI_SURFACE *mrisReadAsciiFile(const char *fname)
{
  MRI_SURFACE *mris;
  char line[STRLEN], *cp;
  int vno, fno, n, nvertices, nfaces, patch, rip;
  FACE *face;
  FILE *fp;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadAsciiFile: could not open file %s", fname));

  patch = 0;
  cp = fgetl(line, STRLEN, fp);
  sscanf(cp, "%d %d\n", &nvertices, &nfaces);
  mris = MRISalloc(nvertices, nfaces);
  mris->type = MRIS_TRIANGULAR_SURFACE;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * v = &mris->vertices[vno];
    float x,y,z;
    fscanf(fp, "%f  %f  %f  %d\n", &x, &y, &z, &rip);
    MRISsetXYZ(mris, vno,x,y,z);
    v->ripflag = rip;
    if (v->ripflag) {
      patch = 1;
    }
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (fno == Gdiag_no) DiagBreak();
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      fscanf(fp, "%d ", &face->v[n]);
      if (face->v[n] < 0 || face->v[n] >= mris->nvertices)
        ErrorExit(ERROR_BADPARM,
                  "%s: face %d vertex %d: %d -- out of range!!! Should be in [0 %d]",
                  Progname,
                  fno,
                  n,
                  face->v[n],
                  mris->nvertices - 1);
      mris->vertices_topology[face->v[n]].num++;
    }
    fscanf(fp, "%d\n", &rip);
    face->ripflag = rip;
  }

  mris->patch = patch;
  if (patch) {
    mris->status = MRIS_PLANE;
  }
  fclose(fp);
  
  return (mris);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static MRI_SURFACE *mrisReadGeoFile(const char *fname)
{
  MRI_SURFACE *mris;
  char line[202], *cp;
  int vno, fno, n, nvertices, nfaces, patch, vertices_per_face, nedges;
  FACE *face;
  FILE *fp;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "mrisReadGeoFile: could not open file %s", fname));

  patch = 0;
  cp = fgetl(line, 100, fp);
  if (sscanf(cp, "%*d %d %d %d\n", &nvertices, &nfaces, &nedges) != 3) {
    fclose(fp);
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "mrisReadGeoFile(%s): could not scan "
                 "dimensions from '%s'",
                 fname,
                 cp));
  }
  vertices_per_face = nedges / nfaces;
  if (vertices_per_face != VERTICES_PER_FACE) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_BADFILE, "mrisReadGeoFile(%s): unsupported vertices/face %d.", fname, vertices_per_face));
  }

  cp = fgetl(line, 200, fp); /* nfaces again */
  mris = MRISalloc(nvertices, nfaces);
  mris->type = MRIS_GEO_TRIANGLE_FILE;
  for (vno = 0; vno < mris->nvertices; vno++) {
    float x,y,z;
    fscanf(fp, "%e %e %e", &x, &y, &z);
    MRISsetXYZ(mris,vno,x,y,z);
    if (ISODD(vno)) {
      fscanf(fp, "\n");
    }
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    int tmp;
    face = &mris->faces[fno];
    for (n = 0; n < VERTICES_PER_FACE - 1; n++) {
      fscanf(fp, "%d ", &face->v[n]);
      face->v[n]--; /* make it 0-based */
      if (face->v[n] < 0) {
        DiagBreak();
      }
      mris->vertices_topology[face->v[n]].num++;
    }
    n = VERTICES_PER_FACE - 1; /* already true - but make it explicit */
    fscanf(fp, "-%d\n", &face->v[n]);
    face->v[n]--; /* make it 0-based */
    mris->vertices_topology[face->v[n]].num++;

    /* swap positions so normal (via cross-product) will point outwards */
    tmp = face->v[1];
    face->v[1] = face->v[2];
    face->v[2] = tmp;
  }

  fclose(fp);
  
  return (mris);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisReadGeoFilePositions(MRI_SURFACE *mris, const char *fname)
{
  char line[202], *cp;
  int vno, nvertices, nfaces, patch, vertices_per_face, nedges;
  FILE *fp;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "mrisReadGeoFile: could not open file %s", fname));

  patch = 0;
  cp = fgetl(line, 100, fp);
  if (sscanf(cp, "%*d %d %d %d\n", &nvertices, &nfaces, &nedges) != 3) {
    fclose(fp);
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,
                 "mrisReadGeoFile(%s): could not scan "
                 "dimensions from '%s'",
                 fname,
                 cp));
  }
  vertices_per_face = nedges / nfaces;
  if (vertices_per_face != VERTICES_PER_FACE) {
    fclose(fp);
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, "mrisReadGeoFile(%s): unsupported vertices/face %d.", fname, vertices_per_face));
  }

  cp = fgetl(line, 200, fp); /* nfaces again */
  for (vno = 0; vno < mris->nvertices; vno++) {
    float x,y,z;
    fscanf(fp, "%e %e %e", &x, &y, &z);
    MRISsetXYZ(mris,vno, x,y,z);
    if (ISODD(vno)) {
      fscanf(fp, "\n");
    }
  }

  fclose(fp);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  ------------------------------------------------------*/
static MRIS* MRISreadOverAlloc_new(const char *fname, double nVFMultiplier);
static MRIS* MRISreadOverAlloc_old(const char *fname, double nVFMultiplier);

MRIS* MRISreadOverAlloc(const char *fname, double nVFMultiplier)
{
  bool useOldBehaviour = false;
  if (useOldBehaviour) {
    switch (copeWithLogicProblem("FREESURFER_fix_MRISreadOverAlloc",
      "was creating triangles with vertices that were the same vno when reading quad files")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      useOldBehaviour = false;
    }
  }
  
  return 
    useOldBehaviour 
    ? MRISreadOverAlloc_old(fname, nVFMultiplier)
    : MRISreadOverAlloc_new(fname, nVFMultiplier);
}

static MRIS* MRISreadOverAlloc_new(const char *fname, double nVFMultiplier)
{
  int const type = MRISfileNameType(fname);  /* using extension to get type */
  
  MRIS* mris    = NULL;
  int   version = -3;
  {
    FILE* fp = NULL;

    // default:

    chklc();                              /* check to make sure license.dat is present */

    if (type == MRIS_ASCII_TRIANGLE_FILE) /* .ASC */
    {
      mris = mrisReadAsciiFile(fname);
      if (!mris) {
        return (NULL);
      }
      version = -3;
    }
    else if (type == MRIS_ICO_FILE) /* .TRI, .ICO */
    {
      mris = ICOreadOverAlloc(fname, nVFMultiplier, 1.0);
      if (!mris) {
        return (NULL);
      }
      return (mris);
      version = -2;
    }
    else if (type == MRIS_GEO_TRIANGLE_FILE) /* .GEO */
    {
      mris = mrisReadGeoFile(fname);
      if (!mris) {
        return (NULL);
      }
      version = -4;
    }
    else if (type == MRIS_STL_FILE) /* .STL */
    {
      mris = mrisReadSTL(fname);
      if (!mris) {
        return (NULL);
      }
      version = -3;
    }
    else if (type == MRIS_VTK_FILE) /* .vtk */
    {
      mris = MRISreadVTK(mris, fname);
      if (!mris) {
        return (NULL);
      }
      version = -3;
    }
    else if (type == MRIS_GIFTI_FILE) /* .gii */
    {
      mris = mrisReadGIFTIfile(fname, NULL);
      if (!mris) {
        return (NULL);
      }
      version = -3; /* Not really sure what is appropriate here */
    }
    else if (type == MRI_MGH_FILE) /* .mgh */
    {
      ErrorExit(ERROR_BADFILE, "ERROR: MRISread: cannot read surface data from file %s!\n", fname);
    }
    else  // default type MRIS_BINARY_QUADRANGLE_FILE ... use magic number
    {
      fp = fopen(fname, "rb");
      if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "MRISread(%s): could not open file", fname));

      int magic = 0;
      int nread = fread3(&magic, fp);
      if (nread != 1) {
        printf("ERROR: reading %s\n", fname);
        printf("Read %d bytes, expected 1\n", nread);
        fclose(fp);
        return (NULL);
      }
      if (magic == QUAD_FILE_MAGIC_NUMBER) {
        version = -1;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
          fprintf(stdout, "new surface file format\n");
        }
      }
      else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) {
        version = -2;
      }
      else if (magic == TRIANGLE_FILE_MAGIC_NUMBER) {
        fclose(fp);
        mris = mrisReadTriangleFile(fname, nVFMultiplier);
        if (!mris) {
          ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n"));
        }
        version = -3;
      }
      else /* no magic number assigned */
      {
        rewind(fp);
        version = 0;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
          printf("surfer: old surface file format\n");
        }
      }
    }

    /* some type of quadrangle file processing */
    if (version >= -2) {
      int nvertices, nquads;
      
      fread3(&nvertices, fp);
      fread3(&nquads, fp); /* # of quadrangles - not triangles */

      if (nvertices <= 0) /* sanity-checks */
        ErrorExit(ERROR_BADFILE,
                  "ERROR: MRISread: file '%s' has %d vertices!\n"
                  "Probably trying to use a scalar data file as a surface!\n",
                  fname,
                  nvertices);
      
      if (nquads > 4 * nvertices) /* sanity-checks */
      {
        fprintf(stderr, "nquads=%d,  nvertices=%d\n", nquads, nvertices);
        ErrorExit(ERROR_BADFILE,
                  "ERROR: MRISread: file '%s' has many more faces than vertices!\n"
                  "Probably trying to use a scalar data file as a surface!\n",
                  fname);
      }

      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stdout, "reading %d vertices and %d faces.\n", nvertices, 2 * nquads);

      mris = MRISoverAlloc(nVFMultiplier * nvertices, nVFMultiplier * 2 * nquads, nvertices, 0);    // don't know yet how many faces there will be
      mris->type = MRIS_BINARY_QUADRANGLE_FILE;

      /* read vertices *************************************************/
      int vno;
      for (vno = 0; vno < nvertices; vno++) {
        VERTEX_TOPOLOGY * const vertext = &mris->vertices_topology[vno];    
        if (version == -1) /* QUAD_FILE_MAGIC_NUMBER */
        {
          int ix,iy,iz;
          fread2(&ix, fp);
          fread2(&iy, fp);
          fread2(&iz, fp);
          MRISsetXYZ(mris,vno,
            ix / 100.0,
            iy / 100.0,
            iz / 100.0);
        }
        else /* version == -2 */ /* NEW_QUAD_FILE_MAGIC_NUMBER */
        {
          float x = freadFloat(fp);
          float y = freadFloat(fp);
          float z = freadFloat(fp);
          MRISsetXYZ(mris,vno, x,y,z);
        }
        if (version == 0) /* old surface format */
        {
          int num;
          fread1(&num, fp); /* # of faces we are part of */
          vertext->num = num;
          vertext->f = (int *)calloc(vertext->num, sizeof(int));
          if (!vertext->f) ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces", vertext->num);
          vertext->n = (uchar *)calloc(vertext->num, sizeof(uchar));
          if (!vertext->n) ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d nbrs", vertext->n);
          
          int n;
          for (n = 0; n < vertext->num; n++) {
            fread3(&vertext->f[n], fp);
          }
        }
        else {
          vertext->num = 0; /* will figure it out */
        }
      }
      
      /* read face vertices *******************************************/
      int fno = 0;
      int quad;
      for (quad = 0; quad < nquads; quad++) {

        cheapAssert(VERTICES_PER_FACE == 3);
        int vertices[4];
        
        int n;
        for (n = 0; n < 4; n++) /* read quandrangular face */ {
          fread3(&vertices[n], fp);
        }

        /* if we're going to be arbitrary,
           we might as well be really arbitrary */
        /*
          NOTE: for this to work properly in the write, the first two
          vertices in the first face (EVEN and ODD) must be 0 and 1.
        */
        int which = WHICH_FACE_SPLIT(vertices[0], vertices[1]);

        /* 1st triangle */
        int va_0, va_1, va_2, vb_0, vb_1, vb_2;
        
        if (EVEN(which)) {
          va_0 = vertices[0];   vb_0 = vertices[2];
          va_1 = vertices[1];   vb_1 = vertices[3];
          va_2 = vertices[3];   vb_2 = vertices[1];
        } else {
          va_0 = vertices[0];   vb_0 = vertices[0];
          va_1 = vertices[1];   vb_1 = vertices[2];
          va_2 = vertices[2];   vb_2 = vertices[3];
        }

        // make faces for the true triangles        
        for (n = 0; n < 2; n++) {
          if (va_0 == va_1 || va_0 == va_2 || va_1 == va_2) continue;   // degenerate
          mris->faces[fno].v[0] = va_0;
          mris->faces[fno].v[1] = va_1;
          mris->faces[fno].v[2] = va_2;
          int m;
          for (m = 0; m < VERTICES_PER_FACE; m++) {
            mris->vertices_topology[mris->faces[fno].v[m]].num++;
          }
          fno++;
          va_0 = vb_0; va_1 = vb_1; va_2 = vb_2;
        }
      }
      cheapAssert(fno <= mris->max_faces);
      MRISgrowNFaces(mris, fno);
      
      mris->useRealRAS = 0;

      // read tags
      {
        long long len;

        int tag;
        while ((tag = TAGreadStart(fp, &len)) != 0) {
          switch (tag) {
            case TAG_GROUP_AVG_SURFACE_AREA:
              mris->group_avg_surface_area = freadFloat(fp);
              fprintf(
                  stdout, "reading group avg surface area %2.0f cm^2 from file\n", mris->group_avg_surface_area / 100.0);
              break;
            case TAG_OLD_SURF_GEOM:
              readVolGeom(fp, &mris->vg);
              break;
            case TAG_OLD_USEREALRAS:
            case TAG_USEREALRAS:
              if (!freadIntEx(&mris->useRealRAS, fp))  // set useRealRAS
              {
                mris->useRealRAS = 0;  // if error, set to default
              }
              break;
            case TAG_CMDLINE:
              if (mris->ncmds > MAX_CMDS)
                ErrorExit(ERROR_NOMEMORY, "mghRead(%s): too many commands (%d) in file", fname, mris->ncmds);
              mris->cmdlines[mris->ncmds] = (char *)calloc(len + 1, sizeof(char));
              fread(mris->cmdlines[mris->ncmds], sizeof(char), len, fp);
              mris->cmdlines[mris->ncmds][len] = 0;
              mris->ncmds++;
              break;
            default:
              TAGskip(fp, tag, (long long)len);
              break;
          }
        }
      }
      fclose(fp);
      fp = NULL;
    }
    /* end of quadrangle file processing */
    /* file is closed now for all types ***********************************/
  }
  
  /* find out if this surface is lh or rh from fname */
  strcpy(mris->fname, fname);
  {
    const char *surf_name;

    surf_name = strrchr(fname, '/');
    if (surf_name == NULL) {
      surf_name = fname;
    }
    else {
      surf_name++; /* past the last slash */
    }
    if (toupper(*surf_name) == 'R') {
      mris->hemisphere = RIGHT_HEMISPHERE;
    }
    else if (toupper(*surf_name) == 'L') {
      mris->hemisphere = LEFT_HEMISPHERE;
    }
    else if (toupper(*surf_name) == 'B') {
      mris->hemisphere = BOTH_HEMISPHERES;
    }
    else {
      mris->hemisphere = NO_HEMISPHERE;
    }
  }

  /***********************************************************************/
  /* build members of mris structure                                     */
  /***********************************************************************/
  if ((version < 0) || type == MRIS_ASCII_TRIANGLE_FILE) {
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      mris->vertices_topology[vno].f = (int *)calloc(mris->vertices_topology[vno].num, sizeof(int));
      if (!mris->vertices_topology[vno].f)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d faces at %dth vertex",
                  fname,
                  vno,
                  mris->vertices_topology[vno].num);

      mris->vertices_topology[vno].n = (uchar *)calloc(mris->vertices_topology[vno].num, sizeof(uchar));
      if (!mris->vertices_topology[vno].n)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d indices at %dth vertex",
                  fname,
                  vno,
                  mris->vertices_topology[vno].num);
      mris->vertices_topology[vno].num = 0;
    }
    
    // This is probably unnecessary, given the mrisCompleteTopology below
    // but I am worried that code won't get them in the same, and hence get equivalent but different results
    //
    int fno;
    for (fno = 0; fno < mris->nfaces; fno++) {
      FACE* face = &mris->faces[fno];
      int n;
      for (n = 0; n < VERTICES_PER_FACE; n++) mris->vertices_topology[face->v[n]].f[mris->vertices_topology[face->v[n]].num++] = fno;
    }
  }

  {
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {

      if (vno == Gdiag_no) {
        DiagBreak();
      }

      mris->vertices[vno].curv     = 0;
      mris->vertices[vno].origarea = -1;
      mris->vertices[vno].border   = 0;

      // This is probably unnecessary, given the mrisCompleteTopology below
      // but I am worried that code won't get them in the same, and hence get equivalent but different results
      //
      int n;
      for (n = 0; n < mris->vertices_topology[vno].num; n++) {
        int m;
        for (m = 0; m < VERTICES_PER_FACE; m++) {
          if (mris->faces[mris->vertices_topology[vno].f[n]].v[m] == vno) {
            mris->vertices_topology[vno].n[n] = m;
          }
        }
      }
    }
  }
  
  mrisCompleteTopology(mris);
  
  mrisCheckVertexFaceTopology(mris);

  mrisComputeSurfaceDimensions(mris);
  mrisComputeVertexDistances(mris);
  MRIScomputeNormals(mris);

  mrisReadTransform(mris, fname);

  mris->radius = MRISaverageRadius(mris);

  MRIScomputeMetricProperties(mris);

  MRISstoreCurrentPositions(mris);

  // Check whether there is an area file for group average
  char tmpstr[2000];
  sprintf(tmpstr, "%s.avg.area.mgh", fname);

  if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
    printf("Trying to read average area %s\n", tmpstr);
  }

  if (fio_FileExistsReadable(tmpstr)) {
    if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
      printf("Reading in average area %s\n", tmpstr);
    }
    MRI* mri = MRIread(tmpstr);
    if (!mri) {
      printf("ERROR: reading in average area %s\n", tmpstr);
      return (NULL);
    }
    MRIScopyMRI(mris, mri, 0, "group_avg_area");
    MRIfree(&mri);
    mris->group_avg_vtxarea_loaded = 1;
  }
  else {
    mris->group_avg_vtxarea_loaded = 0;
  }

  if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
    printf("Average area loaded %d\n", mris->group_avg_vtxarea_loaded);
  }

  return (mris);
}

static MRIS* MRISreadOverAlloc_old(const char *fname, double nVFMultiplier)
{
  MRI_SURFACE *mris = NULL;
  int nquads, nvertices, magic, version, ix, iy, iz, vno, fno, n, m;
  int imnr, imnr0, imnr1, type, vertices[VERTICES_PER_FACE + 1], num;
  float x, y, z, xhi, xlo, yhi, ylo, zhi, zlo;
  FILE *fp = NULL;
  FACE *face;
  int tag, nread;
  char tmpstr[2000];
  MRI *mri;

  // default:
  version = -3;

  chklc();                              /* check to make sure license.dat is present */
  type = MRISfileNameType(fname);       /* using extension to get type */
  if (type == MRIS_ASCII_TRIANGLE_FILE) /* .ASC */
  {
    mris = mrisReadAsciiFile(fname);
    if (!mris) {
      return (NULL);
    }
    version = -3;
  }
  else if (type == MRIS_ICO_FILE) /* .TRI, .ICO */
  {
    mris = ICOreadOverAlloc(fname, nVFMultiplier, 1.0);
    if (!mris) {
      return (NULL);
    }
    return (mris);
    version = -2;
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE) /* .GEO */
  {
    mris = mrisReadGeoFile(fname);
    if (!mris) {
      return (NULL);
    }
    version = -4;
  }
  else if (type == MRIS_STL_FILE) /* .STL */
  {
    mris = mrisReadSTL(fname);
    if (!mris) {
      return (NULL);
    }
    version = -3;
  }
  else if (type == MRIS_VTK_FILE) /* .vtk */
  {
    mris = MRISreadVTK(mris, fname);
    if (!mris) {
      return (NULL);
    }
    version = -3;
  }
  else if (type == MRIS_GIFTI_FILE) /* .gii */
  {
    mris = mrisReadGIFTIfile(fname, NULL);
    if (!mris) {
      return (NULL);
    }
    version = -3; /* Not really sure what is appropriate here */
  }
  else if (type == MRI_MGH_FILE) /* .mgh */
  {
    ErrorExit(ERROR_BADFILE, "ERROR: MRISread: cannot read surface data from file %s!\n", fname);
  }
  else  // default type MRIS_BINARY_QUADRANGLE_FILE ... use magic number
  {
    fp = fopen(fname, "rb");
    if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "MRISread(%s): could not open file", fname));

    magic = 0;
    nread = fread3(&magic, fp);
    if (nread != 1) {
      printf("ERROR: reading %s\n", fname);
      printf("Read %d bytes, expected 1\n", nread);
      fclose(fp);
      return (NULL);
    }
    if (magic == QUAD_FILE_MAGIC_NUMBER) {
      version = -1;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        fprintf(stdout, "new surface file format\n");
      }
    }
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) {
      version = -2;
    }
    else if (magic == TRIANGLE_FILE_MAGIC_NUMBER) {
      fclose(fp);
      mris = mrisReadTriangleFile(fname, nVFMultiplier);
      if (!mris) {
        ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n"));
      }
      version = -3;
    }
    else /* no magic number assigned */
    {
      rewind(fp);
      version = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        printf("surfer: old surface file format\n");
      }
    }
  }
  /* some type of quadrangle file processing */
  if (version >= -2) {
    fread3(&nvertices, fp);
    fread3(&nquads, fp); /* # of qaudrangles - not triangles */

    if (nvertices <= 0) /* sanity-checks */
      ErrorExit(ERROR_BADFILE,
                "ERROR: MRISread: file '%s' has %d vertices!\n"
                "Probably trying to use a scalar data file as a surface!\n",
                fname,
                nvertices);
    if (nquads > 4 * nvertices) /* sanity-checks */
    {
      fprintf(stderr, "nquads=%d,  nvertices=%d\n", nquads, nvertices);
      ErrorExit(ERROR_BADFILE,
                "ERROR: MRISread: file '%s' has many more faces than vertices!\n"
                "Probably trying to use a scalar data file as a surface!\n",
                fname);
    }

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "reading %d vertices and %d faces.\n", nvertices, 2 * nquads);

    mris = MRISoverAlloc(nVFMultiplier * nvertices, nVFMultiplier * 2 * nquads, nvertices, 2 * nquads);
    mris->type = MRIS_BINARY_QUADRANGLE_FILE;

    imnr0 = 1000;
    imnr1 = 0;
    /* read vertices *************************************************/
    for (vno = 0; vno < nvertices; vno++) {
      VERTEX_TOPOLOGY * const vertext = &mris->vertices_topology[vno];    
      VERTEX          * const vertex  = &mris->vertices         [vno];
      if (version == -1) /* QUAD_FILE_MAGIC_NUMBER */
      {
        fread2(&ix, fp);
        fread2(&iy, fp);
        fread2(&iz, fp);
        MRISsetXYZ(mris,vno,
          ix / 100.0,
          iy / 100.0,
          iz / 100.0);
      }
      else /* version == -2 */ /* NEW_QUAD_FILE_MAGIC_NUMBER */
      {
        float x = freadFloat(fp);
        float y = freadFloat(fp);
        float z = freadFloat(fp);
        MRISsetXYZ(mris,vno, x,y,z);
      }
      /* brain-dead code and never used again either */
      imnr = (int)((vertex->y - START_Y) / SLICE_THICKNESS + 0.5);
      if (imnr > imnr1) {
        imnr1 = imnr;
      }
      if (imnr < imnr0) {
        imnr0 = imnr;
      }
      if (version == 0) /* old surface format */
      {
        fread1(&num, fp); /* # of faces we are part of */
        vertext->num = num;
        vertext->f = (int *)calloc(vertext->num, sizeof(int));
        if (!vertext->f) ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d faces", vertext->num);
        vertext->n = (uchar *)calloc(vertext->num, sizeof(uchar));
        if (!vertext->n) ErrorExit(ERROR_NO_MEMORY, "MRISread: could not allocate %d nbrs", vertext->n);
        for (n = 0; n < vertext->num; n++) {
          fread3(&vertext->f[n], fp);
        }
      }
      else {
        vertext->num = 0; /* will figure it out */
      }
    }
    /* read face vertices *******************************************/
    for (fno = 0; fno < mris->nfaces; fno += 2) {
      int which;

      if (fno == 86) {
        DiagBreak();
      }
      for (n = 0; n < 4; n++) /* read quandrangular face */
      {
        fread3(&vertices[n], fp);
        if (vertices[n] == 22) {
          DiagBreak();
        }
      }

      /* if we're going to be arbitrary,
         we might as well be really arbitrary */
      /*
        NOTE: for this to work properly in the write, the first two
        vertices in the first face (EVEN and ODD) must be 0 and 1.
      */
      which = WHICH_FACE_SPLIT(vertices[0], vertices[1]);

      /* 1st triangle */
      if (EVEN(which)) {
        mris->faces[fno].v[0] = vertices[0];
        mris->faces[fno].v[1] = vertices[1];
        mris->faces[fno].v[2] = vertices[3];

        /* 2nd triangle */
        mris->faces[fno + 1].v[0] = vertices[2];
        mris->faces[fno + 1].v[1] = vertices[3];
        mris->faces[fno + 1].v[2] = vertices[1];
      }
      else {
        mris->faces[fno].v[0] = vertices[0];
        mris->faces[fno].v[1] = vertices[1];
        mris->faces[fno].v[2] = vertices[2];

        /* 2nd triangle */
        mris->faces[fno + 1].v[0] = vertices[0];
        mris->faces[fno + 1].v[1] = vertices[2];
        mris->faces[fno + 1].v[2] = vertices[3];
      }
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        mris->vertices_topology[mris->faces[fno    ].v[n]].num++;
        mris->vertices_topology[mris->faces[fno + 1].v[n]].num++;
      }
    }
    mris->useRealRAS = 0;

    // read tags
    {
      long long len;

      while ((tag = TAGreadStart(fp, &len)) != 0) {
        switch (tag) {
          case TAG_GROUP_AVG_SURFACE_AREA:
            mris->group_avg_surface_area = freadFloat(fp);
            fprintf(
                stdout, "reading group avg surface area %2.0f cm^2 from file\n", mris->group_avg_surface_area / 100.0);
            break;
          case TAG_OLD_SURF_GEOM:
            readVolGeom(fp, &mris->vg);
            break;
          case TAG_OLD_USEREALRAS:
          case TAG_USEREALRAS:
            if (!freadIntEx(&mris->useRealRAS, fp))  // set useRealRAS
            {
              mris->useRealRAS = 0;  // if error, set to default
            }
            break;
          case TAG_CMDLINE:
            if (mris->ncmds > MAX_CMDS)
              ErrorExit(ERROR_NOMEMORY, "mghRead(%s): too many commands (%d) in file", fname, mris->ncmds);
            mris->cmdlines[mris->ncmds] = (char *)calloc(len + 1, sizeof(char));
            fread(mris->cmdlines[mris->ncmds], sizeof(char), len, fp);
            mris->cmdlines[mris->ncmds][len] = 0;
            mris->ncmds++;
            break;
          default:
            TAGskip(fp, tag, (long long)len);
            break;
        }
      }
    }
    fclose(fp);
    fp = NULL;
  }
  /* end of quadrangle file processing */
  /* file is closed now for all types ***********************************/

  /* find out if this surface is lh or rh from fname */
  strcpy(mris->fname, fname);
  {
    const char *surf_name;

    surf_name = strrchr(fname, '/');
    if (surf_name == NULL) {
      surf_name = fname;
    }
    else {
      surf_name++; /* past the last slash */
    }
    if (toupper(*surf_name) == 'R') {
      mris->hemisphere = RIGHT_HEMISPHERE;
    }
    else if (toupper(*surf_name) == 'L') {
      mris->hemisphere = LEFT_HEMISPHERE;
    }
    else if (toupper(*surf_name) == 'B') {
      mris->hemisphere = BOTH_HEMISPHERES;
    }
    else {
      mris->hemisphere = NO_HEMISPHERE;
    }
  }

  /***********************************************************************/
  /* build members of mris structure                                     */
  /***********************************************************************/
  if ((version < 0) || type == MRIS_ASCII_TRIANGLE_FILE) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      mris->vertices_topology[vno].f = (int *)calloc(mris->vertices_topology[vno].num, sizeof(int));
      if (!mris->vertices_topology[vno].f)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d faces at %dth vertex",
                  fname,
                  vno,
                  mris->vertices_topology[vno].num);

      mris->vertices_topology[vno].n = (uchar *)calloc(mris->vertices_topology[vno].num, sizeof(uchar));
      if (!mris->vertices_topology[vno].n)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d indices at %dth vertex",
                  fname,
                  vno,
                  mris->vertices_topology[vno].num);
      mris->vertices_topology[vno].num = 0;
    }
    for (fno = 0; fno < mris->nfaces; fno++) {
      face = &mris->faces[fno];
      for (n = 0; n < VERTICES_PER_FACE; n++) mris->vertices_topology[face->v[n]].f[mris->vertices_topology[face->v[n]].num++] = fno;
    }
  }

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    mris->vertices[vno].curv = 0;
    mris->vertices[vno].origarea = -1;
    mris->vertices[vno].border = 0;

    for (n = 0; n < mris->vertices_topology[vno].num; n++) {
      for (m = 0; m < VERTICES_PER_FACE; m++) {
        if (mris->faces[mris->vertices_topology[vno].f[n]].v[m] == vno) {
          mris->vertices_topology[vno].n[n] = m;
        }
      }
    }
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;
    if (x > xhi) {
      xhi = x;
    }
    if (x < xlo) {
      xlo = x;
    }
    if (y > yhi) {
      yhi = y;
    }
    if (y < ylo) {
      ylo = y;
    }
    if (z > zhi) {
      zhi = z;
    }
    if (z < zlo) {
      zlo = z;
    }
  }
  mris->xlo = xlo;
  mris->ylo = ylo;
  mris->zlo = zlo;
  mris->xhi = xhi;
  mris->yhi = yhi;
  mris->zhi = zhi;
  mris->xctr = (xhi + xlo) / 2;
  mris->yctr = (yhi + ylo) / 2;
  mris->zctr = (zhi + zlo) / 2;
  mrisCompleteTopology(mris);
  
  mrisCheckVertexFaceTopology(mris);
  
  MRIScomputeNormals(mris);
  mrisComputeVertexDistances(mris);

  mrisReadTransform(mris, fname);

  mris->radius = MRISaverageRadius(mris);

  MRIScomputeMetricProperties(mris);

  MRISstoreCurrentPositions(mris);

  // Check whether there is an area file for group average
  sprintf(tmpstr, "%s.avg.area.mgh", fname);
  if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
    printf("Trying to read average area %s\n", tmpstr);
  }
  if (fio_FileExistsReadable(tmpstr)) {
    if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
      printf("Reading in average area %s\n", tmpstr);
    }
    mri = MRIread(tmpstr);
    if (!mri) {
      printf("ERROR: reading in average area %s\n", tmpstr);
      return (NULL);
    }
    MRIScopyMRI(mris, mri, 0, "group_avg_area");
    MRIfree(&mri);
    mris->group_avg_vtxarea_loaded = 1;
  }
  else {
    mris->group_avg_vtxarea_loaded = 0;
  }

  if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
    printf("Average area loaded %d\n", mris->group_avg_vtxarea_loaded);
  }

  return (mris);
}

/*-----------------------------------------------------
  MRISfastRead() just calls MRISRead()
  Parameters:
  Returns value:
  Description
  ------------------------------------------------------*/
MRI_SURFACE *MRISfastRead(const char *fname)
{
  return (MRISread(fname));
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRIS * MRISread(const char *fname)
{
  MRIS *mris = MRISreadOverAlloc(fname, 1.0);
  if (mris == NULL) return (NULL);

  // save xyz coordinates space because mris->useRealRAS is changed after conversion
  mris->orig_xyzspace = mris->useRealRAS;

  // convert surface xyz coordinates scanner space to tkregister space 
  // (ie, the native space that the surface xyz coords are generated in)
  if (mris->useRealRAS)
    MRISscanner2Tkr(mris);

  MRISsetNeighborhoodSizeAndDist(mris, 3);  // find nbhds out to 3-nbrs
  MRISresetNeighborhoodSize(mris, 1);       // reset current size to 1-nbrs
  return (mris);
}

int MRISwriteVertexLocations(MRIS *mris, char *fname, int which_vertices)
{
  int retval, i;
  float *coords[3];

  for (i = 0; i < 3; i++) {
    coords[i] = (float *)calloc(mris->nvertices, sizeof(float));
    if (coords[i] == NULL) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d lh coords", Progname, mris->nvertices);
  }

  MRISextractVertexCoords(mris, coords, CURRENT_VERTICES);
  MRISrestoreVertexPositions(mris, which_vertices);
  retval = MRISwrite(mris, fname);
  MRISimportVertexCoords(mris, coords, CURRENT_VERTICES);

  for (i = 0; i < 3; i++) free(coords[i]);

  return (retval);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define USE_NEW_QUAD_FILE 1  // new style stores float instead of int
static int MRISwrite_new(MRI_SURFACE *mris, const char *name);
static int MRISwrite_old(MRI_SURFACE *mris, const char *name);

int MRISwrite(MRIS *mris, const char *name)
{
  bool useOldBehaviour = false;
  if (useOldBehaviour) {
    switch (copeWithLogicProblem("FREESURFER_fix_MRISwrite",
      "was combining non-abutting triangles into a quad when writing quad files")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      useOldBehaviour = false;
    }
  }
  
  return useOldBehaviour
    ? MRISwrite_old(mris, name)
    : MRISwrite_new(mris, name);
}

static bool quadCombine(int quad[4], int vA[3], int vB[3])
{
  if (0) {  // this has been seen to PASS
  
    static bool laterTime;
    if (!laterTime) { laterTime = true;
      fprintf(stdout, "%s:%d testing quadCombine\n", __FILE__, __LINE__); 
      
      // cases where it should be found
#define QUAD_COMBINE_TEST(A, VA0,VA1,VA2,VB0,VB1,VB2)       \
        {   int va[3] = {VA0,VA1,VA2};                      \
            int vb[3] = {VB0,VB1,VB2};                      \
            cheapAssert(A == quadCombine(quad, va, vb));    \
        }

      QUAD_COMBINE_TEST(TRUE, 0,1,2, 2,1,4);
      QUAD_COMBINE_TEST(TRUE, 0,1,2, 4,2,1);
      QUAD_COMBINE_TEST(TRUE, 0,1,2, 1,4,2);

      QUAD_COMBINE_TEST(TRUE, 1,2,0, 2,1,4);
      QUAD_COMBINE_TEST(TRUE, 1,2,0, 4,2,1);
      QUAD_COMBINE_TEST(TRUE, 1,2,0, 1,4,2);

      QUAD_COMBINE_TEST(TRUE, 2,0,1, 2,1,4);
      QUAD_COMBINE_TEST(TRUE, 2,0,1, 4,2,1);
      QUAD_COMBINE_TEST(TRUE, 2,0,1, 1,4,2);
      
      // cases where it should not be found
      QUAD_COMBINE_TEST(FALSE, 0,1,2, 2,1,0);    // two edges
      QUAD_COMBINE_TEST(FALSE, 0,1,2, 4,1,2);    // one edge the same direction
      QUAD_COMBINE_TEST(FALSE, 0,1,2, 1,2,0);    // two edge the same direction
      QUAD_COMBINE_TEST(FALSE, 0,1,2, 3,4,5);    // no vertices
      QUAD_COMBINE_TEST(FALSE, 1,2,0, 3,1,4);    // one vertex

#undef QUAD_COMBINE_TEST
    }
  }
  
  // Combine into a quad using a shared edge in triangles vA0 vA1 vA2 and vB0 vB1 vB2.
  // Don't combine them when this would flip a triangle - the combined edge must be in reverse order ... vA0 vA1 == vB2 vB1
  //
  int qi = 0;

  int vAi = 0;
  for (vAi = 0; vAi < 3; vAi++) {

    // Get the next edge in A
    //  
    int edge0 = vA[vAi], edge1 = vA[(vAi==2)?0:vAi+1];

    // Put the next vertex from A into quad
    //
    if (qi == 4) return false;              // must have shared two edges!
    quad[qi++] = edge0;
    
    // Is there a shared edge in B?
    //
    int vBi;
    for (vBi = 0; vBi < 3; vBi++) {
      if (edge0 != vB[vBi]) continue;       // not shared vertex
      int prev = vB[(vBi>0)?vBi-1:2];
      if (prev != edge1) break;             // not shared edge
      if (qi == 4) return false;            // must have shared two edges!
      quad[qi++] = vB[(vBi<2)?vBi+1:0];     // put the third vertex of B into the quad
    }
  }
  
  return qi == 4;                           // built a valid quad?
}

static int MRISwrite_new(MRI_SURFACE *mris, const char *name)
{
  int k, type;
  float x, y, z;
  FILE *fp;
  char fname[STRLEN];

  chklc();

  strcpy(fname, name);

  type = MRISfileNameType(fname);
  if (type == MRIS_ASCII_TRIANGLE_FILE) {
    return (MRISwriteAscii(mris, fname));
  }
  else if (type == MRIS_VTK_FILE) {
    return (MRISwriteVTK(mris, fname));
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE) {
    return (MRISwriteGeo(mris, fname));
  }
  else if (type == MRIS_ICO_FILE) {
    return MRISwriteICO(mris, fname);
  }
  else if (type == MRIS_STL_FILE) {
    return mrisWriteSTL(mris, fname);
  }
  else if (type == MRIS_GIFTI_FILE) {
    return MRISwriteGIFTI(mris, NIFTI_INTENT_POINTSET, fname, NULL);
  }

  if (mris->type == MRIS_TRIANGULAR_SURFACE) {
    return (MRISwriteTriangularSurface(mris, fname));
  }

  fp = fopen(fname, "w");
  if (fp == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISwrite(%s): can't create file\n", fname));

#if USE_NEW_QUAD_FILE
  fwrite3(NEW_QUAD_FILE_MAGIC_NUMBER, fp);
#else
  fwrite3(QUAD_FILE_MAGIC_NUMBER, fp);
#endif
  fwrite3(mris->nvertices, fp);

  // Below was combining two adjacent faces without checking whether they had an abutting edge!
  // Calculate how many are really needed
  //
  int quadsNeeded = 0;
  for (k = 0; k < mris->nfaces; k++) {
    FACE* f = &mris->faces[k];

    int quad[4]; 
    if ((k+1 < mris->nfaces) && quadCombine(quad, f->v, mris->faces[k+1].v)) {
      k++;
    }
      
    quadsNeeded++;
  }

  fwrite3(quadsNeeded, fp);

  for (k = 0; k < mris->nvertices; k++) {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
#if USE_NEW_QUAD_FILE
    fwriteFloat(x, fp);
    fwriteFloat(y, fp);
    fwriteFloat(z, fp);
#else
    fwrite2((int)(x * 100), fp);
    fwrite2((int)(y * 100), fp);
    fwrite2((int)(z * 100), fp);
#endif
  }

  // This was combining two adjacent faces without checking whether they had an abutting edge!
  // Now it checks, and writes degenerate quads instead if necessary
  //
  int quadsWritten = 0;
  for (k = 0; k < mris->nfaces; k++) {
    FACE* f = &mris->faces[k];

    int quad[4]; 
    if ((k+1 < mris->nfaces) && quadCombine(quad, f->v, mris->faces[k+1].v)) {
      k++;
    } else {
      quad[0] = f->v[0];
      quad[1] = f->v[1];
      quad[2] = f->v[1];
      quad[3] = f->v[2];
    }
      
    fwrite3(quad[0], fp);
    fwrite3(quad[1], fp);
    fwrite3(quad[2], fp);
    fwrite3(quad[3], fp);
    quadsWritten++;
  }
  cheapAssert(quadsNeeded == quadsWritten);

  /* write whether vertex data was using the
     real RAS rather than conformed RAS */
  fwriteInt(TAG_OLD_USEREALRAS, fp);
  fwriteInt(mris->useRealRAS, fp);
  // volume info
  fwriteInt(TAG_OLD_SURF_GEOM, fp);
  writeVolGeom(fp, &mris->vg);

  if (!FZERO(mris->group_avg_surface_area)) {
    long long here;
    printf("writing group avg surface area %2.0f cm^2 into surface file\n", mris->group_avg_surface_area / 100.0);
    TAGwriteStart(fp, TAG_GROUP_AVG_SURFACE_AREA, &here, sizeof(float));
    fwriteFloat(mris->group_avg_surface_area, fp);
    TAGwriteEnd(fp, here);
  }
  // write other tags
  {
    int i;

    for (i = 0; i < mris->ncmds; i++) TAGwrite(fp, TAG_CMDLINE, mris->cmdlines[i], strlen(mris->cmdlines[i]) + 1);
  }
  fclose(fp);
  return (NO_ERROR);
}


static int MRISwrite_old(MRI_SURFACE *mris, const char *name)
{
  int k, type;
  float x, y, z;
  FILE *fp;
  char fname[STRLEN];

  chklc();

  strcpy(fname, name);

  type = MRISfileNameType(fname);
  if (type == MRIS_ASCII_TRIANGLE_FILE) {
    return (MRISwriteAscii(mris, fname));
  }
  else if (type == MRIS_VTK_FILE) {
    return (MRISwriteVTK(mris, fname));
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE) {
    return (MRISwriteGeo(mris, fname));
  }
  else if (type == MRIS_ICO_FILE) {
    return MRISwriteICO(mris, fname);
  }
  else if (type == MRIS_STL_FILE) {
    return mrisWriteSTL(mris, fname);
  }
  else if (type == MRIS_GIFTI_FILE) {
    return MRISwriteGIFTI(mris, NIFTI_INTENT_POINTSET, fname, NULL);
  }

  if (mris->type == MRIS_TRIANGULAR_SURFACE) {
    return (MRISwriteTriangularSurface(mris, fname));
  }

  fp = fopen(fname, "w");
  if (fp == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISwrite(%s): can't create file\n", fname));

#if USE_NEW_QUAD_FILE
  fwrite3(NEW_QUAD_FILE_MAGIC_NUMBER, fp);
#else
  fwrite3(QUAD_FILE_MAGIC_NUMBER, fp);
#endif
  fwrite3(mris->nvertices, fp);
  fwrite3(mris->nfaces / 2, fp); /* # of quadrangles */
  for (k = 0; k < mris->nvertices; k++) {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
#if USE_NEW_QUAD_FILE
    fwriteFloat(x, fp);
    fwriteFloat(y, fp);
    fwriteFloat(z, fp);
#else
    fwrite2((int)(x * 100), fp);
    fwrite2((int)(y * 100), fp);
    fwrite2((int)(z * 100), fp);
#endif
  }
  for (k = 0; k < mris->nfaces; k += 2) {
    int which;
    FACE *f;

    f = &mris->faces[k];
    {
      int n;
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        if ((mris->faces[k].v[n] == 22) || (mris->faces[k + 1].v[n] == 22)) {
          DiagBreak();
        }
      }
    }
    which = WHICH_FACE_SPLIT(f->v[0], f->v[1]);
    if (EVEN(which)) {
      fwrite3(mris->faces[k].v[0], fp);
      fwrite3(mris->faces[k].v[1], fp);
      fwrite3(mris->faces[k + 1].v[0], fp);
      fwrite3(mris->faces[k].v[2], fp);
    }
    else {
      fwrite3(mris->faces[k].v[0], fp);
      fwrite3(mris->faces[k].v[1], fp);
      fwrite3(mris->faces[k].v[2], fp);
      fwrite3(mris->faces[k + 1].v[2], fp);
    }
  }
  /* write whether vertex data was using the
     real RAS rather than conformed RAS */
  fwriteInt(TAG_OLD_USEREALRAS, fp);
  fwriteInt(mris->useRealRAS, fp);
  // volume info
  fwriteInt(TAG_OLD_SURF_GEOM, fp);
  writeVolGeom(fp, &mris->vg);

  if (!FZERO(mris->group_avg_surface_area)) {
    long long here;
    printf("writing group avg surface area %2.0f cm^2 into surface file\n", mris->group_avg_surface_area / 100.0);
    TAGwriteStart(fp, TAG_GROUP_AVG_SURFACE_AREA, &here, sizeof(float));
    fwriteFloat(mris->group_avg_surface_area, fp);
    TAGwriteEnd(fp, here);
  }
  // write other tags
  {
    int i;

    for (i = 0; i < mris->ncmds; i++) TAGwrite(fp, TAG_CMDLINE, mris->cmdlines[i], strlen(mris->cmdlines[i]) + 1);
  }
  fclose(fp);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
// here reads mri/transforms/talairach.xfm (MNI style transform)

MATRIX *getSRASToTalSRAS(LT *lt)
{
  MATRIX *RASFromSRAS = 0;
  MATRIX *sRASFromRAS = 0;
  MATRIX *tmpM = 0;
  MATRIX *SRASToTalSRAS = 0;

  // now calculate transform
  //  conformed -------> surfaceRAS
  //      |                  |
  //      V                  V  RASFromSRAS
  //     src     ------->   RAS
  //      |                  |   xfm
  //      V                  V
  //     tal    ------->   talRAS
  //      |                  |  sRASFromRAS
  //      V                  V
  //  conformed --------> surfaceRAS

  // must convert to RAS first
  RASFromSRAS = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixIdentity(4, RASFromSRAS);
  *MATRIX_RELT(RASFromSRAS, 1, 4) = lt->src.c_r;
  *MATRIX_RELT(RASFromSRAS, 2, 4) = lt->src.c_a;
  *MATRIX_RELT(RASFromSRAS, 3, 4) = lt->src.c_s;

  sRASFromRAS = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixIdentity(4, sRASFromRAS);
  *MATRIX_RELT(sRASFromRAS, 1, 4) = -lt->dst.c_r;
  *MATRIX_RELT(sRASFromRAS, 2, 4) = -lt->dst.c_a;
  *MATRIX_RELT(sRASFromRAS, 3, 4) = -lt->dst.c_s;

  tmpM = MatrixMultiply(lt->m_L, RASFromSRAS, NULL);
  SRASToTalSRAS = MatrixMultiply(sRASFromRAS, tmpM, NULL);

  MatrixFree(&RASFromSRAS);
  MatrixFree(&sRASFromRAS);
  MatrixFree(&tmpM);

  return SRASToTalSRAS;
}


// this function assumes mris_fname is a surface under recon-all output directory structure
int mrisReadTransform(MRIS *mris, const char *mris_fname)
{
  char transform_fname[STRLEN], fpref[300];
  LT *lt = 0;
  MRI *orig = 0;
  struct stat info;
  int rStat;

  // here it is assumed that subjects is set
  FileNamePath(mris_fname, fpref);
  sprintf(transform_fname, "%s/../mri/transforms/talairach.xfm", fpref);
  if (!FileExists(transform_fname)) {
    return (ERROR_NO_FILE);
  }

  if (!(mris->lta = LTAreadEx(transform_fname))) {
    ErrorReturn(ERROR_NO_FILE, (ERROR_NOFILE, "mrisReadTransform: could not read xform file '%s'", transform_fname));
  }
  else {
    if (mris->lta->type != LINEAR_RAS_TO_RAS)
      ErrorExit(ERROR_BADPARM, "the transform is not RAS-TO-RAS.  not supported.");
  }
  
  //////////////////////////////////////////////////////////////////////
  // thus if transform->src is not set, set it to the orig
  lt = &mris->lta->xforms[0];
  // src information
  if (!lt->src.valid) {
    // first try to get it from surface itself
    if (mris->vg.valid) {
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "INFO: found the orig volume info on %s in the surface data.\n", mris->vg.fname);
      lt->src.c_r = mris->vg.c_r;
      lt->src.c_a = mris->vg.c_a;
      lt->src.c_s = mris->vg.c_s;
    }
    else {
      // first try to get it from mri/orig
      sprintf(transform_fname, "%s/../mri/orig", fpref);  // reuse of the buffer
      rStat = stat(transform_fname, &info);
      if (!rStat && S_ISREG(info.st_mode)) {
        orig = MRIreadHeader(transform_fname, -1);
      }
      if (orig) {
        getVolGeom(orig, &lt->src);
        getVolGeom(orig, &mris->vg);
        // add orig volume info in the surface
        MRIfree(&orig);
        orig = 0;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
          fprintf(stderr,
                  "INFO: found the orig volume (mri/orig) "
                  "to get c_(ras) information for src\n");
          fprintf(stderr, "INFO: added info to the surface.\n");
        }
      }
      else {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
          fprintf(stderr,
                  "INFO: cannot find mri/orig volume to "
                  "get c_(ras) information.\n");
          fprintf(stderr,
                  "INFO: transform src volume "
                  "information cannot be found. assume c_(ras) = 0\n");
          fprintf(stderr,
                  "INFO: destination surface "
                  "points may be shifted in the volume.\n");
          fprintf(stderr,
                  "INFO: you should put the "
                  "src info in the transform.\n");
        }
        lt->src.c_r = 0;
        lt->src.c_a = 0;
        lt->src.c_s = 0;
      }
    }
  }
  else  // lt->src.valid == 1
  {
    // verify
    if (mris->vg.valid) {
      if (!FZERO(lt->src.c_r - mris->vg.c_r) || !FZERO(lt->src.c_a - mris->vg.c_a) ||
          !FZERO(lt->src.c_s - mris->vg.c_s)) {
        fprintf(stderr,
                "WARNING: the source volume info "
                "is not consistent between the info contained\n");
        fprintf(stderr,
                "WARNING: in the surface data (%f,%f,%f) and "
                "that of the transform (%f,%f,%f).\n",
                mris->vg.c_r,
                mris->vg.c_a,
                mris->vg.c_s,
                lt->src.c_r,
                lt->src.c_a,
                lt->src.c_s);
      }
    }
  }
  // check dst info
  if (!lt->dst.valid) {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      fprintf(stderr,
              "INFO: transform dst volume information "
              "cannot be found.\n");
      fprintf(stderr,
              "INFO: if the target is MNI average_305, "
              "then you can do 'setenv USE_AVERAGE305 true'\n");
      fprintf(stderr,
              "INFO: otherwise c_(ras) is set to 0. "
              "destination surface points may be shifted.\n");
    }
    if (getenv("USE_AVERAGE305")) {
      lt->dst.c_r = -0.0950;
      lt->dst.c_a = -16.5100;
      lt->dst.c_s = 9.7500;
    }
  }
  // cache the transform
  mris->SRASToTalSRAS_ = getSRASToTalSRAS(lt);
  mris->TalSRASToSRAS_ = MatrixInverse(mris->SRASToTalSRAS_, NULL);

  // mark to make sure it is freed
  mris->free_transform = 1;

  mrisCheckVertexFaceTopology(mris);
  
  return (NO_ERROR);
}

// public
int MRISreadTransform(MRIS *mris, const char *fname) { 
    return mrisReadTransform(mris, fname); 
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadBinaryCurvature(MRI_SURFACE *mris, const char *mris_fname)
{
  char fname[STRLEN], fpref[STRLEN], hemi[20];

  FileNamePath(mris_fname, fpref);
  strcpy(hemi, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
  int req = snprintf(fname, STRLEN, "%s/%s.curv", fpref, hemi); 
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  return (MRISreadCurvatureFile(mris, fname));
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisReadAsciiCurvatureFile(MRI_SURFACE *mris, const char *fname)
{
  FILE *fp;
  int vno;
  char line[STRLEN], *cp;
  VERTEX *v;

  fp = fopen(fname, "r");
  if (!fp)
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "%s could not open file %s.\n", mrisReadAsciiCurvatureFile, fname));
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    cp = fgetl(line, 100, fp);
    if (!cp) {
      break;
    }
    if (sscanf(line, "%*d %*f %*f %*f %f\n", &v->curv) != 1)
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "mrisReadAsciiCurvatureFile(%s): "
                   "could not scan curvature from line '%s'",
                   fname,
                   line));
  }

  fclose(fp);
  return (NO_ERROR);
}

int MRISreadCurvatureFile(MRI_SURFACE *mris, const char *sname)
{
  int k, i, vnum, fnum;
  float curv = 0, curvmin, curvmax;
  FILE *fp;
  const char *cp;
  char path[STRLEN], fname[STRLEN], type;
  int mritype, frame, nv, c, r, s, vno;
  MRI *TempMRI;





  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    if (getenv("FS_POSIX")) {
      // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
      int req = snprintf(fname, STRLEN, "./%s", sname);    
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else {
      cp = strchr(sname, '.');
      FileNamePath(mris->fname, path);
      if (cp && ((strncmp(cp - 2, "lh", 2) == 0) || (strncmp(cp - 2, "rh", 2) == 0))) {
        int req = snprintf(fname, STRLEN, "%s/%s", path, sname); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else /* no hemisphere specified */
      {
        int req = snprintf(fname, STRLEN, "%s/%s.%s", 
			   path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}

      }
    }
  }
  else {
    strcpy(fname, sname); /* path specified explicitly */
  }
  mritype = mri_identify(sname);
  if (mritype == GIFTI_FILE) {
    mris = mrisReadGIFTIfile(fname, mris);
    if (mris) {
      return (NO_ERROR);
    }
    else {
      return (ERROR_BADFILE);
    }
  }
  if (mritype == VTK_FILE) {
    mris = MRISreadVTK(mris, fname);
    return (NO_ERROR);
  }
  if (mritype != MRI_VOLUME_TYPE_UNKNOWN) {
    frame = MRISgetReadFrame();
    TempMRI = MRIreadHeader(fname, mritype);
    if (TempMRI == NULL) {
      return (ERROR_BADFILE);
    }
    if (TempMRI->nframes <= frame) {
      printf("ERROR: attempted to read frame %d from %s\n", frame, fname);
      printf("  but this file only has %d frames.\n", TempMRI->nframes);
      return (ERROR_BADFILE);
    }
    nv = TempMRI->width * TempMRI->height * TempMRI->depth;
    if (nv != mris->nvertices) {
      printf("ERROR: number of vertices in %s does not match surface (%d,%d)\n", sname, nv, mris->nvertices);
      return (1);
    }
    MRIfree(&TempMRI);
    TempMRI = MRIread(fname);
    if (TempMRI == NULL) {
      return (ERROR_BADFILE);
    }
    vno = 0;
    curvmin = 10000.0f;
    curvmax = -10000.0f; /* for compiler warnings */
    for (s = 0; s < TempMRI->depth; s++) {
      for (r = 0; r < TempMRI->height; r++) {
        for (c = 0; c < TempMRI->width; c++) {
          curv = MRIgetVoxVal(TempMRI, c, r, s, frame);
          if (s == 0 && r == 0 && c == 0) {
            curvmin = curvmax = curv;
          }
          if (curv > curvmax) {
            curvmax = curv;
          }
          if (curv < curvmin) {
            curvmin = curv;
          }
          mris->vertices[vno].curv = curv;
          vno++;
        }
      }
    }
    MRIfree(&TempMRI);
    mris->max_curv = curvmax;
    mris->min_curv = curvmin;
    return (NO_ERROR);
  }

  type = MRISfileNameType(fname);
  if (type == MRIS_ASCII_TRIANGLE_FILE) {
    return (mrisReadAsciiCurvatureFile(mris, fname));
  }
  else if (type == MRIS_GIFTI_FILE) {
    mris = mrisReadGIFTIfile(fname, mris);
    if (mris) {
      return (NO_ERROR);
    }
    else {
      return (ERROR_BADFILE);
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "reading curvature file...");
  }

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadCurvature: could not open %s", fname));

  fread3(&vnum, fp);
  if (vnum == NEW_VERSION_MAGIC_NUMBER) {
    fclose(fp);
    return (MRISreadNewCurvatureFile(mris, fname));
  }

  fread3(&fnum, fp);
  if (vnum != mris->nvertices) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadBinaryCurvature: incompatible vertex "
                 "number in file %s",
                 fname));
  }
  curvmin = 10000.0f;
  curvmax = -10000.0f; /* for compiler warnings */
  for (k = 0; k < vnum; k++) {
    fread2(&i, fp);
    curv = i / 100.0;

    if (k == 0) {
      curvmin = curvmax = curv;
    }
    if (curv > curvmax) {
      curvmax = curv;
    }
    if (curv < curvmin) {
      curvmin = curv;
    }
    mris->vertices[k].curv = curv;
  }
  mris->max_curv = curvmax;
  mris->min_curv = curvmin;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "done. min=%2.3f max=%2.3f\n", curvmin, curvmax);
  }
  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
float *MRISreadCurvatureVector(MRI_SURFACE *mris, const char *sname)
{
  const char *cp;
  char path[STRLEN], fname[STRLEN];
  float *cvec = NULL;
  int return_code = ERROR_NONE;

  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.');
    FileNamePath(mris->fname, path);
    if (cp) {
      int req = snprintf(fname, STRLEN, "%s/%s", path, sname); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else {
      /* no hemisphere specified */
      int req = snprintf(fname, STRLEN, "%s/%s.%s",
			 path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }
  else {
    strcpy(fname, sname); /* path specified explcitly */
  }

  /* Try to read an array of values from this file. If we get an
     error, return NULL. */
  return_code = MRISreadCurvatureIntoArray(fname, mris->nvertices, &cvec);
  if (NO_ERROR != return_code) {
    return NULL;
  }

  /* Return the array we read. */
  return cvec;
}

int MRISreadCurvatureIntoArray(const char *sname, int in_array_size, float **out_array)
{
  int k, i, vnum, fnum;
  float *cvec;
  FILE *fp;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "reading curvature file...");
  }

  fp = fopen(sname, "r");
  if (fp == NULL)
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadCurvatureVectorIntoArray(%s): fopen failed", sname));

  fread3(&vnum, fp);
  if (vnum == NEW_VERSION_MAGIC_NUMBER) {
    fclose(fp);
    return (MRISreadNewCurvatureIntoArray(sname, in_array_size, out_array));
  }

  /* If that wasn't the magic number, it should be our number of
     vertices, and should match the array size we were passed. */
  if (vnum != in_array_size) {
    fclose(fp);
    return (ERROR_BADFILE);
  }

  /* This value is ignored. */
  fread3(&fnum, fp);

  /* Allocate vector to size of vnum. */
  cvec = (float *)calloc(in_array_size, sizeof(float));
  if (!cvec) {
    fclose(fp);
    ErrorReturn(ERROR_NOMEMORY, (ERROR_NOMEMORY, "MRISreadCurvatureVectorIntoArray(%s): calloc failed", sname));
  }

  for (k = 0; k < vnum; k++) {
    fread2(&i, fp);
    cvec[k] = i / 100.0;
  }
  fclose(fp);

  /* Return what we read. */
  *out_array = cvec;

  return (ERROR_NONE);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadFloatFile(MRI_SURFACE *mris, const char *sname)
{
  int k, vnum, fnum;
  float f;
  FILE *fp;
  const char *cp;
  char path[STRLEN], fname[STRLEN];

  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.');
    FileNamePath(mris->fname, path);
    if (cp) {
      int req = snprintf(fname, STRLEN, "%s/%s", path, sname);  
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else { /* no hemisphere specified */
      int req = snprintf(fname, STRLEN, "%s/%s.%s",
			 path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname);
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }
  else {
    strcpy(fname, sname); /* path specified explcitly */
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "reading float file...");
  }

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadFloatFile: could not open %s", fname));

  vnum = freadInt(fp);
  fnum = freadInt(fp);
  if (vnum != mris->nvertices) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadFloatFile: incompatible # of vertices "
                 "in file %s",
                 fname));
  }
  if (fnum != mris->nfaces) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadFloatFile: incompatible # of faces "
                 "file %s",
                 fname));
  }
  for (k = 0; k < vnum; k++) {
    f = freadFloat(fp);
    mris->vertices[k].val = f;
  }
  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISreadBinaryAreas(MRI_SURFACE *mris, const char *mris_fname)
{
  int k, vnum, fnum;
  float f;
  FILE *fp;
  char fname[STRLEN], fpref[STRLEN], hemi[20];

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "reading area file...");
  }

  FileNamePath(mris_fname, fpref);
  strcpy(hemi, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
  int req = snprintf(fname, STRLEN, "%s/%s.area", fpref, hemi);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  /*  mris->orig_area = 0.0f ;*/
  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRISreadBinaryAreas: no area file %s\n", fname));
  fread3(&vnum, fp);
  fread3(&fnum, fp);
  if (vnum != mris->nvertices) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadBinaryAreas: incompatible vertex "
                 "number in file %s",
                 fname));
  }

  for (k = 0; k < vnum; k++) {
    f = freadFloat(fp);
    mris->vertices[k].origarea = f;
    /*    mris->orig_area += f;*/
  }
  fclose(fp);

  return (NO_ERROR);
}

SMALL_SURFACE *MRISreadVerticesOnly(char *fname)
{
  SMALL_SURFACE *mriss = NULL;
  int type, magic, version, ix, iy, iz, nquads, nvertices, vno;
  SMALL_VERTEX *vertex;
  FILE *fp;

  type = MRISfileNameType(fname);
  switch (type) {
    case MRIS_ASCII_TRIANGLE_FILE:
    case MRIS_ICO_FILE:
    case MRIS_GEO_TRIANGLE_FILE:
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRISreadVerticesOnly: file type %d not supported", type));
      break; /* not used */
    default:
      break;
  }

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "MRISread(%s): could not open file", fname));

  fread3(&magic, fp);
  if (magic == QUAD_FILE_MAGIC_NUMBER) {
    version = -1;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      fprintf(stdout, "new surface file format\n");
    }
  }
  else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) {
    version = -2;
  }
  else if (magic == TRIANGLE_FILE_MAGIC_NUMBER) {
    fclose(fp);
    mriss = mrisReadTriangleFileVertexPositionsOnly(fname);
    if (!mriss) {
      ErrorReturn(NULL, (Gerror, "mrisReadTriangleFile failed.\n"));
    }
    version = -3;
  }
  else {
    rewind(fp);
    version = 0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      printf("surfer: old surface file format\n");
    }
  }
  if (version >= -2) /* some type of quadrangle file */
  {
    fread3(&nvertices, fp);
    fread3(&nquads, fp); /* # of qaudrangles - not triangles */

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "reading %d vertices and %d faces.\n", nvertices, 2 * nquads);

    mriss = (SMALL_SURFACE *)calloc(1, sizeof(SMALL_SURFACE));
    if (!mriss) ErrorReturn(NULL, (ERROR_NOMEMORY, "MRISreadVerticesOnly: could not allocate surface"));

    mriss->nvertices = nvertices;
    mriss->vertices = (SMALL_VERTEX *)calloc(nvertices, sizeof(SMALL_VERTEX));
    if (!mriss->nvertices) {
      free(mriss);
      ErrorReturn(NULL, (ERROR_NOMEMORY, "MRISreadVerticesOnly: could not allocate surface"));
    }
    for (vno = 0; vno < nvertices; vno++) {
      vertex = &mriss->vertices[vno];
      if (version == -1) {
        fread2(&ix, fp);
        fread2(&iy, fp);
        fread2(&iz, fp);
        vertex->x = ix / 100.0;
        vertex->y = iy / 100.0;
        vertex->z = iz / 100.0;
      }
      else /* version == -2 */
      {
        vertex->x = freadFloat(fp);
        vertex->y = freadFloat(fp);
        vertex->z = freadFloat(fp);
      }
    }
  }

  return (mriss);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteTriangularSurface(MRI_SURFACE *mris, const char *fname)
{
  const char *user = getenv("USER");
  if (!user)  user = getenv("LOGNAME");
  if (!user)  user = "UNKNOWN";

  auto cdt = currentDateTime();
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "writing surface file %s, created by %s on %s.\n", fname, user, cdt.c_str());

  FILE *fp = fopen(fname, "w");
  if (fp == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISwriteTriangularSurface(%s): can't create file\n", fname));
  
  fwrite3(TRIANGLE_FILE_MAGIC_NUMBER, fp);
  fprintf(fp, "created by %s on %s\n\n", user, cdt.c_str());
  fwriteInt(mris->nvertices, fp);
  fwriteInt(mris->nfaces, fp); /* # of triangles */

  for (int k = 0; k < mris->nvertices; k++) {
    fwriteFloat(mris->vertices[k].x, fp);
    fwriteFloat(mris->vertices[k].y, fp);
    fwriteFloat(mris->vertices[k].z, fp);
  }
  for (int k = 0; k < mris->nfaces; k++) {
    for (int n = 0; n < VERTICES_PER_FACE; n++) {
      fwriteInt(mris->faces[k].v[n], fp);
    }
  }
  /* write whether vertex data was using
     the real RAS rather than conformed RAS */
  fwriteInt(TAG_OLD_USEREALRAS, fp);
  fwriteInt(mris->useRealRAS, fp);

  fwriteInt(TAG_OLD_SURF_GEOM, fp);
  writeVolGeom(fp, &mris->vg);

  // write other tags
  if (!FZERO(mris->group_avg_surface_area)) {
    long long here;
    printf("writing group avg surface area %2.0f cm^2 into surface file\n", mris->group_avg_surface_area / 100.0);
    TAGwriteStart(fp, TAG_GROUP_AVG_SURFACE_AREA, &here, sizeof(float));
    fwriteFloat(mris->group_avg_surface_area, fp);
    TAGwriteEnd(fp, here);
  }
  {
    for (int i = 0; i < mris->ncmds; i++) TAGwrite(fp, TAG_CMDLINE, mris->cmdlines[i], strlen(mris->cmdlines[i]) + 1);
  }
  
  fclose(fp);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static SMALL_SURFACE *mrisReadTriangleFileVertexPositionsOnly(const char *fname)
{
  SMALL_VERTEX *v;
  int nvertices, nfaces, magic, vno;
  char line[STRLEN];
  FILE *fp;
  SMALL_SURFACE *mriss;

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "mrisReadTriangleFile(%s): could not open file", fname));

  fread3(&magic, fp);
  fgets(line, 200, fp);
  fscanf(fp, "\n");
  /*  fscanf(fp, "\ncreated by %s on %s\n", user, time_str) ;*/
  nvertices = freadInt(fp);
  nfaces = freadInt(fp);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "surface %s: %d vertices and %d faces.\n", fname, nvertices, nfaces);

  mriss = (SMALL_SURFACE *)calloc(1, sizeof(SMALL_SURFACE));
  if (!mriss) ErrorReturn(NULL, (ERROR_NOMEMORY, "MRISreadVerticesOnly: could not allocate surface"));

  mriss->nvertices = nvertices;
  mriss->vertices = (SMALL_VERTEX *)calloc(nvertices, sizeof(SMALL_VERTEX));
  if (!mriss->nvertices) {
    free(mriss);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRISreadVerticesOnly: could not allocate surface"));
  }
  for (vno = 0; vno < nvertices; vno++) {
    v = &mriss->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    v->x = freadFloat(fp);
    v->y = freadFloat(fp);
    v->z = freadFloat(fp);
    if (fabs(v->x) > 10000 || !std::isfinite(v->x))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d x coordinate %f!", Progname, vno, v->x);
    if (fabs(v->y) > 10000 || !std::isfinite(v->y))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d y coordinate %f!", Progname, vno, v->y);
    if (fabs(v->z) > 10000 || !std::isfinite(v->z))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d z coordinate %f!", Progname, vno, v->z);
  }
  fclose(fp);
  return (mriss);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisReadTriangleFilePositions(MRI_SURFACE *mris, const char *fname)
{
  int nvertices, nfaces, magic, vno;
  char line[STRLEN];
  FILE *fp;

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "mrisReadTriangleFile(%s): could not open file", fname));

  fread3(&magic, fp);
  fgets(line, 200, fp);
  fscanf(fp, "\n");
  /*  fscanf(fp, "\ncreated by %s on %s\n", user, time_str) ;*/
  nvertices = freadInt(fp);
  nfaces = freadInt(fp);

  if (nvertices != mris->nvertices || nfaces != mris->nfaces)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "mrisReadTriangleFile opened %s okay but surface doesn't match %s.  nvertices:%d != "
                 "mris->nvertices:%d || nfaces:%d != mris->nfaces:%d\n",
                 fname,
                 mris->fname ? mris->fname : "NULL",
                 nvertices,
                 mris->nvertices,
                 nfaces,
                 mris->nfaces));

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "surface %s: %d vertices and %d faces.\n", fname, nvertices, nfaces);

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  for (vno = 0; vno < nvertices; vno++) {
    float x = freadFloat(fp);
    float y = freadFloat(fp);
    float z = freadFloat(fp);
    MRISsetXYZ(mris, vno, x,y,z);
  }

  fclose(fp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static MRI_SURFACE *mrisReadTriangleFile(const char *fname, double nVFMultiplier)
{
  FACE *f;
  int nvertices, nfaces, magic, vno, fno, n;
  char line[STRLEN];
  FILE *fp;
  int tag;

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "mrisReadTriangleFile(%s): could not open file", fname));

  fread3(&magic, fp);
  fgets(line, 200, fp);
  fscanf(fp, "\n");
  /*  fscanf(fp, "\ncreated by %s on %s\n", user, time_str) ;*/
  nvertices = freadInt(fp);
  nfaces    = freadInt(fp);
  
  if (nvertices < 0 || nfaces < 0) {
    fflush(stdout);
    fflush(stderr);
    fprintf(stderr, "%s:%d %s freadInt returned nvertices:%d \n", __FILE__,__LINE__,fname,nvertices);
    fprintf(stderr, "%s:%d %s freadInt returned nfaces:%d \n", __FILE__,__LINE__,fname,nfaces);
    exit(1);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "surface %s: %d vertices and %d faces.\n", fname, nvertices, nfaces);

  MRIS * mris = MRISoverAlloc(nVFMultiplier * nvertices, nVFMultiplier * nfaces, nvertices, nfaces);
  mris->type = MRIS_TRIANGULAR_SURFACE;

  for (vno = 0; vno < nvertices; vno++) {
    if (vno % 100 == 0) exec_progress_callback(vno, nvertices, 0, 1);
    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
    VERTEX          * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    float x = freadFloat(fp);
    float y = freadFloat(fp);
    float z = freadFloat(fp);
    
    MRISsetXYZ(mris,vno, x, y, z);

    vt->num = 0; /* will figure it out */
    if (fabs(v->x) > 10000 || !std::isfinite(v->x))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d x coordinate %f!", Progname, vno, v->x);
    if (fabs(v->y) > 10000 || !std::isfinite(v->y))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d y coordinate %f!", Progname, vno, v->y);
    if (fabs(v->z) > 10000 || !std::isfinite(v->z))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d z coordinate %f!", Progname, vno, v->z);
  }

  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      f->v[n] = freadInt(fp);
      if (f->v[n] >= mris->nvertices || f->v[n] < 0)
        ErrorExit(ERROR_BADFILE, "f[%d]->v[%d] = %d - out of range!\n", fno, n, f->v[n]);
    }

    for (n = 0; n < VERTICES_PER_FACE; n++) {
      mris->vertices_topology[mris->faces[fno].v[n]].num++;
    }
  }
  // new addition
  mris->useRealRAS = 0;

  // read tags
  {
    long long len;

    while ((tag = TAGreadStart(fp, &len)) != 0) {
      switch (tag) {
        case TAG_GROUP_AVG_SURFACE_AREA:
          mris->group_avg_surface_area = freadFloat(fp);
          break;
        case TAG_OLD_SURF_GEOM:
          readVolGeom(fp, &mris->vg);
          break;
        case TAG_OLD_USEREALRAS:
          if (!freadIntEx(&mris->useRealRAS, fp))  // set useRealRAS
          {
            mris->useRealRAS = 0;  // if error, set to default
          }
          break;
        case TAG_CMDLINE:
          if (mris->ncmds > MAX_CMDS)
            ErrorExit(ERROR_NOMEMORY, "MRISread(%s): too many commands (%d) in file", fname, mris->ncmds);
          mris->cmdlines[mris->ncmds] = (char *)calloc(len + 1, sizeof(char));
          if (mris->cmdlines[mris->ncmds] == NULL)
            ErrorExit(ERROR_NOMEMORY, "MRISread(%s): could not allocate %d byte cmdline", fname, len);
          mris->cmdlines[mris->ncmds][len] = 0;
          fread(mris->cmdlines[mris->ncmds], sizeof(char), len, fp);
          mris->ncmds++;
          break;
        default:
          TAGskip(fp, tag, (long long)len);
          break;
      }
    }
  }

  fclose(fp);

  // IT IS NOT YET COMPLETE mrisCheckVertexFaceTopology(mris);

  return (mris);
}
/*!
\fn int MRISbuildFileName_read(MRI_SURFACE *mris, const char *sname, char *fname)
\brief This function "builds" a file path (fname) for reading a surface.
\param mris  - input MRI_SURFACE struct
\param sname - input surface name, it is in one of these three formats: <lh|rh.surfacename>, <surfacename>, or <full path to surfacename>
\param fname - output surface name. 
\description This function should be used to contruct surface read path only. 
The surface read path is computed as following:
1. If sname is <full path to surfacename> (determined by a forward slash "/" in sname), copy sname to fname, and return;
2. if hemisphere part is missing in the surface name (sname is in <surfacename> format), the prefix is taken from mris->hemisphere;
3. if environment variable FS_POSIX is set, use cwd as the path, otherwise, use where mris->fname is in.
*/
int MRISbuildFileName_read(MRI_SURFACE *mris, const char *sname, char *fname)
{
  char path[STRLEN];
  const char *slash, *dot;

  slash = strchr(sname, '/');
  if(!slash)   {
    /* no path - use same one as mris was read from */
    dot = strchr(sname, '.');
    FileNamePath(mris->fname, path);
    if (dot && (*(dot - 1) == 'h') && (*(dot - 2) == 'l' || *(dot - 2) == 'r')) {
      if (getenv("FS_POSIX")) {
        // PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
        int req = snprintf(fname, STRLEN, "./%s", sname); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
        // PW 2017/05/15: Legacy behaviour
        int req = snprintf(fname, STRLEN, "%s/%s", path, sname); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
    }
    else {/* no hemisphere specified */
      const char *hemistr;
      if(mris->hemisphere == LEFT_HEMISPHERE)       hemistr = "lh.";
      else if(mris->hemisphere == RIGHT_HEMISPHERE) hemistr = "rh.";
      else if(mris->hemisphere == BOTH_HEMISPHERES)  hemistr = "both.";
      else hemistr = "";
      if (getenv("FS_POSIX")) {
	// PW 2017/05/15: If FS_POSIX is set, write to cwd (as per POSIX:4.11)
	int req = snprintf(fname, STRLEN,"./%s%s",hemistr,sname);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      else {
	// PW 2017/05/15: Legacy behaviour
	sprintf(fname,"%s/%s%s",path,hemistr,sname);
      }
    }
  }
  else {
    strcpy(fname, sname); /* path specified explicitly */
  }
  return (NO_ERROR);
}

int MRISwriteDecimation(MRI_SURFACE *mris, char *fname)
{
  int k;
  FILE *fptr;

  fptr = fopen(fname, "w");
  if (fptr == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISwriteDecimation: could not create %s", fname));
  fputc('\0', fptr);
  fwriteInt(mris->nvertices, fptr);
  for (k = 0; k < mris->nvertices; k++) {
    if (mris->vertices[k].d == 0) {
      fputc('\1', fptr);
    }
    else {
      fputc('\0', fptr);
    }
  }
  fclose(fptr);
  return (NO_ERROR);
}
int MRISreadDecimation(MRI_SURFACE *mris, char *fname)
{
  int k, d, ndec;
  char c;
  FILE *fptr;

  ndec = 0;
  for (k = 0; k < mris->nvertices; k++) {
    mris->vertices[k].undefval = TRUE;
    mris->vertices[k].fixedval = FALSE;
  }
  fptr = fopen(fname, "r");
  if (fptr == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadDecimation: could not create %s", fname));
  c = fgetc(fptr);
  if (c == '#') {
    fscanf(fptr, "%*s");
    fscanf(fptr, "%d", &d);
    if (d != mris->nvertices)
      ErrorReturn(0, (ERROR_BADFILE, "%s: decimation file %s has wrong # of vertices\n", Progname, fname, d));
    for (k = 0; k < mris->nvertices; k++) {
      fscanf(fptr, "%d", &d);
      if (d != 0) {
        mris->vertices[k].d = 0;
        mris->vertices[k].fixedval = TRUE;
        mris->vertices[k].undefval = FALSE;
        ndec++;
      }
    }
  }
  else {
    d = freadInt(fptr);
    if (d != mris->nvertices)
      ErrorReturn(0, (ERROR_BADFILE, "%s: decimation file %s has wrong # of vertices\n", Progname, fname, d));
    for (k = 0; k < mris->nvertices; k++) {
      c = fgetc(fptr);
      if (c != '\0') {
        mris->vertices[k].d = 0;
        mris->vertices[k].fixedval = TRUE;
        mris->vertices[k].undefval = FALSE;
        ndec++;
      }
    }
  }
  fclose(fptr);
  return (ndec);
}

int MRISreadNewCurvatureFile(MRI_SURFACE *mris, const char *sname)
{
  int k, vnum, fnum, vals_per_vertex;
  float curv, curvmin, curvmax;
  FILE *fp;
  const char *cp;
  char path[STRLEN], fname[STRLEN];

  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.');
    FileNamePath(mris->fname, path);
    if (cp) {
      int req = snprintf(fname, STRLEN, "%s/%s", path, sname); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else {
      /* no hemisphere specified */
      int req = snprintf(fname, STRLEN, "%s/%s.%s",
			 path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }
  else {
    strcpy(fname, sname); /* path specified explcitly */
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "reading curvature file...");
  }

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "MRISreadBinaryCurvature: could not open %s", fname));

  fread3(&vnum, fp);
  if (vnum != NEW_VERSION_MAGIC_NUMBER) {
    fclose(fp);
    return (MRISreadCurvatureFile(mris, fname));
  }

  vnum = freadInt(fp);
  fnum = freadInt(fp);
  if (vnum != mris->nvertices) {
    fclose(fp);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "MRISreadNewCurvature: incompatible vertex "
                 "number in file %s",
                 fname));
  }
  vals_per_vertex = freadInt(fp);
  if (vals_per_vertex != 1) {
    fclose(fp);
    ErrorReturn(
        ERROR_NOFILE,
        (ERROR_NOFILE, "MRISreadNewCurvature(%s): vals/vertex %d unsupported (must be 1) ", fname, vals_per_vertex));
  }

  curvmin = 10000.0f;
  curvmax = -10000.0f; /* for compiler warnings */
  for (k = 0; k < vnum; k++) {
    curv = freadFloat(fp);
    if (k == 0) {
      curvmin = curvmax = curv;
    }
    if (curv > curvmax) {
      curvmax = curv;
    }
    if (curv < curvmin) {
      curvmin = curv;
    }
    mris->vertices[k].curv = curv;
  }
  mris->max_curv = curvmax;
  mris->min_curv = curvmin;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "done. min=%2.3f max=%2.3f\n", curvmin, curvmax);
  }
  fclose(fp);
  return (NO_ERROR);
}
float *MRISreadNewCurvatureVector(MRI_SURFACE *mris, const char *sname)
{
  const char *cp;
  char path[STRLEN], fname[STRLEN];
  float *cvec = NULL;
  int return_code = ERROR_NONE;

  cp = strchr(sname, '/');
  if (!cp) /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.');
    FileNamePath(mris->fname, path);
    if (cp) {
      int req = snprintf(fname, STRLEN, "%s/%s", path, sname);  
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    else { /* no hemisphere specified */
      int req = snprintf(fname, STRLEN, "%s/%s.%s",
			 path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname);  
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
  }
  else {
    strcpy(fname, sname); /* path specified explcitly */
  }

  /* Try to read an array of values from this file. If we get an
     error, return NULL. */
  return_code = MRISreadNewCurvatureIntoArray(fname, mris->nvertices, &cvec);
  if (NO_ERROR != return_code) {
    return NULL;
  }

  /* Return the array we read. */
  return cvec;
}

int MRISreadNewCurvatureIntoArray(const char *sname, int in_array_size, float **out_array)
{
  int k, vnum, fnum;
  float *cvec;
  FILE *fp;
  int vals_per_vertex;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "reading curvature file...");
  }

  fp = fopen(sname, "r");
  if (fp == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadNewCurvatureIntoArray(%s): fopen failed"));

  /* This is the version number. */
  fread3(&vnum, fp);
  if (vnum != NEW_VERSION_MAGIC_NUMBER) {
    fclose(fp);
    return (MRISreadCurvatureIntoArray(sname, in_array_size, out_array));
  }

  /* Read number of vertices and faces. */
  vnum = freadInt(fp);
  fnum = freadInt(fp);

  /* Make sure the number of vertices mathces what we're expecting. */
  if (vnum != in_array_size) {
    fclose(fp);
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,
                 "MRISreadNewCurvatureIntoArray(%s): number of vertices (%d) doesn't match what was expected (%d)",
                 sname,
                 vnum,
                 in_array_size));
  }

  /* Number of values per vertex. */
  vals_per_vertex = freadInt(fp);
  if (vals_per_vertex != 1) {
    fclose(fp);
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "MRISreadNewCurvatureIntoArray(%s): vals_per_vertex was not 1", sname));
  }

  /* Allocate to size of vnum.*/
  cvec = (float *)calloc(in_array_size, sizeof(float));
  if (!cvec) ErrorExit(ERROR_NOMEMORY, "MRISreadNewCurvatureVector(%s): calloc failed", sname);

  /* Read in values. */
  for (k = 0; k < vnum; k++) {
    cvec[k] = freadFloat(fp);
  }
  fclose(fp);

  /* Return what we read. */
  *out_array = cvec;

  return (ERROR_NONE);
}


int MRISwriteCropped(MRI_SURFACE *mris, const char *fname)
{
  float *vals;

  vals = (float *)calloc(mris->nvertices, sizeof(*vals));
  MRISexportValVector(mris, vals);

  MRIScopyFromCropped(mris, VERTEX_VALS);
  MRISwriteValues(mris, fname);
  MRISimportValVector(mris, vals);

  free(vals);
  return (NO_ERROR);
}


int MRISwriteMarked(MRI_SURFACE *mris, const char *sname)
{
  float *curv_save;

  curv_save = (float *)calloc(mris->nvertices, sizeof(float));
  if (!curv_save) ErrorExit(ERROR_NOMEMORY, "MRISwriteMarked: could not alloc %d vertex curv storage", mris->nvertices);

  MRISextractCurvatureVector(mris, curv_save);
  MRISmarkedToCurv(mris);
  MRISwriteCurvature(mris, sname);
  MRISimportCurvatureVector(mris, curv_save);
  free(curv_save);
  return (NO_ERROR);
}

int MRISreadMarked(MRI_SURFACE *mris, const char *sname)
{
  float *curv_save;

  curv_save = (float *)calloc(mris->nvertices, sizeof(float));
  if (!curv_save) ErrorExit(ERROR_NOMEMORY, "MRISwriteMarked: could not alloc %d vertex curv storage", mris->nvertices);

  MRISextractCurvatureVector(mris, curv_save);
  if (MRISreadCurvatureFile(mris, sname) != NO_ERROR) {
    return (Gerror);
  }
  MRIScurvToMarked(mris);
  MRISimportCurvatureVector(mris, curv_save);
  free(curv_save);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwriteArea(MRI_SURFACE *mris, const char *sname)
{
  float *curv_save;

  curv_save = (float *)calloc(mris->nvertices, sizeof(float));
  if (!curv_save) ErrorExit(ERROR_NOMEMORY, "MRISwriteArea: could not alloc %d vertex curv storage", mris->nvertices);

  MRISextractCurvatureVector(mris, curv_save);
  MRISareaToCurv(mris);
  MRISwriteCurvature(mris, sname);
  MRISimportCurvatureVector(mris, curv_save);
  free(curv_save);
  return (NO_ERROR);
}


MRI *MRISreadParameterizationToSurface(MRI_SURFACE *mris, char *fname)
{
  // ideally, canonical vertices have been saved, but there is no good way to check this
  // besides probing whether a vertex c-xyz is nonzero
  VERTEX *v = &mris->vertices[0];
  bool has_canonical = !(FZERO(v->cx) && FZERO(v->cy) && FZERO(v->cz));

  // if the canonical vertices have been saved, let's use those - otherwise let's hope that
  // the surface's current vertices are already spherical
  if (has_canonical) {
    MRISsaveVertexPositions(mris, TMP_VERTICES);
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
    MRIScomputeMetricProperties(mris);
  }

  MRI_SP *mrisp = MRISPread(fname) ;
  if (mrisp == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "MRISreadParameterizationToSurface: could not open file %s",fname));

  int nframes = mrisp->Ip->num_frame;
  MRI *mri = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, nframes) ;

  for (int frame = 0 ; frame < nframes ; frame++)
  {
    MRISfromParameterizationBarycentric(mrisp, mris, frame) ;
    for (int vno = 0 ; vno < mris->nvertices ; vno++)
      MRIsetVoxVal(mri, vno, 0, 0, frame, mris->vertices[vno].curv) ;
  }

  // restore the previous vertices if we modified them
  if (has_canonical) MRISrestoreVertexPositions(mris, TMP_VERTICES);

  MRISPfree(&mrisp) ;
  MRIScomputeMetricProperties(mris) ;
  return(mri) ;
}


int MatlabPlotFace(FILE *fp, MRIS *surf, int faceno, char color, double NormLen)
{
  FACE *face = &(surf->faces[faceno]);
  int n,m,vno1,vno2;
  VERTEX *v1,*v2;
  double cx=0,cy=0,cz=0;
  int fontsize=7;
  char normcolor = 'k';

  printf("%% Face %d\n",faceno);
  fprintf(fp,"plot3(");
  for(n=0; n<3; n++){
    m = n + 1;
    if(m>2) m = 0;
    vno1 = face->v[n];
    vno2 = face->v[m];
    v1 = &(surf->vertices[vno1]);
    v2 = &(surf->vertices[vno2]);
    fprintf(fp,"[%6.4f,%6.4f],[%6.4f,%6.4f],[%6.4f,%6.4f],'%c' ",
	    v1->x,v2->x,v1->y,v2->y,v1->z,v2->z,color);
    if(n!=2) fprintf(fp,",");
    cx += v1->x;
    cy += v1->y;
    cz += v1->z;
  }
  cx /= 3;
  cy /= 3;
  cz /= 3;
  fprintf(fp,");\n");
  fprintf(fp,"h=text(%6.4f,%6.4f,%6.4f,'F%d');\n",cx,cy,cz,faceno);
  fprintf(fp,"set(h,'fontsize',%d,'color','%c');\n",fontsize,color);
  if(fabs(NormLen) > 0){
    float snorm[3];
    mrisNormalFace(surf, faceno, 0, snorm);
    //printf("%g %g %g\n",snorm[0],snorm[1],snorm[2]);
    fprintf(fp,"hold on;\n");
    fprintf(fp,"plot3([%6.4f],[%6.4f],[%6.4f],'%c*',",cx,cy,cz,color);
    fprintf(fp,"[%6.4f,%6.4f],[%6.4f,%6.4f],[%6.4f,%6.4f],'%c'",
	    cx,cx+NormLen*snorm[0],cy,cy+NormLen*snorm[1],cz,cz+NormLen*snorm[2],
	    normcolor);
    fprintf(fp,");\n");
    fprintf(fp,"hold off;\n");
  }
  fflush(fp);
  return(0);
}

int MatlabPlotVertex(FILE *fp, MRIS *surf, int vno, char color, double NormLen)
{
  VERTEX *v = &(surf->vertices[vno]);
  int fontsize=7;
  char normcolor = 'k';

  printf("%% Vertex %d\n",vno);
  fprintf(fp,"plot3([%6.4f],[%6.4f],[%6.4f],'%c*')\n",v->x,v->y,v->z,color);
  fprintf(fp,"h=text(%6.4f,%6.4f,%6.4f,'V%d');\n",v->x,v->y,v->z,vno);
  fprintf(fp,"set(h,'fontsize',%d,'color','%c');\n",fontsize,color);
  if(fabs(NormLen) > 0){
    fprintf(fp,"hold on;\n");
    fprintf(fp,"plot3([%6.4f,%6.4f],[%6.4f,%6.4f],[%6.4f,%6.4f],'%c')\n",
	    v->x,v->x+NormLen*v->nx, v->y,v->y+NormLen*v->ny, v->z,v->z+NormLen*v->nz, normcolor);
    fprintf(fp,"hold off;\n");
  }
  fflush(fp);
  return(0);
}

int MatlabPlotVertexNbhd(FILE *fp, MRIS *surf, int cvno, int nhops, char color, double NormLen)
{

  int nthhop, nnbrs, nthnbr, nbrvtxno, faceno,nthface;
  SURFHOPLIST *shl;
  shl = SetSurfHopList(cvno, surf, nhops);

  for(nthhop = 0; nthhop < nhops; nthhop++) {
    nnbrs = shl->nperhop[nthhop];
    // loop through the neighbors nthhop links away
    for(nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      nbrvtxno = shl->vtxlist[nthhop][nthnbr];
      fprintf(fp,"hold on;\n");
      MatlabPlotVertex(fp, surf, nbrvtxno, 'g', NormLen);
      if(nthhop >= nhops-1) continue;
      VERTEX_TOPOLOGY *vt = &(surf->vertices_topology[nbrvtxno]);
      for(nthface=0; nthface <  vt->num; nthface++){
	faceno = vt->f[nthface];
	fprintf(fp,"hold on;\n");
	MatlabPlotFace(fp, surf, faceno, 'b', NormLen);
	//fprintf(stderr,"%2d %3d %d %6d %6d\n",nthhop,nthnbr,nthface,nbrvtxno,faceno);
      }

    } /* end loop over hop neighborhood */
  }   /* end loop over hop */

  // Have to do this at the end
  fprintf(fp,"hold on;\n");
  MatlabPlotVertex(fp, surf, cvno, 'r', NormLen);
  fprintf(fp,"title('Vertex %d');\n",cvno);
  fprintf(fp,"hold off;\n");

  SurfHopListFree(&shl);
  return(0);
}

/*!
  \fn int MRISwriteField(MRIS *surf, char **fields, int nfields, char *outname)
  \brief Converts data in the given fields into "voxel" in an MRI structure
  using MRIcopyMRIS() and writes to the given volume-style (eg, mgz) output file. 
  Field names must be known to MRIcopyMRIS(). 
*/
int MRISwriteField(MRIS *surf, const char **fields, int nfields, const char *outname)
{
  int n;
  MRI *mri=NULL;
  for(n=nfields-1; n >= 0; n--){
    mri = MRIcopyMRIS(mri, surf, n, fields[n]);
    if(mri==NULL) exit(1);
  }
  int err = MRIwrite(mri,outname);
  MRIfree(&mri);
  return(err);
}


