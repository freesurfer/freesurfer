#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "diag.h"
#include "matrix.h"
#include "stats.h"
#include "const.h"
#include "proto.h"
#include "utils.h"

#define REG_ROWS      4
#define REG_COLS      4
#define STRUCT_DIM    256


/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
fMRI_REG *
StatReadRegistration(char *fname)
{
  int        row, col ;
  FILE       *fp ;
  fMRI_REG   *reg ;
  char       line[MAX_LINE_LEN] ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, 
                       "StatReadRegistration: could not open %s", fname)) ;
  reg = (fMRI_REG *)calloc(1, sizeof(fMRI_REG)) ;
  reg->fmri2mri = MatrixAlloc(REG_ROWS, REG_COLS, MATRIX_REAL) ;
  fgetl(reg->name, MAX_LINE_LEN-1, fp) ;

  fgetl(line, MAX_LINE_LEN-1, fp) ;
  if (sscanf(line, "%f", &reg->in_plane_res) != 1)
    ErrorReturn(reg, (ERROR_BADFILE, "StatReadRegistration: could not scan "
                      "in plane resolution from %s", line)) ;
  fgetl(line, MAX_LINE_LEN-1, fp) ;
  if (sscanf(line, "%f", &reg->slice_thickness) != 1)
    ErrorReturn(reg, (ERROR_BADFILE, "StatReadRegistration: could not scan "
                      "slice thickness from %s", line)) ;
  fgetl(line, MAX_LINE_LEN-1, fp) ;
  if (sscanf(line, "%f", &reg->brightness_scale) != 1)
    ErrorReturn(reg, (ERROR_BADFILE, "StatReadRegistration: could not scan "
                      "brightness scale from %s", line)) ;
  
  for (row = 1 ; row <= REG_ROWS ; row++)
  {
    for (col = 1 ; col <= REG_COLS ; col++)
    {
      if (fscanf(fp, "%f  ", &reg->fmri2mri->rptr[row][col]) != 1)
        ErrorReturn(NULL,
                    (ERROR_BADFILE, 
                     "StatReadRegistration: could not scan element (%d, %d)",
                     row, col)) ;
    }
    fscanf(fp, "\n") ;
  }
  fclose(fp) ;
  reg->mri2fmri = MatrixInverse(reg->fmri2mri, NULL) ;
  if (!reg->mri2fmri)
    ErrorExit(ERROR_BADPARM, "StatReadReg: singular registration matrix") ;

  return(reg) ;
}

/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
int
StatFreeRegistration(fMRI_REG **preg)
{
  fMRI_REG   *reg ;

  reg = *preg ;
  MatrixFree(&reg->fmri2mri) ;
  free(reg) ;
  return(NO_ERROR) ;
}

/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
SV *
StatReadVolume(char *prefix)
{
  char         path[100], fname[100], line[MAX_LINE_LEN], *cp, subjects[100] ;
  STAT_VOLUME  *sv ;
  FILE         *fp ;
  int          dof_mean, dof_sigma, event_number, slice_number,
               width, height, nframes, nslices, t, event, nitems, x, y, z ;
  float        *buf ;

  FileNamePath(prefix, path) ;
  sv = (SV *)calloc(1, sizeof(SV)) ;
  if (!sv)
    ErrorExit(ERROR_NOMEMORY, "StatReadVolume(%s): could not allocate sv",
              prefix) ;

  /* read in register.dat */
  sprintf(fname, "%s/register.dat", path) ;
  sv->reg = StatReadRegistration(fname) ;

  /* read in the Talairach transform file */
  cp = getenv("SUBJECTS_DIR") ;
  if (cp)
    strcpy(subjects, cp) ;
  else
    strcpy(subjects, "~inverse/subjects") ;
  sprintf(fname, "%s/%s/mri/transforms/talairach.xfm",subjects, sv->reg->name);

  /*  fprintf(stderr, "reading transform from %s...", fname) ;*/
  if (input_transform_file(fname, &sv->transform) != OK)
    ErrorPrintf(ERROR_NO_FILE, 
                "%s: could not read xform file '%s'\n", Progname, fname) ;
  sv->linear_transform = get_linear_transform_ptr(&sv->transform) ;
  sv->inverse_linear_transform = 
    get_inverse_linear_transform_ptr(&sv->transform) ;
  /*  fprintf(stderr, "done.\n") ;*/

  /* read the global header file */
  sprintf(fname, "%s.dat", prefix) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, 
                       "StatReadVolume: could not open dat file %s", fname)) ;
  fgetl(line, MAX_LINE_LEN-1, fp) ;
  sscanf(line, "%*s %f", &sv->tr) ;
  fgetl(line, MAX_LINE_LEN-1, fp) ;
  sscanf(line, "%*s %f", &sv->timewindow) ;
  fgetl(line, MAX_LINE_LEN-1, fp) ;
  sscanf(line, "%*s %f", &sv->prestim) ;
  fgetl(line, MAX_LINE_LEN-1, fp) ;
  sscanf(line, "%*s %d", &sv->nevents) ;
  fgetl(line, MAX_LINE_LEN-1, fp) ;
  sscanf(line, "%*s %d", &sv->time_per_event) ;
  fclose(fp) ;

  /* now read in the dof file */
  sprintf(fname, "%s_000.dof", prefix) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, 
                       "StatReadVolume: could not open dof file %s", fname)) ;
  while ((cp = fgetl(line, MAX_LINE_LEN-1, fp)) != NULL)
  {
    sscanf(cp, "%d %d %d", &event_number, &dof_mean, &dof_sigma) ;
    sv->mean_dofs[event_number] = (float)dof_mean ;
    sv->std_dofs[event_number] = (float)dof_sigma ;
  }
  fclose(fp) ;

  /* count # of slices */
  slice_number = 0 ;
  do
  {
    sprintf(fname, "%s_%3.3d.bfloat", prefix, slice_number) ;
    fp = fopen(fname, "r") ;
    if (fp)   /* this is a valid slice */
    {
      fclose(fp) ;
      slice_number++ ;
    }
  } while (fp) ;

  sv->nslices = nslices = slice_number ;

  nframes = 2 * sv->time_per_event * sv->nevents ;
  fprintf(stderr, "reading %d slices, %d images per slice\n",nslices, nframes);
  /* first read .hdr file to get slice dimensions */
  sprintf(fname, "%s_%3.3d.hdr", prefix, 0) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "StatReadVolume: could not open %s",fname) ;
  fscanf(fp, "%d %d %d", &width, &height, &nframes) ;
  fclose(fp) ;

  buf = (float *)calloc(width*height, sizeof(float)) ;
  if (!buf)
    ErrorExit(ERROR_NO_MEMORY, "StatReadVolume: could not allocate buffer") ;

  nframes = sv->time_per_event ;  /* one frame per time point */

  StatAllocVolume(sv, sv->nevents, width,height,nslices,sv->time_per_event,0);

  /* read in the actual data */
  nitems = width * height ;
  for (z = 0 ; z < nslices ; z++)
  {
    /* first read .hdr file to get slice dimensions */
    sprintf(fname, "%s_%3.3d.hdr", prefix, z) ;
    fp = fopen(fname, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "StatReadVolume: could not open %s",fname) ;
    fscanf(fp, "%d %d %d", &width, &height, &nframes) ;
    sv->slice_width = width ; sv->slice_height = height ;
    fclose(fp) ;

    /* now read .bfloat file to parse actual data ANDERS */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z) ;
    fp = fopen(fname, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "StatReadVolume: could not open %s",fname) ;

    for (event = 0 ; event < sv->nevents ; event++)
    {
      /* read slice of means for each time point */
      for (t = 0 ; t < sv->time_per_event ; t++)
      {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(sv, (ERROR_BADFILE, "StatReadVolume: could not read "
                           "%dth slice",z)) ;
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            MRIFseq_vox(sv->mri_avgs[event],x,y,z,t) = buf[y*width+x] ;
          }
        }
      }

      /* read slice of standard deviations for each time point */
      for (t = 0 ; t < sv->time_per_event ; t++)
      {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(sv, (ERROR_BADFILE, "StatReadVolume: could not read "
                           "%dth slice",z)) ;
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            MRIFseq_vox(sv->mri_stds[event],x,y,z,t) = buf[y*width+x] ;
          }
        }
      }
    }

    fclose(fp) ;
  }
  for (event = 0 ; event < sv->nevents ; event++)
  {
    MRIsetResolution(sv->mri_avgs[event], sv->reg->in_plane_res, 
                     sv->reg->in_plane_res, sv->reg->slice_thickness) ;
    MRIsetResolution(sv->mri_stds[event], sv->reg->in_plane_res, 
                     sv->reg->in_plane_res, sv->reg->slice_thickness) ;
    MRIsetTransform(sv->mri_avgs[event], &sv->transform) ;
    MRIsetTransform(sv->mri_stds[event], &sv->transform) ;
  }

  free(buf) ;
  return(sv) ;
}

/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
int
StatFree(SV **psv)
{
  SV  *sv ;
  int event, width, height, nslices ;

  sv = *psv ;
  *psv = NULL ;

  width = sv->slice_width ; height = sv->slice_height ; nslices = sv->nslices ;
  for (event = 0 ; event < sv->nevents ; event++)
  {
    MRIfree(&sv->mri_avgs[event]) ;
    MRIfree(&sv->mri_stds[event]) ;
    if (sv->mri_avg_dofs[event])
      MRIfree(&sv->mri_avg_dofs[event]) ;
    if (sv->mri_std_dofs[event])
      MRIfree(&sv->mri_std_dofs[event]) ;
  }

  delete_general_transform(&sv->transform) ;
  StatFreeRegistration(&sv->reg) ;
  free(sv) ;

  return(NO_ERROR) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
#if 0
MRI *
StatVolumeToTalairach(SV *sv, MRI *mri, int resolution)
{
  if (!mri)
  {
    int dim ;

    dim = nint((float)STRUCT_DIM / (float)resolution) ;
    mri = MRIallocSequence(dim, dim, dim, MRI_FLOAT, ) ;
  }

  return(mri) ;
}
#endif
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
STAT_VOLUME *
StatAllocVolume(SV *sv, int nevents, int width, int height, 
                int nslices, int time_points, int track_dofs)
{
  int    event ;
  
  if (!sv)
  {
    sv = (SV *)calloc(1, sizeof(SV)) ;
    if (!sv)
      ErrorExit(ERROR_NOMEMORY, "StatAllocVolume: could not allocate sv") ;

    sv->reg = (fMRI_REG *)calloc(1, sizeof(fMRI_REG)) ;
    sv->reg->in_plane_res = sv->reg->slice_thickness = 1.0f ;
    strcpy(sv->reg->name, "none") ;
    sv->reg->fmri2mri = MatrixIdentity(4, NULL) ;
    sv->reg->mri2fmri = MatrixIdentity(4, NULL) ;
    sv->time_per_event = time_points ;
    sv->nevents = nevents ;
  }

  for (event = 0 ; event < sv->nevents ; event++)
  {
    sv->mri_avgs[event] = 
      MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points) ;
    if (!sv->mri_avgs[event])
      ErrorExit(ERROR_NO_MEMORY, "StatAllocVolume: could not allocate volume");
    sv->mri_stds[event] = 
      MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points) ;
    if (!sv->mri_stds[event])
      ErrorExit(ERROR_NO_MEMORY, "StatAllocVolume: could not allocate volume");
    if (track_dofs)
    {
      sv->mri_avg_dofs[event] = 
        MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points) ;
      if (!sv->mri_avg_dofs[event])
        ErrorExit(ERROR_NO_MEMORY, 
                  "StatAllocVolume: could not allocate volume");
      sv->mri_std_dofs[event] = 
        MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points) ;
      if (!sv->mri_std_dofs[event])
        ErrorExit(ERROR_NO_MEMORY, 
                  "StatAllocVolume: could not allocate volume");
    }
  }
  return(sv) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
SV *
StatAllocTalairachVolume(SV *sv, float fov, float resolution)
{
  SV     *sv_tal ;
  int    width, height, depth, event ;

  width = height = depth = nint(fov / resolution) ;
  sv_tal = StatAllocVolume(NULL, sv->nevents, width, height, depth, 
                           sv->time_per_event,1);

  strcpy(sv_tal->reg->name, "talairach") ;
  sv_tal->nslices = depth ;
  sv_tal->slice_width = width ;
  sv_tal->slice_height = height ;

  for (event = 0 ; event < sv->nevents ; event++)
  {
    MRIsetResolution(sv_tal->mri_avgs[event],resolution,resolution,resolution);
    MRIsetResolution(sv_tal->mri_stds[event],resolution,resolution,resolution);
    MRIsetResolution(sv_tal->mri_avg_dofs[event],
                     resolution,resolution,resolution);
    MRIsetResolution(sv_tal->mri_std_dofs[event],
                     resolution,resolution,resolution);
  }

  sv_tal->timewindow = sv->timewindow ;
  sv_tal->prestim = sv->prestim ;
  sv_tal->tr = sv->tr ;
  sv_tal->timewindow = sv->timewindow ;
  return(sv_tal) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
int
StatAccumulateTalairachVolume(SV *sv_tal, SV *sv)
{
  int    x, y, z, width, height, depth, event, t, xv, yv, zv, swidth, sheight,
         sdepth ;
  Real   xf, yf, zf ;
  float  mean, tal_mean, std, tal_std, tal_dof, dof, xoff, yoff, zoff, sxoff, 
         syoff, szoff ;
  VECTOR *v_struct, *v_func ;
  MRI    *mri_avg, *mri_std ;

  v_func = VectorAlloc(4, MATRIX_REAL) ;
  v_struct = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_func, 4) = VECTOR_ELT(v_struct, 4) = 1.0f ;

  width = sv_tal->mri_avgs[0]->width ;
  height = sv_tal->mri_avgs[0]->height ;
  depth = sv_tal->mri_avgs[0]->depth ;

  swidth = sv->mri_avgs[0]->width ;
  sheight = sv->mri_avgs[0]->height ;
  sdepth = sv->mri_avgs[0]->depth ;

  xoff = (float)(sv_tal->mri_avgs[0]->width-1)/2.0f ;
  yoff = (float)(sv_tal->mri_avgs[0]->height-1)/2.0f ;
  zoff = (float)(sv_tal->mri_avgs[0]->depth-1)/2.0f ;

  sxoff = (float)(sv->slice_width-1)/2.0f ;
  syoff = (float)(sv->slice_height-1)/2.0f ;
  szoff = (float)(sv->nslices-1)/2.0f ;

  for (event = 0 ; event < sv_tal->nevents ; event++)
  {
    sv_tal->mean_dofs[event] += sv->mean_dofs[event] ;
    sv_tal->std_dofs[event] += sv->std_dofs[event] ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "\rprocessing event %d of %d....", 
              event+1, sv_tal->nevents) ;
    mri_avg = sv_tal->mri_avgs[event] ;
    mri_std = sv_tal->mri_stds[event] ;
    mri_avg->linear_transform = sv->mri_avgs[event]->linear_transform;
    mri_avg->inverse_linear_transform = 
      sv->mri_avgs[event]->inverse_linear_transform;

    for (t = 0 ; t < sv_tal->time_per_event ; t++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            /* convert to subject's world coordinates (in mm) */
            MRItalairachVoxelToWorld(mri_avg, (Real)x, (Real)y, (Real)z, 
                                     &xf, &yf, &zf) ;

            /* transform it into functional space */
            VECTOR3_LOAD(v_struct, xf, yf, zf) ;
            MatrixMultiply(sv->reg->mri2fmri, v_struct, v_func) ;

            /* transform it into a voxel coordinate */
            MRIworldToVoxel(sv->mri_avgs[event], (Real)V3_X(v_func), 
                            (Real)V3_Y(v_func), (Real)V3_Z(v_func), 
                            &xf, &yf, &zf) ;

            xv = nint(xf) ; yv = nint(yf) ; zv = nint(zf) ;
            if (xv >= 0 && xv < swidth &&
                yv >= 0 && yv < sheight &&
                zv >= 0 && zv < sdepth)
            {
              /* update means */
              tal_mean = MRIFseq_vox(sv_tal->mri_avgs[event], x, y, z, t) ;
              mean = MRIFseq_vox(sv->mri_avgs[event], xv, yv, zv, t) ;
              dof = sv->mean_dofs[event] ;
              tal_dof = MRIFseq_vox(sv_tal->mri_avg_dofs[event], x, y, z, t) ;
              tal_mean = (tal_mean * tal_dof + mean * dof) / (tal_dof + dof) ;
              tal_dof += dof ;
              MRIFseq_vox(sv_tal->mri_avg_dofs[event], x, y, z, t) = tal_dof ;
              MRIFseq_vox(sv_tal->mri_avgs[event], x, y, z, t) = tal_mean ;

              /* update stds */
              tal_std = MRIFseq_vox(sv_tal->mri_stds[event], x, y, z, t) ;
              std = MRIFseq_vox(sv->mri_stds[event], xv, yv, zv, t) ;
              dof = sv->std_dofs[event] ;
              tal_dof = MRIFseq_vox(sv_tal->mri_std_dofs[event], x, y, z, t) ;
              tal_std = (tal_std * tal_dof + std * dof) / (tal_dof + dof) ;
              tal_dof += dof ;
              MRIFseq_vox(sv_tal->mri_std_dofs[event], x, y, z, t) = tal_dof ;
              MRIFseq_vox(sv_tal->mri_stds[event], x, y, z, t) = tal_std ;
            }
          }
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "done.\n") ;
  VectorFree(&v_func) ;
  VectorFree(&v_struct) ;
  return(NO_ERROR) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
int
StatWriteVolume(SV *sv, char *prefix)
{
  char         path[100], fname[100] ;
  FILE         *fp ;
  int          event_number, width, height, nslices, t, 
               event, nitems, x, y, z ;
  float        *buf ;

  width = sv->slice_width ; height = sv->slice_height ; nslices = sv->nslices ;
  FileNamePath(prefix, path) ;
  sprintf(fname, "%s/register.dat", path) ;
  StatWriteRegistration(sv->reg, fname) ;

  /* write the global header file */
  sprintf(fname, "%s.dat", prefix) ;
  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, 
                       "StatWriteVolume: could not open dat file %s", fname)) ;
  fprintf(fp, "tr %f\n", sv->tr) ;
  fprintf(fp, "timewindow %f\n", sv->timewindow) ;
  fprintf(fp, "prestim %f\n", sv->prestim) ;
  fprintf(fp, "nbins %d\n", sv->nevents) ;
  fprintf(fp, "perevent %d\n", sv->time_per_event) ;
  fclose(fp) ;

  buf = (float *)calloc(sv->slice_width*sv->slice_height, sizeof(float)) ;
  if (!buf)
    ErrorExit(ERROR_NO_MEMORY, "StatReadVolume: could not allocate buffer") ;

  /* write out the dof files, the .hdr files and the .bfloat data files */
  nitems = width * height ;
  for (z = 0 ; z < sv->nslices ; z++)
  {
    /* write out .dof file */
    sprintf(fname, "%s_%3.3d.dof", prefix, z) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorReturn(NULL, (ERROR_NOFILE, "StatWriteVolume: could not open "
                         "dof file %s",fname));
    for (event_number = 0 ; event_number < sv->nevents ; event_number++)
      fprintf(fp, "%d %2.0f %2.0f\n", event_number, 
              sv->mean_dofs[event_number], sv->std_dofs[event_number]) ;
    fclose(fp) ;

    /* write out .hdr file */
    sprintf(fname, "%s_%3.3d.hdr", prefix, z) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "StatWriteVolume: could not open %s",fname) ;
    fprintf(fp, "%d %d %d 0\n", sv->slice_width, sv->slice_height, 
            sv->time_per_event*sv->nevents*2) ;
    fclose(fp) ;

    /* now write out actual data into .bfloat files */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "StatReadVolume: could not open %s",fname) ;

    for (event = 0 ; event < sv->nevents ; event++)
    {
      /* write slice of means for each time point */
      for (t = 0 ; t < sv->time_per_event ; t++)
      {
        for (y = 0 ; y < sv->slice_height ; y++)
        {
          for (x = 0 ; x < sv->slice_width ; x++)
          {
            buf[y*width+x] = MRIFseq_vox(sv->mri_avgs[event],x,y,z,t) ;
          }
        }
        if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(ERROR_BADFILE, 
                      (ERROR_BADFILE, "StatWriteVolume: could not read "
                       "%dth slice",z)) ;
      }

      /* write slice of standard deviations for each time point */
      for (t = 0 ; t < sv->time_per_event ; t++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            buf[y*width+x] = MRIFseq_vox(sv->mri_stds[event],x,y,z,t) ;
          }
        }
        if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(ERROR_BADFILE, 
                      (ERROR_BADFILE, "StatWriteVolume: could not read "
                       "%dth slice", z)) ;
      }
    }

    fclose(fp) ;

    if (sv->mri_avg_dofs[event])
    {
      /* now write out dofs on a per-voxel basis */
      sprintf(fname, "%s_dof_%3.3d.bfloat", prefix, z) ;
      fp = fopen(fname, "w") ;
      if (!fp)
        ErrorExit(ERROR_NOFILE, "StatReadVolume: could not open %s",fname) ;
      
      for (event = 0 ; event < sv->nevents ; event++)
      {
        /* write slice of means for each time point */
        for (t = 0 ; t < sv->time_per_event ; t++)
        {
          for (y = 0 ; y < sv->slice_height ; y++)
          {
            for (x = 0 ; x < sv->slice_width ; x++)
              buf[y*width+x] = MRIFseq_vox(sv->mri_avg_dofs[event],x,y,z,t) ;
          }
          if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(ERROR_BADFILE, 
                        (ERROR_BADFILE, "StatWriteVolume: could not read "
                         "%dth slice",z)) ;
        }
        
        /* write slice of standard deviations for each time point */
        for (t = 0 ; t < sv->time_per_event ; t++)
        {
          for (y = 0 ; y < height ; y++)
          {
            for (x = 0 ; x < width ; x++)
              buf[y*width+x] = MRIFseq_vox(sv->mri_std_dofs[event],x,y,z,t) ;
          }
          if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(ERROR_BADFILE, 
                        (ERROR_BADFILE, "StatWriteVolume: could not read "
                         "%dth slice", z)) ;
        }
      }
      
      fclose(fp) ;
    }
  }
  free(buf) ;
  return(NO_ERROR) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
int
StatWriteRegistration(fMRI_REG *reg, char *fname)
{
  FILE  *fp ;
  int   row, col ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NOFILE, 
                       "StatWriteRegistration: could not open %s", fname)) ;

  fprintf(fp, "%f\n", reg->in_plane_res) ;
  fprintf(fp, "%f\n", reg->slice_thickness) ;
  fprintf(fp, "%f\n", reg->brightness_scale) ;
  for (row = 1 ; row <= REG_ROWS ; row++)
  {
    for (col = 1 ; col <= REG_COLS ; col++)
      fprintf(fp, "%f  ", reg->fmri2mri->rptr[row][col]) ;

    fprintf(fp, "\n") ;
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

