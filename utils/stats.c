#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "error.h"
#include "diag.h"
#include "matrix.h"
#include "stats.h"
#include "const.h"
#include "proto.h"
#include "utils.h"
#include "machine.h"
#include "mrinorm.h"

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
  reg->mri2fmri = MatrixAlloc(REG_ROWS, REG_COLS, MATRIX_REAL) ;
  fgetl(reg->name, MAX_LINE_LEN-1, fp) ;

  if (!stricmp(reg->name, "margaret"))
    DiagBreak() ;

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
      if (fscanf(fp, "%f  ", &reg->mri2fmri->rptr[row][col]) != 1)
        ErrorReturn(NULL,
                    (ERROR_BADFILE, 
                  "StatReadRegistration(%s): could not scan element (%d, %d)",
                     fname, row, col)) ;
    }
    fscanf(fp, "\n") ;
  }
  fclose(fp) ;
  reg->fmri2mri = MatrixInverse(reg->mri2fmri, NULL) ;
  if (!reg->fmri2mri)
    ErrorExit(ERROR_BADPARM, "StatReadReg(%s): singular registration matrix",
              fname) ;

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
  MatrixFree(&reg->mri2fmri) ;
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
  char         path[STRLEN], fname[STRLEN], line[MAX_LINE_LEN], *cp ;
  STAT_VOLUME  *sv ;
  FILE         *fp ;
  int          dof_mean, dof_sigma, event_number, slice_number, which_alloc,
               width, height, nframes, nslices, t, event, nitems, x, y, z ;
  float        *buf, fval ;
  int DatVersion, DOF;
  float TER;
  char  *regfile = NULL;

  FileNamePath(prefix, path) ;
  sv = (SV *)calloc(1, sizeof(SV)) ;
  if (!sv)
    ErrorExit(ERROR_NOMEMORY, "StatReadVolume(%s): could not allocate sv",
              prefix) ;

  /* read in register.dat */
  if(regfile != NULL)  sprintf(fname, "%s", regfile);
  else                 sprintf(fname, "%s/register.dat", path) ;

  sv->reg = StatReadRegistration(fname) ;
  if (!sv->reg)
    return(NULL) ;

  /* read the selavg/selxavg dat file, if it exists */
  sprintf(fname, "%s.dat", prefix) ;
  fp = fopen(fname, "r") ;
  which_alloc = ALLOC_MEANS ;
  if (fp)  /* means there are time points and means and sigmas */
  {
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
    if(fscanf(fp,"%*s %d",&DatVersion) != EOF){
      fscanf(fp,"%*s %f",&TER);
      fscanf(fp,"%*s %d",&DOF);
      printf("SelXAvg Format TER = %g, DOF =  %d\n",TER,DOF);
      sv->voltype = 2;
      fprintf(stderr,"INFO: detected volume %s as type selxavg\n",prefix);
    }
    else{
      sv->voltype = 1;
      fprintf(stderr,"INFO: detected volume %s as type selavg\n",prefix);
    }
    fclose(fp) ;
    which_alloc |= ALLOC_STDS ;
  }
  else
  {
    /*fprintf(stderr,"WARNING: %s: StatReadVolume():\n",Progname);
    fprintf(stderr,"%s does not exist\n",fname);*/
    fprintf(stderr,"INFO: detected volume %s as type raw\n",prefix);
    sv->nevents = 1 ;
    sv->time_per_event = 0 ;  /* will be filled in later by .hdr file */
    sv->voltype = 0;
  }

  if(sv->nevents > MAX_EVENTS){
    fprintf(stderr,"ERROR: %s, StatReadVolume():\n",Progname);
    fprintf(stderr,"Number of events (%d) exceeds maximum allowed (%d)\n",
      sv->nevents, MAX_EVENTS);
    exit(1);
  }

  if(sv->voltype == 1){
    /* read in the dof file */
    sprintf(fname, "%s_000.dof", prefix) ;
    fp = fopen(fname, "r") ;
    if (fp){
      while ((cp = fgetl(line, MAX_LINE_LEN-1, fp)) != NULL){
  sscanf(cp, "%d %d %d", &event_number, &dof_mean, &dof_sigma) ;
  sv->mean_dofs[event_number] = (float)dof_mean ;
  sv->std_dofs[event_number] = (float)dof_sigma ;
      }
      fclose(fp) ;
    }
    else{
      fprintf(stderr,"WARNING: %s: StatReadVolume():\n",Progname);
      fprintf(stderr,"%s does not exist\n",fname);
    }
  }
  else{
    if(sv->voltype == 0) DOF = 1; /* for raw type */
    for(event_number = 0; event_number < sv->nevents; event_number++){
      sv->mean_dofs[event_number] = (float)DOF+1;
      sv->std_dofs[event_number]  = (float)DOF;
    }
  }

  /* count # of slices */
  slice_number = 0 ;
  do
  {
    sprintf(fname, "%s_%3.3d.hdr", prefix, slice_number) ;
    fp = fopen(fname, "r") ;
    if (fp)   /* this is a valid slice */
    {
      fclose(fp) ;
      slice_number++ ;
    }
  } while (fp) ;

  sv->nslices = nslices = slice_number ;

  /* first read .hdr file to get slice dimensions */
  sprintf(fname, "%s_%3.3d.hdr", prefix, 0) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "StatReadVolume: could not open %s",fname) ;
  fscanf(fp, "%d %d %d", &width, &height, &nframes) ;
  if (!sv->time_per_event)
    sv->time_per_event = nframes ;  /* no global .dat file */
  fclose(fp) ;

  nframes = sv->time_per_event * sv->nevents ;
  if (which_alloc & ALLOC_STDS)
    nframes *= 2.0f ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "reading %d slices, %d images per slice\n",
            nslices, nframes);

  buf = (float *)calloc(width*height, sizeof(float)) ;
  if (!buf)
    ErrorExit(ERROR_NO_MEMORY, "StatReadVolume: could not allocate buffer") ;

  nframes = sv->time_per_event ;  /* one frame per time point */

  StatAllocVolume(sv, sv->nevents, width,height,nslices,sv->time_per_event,
                  which_alloc);

#if 0
  for(w=0;w<width;w++){
    for(h=0;h<height;h++){
      for(s=0;s<nslices;s++){
  for(t=0;t<sn->time_per_event;t++){

  }
      }
    }
  }
#endif

  /* read it after nevents */
  if (stricmp(sv->reg->name, "talairach") && 
      stricmp(sv->reg->name, "spherical"))
    StatReadTransform(sv, sv->reg->name) ;

  /* read in the actual data */
  nitems = width * height ;
  for (z = 0 ; z < nslices ; z++){

    /* first read .hdr file to get slice dimensions */
    sprintf(fname, "%s_%3.3d.hdr", prefix, z) ;
    fp = fopen(fname, "r") ;
    if (!fp)
      ErrorReturn(NULL,
                  (ERROR_NOFILE, "StatReadVolume: could not open %s",fname)) ;
    fscanf(fp, "%d %d %d", &width, &height, &nframes) ;
    sv->slice_width = width ; sv->slice_height = height ;
    fclose(fp) ;

    /* now read .bfloat file to parse actual data */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z) ;
    fp = fopen(fname, "r") ;
    if (!fp)
    {
      sprintf(fname, "%s_%3.3d.bshort", prefix, z) ;
      if (!fp)
      {
        ErrorReturn(NULL,
                    (ERROR_NOFILE, "StatReadVolume: could not open %s",fname));
      }
      else
        fprintf(stderr,"ERROR: %s does not support bshort volumes\n",Progname);
    }

    for (event = 0 ; event < sv->nevents ; event++)
    {
      /* read slice of means for each time point */
      for (t = 0 ; t < sv->time_per_event ; t++)
      {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(NULL, (ERROR_BADFILE, "StatReadVolume: could not read "
                           "%dth slice",z)) ;
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            fval = buf[y*width+x] ;
#ifdef Linux
            fval = swapFloat(fval) ;
#endif
            MRIFseq_vox(sv->mri_avgs[event],x,y,z,t) = fval ;
          }
        }
      }

      /* read slice of standard deviations for each time point */
      if (sv->mri_stds[event]) for (t = 0 ; t < sv->time_per_event ; t++)
      {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(NULL, (ERROR_BADFILE, "StatReadVolume: could not read "
                           "%dth slice",z)) ;
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            fval = buf[y*width+x] ;
#ifdef Linux
            fval = swapFloat(fval) ;
#endif
            MRIFseq_vox(sv->mri_stds[event],x,y,z,t) = fval ;
          }
        }
      }
    }

    fclose(fp) ;


    if (!sv->mri_avg_dofs[event])
      continue ;

    /* if dof files exist for each slice, read them in */
    sprintf(fname, "%s_dof_%3.3d.bfloat", prefix, z) ;
    fp = fopen(fname, "r") ;
    if (!fp)
      ErrorReturn(NULL,
                  (ERROR_NOFILE, "StatReadVolume: could not open %s",fname)) ;

    for (event = 0 ; event < sv->nevents ; event++)
    {
      /* read slice of mean dofs for each time point */
      for (t = 0 ; t < sv->time_per_event ; t++)
      {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(NULL, (ERROR_BADFILE, "StatReadVolume: could not read "
                           "%dth slice",z)) ;
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            fval = buf[y*width+x] ;
#ifdef Linux
            fval = swapFloat(fval) ;
#endif
            MRIFseq_vox(sv->mri_avg_dofs[event],x,y,z,t) = fval ;
          }
        }
      }

      /* read slice of standard deviation dofs for each time point */
      if (sv->mri_stds[event]) for (t = 0 ; t < sv->time_per_event ; t++)
      {
        if (fread(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(NULL, (ERROR_BADFILE, "StatReadVolume: could not read "
                           "%dth slice",z)) ;
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            fval = buf[y*width+x] ;
#ifdef Linux
            fval = swapFloat(fval) ;
#endif
            MRIFseq_vox(sv->mri_std_dofs[event],x,y,z,t) = fval ;
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
    if (sv->mri_stds[event])
      MRIsetResolution(sv->mri_stds[event], sv->reg->in_plane_res, 
                       sv->reg->in_plane_res, sv->reg->slice_thickness) ;
      
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
    if(sv->voltype != 0)
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
                int nslices, int time_points, int which_alloc)
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
    if (which_alloc & ALLOC_STDS)
    {
      sv->mri_stds[event] = 
        MRIallocSequence(width, height, nslices, MRI_FLOAT, time_points) ;
      if (!sv->mri_stds[event])
        ErrorExit(ERROR_NO_MEMORY, "StatAllocVolume: could not allocate "
                  "volume");
    }
    if (which_alloc & ALLOC_DOFS)
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
StatAllocStructuralVolume(SV *sv, float fov, float resolution, char *name)
{
  SV     *sv_tal ;
  int    width, height, depth, event ;

  width = height = depth = nint(fov / resolution) ;
  sv_tal = StatAllocVolume(NULL, sv->nevents, width, height, depth, 
                     sv->time_per_event,ALLOC_MEANS|ALLOC_STDS|ALLOC_DOFS);

  sv_tal->voltype = sv->voltype;
  strcpy(sv_tal->reg->name, name) ;
  sv_tal->nslices = depth ;
  sv_tal->slice_width = width ;
  sv_tal->slice_height = height ;
  sv_tal->reg->slice_thickness = sv_tal->reg->in_plane_res = resolution ;
  sv_tal->reg->brightness_scale = 2.0f ;

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
StatAccumulateSurfaceVolume(SV *sv_surf, SV *sv, MRI_SURFACE *mris)
{
  int    x, y, z, width, height, depth, event, t, xv, yv, zv, swidth, sheight,
         sdepth, vno ;
  Real   xf, yf, zf, xr, yr, zr ;
  float  mean, surf_mean, std, surf_std, surf_dof, dof, xoff, yoff, zoff, 
         sxoff, syoff, szoff, xs, ys, zs ;
  VECTOR *v_struct, *v_func ;
  MRI    *mri_avg, *mri_std, *mri_ctrl ;
  VERTEX *vertex ;

  v_func = VectorAlloc(4, MATRIX_REAL) ;
  v_struct = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_func, 4) = VECTOR_ELT(v_struct, 4) = 1.0f ;

  width = sv_surf->mri_avgs[0]->width ;
  height = sv_surf->mri_avgs[0]->height ;
  depth = sv_surf->mri_avgs[0]->depth ;

  swidth = sv->mri_avgs[0]->width ;
  sheight = sv->mri_avgs[0]->height ;
  sdepth = sv->mri_avgs[0]->depth ;

  xoff = (float)(width-1)/2.0f ;
  yoff = (float)(height-1)/2.0f ;
  zoff = (float)(depth-1)/2.0f ;

  sxoff = (float)(sv->slice_width-1)/2.0f ;
  syoff = (float)(sv->slice_height-1)/2.0f ;
  szoff = (float)(sv->nslices-1)/2.0f ;

  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR) ;
  mri_avg = MRIallocSequence(width, height, depth, MRI_FLOAT, 2) ;
  MRIcopyHeader(sv_surf->mri_avgs[0], mri_avg) ;
  mri_std = MRIallocSequence(width, height, depth, MRI_FLOAT, 2) ;
  MRIcopyHeader(sv_surf->mri_stds[0], mri_std) ;

  if (sv_surf->nevents != sv->nevents)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, 
                 "StatAccumulateTalairachVolume: inconsistent # of events "
                 "(%d vs %d)", sv_surf->nevents, sv->nevents)) ;

  for (event = 0 ; event < sv_surf->nevents ; event++)
  {
    sv_surf->mean_dofs[event] += sv->mean_dofs[event] ;
    sv_surf->std_dofs[event] += sv->std_dofs[event] ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "\rprocessing event %d of %d....", 
              event+1, sv_surf->nevents) ;

    for (t = 0 ; t < sv_surf->time_per_event ; t++)
    {
      MRIclear(mri_avg) ; MRIclear(mri_std) ; MRIclear(mri_ctrl) ;
/* 
   sample from the surface, through the strucural space, into the functional 
   one. 
 */

/* 
   first build average volume for this subjects. This is to avoid sampling
   from the volume onto the surface which would take forever since it
   would require searching the entire surface for the nearest point for
   every point in the volume.
   */
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        vertex = &mris->vertices[vno] ;

/* 
   read the functional data from the coordinates in 'folded' space, then
   write them into the structural volume using the canonical coordinates.
*/
        xs = vertex->x ; ys = vertex->y ; zs = vertex->z ;

        /* transform it into functional space */
        VECTOR3_LOAD(v_struct, xs, ys, zs) ;
        MatrixMultiply(sv->reg->mri2fmri, v_struct, v_func) ;

        /* transform it into a (functional) voxel coordinate */
        MRIworldToVoxel(sv->mri_avgs[event], (Real)V3_X(v_func), 
                        (Real)V3_Y(v_func), (Real)V3_Z(v_func), &xf, &yf, &zf);
        xv = nint(xf) ; yv = nint(yf) ; zv = nint(zf) ;

        if (xv >= 0 && xv < swidth &&
            yv >= 0 && yv < sheight &&
            zv >= 0 && zv < sdepth)
        {
          /* convert from canonical surface coordinate to voxel coordinate */
          xs = vertex->cx ; ys = vertex->cy ; zs = vertex->cz ;
          MRIworldToVoxel(sv_surf->mri_avgs[event], xs, ys, zs, &xr, &yr, &zr);
          x = nint(xr) ; y = nint(yr) ; z = nint(zr) ;
          if ((vno == 161) ||  (x == 15 && y == 20 && z == 3 && t == 0))
            DiagBreak() ;

          if ((x < 0) || (y < 0) || (z < 0) || (z >= mri_avg->depth) || 
              (y >= mri_avg->height) || (x >= mri_avg->width))
            continue ;
              
/* 
   at this point (xv,yv,zv) is a functional coordinate, and (x,y,z) is
   the corresponding structural coordinate.
 */
          /* update means */
          surf_mean = MRIFseq_vox(mri_avg, x, y, z, 0) ;
          surf_dof  = MRIFseq_vox(mri_avg, x, y, z, 1) ;
          mean = MRIFseq_vox(sv->mri_avgs[event], xv, yv, zv, t) ;
          surf_mean = (surf_mean * surf_dof + mean) / (surf_dof+1) ;
          MRIFseq_vox(mri_avg, x, y, z, 0) = surf_mean ;
          MRIFseq_vox(mri_avg, x, y, z, 1) = ++surf_dof ;
          MRIvox(mri_ctrl, x, y, z) = CONTROL_MARKED ;

#if 0         
          if (x == 15 && y == 20 && z == 3 && t == 0)
            fprintf(stderr, "adding mean %2.3f, avg = %2.3f\n",mean,surf_mean);
#endif

          /* update stds */
          surf_std = MRIFseq_vox(mri_std, x, y, z, 0) ;
          surf_dof = MRIFseq_vox(mri_std, x, y, z, 1) ;
          std = MRIFseq_vox(sv->mri_stds[event], xv, yv, zv, t) ;
          
          /* work with variances so things are linear */
          surf_std *= surf_std ; std *= std ;
          surf_std = sqrt((surf_std * surf_dof + std)/(surf_dof+1));
          MRIFseq_vox(mri_std, x, y, z, 0) = surf_std ;
          MRIFseq_vox(mri_std, x, y, z, 1) = ++surf_dof ;
        }
      }
      if (getenv("SOAP_STATS"))
      {
        int niter = atoi(getenv("SOAP_STATS")) ;

        fprintf(stderr, "performing soap bubble for %d iterations...\n", niter) ;
        /*        MRIsoapBubbleExpand(mri_avg, mri_ctrl, mri_avg, 1) ;*/
        MRIbuildVoronoiDiagram(mri_avg, mri_ctrl, mri_avg) ;
        MRIsoapBubble(mri_avg, mri_ctrl, mri_avg, niter) ;
        MRIbuildVoronoiDiagram(mri_avg, mri_ctrl, mri_avg) ;
        /*        MRIsoapBubbleExpand(mri_std, mri_ctrl, mri_std, 1) ;*/
        MRIsoapBubble(mri_std, mri_ctrl, mri_std, niter) ;
      }

/* 
   now use the intermediate structural volumes to update the
   cross-subject average volume.
   */
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            if (x == 15 && y == 20 && z == 3 && t == 0)
              DiagBreak() ;

            /* update means */
            surf_mean = MRIFseq_vox(sv_surf->mri_avgs[event], x, y, z, t) ;
            surf_dof = MRIFseq_vox(sv_surf->mri_avg_dofs[event], x, y, z, t) ;
            mean = MRIFvox(mri_avg, x, y, z) ;
            dof = sv->mean_dofs[event] ;
            surf_mean = (surf_mean * surf_dof + mean * dof) / (surf_dof + dof);
            surf_dof += dof ;
            MRIFseq_vox(sv_surf->mri_avg_dofs[event], x, y, z, t) = surf_dof ;
            MRIFseq_vox(sv_surf->mri_avgs[event], x, y, z, t) = surf_mean ;
          
            /* update stds */
            surf_std = MRIFseq_vox(sv_surf->mri_stds[event], x, y, z, t) ;
            surf_dof = MRIFseq_vox(sv_surf->mri_std_dofs[event], x, y, z, t) ;
            std = MRIFvox(mri_std, x, y, z) ;
            dof = sv->std_dofs[event] ;
          
            /* work with variances so things are linear */
            surf_std *= surf_std ; std *= std ;
            surf_std = sqrt((surf_std * surf_dof + std * dof)/(surf_dof+dof));
            surf_dof += dof ;
            MRIFseq_vox(sv_surf->mri_std_dofs[event], x, y, z, t) = surf_dof ;
            MRIFseq_vox(sv_surf->mri_stds[event], x, y, z, t) = surf_std ;
          }
        }
      }
    }
  }

  MRIfree(&mri_std) ; MRIfree(&mri_avg) ;

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
StatAccumulateTalairachVolume(SV *sv_tal, SV *sv)
{
  int    x, y, z, width, height, depth, event, t, xv, yv, zv, swidth, sheight,
         sdepth ;
  Real   xf, yf, zf ;
  float  mean, tal_mean, std, tal_std, tal_dof, dof, xoff, yoff, zoff, sxoff, 
         syoff, szoff ;
  VECTOR *v_struct, *v_func ;
  MRI    *mri_avg, *mri_std ;

  if(!sv){
    fprintf(stderr,"ERROR: %s: StatAccumulateTalairachVolume():\n",Progname);
    fprintf(stderr,"Input stat volume is null\n");
    fprintf(stderr,"%s: %d\n",__FILE__,__LINE__);
    exit(1);
  }

  v_func   = VectorAlloc(4, MATRIX_REAL) ;
  v_struct = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_func, 4) = VECTOR_ELT(v_struct, 4) = 1.0f ;

  width  = sv_tal->mri_avgs[0]->width ;
  height = sv_tal->mri_avgs[0]->height ;
  depth  = sv_tal->mri_avgs[0]->depth ;

  swidth  = sv->mri_avgs[0]->width ;
  sheight = sv->mri_avgs[0]->height ;
  sdepth  = sv->mri_avgs[0]->depth ;

  xoff = (float)(sv_tal->mri_avgs[0]->width-1)/2.0f ;
  yoff = (float)(sv_tal->mri_avgs[0]->height-1)/2.0f ;
  zoff = (float)(sv_tal->mri_avgs[0]->depth-1)/2.0f ;

  sxoff = (float)(sv->slice_width-1)/2.0f ;
  syoff = (float)(sv->slice_height-1)/2.0f ;
  szoff = (float)(sv->nslices-1)/2.0f ;

  if (sv_tal->nevents != sv->nevents)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, 
                 "StatAccumulateTalairachVolume: inconsistent # of events "
                 "(%d vs %d)", sv_tal->nevents, sv->nevents)) ;

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
            if (x == 35 && y == 35 && z == 6)
              DiagBreak() ;

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
              if (xv == 30 && yv == 51 && zv == 2 && t == 2)
                DiagBreak() ;

              /* update means */
              tal_mean = MRIFseq_vox(sv_tal->mri_avgs[event], x, y, z, t) ;
              mean = MRIFseq_vox(sv->mri_avgs[event], xv, yv, zv, t) ;
              dof = sv->mean_dofs[event] ;
              tal_dof = MRIFseq_vox(sv_tal->mri_avg_dofs[event], x, y, z, t) ;
              tal_mean = (tal_mean * tal_dof + mean * dof) / (tal_dof + dof) ;
              tal_dof += dof ;
              MRIFseq_vox(sv_tal->mri_avg_dofs[event], x, y, z, t) = tal_dof ;
              MRIFseq_vox(sv_tal->mri_avgs[event], x, y, z, t) = tal_mean ;

        if(sv->voltype != 0) {
    tal_std = MRIFseq_vox(sv_tal->mri_stds[event], x, y, z, t) ;
    std = MRIFseq_vox(sv->mri_stds[event], xv, yv, zv, t) ;
    dof = sv->std_dofs[event] ;
    tal_dof = MRIFseq_vox(sv_tal->mri_std_dofs[event], x, y, z, t) ;
    
    /* work with variances so things are linear */
    tal_std *= tal_std ; std *= std ;
    tal_std = sqrt((tal_std * tal_dof + std * dof)/(tal_dof + dof));
    tal_dof += dof ;
    MRIFseq_vox(sv_tal->mri_std_dofs[event], x, y, z, t) = tal_dof ;
    MRIFseq_vox(sv_tal->mri_stds[event], x, y, z, t) = tal_std ;
        }

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
  char         path[STRLEN], fname[STRLEN] ;
  FILE         *fp ;
  int          event_number, width, height, nslices, t, 
               event, nitems, x, y, z, nframes ;
  float        *buf, fval ;

  width = sv->slice_width ; height = sv->slice_height ; nslices = sv->nslices ;
  FileNamePath(prefix, path) ;
  sprintf(fname, "%s/register.dat", path) ;
  StatWriteRegistration(sv->reg, fname) ;

  if(sv->voltype != 0) {
    /* write the global header file */
    sprintf(fname, "%s.dat", prefix) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, 
         "StatWriteVolume: could not open dat file %s", fname)) ;
    fprintf(fp, "tr %f\n", sv->tr) ;
    fprintf(fp, "timewindow %f\n", sv->timewindow) ;
    fprintf(fp, "prestim %f\n", sv->prestim) ;
    fprintf(fp, "nbins %d\n", sv->nevents) ;
    fprintf(fp, "perevent %d\n", sv->time_per_event) ;
    fclose(fp) ;
  }

  buf = (float *)calloc(sv->slice_width*sv->slice_height, sizeof(float)) ;
  if (!buf)
    ErrorExit(ERROR_NO_MEMORY, "StatWriteVolume: could not allocate buffer") ;

  if(sv->voltype == 0)  nframes =   sv->time_per_event*sv->nevents;
  else                  nframes = 2*sv->time_per_event*sv->nevents;

  /* write out the dof files, the .hdr files and the .bfloat data files */
  nitems = width * height ;
  for (z = 0 ; z < sv->nslices ; z++)
  {
    if(sv->voltype == 1) {
      /* write out .dof file */
      sprintf(fname, "%s_%3.3d.dof", prefix, z) ;
      fp = fopen(fname, "w") ;
      if (!fp)
  ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,"StatWriteVolume: could not open "
          "dof file %s",fname));
      for (event_number = 0 ; event_number < sv->nevents ; event_number++)
  fprintf(fp, "%d %2.0f %2.0f\n", event_number, 
    sv->mean_dofs[event_number], sv->std_dofs[event_number]) ;
      fclose(fp) ;
    }

    /* write out .hdr file */
    sprintf(fname, "%s_%3.3d.hdr", prefix, z) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "StatWriteVolume: could not open %s",fname) ;
    fprintf(fp, "%d %d %d 0\n", sv->slice_width, sv->slice_height, nframes);
    fclose(fp) ;

    /* now write out actual data into .bfloat files */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "StatWriteVolume: could not open %s",fname) ;

    for (event = 0 ; event < sv->nevents ; event++) {
      /* write slice of means for each time point */
      for (t = 0 ; t < sv->time_per_event ; t++)
      {
        for (y = 0 ; y < sv->slice_height ; y++)
        {
          for (x = 0 ; x < sv->slice_width ; x++)
          {
            fval = MRIFseq_vox(sv->mri_avgs[event],x,y,z,t) ;
#ifdef Linux
            fval = swapFloat(fval) ;
#endif
            buf[y*width+x] = fval ;
          }
        }
        if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
          ErrorReturn(ERROR_BADFILE, 
                      (ERROR_BADFILE, "StatWriteVolume: could not write "
                       "%dth slice",z)) ;
      }

      if(sv->voltype != 0){
  /* write slice of standard deviations for each time point */
  for (t = 0 ; t < sv->time_per_event ; t++){
    for (y = 0 ; y < height ; y++){
      for (x = 0 ; x < width ; x++){
        fval = MRIFseq_vox(sv->mri_stds[event],x,y,z,t) ;
#ifdef Linux
        fval = swapFloat(fval) ;
#endif
        buf[y*width+x] = fval ;
      }
    }
    if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
      ErrorReturn(ERROR_BADFILE, 
      (ERROR_BADFILE, "StatWriteVolume: could not write "
       "%dth slice", z)) ;
  }
      }
    }

    fclose(fp) ;

    if (sv->mri_avg_dofs[0]) /* keeping track of dofs on a per-voxel basis */
    {
#if 0 /* don't write out _dof_ file */
      /* now write out dofs on a per-voxel basis */
      sprintf(fname, "%s_dof_%3.3d.bfloat", prefix, z) ;
      fp = fopen(fname, "w") ;
      if (!fp)
        ErrorExit(ERROR_NOFILE, "StatWriteVolume: could not open %s",fname) ;
      
      for (event = 0 ; event < sv->nevents ; event++)
      {
        /* write slice of means for each time point */
        for (t = 0 ; t < sv->time_per_event ; t++)
        {
          for (y = 0 ; y < sv->slice_height ; y++)
          {
            for (x = 0 ; x < sv->slice_width ; x++)
            {
              fval = MRIFseq_vox(sv->mri_avg_dofs[event],x,y,z,t) ;
#ifdef Linux
              fval = swapFloat(fval) ;
#endif
              buf[y*width+x] = fval ;
            }
          }
          if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(ERROR_BADFILE, 
                        (ERROR_BADFILE, "StatWriteVolume: could not write "
                         "%dth slice",z)) ;
        }
        
        /* write slice of standard deviations for each time point */
        for (t = 0 ; t < sv->time_per_event ; t++)
        {
          for (y = 0 ; y < height ; y++)
          {
            for (x = 0 ; x < width ; x++)
            {
              fval = MRIFseq_vox(sv->mri_std_dofs[event],x,y,z,t) ;
#ifdef Linux
              fval = swapFloat(fval) ;
#endif
              buf[y*width+x] = fval ;
            }
          }
          if (fwrite(buf, sizeof(float), nitems, fp) != nitems)
            ErrorReturn(ERROR_BADFILE, 
                        (ERROR_BADFILE, "StatWriteVolume: could not write "
                         "%dth slice", z)) ;
        }
      }
      
      fclose(fp) ;
#endif      
    }
    }

  /* now remove old files which may have had more slices */
  for ( ; z < sv->nslices*4 ; z++) {
    /* write out .dof file */
    sprintf(fname, "%s_%3.3d.bfloat", prefix, z) ;
    unlink(fname) ;
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
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, 
                       "StatWriteRegistration: could not open %s", fname)) ;
  fprintf(fp, "%s\n", reg->name) ;
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
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
int
StatReadTransform(STAT_VOLUME *sv, char *name)
{
  char  *cp, subjects[STRLEN], fname[STRLEN] ;
  int   event ;

  /* read in the Talairach transform file */
  cp = getenv("SUBJECTS_DIR") ;
  if (cp)
    strcpy(subjects, cp) ;
  else
    strcpy(subjects, "~inverse/subjects") ;
  sprintf(fname,"%s/%s/mri/transforms/talairach.xfm",subjects, name);
  
  if (input_transform_file(fname, &sv->transform) != OK)
    ErrorPrintf(ERROR_NO_FILE, 
                "%s: could not read xform file '%s'\n", Progname, fname) ;
  sv->linear_transform = get_linear_transform_ptr(&sv->transform) ;
  sv->inverse_linear_transform = 
    get_inverse_linear_transform_ptr(&sv->transform) ;
  
  for (event = 0 ; event < sv->nevents ; event++)
  {
    MRIsetTransform(sv->mri_avgs[event], &sv->transform) ;
    if (sv->mri_stds[event])
      MRIsetTransform(sv->mri_stds[event], &sv->transform) ;
  }
  return(NO_ERROR) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
int       
StatVolumeExists(char *prefix)
{
  char   fname[STRLEN] ;
  FILE   *fp ;

  sprintf(fname, "%s_%3.3d.bfloat", prefix, 0) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    return(0) ;
  fclose(fp) ;
  return(1) ;
}
