/**
 * @brief FS interface to Diffusion Toolkit and TrackVis data.
 *
 * DTk http://trackvis.org/docs/?subsect=fileformat
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double round(double x);
#include "diag.h"
#include "error.h"
#include "matrix.h"
#include "mri.h"
#include "mri2.h"
#include "mri_identify.h"
#include "proto.h"
#include "utils.h"

#include "dtk.fs.h"

#ifdef X
#undef X
#endif

double trkvisOffset = 0.5;

/*----------------------------------------------------------------*/
DTK_TRACK_SET *DTKloadTrackSet(char *trkfile, char *mrifile)
{
  DTK_TRACK_SET *dtkset;
  FILE *fp;
  size_t nread;
  int nthtrk;
  std::string stem;
  char *fname;
  MRI *mri = NULL;
  float dc, dr, ds;

  // Open the track file
  fp = fopen(trkfile, "r");
  if (fp == NULL) {
    printf("ERROR: DTKloadTrackSet(): cannot open %s\n", trkfile);
    return (NULL);
  }

  // Read in template MRI to get geometry
  if (mrifile == NULL) {
    stem = trkfile;
    stem.erase(stem.size()-4);
    fname = IDnameFromStem(stem.c_str());
    if (fname == NULL) {
      printf("WARNING: DTKloadTrackSet(): cannot find matching MRI file for %s\n", trkfile);
    }
  }
  else
    fname = mrifile;
  if (fname) {
    printf("DTKloadTrackSet(): reading geometry from %s\n", fname);
    mri = MRIreadHeader(fname, MRI_VOLUME_TYPE_UNKNOWN);
    if (mri == NULL) {
      fclose(fp);
      return (NULL);
    }
  }

  dtkset = (DTK_TRACK_SET *)calloc(sizeof(DTK_TRACK_SET), 1);
  dtkset->mri_template = mri;
  dtkset->hdr = (DTK_HDR *)calloc(sizeof(DTK_HDR), 1);

  // This should be 1000 bytes
  nread = fread(dtkset->hdr, sizeof(DTK_HDR), 1, fp);
  if (nread != 1) {
    printf("ERROR: DTKloadTrackSet(): Read only %d items\n", (int)nread);
    fclose(fp);
    return (NULL);
  }
  if (dtkset->hdr->hdr_size != 1000) {
    printf("ERROR: DTKloadTrackSet(): %s may need byte swapping\n", trkfile);
    fclose(fp);
    return (NULL);
  }
  if (dtkset->hdr->n_scalars > DTK_N_SCALARS_MAX) {
    printf(
        "ERROR: DTKloadTrackSet(): n_scalars=%d > DTK_N_SCALARS_MAX=%d\n", dtkset->hdr->n_scalars, DTK_N_SCALARS_MAX);
    fclose(fp);
    return (NULL);
  }

  DTKprintHeader(stdout, dtkset->hdr);
  printf("------------------------\n");

  dc = dtkset->hdr->voxel_size[0];
  dr = dtkset->hdr->voxel_size[1];
  ds = dtkset->hdr->voxel_size[2];

  printf("Allocating %d tracks\n", dtkset->hdr->n_count);
  dtkset->trk = (DTK_TRACK **)calloc(sizeof(DTK_TRACK *), dtkset->hdr->n_count);
  printf("Loading %d tracks\n", dtkset->hdr->n_count);
  for (nthtrk = 0; nthtrk < dtkset->hdr->n_count; nthtrk++) {
    dtkset->trk[nthtrk] = DTKreadTrack(fp, dtkset->hdr->n_scalars, dtkset->hdr->n_properties, dc, dr, ds);
    // if(nthtrk % 1000 == 0) printf("%5d %d\n",nthtrk,dtkset->trk[nthtrk]->npoints);
    if (dtkset->trk[nthtrk] == NULL) {
      printf("ERROR: DTKloadTrackSet(): loading track %d\n", nthtrk);
      fclose(fp);
      return (NULL);
    }
  }
  fclose(fp);
  DTKprintTrack(stdout, dtkset->trk[0]);
  printf("Done loading\n");
  return (dtkset);
}
/* --------------------------------------------- */
int DTKwriteTrackSet(char *trkfile, DTK_TRACK_SET *dtkset)
{
  FILE *fp;
  int nthtrk;
  float dc, dr, ds;

  dc = dtkset->hdr->voxel_size[0];
  dr = dtkset->hdr->voxel_size[1];
  ds = dtkset->hdr->voxel_size[2];

  fp = fopen(trkfile, "w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", trkfile);
    return (1);
  }
  fwrite(dtkset->hdr, sizeof(DTK_HDR), 1, fp);
  for (nthtrk = 0; nthtrk < dtkset->hdr->n_count; nthtrk++) DTKwriteTrack(fp, dtkset->trk[nthtrk], dc, dr, ds);
  fclose(fp);

  return (0);
}
/*---------------------------------------------------------------------*/
DTK_TRACK *DTKreadTrack(FILE *fp, int n_scalars, int n_properties, float dc, float dr, float ds)
{
  int npoints, n, ns;
  DTK_TRACK *trk;
  float x, y, z;

  if (fread(&npoints, sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "Could not read value(s)\n");
  }
  trk = DTKallocTrack(npoints, n_scalars, n_properties);
  if (trk == NULL) {
    return (NULL);
  }
  for (n = 0; n < trk->npoints; n++) {
    if (fread(&x, sizeof(float), 1, fp) != 1) {
      fprintf(stderr, "Could not read value(s)\n");
    }
    if (fread(&y, sizeof(float), 1, fp) != 1) {
      fprintf(stderr, "Could not read value(s)\n");
    }
    if (fread(&z, sizeof(float), 1, fp) != 1) {
      fprintf(stderr, "Could not read value(s)\n");
    }
    trk->c[n] = x / dc - trkvisOffset;
    trk->r[n] = y / dr - trkvisOffset;
    trk->s[n] = z / ds - trkvisOffset;
    for (ns = 0; ns < trk->n_scalars; ns++) {
      if (fread(&(trk->scalars[ns][n]), sizeof(float), 1, fp) != 1) {
        fprintf(stderr, "Could not read value(s)\n");
      }
    }
  }
  for (n = 0; n < trk->n_properties; n++) {
    if (fread(&(trk->properties[n]), sizeof(float), 1, fp) != 1) {
      fprintf(stderr, "Could not read value(s)\n");
    }
  }

  return (trk);
}
/*-------------------------------------------------------------*/
int DTKwriteTrack(FILE *fp, DTK_TRACK *trk, float dc, float dr, float ds)
{
  int n, ns;
  float x, y, z;

  fwrite(&(trk->npoints), sizeof(int), 1, fp);
  for (n = 0; n < trk->npoints; n++) {
    x = (trk->c[n] + trkvisOffset) * dc;
    y = (trk->r[n] + trkvisOffset) * dr;
    z = (trk->s[n] + trkvisOffset) * ds;
    fwrite(&x, sizeof(float), 1, fp);
    fwrite(&y, sizeof(float), 1, fp);
    fwrite(&z, sizeof(float), 1, fp);
    for (ns = 0; ns < trk->n_scalars; ns++) fwrite(&(trk->scalars[ns][n]), sizeof(float), 1, fp);
  }
  for (n = 0; n < trk->n_properties; n++) fwrite(&(trk->properties[n]), sizeof(float), 1, fp);

  return (0);
}
/*-----------------------------------------------------------*/
int DTKprintTrack(FILE *fp, DTK_TRACK *trk)
{
  int n, ns;
  for (n = 0; n < trk->npoints; n++) {
    fprintf(fp, "%d %f %f %f ", n, trk->c[n], trk->r[n], trk->s[n]);
    for (ns = 0; ns < trk->n_scalars; ns++) fprintf(fp, "%f ", trk->scalars[ns][n]);
    fprintf(fp, "\n");
  }
  return (0);
}
/*-------------------------------------------------------------------*/
DTK_TRACK *DTKcopyTrack(DTK_TRACK *trk)
{
  DTK_TRACK *newtrk;
  int n, ns;
  newtrk = DTKallocTrack(trk->npoints, trk->n_scalars, trk->n_properties);
  // if(trklabel)

  for (n = 0; n < trk->npoints; n++) {
    newtrk->c[n] = trk->c[n];
    newtrk->r[n] = trk->r[n];
    newtrk->s[n] = trk->s[n];
    // if(trklabel)
    for (ns = 0; ns < trk->n_scalars; ns++) newtrk->scalars[ns][n] = trk->scalars[ns][n];
  }
  for (n = 0; n < trk->n_properties; n++) newtrk->properties[n] = trk->properties[n];
  return (newtrk);
}
/*-------------------------------------------------------------------*/
DTK_TRACK *DTKallocTrack(int npoints, int n_scalars, int n_properties)
{
  DTK_TRACK *trk;
  int ns;

  if (n_scalars > DTK_N_SCALARS_MAX) {
    printf("ERROR: DTKallocTrack(): n_scalars=%d > DTK_N_SCALARS_MAX=%d\n", n_scalars, DTK_N_SCALARS_MAX);
    return (NULL);
  }

  trk = (DTK_TRACK *)calloc(sizeof(DTK_TRACK), 1);
  trk->npoints = npoints;
  trk->n_scalars = n_scalars;
  trk->n_properties = n_properties;

  trk->c = (float *)calloc(sizeof(float), trk->npoints);
  trk->r = (float *)calloc(sizeof(float), trk->npoints);
  trk->s = (float *)calloc(sizeof(float), trk->npoints);

  for (ns = 0; ns < n_scalars; ns++) trk->scalars[ns] = (float *)calloc(sizeof(float), trk->npoints);

  trk->properties = (float *)calloc(sizeof(float), trk->n_properties);
  return (trk);
}
/*-------------------------------------------------------------*/
MRI *DTKmapEndPoints(DTK_TRACK_SET *dtkset)
{
  MRI *mri, *t;
  int ntrk, c, r, s, nlast, n;
  DTK_TRACK *trk;

  t = dtkset->mri_template;
  mri = MRIalloc(t->width, t->height, t->depth, MRI_INT);
  MRIcopyHeader(t, mri);

  for (ntrk = 0; ntrk < dtkset->hdr->n_count; ntrk++) {
    trk = dtkset->trk[ntrk];
    // First point
    c = nint(trk->c[0]);
    r = nint(trk->r[0]);
    s = nint(trk->s[0]);
    if (c < 0 || c >= t->width) continue;
    if (r < 0 || r >= t->height) continue;
    if (s < 0 || s >= t->depth) continue;
    n = MRIgetVoxVal(mri, c, r, s, 0);
    MRIsetVoxVal(mri, c, r, s, 0, n + 1);
    // Last point
    nlast = trk->npoints - 1;
    c = nint(trk->c[nlast]);
    r = nint(trk->r[nlast]);
    s = nint(trk->s[nlast]);
    if (c < 0 || c >= t->width) continue;
    if (r < 0 || r >= t->height) continue;
    if (s < 0 || s >= t->depth) continue;
    n = MRIgetVoxVal(mri, c, r, s, 0);
    MRIsetVoxVal(mri, c, r, s, 0, n + 1);
  }
  return (mri);
}

/*-------------------------------------------------------------*/
int DTKlabelTracks(DTK_TRACK_SET *trkset, MRI *seg)
{
  int nthtrk, c, r, s, n;
  DTK_TRACK *trk;
  int ndiv = 2000, nthdiv;
  double dc, dr, ds;

  // printf("Labeling tracks\n");
  for (nthtrk = 0; nthtrk < trkset->hdr->n_count; nthtrk++) {
    trk = trkset->trk[nthtrk];
    trk->label = (int *)calloc(sizeof(int), trk->npoints);
    for (n = 0; n < trk->npoints - 1; n++) {
      dc = (trk->c[n + 1] - trk->c[n]) / ndiv;
      dr = (trk->r[n + 1] - trk->r[n]) / ndiv;
      ds = (trk->s[n + 1] - trk->s[n]) / ndiv;
      for (nthdiv = 0; nthdiv < ndiv; nthdiv++) {
        c = (int)round(trk->c[n] + dc * nthdiv);
        r = (int)round(trk->r[n] + dr * nthdiv);
        s = (int)round(trk->s[n] + ds * nthdiv);
        if (c < 0 || c >= seg->width) continue;
        if (r < 0 || r >= seg->height) continue;
        if (s < 0 || s >= seg->depth) continue;
        if (trk->label[n]) continue;
        trk->label[n] = MRIgetVoxVal(seg, c, r, s, 0);
        // printf("%d %d  %d %d %d  %4.1lf %4.1lf %4.1lf  %d\n",nthtrk,n,c,r,s,
        //     trk->c[n],trk->r[n],trk->s[n],trk->label[n]);
      }
    }
    // printf("Done labeling tracks\n");
  }
  return (0);
}

/*-------------------------------------------------------------*/
int DTKlabelTracks0(DTK_TRACK_SET *trkset, MRI *seg)
{
  int nthtrk, c, r, s, n;
  DTK_TRACK *trk;
  // printf("Labeling tracks\n");
  for (nthtrk = 0; nthtrk < trkset->hdr->n_count; nthtrk++) {
    trk = trkset->trk[nthtrk];
    trk->label = (int *)calloc(sizeof(int), trk->npoints);
    for (n = 0; n < trk->npoints; n++) {
      c = (int)round(trk->c[n]);
      r = (int)round(trk->r[n]);
      s = (int)round(trk->s[n]);
      if (c < 0 || c >= seg->width) continue;
      if (r < 0 || r >= seg->height) continue;
      if (s < 0 || s >= seg->depth) continue;
      trk->label[n] = MRIgetVoxVal(seg, c, r, s, 0);
      // printf("%d %d  %d %d %d  %4.1lf %4.1lf %4.1lf  %d\n",nthtrk,n,c,r,s,
      //     trk->c[n],trk->r[n],trk->s[n],trk->label[n]);
    }
    // printf("Done labeling tracks\n");
  }
  return (0);
}

/*--------------------------------------------------------*/
// must first DTKlabelTracks();
DTK_TRACK_SET *DTKextractCC(DTK_TRACK_SET *trkset)
{
  int n, nthtrk, id, isCC;
  DTK_TRACK *trk;
  DTK_TRACK_SET *newset;

  newset = (DTK_TRACK_SET *)calloc(sizeof(DTK_TRACK_SET), 1);
  newset->mri_template = trkset->mri_template;
  newset->hdr = (DTK_HDR *)calloc(sizeof(DTK_HDR), 1);
  memcpy(newset->hdr, trkset->hdr, sizeof(DTK_HDR));
  newset->hdr->n_count = 0;

  for (nthtrk = 0; nthtrk < trkset->hdr->n_count; nthtrk++) {
    trk = trkset->trk[nthtrk];
    isCC = 0;
    for (n = 0; n < trk->npoints; n++) {
      id = trk->label[n];
      if (id == 251 || id == 252 || id == 253 || id == 254 || id == 255) {
        isCC = 1;
        break;
      }
    }
    if (isCC) {
      newset->trk = (DTK_TRACK **)realloc(newset->trk, sizeof(DTK_TRACK *) * (newset->hdr->n_count + 1));
      newset->trk[newset->hdr->n_count] = DTKcopyTrack(trk);
      newset->hdr->n_count++;
    }
  }
  return (newset);
}

/*---------------------------------------------------------------*/
DTK_TRACK_SET *DTKextractSeg(DTK_TRACK_SET *trkset, int segid)
{
  int n, nthtrk, id, isSeg;
  DTK_TRACK *trk;
  DTK_TRACK_SET *newset;

  newset = (DTK_TRACK_SET *)calloc(sizeof(DTK_TRACK_SET), 1);
  newset->mri_template = trkset->mri_template;
  newset->hdr = (DTK_HDR *)calloc(sizeof(DTK_HDR), 1);
  memcpy(newset->hdr, trkset->hdr, sizeof(DTK_HDR));
  newset->hdr->n_count = 0;

  for (nthtrk = 0; nthtrk < trkset->hdr->n_count; nthtrk++) {
    trk = trkset->trk[nthtrk];
    isSeg = 0;
    for (n = 0; n < trk->npoints; n++) {
      id = trk->label[n];
      if (id == segid) {
        isSeg = 1;
        break;
      }
    }
    if (isSeg) {
      // printf("nthtrk=%d isSeg=%d %d %d\n",nthtrk,isSeg,id,segid);
      newset->trk = (DTK_TRACK **)realloc(newset->trk, sizeof(DTK_TRACK *) * (newset->hdr->n_count + 1));
      newset->trk[newset->hdr->n_count] = DTKcopyTrack(trk);
      newset->hdr->n_count++;
    }
  }
  printf("Found %d \n", newset->hdr->n_count);
  return (newset);
}

/*---------------------------------------------------------------*/
DTK_TRACK_SET *DTKextractSegEndPoints(DTK_TRACK_SET *trkset, int segid)
{
  int nthtrk;
  DTK_TRACK *trk;
  DTK_TRACK_SET *newset;

  newset = (DTK_TRACK_SET *)calloc(sizeof(DTK_TRACK_SET), 1);
  newset->mri_template = trkset->mri_template;
  newset->hdr = (DTK_HDR *)calloc(sizeof(DTK_HDR), 1);
  memcpy(newset->hdr, trkset->hdr, sizeof(DTK_HDR));
  newset->hdr->n_count = 0;

  for (nthtrk = 0; nthtrk < trkset->hdr->n_count; nthtrk++) {
    trk = trkset->trk[nthtrk];

    if (trk->label[0] == segid || trk->label[trk->npoints - 1] == segid) {
      // printf("nthtrk=%d isSeg=%d %d %d\n",nthtrk,isSeg,id,segid);
      newset->trk = (DTK_TRACK **)realloc(newset->trk, sizeof(DTK_TRACK *) * (newset->hdr->n_count + 1));
      newset->trk[newset->hdr->n_count] = DTKcopyTrack(trk);
      newset->hdr->n_count++;
    }
  }
  printf("Found %d \n", newset->hdr->n_count);
  return (newset);
}

/*-------------------------------------------------------------*/
int DTKprintHeader(FILE *fp, DTK_HDR *trkh)
{
  fprintf(fp, "id_string %s\n", trkh->id_string);
  fprintf(fp, "dim %d %d %d\n", trkh->dim[0], trkh->dim[1], trkh->dim[2]);
  fprintf(fp, "voxel_size  %f %f %f\n", trkh->voxel_size[0], trkh->voxel_size[1], trkh->voxel_size[2]);
  fprintf(fp, "origin %f %f %f\n", trkh->origin[0], trkh->origin[1], trkh->origin[2]);
  fprintf(fp, "n_scalars %d\n", trkh->n_scalars);
  fprintf(fp, "n_properties %d\n", trkh->n_properties);
  fprintf(fp, "voxel_order %s\n", trkh->voxel_order);
  fprintf(fp, "n_count %d\n", trkh->n_count);
  fprintf(fp, "version %d\n", trkh->version);
  fprintf(fp, "hdr_size %d\n", trkh->hdr_size);
  return (0);
}
LABEL *DTKtrack2Label(DTK_TRACK *trk)
{
  LABEL *trklabel;
  int n;
  MRI *mri;
  MATRIX *M, *crs1, *xyz1;

  mri = MRIread("dsi.nii");
  M = MRIxfmCRS2XYZtkreg(mri);
  MRIfree(&mri);

  crs1 = MatrixAlloc(4, 1, MATRIX_REAL);
  crs1->rptr[4][1] = 1;
  xyz1 = MatrixAlloc(4, 1, MATRIX_REAL);

  trklabel = LabelAlloc(trk->npoints, NULL, NULL);
  trklabel->n_points = trk->npoints;
  for (n = 0; n < trk->npoints; n++) {
    crs1->rptr[1][1] = trk->c[n];
    crs1->rptr[2][1] = trk->r[n];
    crs1->rptr[3][1] = trk->s[n];
    xyz1 = MatrixMultiply(M, crs1, xyz1);
    trklabel->lv[n].x = xyz1->rptr[1][1];
    trklabel->lv[n].y = xyz1->rptr[2][1];
    trklabel->lv[n].z = xyz1->rptr[3][1];
  }

  MatrixFree(&M);
  MatrixFree(&crs1);
  MatrixFree(&xyz1);

  return (trklabel);
}

/*--------------------------------------------------------------
  MRI *DTKmapTrackNos() - first tp is the number of tracks that pass
  through the voxel, the remainder are the track numbers of the
  passing tracks.
  --------------------------------------------------------------*/
MRI *DTKmapTrackNos(DTK_TRACK_SET *trkset)
{
  MRI *mri, *t;
  int ntrk, c, r, s, n, nhits, nmaxhits = 10000;
  DTK_TRACK *trk;

  t = trkset->mri_template;
  mri = MRIallocSequence(t->width, t->height, t->depth, MRI_INT, nmaxhits);
  if (mri == NULL) return (NULL);
  MRIcopyHeader(t, mri);

  printf("Mapping Track Numbers\n");
  // Go thru each track
  for (ntrk = 0; ntrk < trkset->hdr->n_count; ntrk++) {
    trk = trkset->trk[ntrk];
    // Go thru each point in the track
    for (n = 0; n < trk->npoints; n++) {
      // CRS of voxel for this point
      c = (int)round(trk->c[n]);
      r = (int)round(trk->r[n]);
      s = (int)round(trk->s[n]);
      if (c < 0 || c >= t->width) continue;
      if (r < 0 || r >= t->height) continue;
      if (s < 0 || s >= t->depth) continue;
      // Current number of tracks passing through this voxel
      nhits = MRIgetVoxVal(mri, c, r, s, 0);
      // Track number for this track
      MRIsetVoxVal(mri, c, r, s, nhits, ntrk);
      // Update number of tracks passing through this voxel
      nhits++;
      if (nhits >= nmaxhits) {
        printf("ERROR: nhits exceeds %d   %d %d %d\n", nmaxhits, c, r, s);
        return (NULL);
      }
      MRIsetVoxVal(mri, c, r, s, 0, nhits);
    }
  }
  printf("Done mapping Track Numbers\n");
  return (mri);
}

/*----------------------------------------------------------------*/
MRI *DTKsegROI(DTK_TRACK_SET *trkset, MRI *seg, int segid)
{
  MRI *mri, *tracknos;
  int c, r, s, id, nthid, idlist[1000], ntracks, id0, n0, trackno, n;
  int term, cT = 0, rT = 0, sT = 0, idT = 0, pct = 0;
  DTK_TRACK *trk;

  tracknos = DTKmapTrackNos(trkset);
  if (tracknos == NULL) return (NULL);

  mri = MRIallocSequence(seg->width, seg->height, seg->depth, MRI_INT, 2);
  MRIcopyHeader(seg, mri);

  for (c = 0; c < seg->width; c++) {
    for (r = 0; r < seg->height; r++) {
      for (s = 0; s < seg->depth; s++) {
        MRIsetVoxVal(mri, c, r, s, 0, 0);  // make sure it starts at 0
        id = MRIgetVoxVal(seg, c, r, s, 0);
        // Is this CRS a voxel in the seed region?
        if (id != segid) {
          if (id == 2 || id == 41) continue;
          if (id == 11 || id == 12 || id == 13) continue;
          if (id == 50 || id == 51 || id == 52) continue;
          if (id == 4 || id == 5 || id == 13 || id == 28) continue;
          if (id == 43 || id == 44 || id == 15 || id == 60) continue;
          if (id == 31 || id == 63) continue;
          if (id == 26 || id == 58) continue;
          // Just set id to target seg id
          MRIsetVoxVal(mri, c, r, s, 0, id);
          continue;
        }
        // Get number of tracks through this CRS
        ntracks = MRIgetVoxVal(tracknos, c, r, s, 0);
        if (ntracks == 0) continue;
        // Go through each track and find where the end
        // points land. Get list of all end point ids
        nthid = 0;
        for (n = 0; n < ntracks; n++) {
          trackno = MRIgetVoxVal(tracknos, c, r, s, n + 1);
          trk = trkset->trk[trackno];
          for (term = 0; term <= 1; term++) {
            if (term == 0) {
              // "Start" of the track
              cT = (int)round(trk->c[0]);
              rT = (int)round(trk->r[0]);
              sT = (int)round(trk->s[0]);
            }
            if (term == 1) {
              // "End" of the track
              cT = (int)round(trk->c[trk->npoints - 1]);
              rT = (int)round(trk->r[trk->npoints - 1]);
              sT = (int)round(trk->s[trk->npoints - 1]);
            }
            if (cT < 0 || cT >= seg->width) continue;
            if (rT < 0 || rT >= seg->height) continue;
            if (sT < 0 || sT >= seg->depth) continue;
            // Seg Id at this termination
            idT = MRIgetVoxVal(seg, cT, rT, sT, 0);
            if (idT == segid) continue;  // Dont map to self
            if ((idT >= 1000 && idT <= 1035) || (idT >= 2000 && idT <= 2035)) {
              idlist[nthid] = idT;
              nthid++;
            }
          }
        }
        if (nthid != 0) {
          // Get the most frequently occuring label
          id0 = most_frequent_int_list(idlist, nthid, &n0);
          pct = nint(1000.0 * n0 / nthid);  // 10x percent
        }
        else
          id0 = 0;  // nothing appropriate

        MRIsetVoxVal(mri, c, r, s, 0, id0);
        MRIsetVoxVal(mri, c, r, s, 1, pct);

      }  // C
    }    // R
  }      // S

  MRIfree(&tracknos);
  return (mri);
}
