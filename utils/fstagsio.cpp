#include <string.h>
#include <stdlib.h>

#include "fstagsio.h"
#include "tags.h"
#include "matrix.h"
#include "colortab.h"
#include "mri.h"
#include "log.h"

/* The following TAGs are in length-less format for mgz output:
 *   TAG_OLD_COLORTABLE
 *   TAG_GCAMORPH_GEOM
 *   TAG_GCAMORPH_TYPE
 *   TAG_GCAMORPH_LABELS
 *
 * If niftiheaderext = true, all TAGs are in this format:
 *   tagid + sizeof(long long) + len(tagdata)
 */

// constructor
FStagsIO::FStagsIO(znzFile infp, bool niftiheaderext0)
{
  fp = infp;
  niftiheaderext = niftiheaderext0;
}


// destructor
FStagsIO::~FStagsIO()
{
  fp = NULL;
}


// the following getlen_*() methods return TAG length as following:
//
// if addtaglength == true,
//   if niftiheaderext == false,
//       TAG w/o  data-length: tagid + len(tagdata)
//       TAG w/ a data-length: tagid + sizeof(long long) + len(tagdata)
//    else (niftiheaderext == true)
//       TAG w/ a data-length: tagid + sizeof(long long) + len(tagdata)
// otherwise,
//   return len(tagdata)
//
// notes:
//   1. when getlen_*() are called from write_*(), we are getting len(tagdata) only,
//        set 'addtaglength = false'
//      * len(tagdata) returned could be different depending on niftiheaderext
//      * getlen_matrix() and getlen_mri_frames() return different len(tagdata) for
//            .mgz/.mgh and fs nifti header extension
//   2. when getlen_*() are called to calculate the total length of nifti header extension,
//      we would like to have a data-length field for all TAGs,
//        set 'addtaglength = true' (default), 'niftiheaderext = true'
//        
long long FStagsIO::getlen_tag(int tag, long long len, bool niftiheaderext, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
  
    if (niftiheaderext ||
        (tag != TAG_OLD_COLORTABLE && tag != TAG_GCAMORPH_GEOM &&
	 tag != TAG_GCAMORPH_TYPE && tag != TAG_GCAMORPH_LABELS))
      dlen += sizeof(long long);
  }

  dlen += len;

  return dlen;
}


// data length is different depends on niftiheaderext
long long FStagsIO::getlen_matrix(bool niftiheaderext, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }

  if (niftiheaderext)
    dlen += 16 * sizeof(float);
  else
    dlen += MATRIX_STRLEN;  // MATRIX_STRLEN = (4 * 4 * 100)
  
  return dlen;
}


long long FStagsIO::getlen_old_colortable(COLOR_TABLE *ctab, bool niftiheaderext, bool addtaglength)
{
  long long dlen = 0;

  if (addtaglength)
  {
    dlen += 4;

    if (niftiheaderext)
      dlen += sizeof(long long);
  }

  // needs to be consistent with znzCTABwriteIntoBinary()
  int version = ctab->version;
  if (version == 2)
    dlen += sizeof(int); // version (not in v1)
  dlen += sizeof(int); // nentries
  dlen += sizeof(int); // len(fname)
  dlen += strlen(ctab->fname) + 1;
  if (version == 2)
    dlen += sizeof(int); // num_entries (not in v1)
  for (int structure = 0; structure < ctab->nentries; structure++) {
    if (NULL != ctab->entries[structure]) {
      if (version == 2)
	dlen += sizeof(int); // structure (not in v1)
      dlen += sizeof(int); // len(structure-name)+1
      dlen += strlen(ctab->entries[structure]->name) + 1;
      dlen += sizeof(int); // ri
      dlen += sizeof(int); // gi
      dlen += sizeof(int); // bi
      dlen += sizeof(int); // t-ai
    }
  }
  
  return dlen;
}


// data length is different depends on niftiheaderext
// .mgz/.mgh outputs '10 * mri->nframes * sizeof(MRI_FRAME)' bytes to disk
// Freesurfer nifti header extension only outputs fields label, name, and thresh
// see comments in mri.h about struct MRI_FRAME
long long FStagsIO::getlen_mri_frames(MRI *mri, bool niftiheaderext, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }

  if (niftiheaderext)
  {
    for (int fno = 0; fno < mri->nframes; fno++) {
      MRI_FRAME *frame = &mri->frames[fno];
      dlen += sizeof(frame->label);
      dlen += sizeof(frame->name);
      dlen += sizeof(frame->thresh);
    }
  }
  else
    dlen += 10 * mri->nframes * sizeof(MRI_FRAME);

  return dlen;
}


// return different length depends on niftiheaderext
long long FStagsIO::getlen_gcamorph_geom(const char *source_fname, const char *target_fname, bool niftiheaderext, bool addtaglength, bool shearless)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
  
    if (niftiheaderext)
      dlen += sizeof(long long);
  }

  // this needs to be consistent with write_gcamorph_geom()/VOL_GEOM.write()
  if (!niftiheaderext)
  {
    int geom_len = 4 * sizeof(int) + 15 * sizeof(float) + 512;
    dlen += 2 * geom_len;
  }
  else
  {
    int geom_len = 4 * sizeof(int) + 15 * sizeof(float) + sizeof(int);
    geom_len *= 2;
    geom_len += strlen(source_fname) + strlen(target_fname);
    dlen += geom_len;
  }

  if (!shearless)
    dlen += 6 * sizeof(float);  // 3 floats for each geom

  return dlen;
}


// needs to be consistent with write_gcamorph_meta()
long long FStagsIO::getlen_gcamorph_meta(bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }
  dlen += sizeof(int) + sizeof(int) + sizeof(float);
  
  return dlen;
}


// TAG_GCAMORPH_LABELS is in legnth-less format for .mgz;
long long FStagsIO::getlen_gcamorph_labels(int x, int y, int z, int len, bool niftiheaderext, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
  
    if (niftiheaderext)
      dlen += sizeof(long long);
  }
  
  dlen += (x * y * z * len);

  return dlen;
}


long long FStagsIO::getlen_dof(int dof, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }
  dlen += sizeof(dof);
  
  return dlen;  
}


long long FStagsIO::getlen_scan_parameters(MRI *mri, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }

  // this needs to be consistent with write_scan_parameters()
  dlen += sizeof(mri->te);
  dlen += sizeof(mri->ti);
  dlen += sizeof(mri->flip_angle);
  dlen += sizeof(mri->FieldStrength);
  if (mri->pedir)
    dlen += strlen(mri->pedir) + 1;
  else
    dlen += strlen("UNKNOWN");

  return dlen;  
}


long long FStagsIO::getlen_ras_xform(VOL_GEOM *vg, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }

  // this needs to be consistent with write_ras_xform()
  dlen += sizeof(vg->x_r); dlen += sizeof(vg->x_a); dlen += sizeof(vg->x_s);
  dlen += sizeof(vg->y_r); dlen += sizeof(vg->y_a); dlen += sizeof(vg->y_s);
  dlen += sizeof(vg->z_r); dlen += sizeof(vg->z_a); dlen += sizeof(vg->z_s);
  dlen += sizeof(vg->c_r); dlen += sizeof(vg->c_a); dlen += sizeof(vg->c_s);

  return dlen;  
}


long long FStagsIO::getlen_ras_xform_extra(VOL_GEOM *vg, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }

  // this needs to be consistent with write_ras_xform_extra()
  dlen += sizeof(vg->s_r);   dlen += sizeof(vg->s_a);    dlen += sizeof(vg->s_s);
  dlen += sizeof(vg->width); dlen += sizeof(vg->height); dlen += sizeof(vg->depth);
  dlen += sizeof(vg->xsize); dlen += sizeof(vg->ysize);  dlen += sizeof(vg->zsize);  
  
  return dlen;  
}

// this is for nifti1 header extension only
//   TAG_END_NIIHDREXTENSION data-length=1 '*'
// needs to be consistent with FStagsIO::write_endtag()
long long FStagsIO::getlen_endtag(bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }

  dlen += 1; // extra char '*'

  return dlen;  
}


// tags.cpp::znzTAGwrite()
// 
// output TAG in the following format:
//   TAG w/o  length: tagid + tagdata
//   TAG w/ a length: tagid + len(tagdata) + tagdata
int FStagsIO::write_tag(int tag, void *data, long long dlen)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(tag, fp);
  
  // ??? todo: check if tag is length-less
  if (niftiheaderext ||
      (tag != TAG_OLD_COLORTABLE && tag != TAG_GCAMORPH_GEOM &&
       tag != TAG_GCAMORPH_TYPE && tag != TAG_GCAMORPH_LABELS))
    znzwriteLong(dlen, fp);
  
  znzwrite(data, sizeof(char), dlen, fp);

  long long fend = znztell(fp);
    
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    printf("[DEBUG] FStagsIO::write_tag() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld) (%-6lld)\n", tag, fend-fstart, fstart, fend, dlen);
  
  return NO_ERROR;
}


// tags.cpp::znzWriteMatrix()
int FStagsIO::write_matrix(MATRIX *M, int tag)
{
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    printf("[DEBUG] FStagsIO::write_matrix()\n");
    MatrixPrint(stdout, M);
  }
  
  if (niftiheaderext)
    return __write_matrix_niftiheaderext(M, tag);
  
  long long dlen = MATRIX_STRLEN;
  char matbuf[dlen];

  bzero(matbuf, dlen);
  sprintf(matbuf,
          "Matrix %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf",
          M->rptr[1][1],
          M->rptr[1][2],
          M->rptr[1][3],
          M->rptr[1][4],
          M->rptr[2][1],
          M->rptr[2][2],
          M->rptr[2][3],
          M->rptr[2][4],
          M->rptr[3][1],
          M->rptr[3][2],
          M->rptr[3][3],
          M->rptr[3][4],
          M->rptr[4][1],
          M->rptr[4][2],
          M->rptr[4][3],
          M->rptr[4][4]);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    printf("[DEBUG] FStagsIO::write_matrix() TAG = %-4d, len = %-6lld\n", tag, dlen);
    printf("[DEBUG] FStagsIO::write_matrix() %s\n", matbuf);
  }
    
  return write_tag(tag, matbuf, dlen);
}


// write binary colortable
int FStagsIO::write_old_colortable(COLOR_TABLE *ctab)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_OLD_COLORTABLE, fp);
  if (niftiheaderext)
  {
    long long dlen = getlen_old_colortable(ctab, niftiheaderext, false);
    znzwriteLong(dlen, fp);
  }
  
  znzCTABwriteIntoBinary(ctab, fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_old_colortable() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_OLD_COLORTABLE, fend-fstart, fstart, fend);
  }
    
  return NO_ERROR;
}


// mriio.cpp::znzTAGwriteMRIframes()
int FStagsIO::write_mri_frames(MRI *mri)
{
  if (niftiheaderext)
    return __write_mri_frames_niftiheaderext(mri);

  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  // write some extra space so that we have enough room (can't seek in zz files)
  long long len = getlen_mri_frames(mri, niftiheaderext, false);  //10 * mri->nframes * sizeof(MRI_FRAME);

  znzwriteInt(TAG_MRI_FRAME, fp);
  znzwriteLong(len, fp);
  
  long long here = znztell(fp);
  for (int fno = 0; fno < mri->nframes; fno++) {
    MRI_FRAME *frame = &mri->frames[fno];
    znzwriteInt(frame->type, fp);
    znzwriteFloat(frame->TE, fp);
    znzwriteFloat(frame->TR, fp);
    znzwriteFloat(frame->flip, fp);
    znzwriteFloat(frame->TI, fp);
    znzwriteFloat(frame->TD, fp);
    znzwriteFloat(frame->TM, fp);
    znzwriteInt(frame->sequence_type, fp);
    znzwriteFloat(frame->echo_spacing, fp);
    znzwriteFloat(frame->echo_train_len, fp);
    for (int i = 0; i < 3; i++) znzwriteFloat(frame->read_dir[i], fp);
    for (int i = 0; i < 3; i++) znzwriteFloat(frame->pe_dir[i], fp);
    for (int i = 0; i < 3; i++) znzwriteFloat(frame->slice_dir[i], fp);
    znzwriteInt(frame->label, fp);
    znzwrite(frame->name, sizeof(char), STRLEN, fp);
    znzwriteInt(frame->dof, fp);

    if (frame->m_ras2vox && frame->m_ras2vox->rows > 0) {
      // znzwriteMatrix(fp, frame->m_ras2vox, 0);
      write_matrix(frame->m_ras2vox, 0);
    }
    else {
      MATRIX *m = MatrixAlloc(4, 4, MATRIX_REAL);
      // znzwriteMatrix(fp, m, 0);
      write_matrix(m, 0);      
      MatrixFree(&m);
    }

    znzwriteFloat(frame->thresh, fp);
    znzwriteInt(frame->units, fp);
    if (frame->type == FRAME_TYPE_DIFFUSION_AUGMENTED)  // also store diffusion info
    {
      znzwriteDouble(frame->DX, fp);
      znzwriteDouble(frame->DY, fp);
      znzwriteDouble(frame->DZ, fp);

      znzwriteDouble(frame->DR, fp);
      znzwriteDouble(frame->DP, fp);
      znzwriteDouble(frame->DS, fp);
      znzwriteDouble(frame->bvalue, fp);
      znzwriteDouble(frame->TM, fp);

      znzwriteLong(frame->D1_ramp, fp);
      znzwriteLong(frame->D1_flat, fp);
      znzwriteDouble(frame->D1_amp, fp);

      znzwriteLong(frame->D2_ramp, fp);
      znzwriteLong(frame->D2_flat, fp);
      znzwriteDouble(frame->D2_amp, fp);

      znzwriteLong(frame->D3_ramp, fp);
      znzwriteLong(frame->D3_flat, fp);
      znzwriteDouble(frame->D3_amp, fp);

      znzwriteLong(frame->D4_ramp, fp);
      znzwriteLong(frame->D4_flat, fp);
      znzwriteDouble(frame->D4_amp, fp);
    }
  }

  long long fend = znztell(fp);
  len -= (fend - here);  // unused space
  if (len > 0) {
    char *buf = (char *)calloc(len, sizeof(char));
    znzwrite(buf, len, sizeof(char), fp);
    free(buf);
  }
  //znzTAGwriteEnd(fp, fend);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_mri_frames() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_MRI_FRAME, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


// write TAG_GCAMORPH_GEOM/TAG_GCAMORPH_GEOM_PLUSSHEAR
// TAG_GCAMORPH_GEOM is in length-less format if niftiheaderext = false
// TAG_GCAMORPH_GEOM_PLUSSHEAR has a length (shearless = false || niftiheader = true)
int FStagsIO::write_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target, bool shearless)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  int tag = TAG_GCAMORPH_GEOM;
  if (!shearless)
    tag = TAG_GCAMORPH_GEOM_PLUSSHEAR;
  znzwriteInt(tag, fp);
  
  if (niftiheaderext || !shearless)
  {
    long long dlen = getlen_gcamorph_geom(source->fname, target->fname, niftiheaderext, false, shearless);
    znzwriteLong(dlen, fp);
  }

  VOL_GEOM src_geom = *source;
  VOL_GEOM trg_geom = *target;
  if (shearless)
  {
    src_geom.shearless_components();
    trg_geom.shearless_components();
  }
  src_geom.write(fp, niftiheaderext, shearless);
  trg_geom.write(fp, niftiheaderext, shearless);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    source->vgprint();
    target->vgprint();
    
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_gcamorph_geom() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", tag, fend-fstart, fstart, fend);

    src_geom.vgprint();
    trg_geom.vgprint();
  }
  
  return NO_ERROR;
}


// TAG_GCAMORPH_META
int FStagsIO::write_gcamorph_meta(int warpFieldFormat, int gcamorphSpacing, double gcamorphExp_k)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_GCAMORPH_META, fp);

  long long dlen = getlen_gcamorph_meta(false);
  znzwriteLong(dlen, fp);
  znzwriteInt(warpFieldFormat, fp);
  znzwriteInt(gcamorphSpacing, fp);
  znzwriteFloat(gcamorphExp_k, fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_gcamorph_meta() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_GCAMORPH_META, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


// TAG_GCAMORPH_LABELS is in length-less format
int FStagsIO::write_gcamorph_labels(int x0, int y0, int z0, int ***gcamorphLabel)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_GCAMORPH_LABELS, fp);
  
  if (niftiheaderext)
  {
    long long dlen = getlen_gcamorph_labels(x0, y0, z0, sizeof(int), niftiheaderext, false);
    znzwriteLong(dlen, fp);
  }
  
  for (int x = 0; x < x0; x++)
    for (int y = 0; y < y0; y++)
      for (int z = 0; z < z0; z++)
        znzwriteInt(gcamorphLabel[x][y][z], fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_gcamorph_labels() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_GCAMORPH_LABELS, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


// write TAG_DOF (nifti header extension only)
int FStagsIO::write_dof(int dof)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_DOF, fp);

  long long dlen = getlen_dof(dof, false);
  znzwriteLong(dlen, fp);
  znzwriteInt(dof, fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_dof() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_DOF, fend-fstart, fstart, fend);
  }

  return NO_ERROR;  
}


// write TAG_SCAN_PARAMETERS (nifti header extension only)
int FStagsIO::write_scan_parameters(MRI *mri)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_SCAN_PARAMETERS, fp);

  long long dlen = getlen_scan_parameters(mri, false);
  znzwriteLong(dlen, fp);
  znzwriteFloat(mri->te, fp);
  znzwriteFloat(mri->ti, fp);
  znzwriteDouble(mri->flip_angle, fp);
  // ??? todo: check how is fov calculated ???
  // skip fov, it can be calculated from other parameters
  // znzwriteFloat(mri->fov, fp);
  znzwriteFloat(mri->FieldStrength, fp);
  if (mri->pedir)
    znzwrite(mri->pedir, sizeof(char), strlen(mri->pedir) + 1, fp);
  else
    znzwrite((void *)"UNKNOWN", sizeof(char), strlen("UNKNOWN"), fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_scan_parameters() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_SCAN_PARAMETERS, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;  
}

  
// write TAG_RAS_XFORM (nifti header extension only)
int FStagsIO::write_ras_xform(VOL_GEOM *vg)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);

  znzwriteInt(TAG_RAS_XFORM, fp);

  long long dlen = getlen_ras_xform(vg, false);
  znzwriteLong(dlen, fp);
  znzwriteFloat(vg->x_r, fp); znzwriteFloat(vg->x_a, fp); znzwriteFloat(vg->x_s, fp);
  znzwriteFloat(vg->y_r, fp); znzwriteFloat(vg->y_a, fp); znzwriteFloat(vg->y_s, fp);
  znzwriteFloat(vg->z_r, fp); znzwriteFloat(vg->z_a, fp); znzwriteFloat(vg->z_s, fp);
  znzwriteFloat(vg->c_r, fp); znzwriteFloat(vg->c_a, fp); znzwriteFloat(vg->c_s, fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    vg->geomprint("[DEBUG] FStagsIO::write_ras_xform() ras xform info:\n");

    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_ras_xform() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_RAS_XFORM, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;  
}


// write TAG_RAS_XFORM_EXTRA (nifti header extension only)
int FStagsIO::write_ras_xform_extra(VOL_GEOM *vg)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);

  znzwriteInt(TAG_RAS_XFORM_EXTRA, fp);

  long long dlen = getlen_ras_xform_extra(vg, false);
  znzwriteLong(dlen, fp);
  znzwriteFloat(vg->s_r, fp);   znzwriteFloat(vg->s_a, fp);   znzwriteFloat(vg->s_s, fp);
  znzwriteInt(vg->width, fp);   znzwriteInt(vg->height, fp);  znzwriteInt(vg->depth, fp);
  znzwriteFloat(vg->xsize, fp); znzwriteFloat(vg->ysize, fp); znzwriteFloat(vg->zsize, fp);  

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    vg->geomprint("[DEBUG] FStagsIO::write_ras_xform_extra() ras xform info:\n");

    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_ras_xform_extra() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_RAS_XFORM_EXTRA, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;  
}

/* write TAG_END_NIIHDREXTENSION (nifti header extension only)
 * this needs to be the last tag.
 *
 * write TAG_END_NIIHDREXTENSION at the end of extension data to avoid the data to be truncated:
 *   TAG_END_NIIHDREXTENSION (-1)  data-length (1) '*'
 *
 * If the extension data has trailing null characters or zeros at the end,
 * nibabel.nifti1.Nifti1Extension.get_content() will truncate the data.
 * See https://github.com/nipy/nibabel/blob/master/nibabel/nifti1.py#L629C1-L630C1,
 * line 629:  'evalue = evalue.rstrip(b'\x00')'
 */
int FStagsIO::write_endtag()
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_END_NIIHDREXTENSION, fp);

  long long dlen = getlen_endtag(false);
  znzwriteLong(dlen, fp);

  char endchar = '*';
  znzwrite(&endchar, sizeof(char), dlen, fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::write_endtag() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_END_NIIHDREXTENSION, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;  
}


// tags.cpp::znzTAGreadStart()
int FStagsIO::read_tagid_len(long long *plen, int tagwithzerolen)
{
  int tag;

  if (znzeof(fp)) return (0);
  tag = znzreadInt(fp);
  if (znzeof(fp)) return (0);

  // for nifti header extension, there is a data-length for all TAGs
  if (niftiheaderext)
  {
    *plen = znzreadLong(fp);  // read data-length for the tagid
    return tag;
  }
  
  if (tagwithzerolen && tagwithzerolen == tag)
  {
    /* This is to handle following situation: 
     * TAG_MGH_XFORM is used in both mgz and m3z, but in different format.
     * in mgz, data-length is output after TAG_MGH_XFORM
     * in m3z, no data-length is output after TAG_MGH_XFORM
     *
     * in __m3zRead(), the function is called as znzTAGreadStart(file, &len, TAG_MGH_XFORM);
     */
    *plen = 0;
    return tag;
  }
  
  switch (tag) {
    case TAG_OLD_MGH_XFORM:  // the tag doesn't look like being used
      *plen = (long long)znzreadInt(fp); /* sorry - backwards compatibility
                                            with Tosa's stuff */
      *plen = *plen - 1;                 // doesn't include null
      break;
    case TAG_OLD_SURF_GEOM:  // these don't take lengths at all
    case TAG_OLD_USEREALRAS:
    case TAG_OLD_COLORTABLE:
    case TAG_GCAMORPH_GEOM:
    case TAG_GCAMORPH_TYPE:
    case TAG_GCAMORPH_LABELS:
      *plen = 0;  // these tags have no data-length output after tagid
      break;
    default:
      *plen = znzreadLong(fp);  // read data-length for the tagid
  }

  return (tag);  
}


// read len of bytes at current file position into databuf
int FStagsIO::read_data(void *databuf, long long len)
{
  znzread(databuf, sizeof(char), len, fp);
  return  NO_ERROR;
}


// tags.cpp::znzReadAutoAlignMatrix()
MATRIX* FStagsIO::read_matrix()
{
  MATRIX *M = NULL;
  if (niftiheaderext)
    M = __read_matrix_niftiheaderext();
  else
  {
    char buf[MATRIX_STRLEN];

    /* no fscanf equivalent in zlib!! have to hack it */
    znzread(buf, sizeof(unsigned char), MATRIX_STRLEN, fp);
    M = MatrixAlloc(4, 4, MATRIX_REAL);
    char ch[100];
    sscanf(buf,
           "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
             ch,
           &(M->rptr[1][1]),
           &(M->rptr[1][2]),
           &(M->rptr[1][3]),
           &(M->rptr[1][4]),
           &(M->rptr[2][1]),
           &(M->rptr[2][2]),
           &(M->rptr[2][3]),
           &(M->rptr[2][4]),
           &(M->rptr[3][1]),
           &(M->rptr[3][2]),
           &(M->rptr[3][3]),
           &(M->rptr[3][4]),
           &(M->rptr[4][1]),
           &(M->rptr[4][2]),
           &(M->rptr[4][3]),
           &(M->rptr[4][4]));
  }

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    printf("[DEBUG] FStagsIO::read_matrix()\n");
    MatrixPrint(stdout, M);
  }
    
  return (M);
}


// read binary colortable
COLOR_TABLE* FStagsIO::read_old_colortable()
{
  return znzCTABreadFromBinary(fp);
}


// mriio.cpp::znzTAGreadMRIframes()
int FStagsIO::read_mri_frames(MRI *mri, long long len)
{
  if (niftiheaderext)
    return __read_mri_frames_niftiheaderext(mri, len);
  
  long long fstart = znztell(fp);
  for (int fno = 0; fno < mri->nframes; fno++) {
    MRI_FRAME *frame = &mri->frames[fno];
    frame->type = znzreadInt(fp);
    frame->TE = znzreadFloat(fp);
    frame->TR = znzreadFloat(fp);
    frame->flip = znzreadFloat(fp);
    frame->TI = znzreadFloat(fp);
    frame->TD = znzreadFloat(fp);
    frame->TM = znzreadFloat(fp);
    frame->sequence_type = znzreadInt(fp);
    frame->echo_spacing = znzreadFloat(fp);
    frame->echo_train_len = znzreadFloat(fp);
    for (int i = 0; i < 3; i++) frame->read_dir[i] = znzreadFloat(fp);
    for (int i = 0; i < 3; i++) frame->pe_dir[i] = znzreadFloat(fp);
    for (int i = 0; i < 3; i++) frame->slice_dir[i] = znzreadFloat(fp);
    frame->label = znzreadInt(fp);
    znzread(frame->name, sizeof(char), STRLEN, fp);
    frame->dof = znzreadInt(fp);

    // the embedded matrix has tag and data-length 
    long long matlen = 0;
    int mattag = read_tagid_len(&matlen);
    frame->m_ras2vox = read_matrix();  // znzReadMatrix(fp);
    if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    {
      printf("[DEBUG] FStagsIO::read_mri_frame() TAG = %-4d, len = %-6lld\n", mattag, matlen);
      MatrixPrint(stdout, frame->m_ras2vox);
    }    

    frame->thresh = znzreadFloat(fp);
    frame->units = znzreadInt(fp);
    if (frame->type == FRAME_TYPE_DIFFUSION_AUGMENTED) {
      frame->DX = znzreadDouble(fp);
      frame->DY = znzreadDouble(fp);
      frame->DZ = znzreadDouble(fp);

      frame->DR = znzreadDouble(fp);
      frame->DP = znzreadDouble(fp);
      frame->DS = znzreadDouble(fp);
      frame->bvalue = znzreadDouble(fp);
      frame->TM = znzreadDouble(fp);

      frame->D1_ramp = znzreadLong(fp);
      frame->D1_flat = znzreadLong(fp);
      frame->D1_amp = znzreadDouble(fp);

      frame->D2_ramp = znzreadLong(fp);
      frame->D2_flat = znzreadLong(fp);
      frame->D2_amp = znzreadDouble(fp);

      frame->D3_ramp = znzreadLong(fp);
      frame->D3_flat = znzreadLong(fp);
      frame->D3_amp = znzreadDouble(fp);

      frame->D4_ramp = znzreadLong(fp);
      frame->D4_flat = znzreadLong(fp);
      frame->D4_amp = znzreadDouble(fp);
    }
  }

  long long fend = znztell(fp);
  len -= (fend - fstart);
  if (len > 0) {
    // write_mri_frames() outputs more than it is needed to disk
    // skip any extra bytes
    if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
      printf("[DEBUG] read_mri_frames() TAG = %-4d, bytes_read = %-6lld (%-6lld - %-6lld), skip extra bytes %lld\n", TAG_MRI_FRAME, fend-fstart, fstart, fend, len);
	
    char *buf = (char *)calloc(len, sizeof(char));
    znzread(buf, len, sizeof(char), fp);
    free(buf);
  }

  return NO_ERROR;  
}


// read TAG_GCAMORPH_GEOM/TAG_GCAMORPH_GEOM_PLUSSHEAR data
int FStagsIO::read_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target, bool shearless)
{
  source->read(fp, niftiheaderext, shearless);
  target->read(fp, niftiheaderext, shearless);

  return NO_ERROR;
}


// read TAG_GCAMORPH_META data
int FStagsIO::read_gcamorph_meta(int *warpFieldFormat, int *gcamorphSpacing, double *gcamorphExp_k)
{
  *warpFieldFormat = znzreadInt(fp);
  *gcamorphSpacing = znzreadInt(fp);
  *gcamorphExp_k   = znzreadFloat(fp);

  return NO_ERROR;
}


// read TAG_GCAMORPH_LABELS data
int FStagsIO::read_gcamorph_labels(int x0, int y0, int z0, int ***gcamorphLabel)
{
  for (int x = 0; x < x0; x++)
    for (int y = 0; y < y0; y++)
      for (int z = 0; z < z0; z++)
        gcamorphLabel[x][y][z] = znzreadInt(fp);
  
  return NO_ERROR;
}


// read TAG_DOF (nifti header extension only)
int FStagsIO::read_dof(int *dof)
{
  *dof = znzreadInt(fp);

  return NO_ERROR;  
}


// read TAG_SCAN_PARAMETERS (nifti header extension only)
int FStagsIO::read_scan_parameters(MRI *mri, long long dlen)
{
  mri->te = znzreadFloat(fp);
  dlen -= sizeof(mri->te);
  
  mri->ti = znzreadFloat(fp);
  dlen -= sizeof(mri->ti);
  
  mri->flip_angle = znzreadDouble(fp);
  dlen -= sizeof(mri->flip_angle);
  
  // ??? todo: check how is fov calculated ???
  // skip fov, it can be calculated from other parameters
  // znzwriteFloat(mri->fov, fp);
  mri->FieldStrength = znzreadFloat(fp);
  dlen -= sizeof(mri->FieldStrength);
  
  mri->pedir = (char *)calloc(dlen + 1, sizeof(char));
  long long bytesread = znzread(mri->pedir, sizeof(char), dlen, fp);
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    printf("[DEBUG] read_scan_parameters(): bytesread=%lld, dlen=%lld\n", bytesread, dlen);

  return NO_ERROR;  
}


// read TAG_RAS_XFORM (nifti header extension only)
int FStagsIO::read_ras_xform(VOL_GEOM *vg)
{
  vg->x_r = znzreadFloat(fp); vg->x_a = znzreadFloat(fp); vg->x_s = znzreadFloat(fp);
  vg->y_r = znzreadFloat(fp); vg->y_a = znzreadFloat(fp); vg->y_s = znzreadFloat(fp);
  vg->z_r = znzreadFloat(fp); vg->z_a = znzreadFloat(fp); vg->z_s = znzreadFloat(fp);
  vg->c_r = znzreadFloat(fp); vg->c_a = znzreadFloat(fp); vg->c_s = znzreadFloat(fp);
  
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    vg->geomprint("[DEBUG] FStagsIO::read_ras_xform(): from FS header extension:\n");

  return NO_ERROR;  
}


// read TAG_RAS_XFORM_EXTRA (nifti header extension only)
int FStagsIO::read_ras_xform_extra(VOL_GEOM *vg)
{
  vg->s_r   = znzreadFloat(fp); vg->s_a    = znzreadFloat(fp); vg->s_s   = znzreadFloat(fp);
  vg->width = znzreadInt(fp);   vg->height = znzreadInt(fp);   vg->depth = znzreadInt(fp);
  vg->xsize = znzreadFloat(fp); vg->ysize  = znzreadFloat(fp); vg->zsize = znzreadFloat(fp);
  
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    vg->geomprint("[DEBUG] FStagsIO::read_ras_xform_extra(): from FS header extension:\n");

  return NO_ERROR;  
}


// tags.cpp::znzTAGskip()
int FStagsIO::skip_tag(int tag, long long len)
{
#if 1
  unsigned char *buf;
  int ret;

  buf = (unsigned char *)calloc(len, sizeof(unsigned char));
  if (buf == NULL) ErrorExit(ERROR_NOMEMORY, "FStagsIO::skip_tag(): tag=%-4d, failed to calloc %u bytes!\n", tag, len);
  ret = znzread(buf, sizeof(unsigned char), len, fp);
  free(buf);
  return (ret);
#else
  return (znzseek(fp, len, SEEK_CUR));  // doesn't work for gzipped files
#endif  
}


int FStagsIO::__write_matrix_niftiheaderext(MATRIX *M, int tag)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);

  if (tag > 0)
  {
    znzwriteInt(tag, fp);

    long long dlen = getlen_matrix(niftiheaderext, false);
    znzwriteLong(dlen, fp);
  }
  
  znzwriteFloat(M->rptr[1][1], fp);
  znzwriteFloat(M->rptr[1][2], fp);
  znzwriteFloat(M->rptr[1][3], fp);
  znzwriteFloat(M->rptr[1][4], fp);
  znzwriteFloat(M->rptr[2][1], fp);
  znzwriteFloat(M->rptr[2][2], fp);
  znzwriteFloat(M->rptr[2][3], fp);
  znzwriteFloat(M->rptr[2][4], fp);
  znzwriteFloat(M->rptr[3][1], fp);
  znzwriteFloat(M->rptr[3][2], fp);
  znzwriteFloat(M->rptr[3][3], fp);
  znzwriteFloat(M->rptr[3][4], fp);
  znzwriteFloat(M->rptr[4][1], fp);
  znzwriteFloat(M->rptr[4][2], fp);
  znzwriteFloat(M->rptr[4][3], fp);
  znzwriteFloat(M->rptr[4][4], fp);

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::__write_matrix_niftiheaderext() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", tag, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


MATRIX* FStagsIO::__read_matrix_niftiheaderext()
{
  MATRIX *M = MatrixAlloc(4, 4, MATRIX_REAL);

  M->rptr[1][1] = znzreadFloat(fp);
  M->rptr[1][2] = znzreadFloat(fp);
  M->rptr[1][3] = znzreadFloat(fp);
  M->rptr[1][4] = znzreadFloat(fp);
  M->rptr[2][1] = znzreadFloat(fp);
  M->rptr[2][2] = znzreadFloat(fp);
  M->rptr[2][3] = znzreadFloat(fp);
  M->rptr[2][4] = znzreadFloat(fp);
  M->rptr[3][1] = znzreadFloat(fp);
  M->rptr[3][2] = znzreadFloat(fp);
  M->rptr[3][3] = znzreadFloat(fp);
  M->rptr[3][4] = znzreadFloat(fp);
  M->rptr[4][1] = znzreadFloat(fp);
  M->rptr[4][2] = znzreadFloat(fp);
  M->rptr[4][3] = znzreadFloat(fp);
  M->rptr[4][4] = znzreadFloat(fp);

  return (M);
}


int FStagsIO::__write_mri_frames_niftiheaderext(MRI *mri)
{
  long long fstart = 0;
  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_MRI_FRAME, fp);

  long long dlen = getlen_mri_frames(mri, niftiheaderext, false);
  znzwriteLong(dlen, fp);
  
  for (int fno = 0; fno < mri->nframes; fno++) {
    MRI_FRAME *frame = &mri->frames[fno];
    znzwriteInt(frame->label, fp);
    znzwrite(frame->name, sizeof(char), STRLEN, fp);
    znzwriteFloat(frame->thresh, fp);
  }

  if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] FStagsIO::__write_mri_frames_niftiheaderext() TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_MRI_FRAME, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;  
}


int FStagsIO::__read_mri_frames_niftiheaderext(MRI *mri, long long len)
{
  long long fstart = znztell(fp);
  
  for (int fno = 0; fno < mri->nframes; fno++) {
    MRI_FRAME *frame = &mri->frames[fno];
    frame->label = znzreadInt(fp);
    znzread(frame->name, sizeof(char), STRLEN, fp);
    frame->thresh = znzreadFloat(fp);
  }

  long long fend = znztell(fp);
  len -= (fend - fstart);
  if (len > 0) {
    // previous version wrote more data under this TAG
    // if len > 0, it is reading file generated with previous version write
    // this is to skip those extra data, but the data read in are also wrong
    if (getenv("NII_HEADER_EXT_DEBUG") != NULL)
      printf("[DEBUG] __read_mri_frames_niftiheaderext() TAG = %-4d, bytes_read = %-6lld (%-6lld - %-6lld), skip extra bytes %lld\n", TAG_MRI_FRAME, fend-fstart, fstart, fend, len);
  
    char *buf = (char *)calloc(len, sizeof(char));
    znzread(buf, len, sizeof(char), fp);
    free(buf);
  }

  return NO_ERROR;
}
