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
//        set 'addtaglength = false', niftiheaderext is ignored
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
        (tag != TAG_OLD_COLORTABLE &&
         tag != TAG_GCAMORPH_GEOM && tag != TAG_GCAMORPH_TYPE && tag != TAG_GCAMORPH_LABELS))
      dlen += sizeof(long long);
  }

  dlen += len;

  return dlen;
}


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


// ??? can we get rid of TAG_MRI_FRAME
// ??? it writes '10 * mri->nframes * sizeof(MRI_FRAME)' bytes to disk
// ??? only thing set seems to be the identity matrix m_ras2vox
long long FStagsIO::getlen_mri_frames(MRI *mri, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }
  dlen += 10 * mri->nframes * sizeof(MRI_FRAME);

  return dlen;
}


long long FStagsIO::getlen_gcamorph_geom(bool niftiheaderext, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
  
    if (niftiheaderext)
      dlen += sizeof(long long);
  }

  // this needs to be consistent with write_gcamorph_geom()/VOL_GEOM.write()
  int geom_len = 4 * sizeof(int) + 5 * sizeof(float) + 512;
  dlen += 2 * geom_len;

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


long long FStagsIO::getlen_ras_xform(MRI *mri, bool addtaglength)
{
  long long dlen = 0;
  if (addtaglength)
  {
    dlen += 4;
    dlen += sizeof(long long);
  }

  // this needs to be consistent with write_scan_parameters()
  dlen += sizeof(mri->x_r); dlen += sizeof(mri->x_a); dlen += sizeof(mri->x_s);
  dlen += sizeof(mri->y_r); dlen += sizeof(mri->y_a); dlen += sizeof(mri->y_s);
  dlen += sizeof(mri->z_r); dlen += sizeof(mri->z_a); dlen += sizeof(mri->z_s);
  dlen += sizeof(mri->c_r); dlen += sizeof(mri->c_a); dlen += sizeof(mri->c_s);

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
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  znzwriteInt(tag, fp);
  
  // ??? todo: check if tag is length-less
  if (niftiheaderext ||
      (tag != TAG_OLD_COLORTABLE &&
       tag != TAG_GCAMORPH_GEOM && tag != TAG_GCAMORPH_TYPE && tag != TAG_GCAMORPH_LABELS))
    znzwriteLong(dlen, fp);
  
  znzwrite(data, sizeof(char), dlen, fp);

  long long fend = znztell(fp);
    
  if (Gdiag & DIAG_INFO)
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld) (%-6lld)\n", tag, fend-fstart, fstart, fend, dlen);
  
  return NO_ERROR;
}


// tags.cpp::znzWriteMatrix()
int FStagsIO::write_matrix(MATRIX *M, int tag)
{
  if (Gdiag & DIAG_INFO)
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

  if (Gdiag & DIAG_INFO)
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
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_OLD_COLORTABLE, fp);
  if (niftiheaderext)
  {
    long long dlen = getlen_old_colortable(ctab, false, false);
    znzwriteLong(dlen, fp);
  }
  
  znzCTABwriteIntoBinary(ctab, fp);

  if (Gdiag & DIAG_INFO)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_OLD_COLORTABLE, fend-fstart, fstart, fend);
  }
    
  return NO_ERROR;
}


// mriio.cpp::znzTAGwriteMRIframes()
int FStagsIO::write_mri_frames(MRI *mri)
{
  long long fstart = 0;
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  // write some extra space so that we have enough room (can't seek in zz files)
  long long len = 10 * mri->nframes * sizeof(MRI_FRAME);

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

  if (Gdiag & DIAG_INFO)
  {
    fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_MRI_FRAME, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


// TAG_GCAMORPH_GEOM is in length-less format
int FStagsIO::write_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target)
{
  long long fstart = 0;
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_GCAMORPH_GEOM, fp);
  
  if (niftiheaderext)
  {
    long long dlen = getlen_gcamorph_geom(false, false);
    znzwriteLong(dlen, fp);
  }
  
  source->write(fp);
  target->write(fp);

  if (Gdiag & DIAG_INFO)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_GCAMORPH_GEOM, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


// TAG_GCAMORPH_META
int FStagsIO::write_gcamorph_meta(int warpFieldFormat, int gcamorphSpacing, double gcamorphExp_k)
{
  long long fstart = 0;
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_GCAMORPH_META, fp);

  long long dlen = getlen_gcamorph_meta(false);
  znzwriteLong(dlen, fp);
  znzwriteInt(warpFieldFormat, fp);
  znzwriteInt(gcamorphSpacing, fp);
  znzwriteFloat(gcamorphExp_k, fp);

  if (Gdiag & DIAG_INFO)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_GCAMORPH_META, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


// TAG_GCAMORPH_LABELS is in length-less format
int FStagsIO::write_gcamorph_labels(int x0, int y0, int z0, int ***gcamorphLabel)
{
  long long fstart = 0;
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_GCAMORPH_LABELS, fp);
  
  if (niftiheaderext)
  {
    long long dlen = getlen_gcamorph_labels(x0, y0, z0, sizeof(int), false, false);
    znzwriteLong(dlen, fp);
  }
  
  for (int x = 0; x < x0; x++)
    for (int y = 0; y < y0; y++)
      for (int z = 0; z < z0; z++)
        znzwriteInt(gcamorphLabel[x][y][z], fp);

  if (Gdiag & DIAG_INFO)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_GCAMORPH_LABELS, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;
}


// write TAG_DOF (nifti header extension only)
int FStagsIO::write_dof(int dof)
{
  long long fstart = 0;
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_DOF, fp);

  long long dlen = getlen_dof(dof, false);
  znzwriteLong(dlen, fp);
  znzwriteInt(dof, fp);

  if (Gdiag & DIAG_INFO)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_DOF, fend-fstart, fstart, fend);
  }

  return NO_ERROR;  
}


// write TAG_SCAN_PARAMETERS (nifti header extension only)
int FStagsIO::write_scan_parameters(MRI *mri)
{
  long long fstart = 0;
  if (Gdiag & DIAG_INFO)
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

  if (Gdiag & DIAG_INFO)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_SCAN_PARAMETERS, fend-fstart, fstart, fend);
  }
  
  return NO_ERROR;  
}

  
// write TAG_RAS_XFORM (nifti header extension only)
int FStagsIO::write_ras_xform(MRI *mri)
{
  long long fstart = 0;
  if (Gdiag & DIAG_INFO)
    fstart = znztell(fp);
  
  znzwriteInt(TAG_RAS_XFORM, fp);

  long long dlen = getlen_ras_xform(mri, false);
  znzwriteLong(dlen, fp);
  znzwriteFloat(mri->x_r, fp); znzwriteFloat(mri->x_a, fp); znzwriteFloat(mri->x_s, fp);
  znzwriteFloat(mri->y_r, fp); znzwriteFloat(mri->y_a, fp); znzwriteFloat(mri->y_s, fp);
  znzwriteFloat(mri->z_r, fp); znzwriteFloat(mri->z_a, fp); znzwriteFloat(mri->z_s, fp);
  znzwriteFloat(mri->c_r, fp); znzwriteFloat(mri->c_a, fp); znzwriteFloat(mri->c_s, fp);

  if (Gdiag & DIAG_INFO)
  {
    long long fend = znztell(fp);
    printf("[DEBUG] TAG = %-4d, dlen = %-6lld (%-6lld - %-6lld)\n", TAG_RAS_XFORM, fend-fstart, fstart, fend);
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

  if (Gdiag & DIAG_INFO)
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
int FStagsIO::read_mri_frames(MRI *mri, long len)
{
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
    if (Gdiag & DIAG_INFO)
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
    char *buf = (char *)calloc(len, sizeof(char));
    znzread(buf, len, sizeof(char), fp);
    free(buf);
  }

  return NO_ERROR;  
}


// read TAG_GCAMORPH_GEOM data
int FStagsIO::read_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target)
{
  source->read(fp);
  target->read(fp);

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


// read TAG_INTENT_ENCODED_VERSION (nifti header extension only)
int FStagsIO::read_intent_encoded_version(int *version)
{
  *version = znzreadInt(fp);
  
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
  znzread(mri->pedir, sizeof(char), dlen, fp);

  return NO_ERROR;  
}


// read TAG_RAS_XFORM (nifti header extension only)
int FStagsIO::read_ras_xform(MRI *mri)
{
  mri->x_r = znzreadFloat(fp); mri->x_a = znzreadFloat(fp); mri->x_s = znzreadFloat(fp);
  mri->y_r = znzreadFloat(fp); mri->y_a = znzreadFloat(fp); mri->y_s = znzreadFloat(fp);
  mri->z_r = znzreadFloat(fp); mri->z_a = znzreadFloat(fp); mri->z_s = znzreadFloat(fp);
  mri->c_r = znzreadFloat(fp); mri->c_a = znzreadFloat(fp); mri->c_s = znzreadFloat(fp);
  
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
  if (Gdiag & DIAG_INFO)
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

  if (Gdiag & DIAG_INFO)
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
