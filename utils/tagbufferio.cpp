#include <string.h>
#include <stdlib.h>

#include "tagbufferio.h"
#include "fsbufferio.h"
#include "tags.h"
#include "matrix.h"
#include "colortab.h"
#include "mri.h"

TAGbufferIO::TAGbufferIO(int niftiheaderextension)
{
  tagbuflen = 0;
  
  int tagbufsize = 16 * 1024;
  fsbufio = new FSbufferIO(tagbufsize);
}


TAGbufferIO::~TAGbufferIO()
{
  delete fsbufio;
  
  tagbuflen = 0;
  fsbufio = NULL;
}


unsigned char *TAGbufferIO::getbuffer()
{
  return fsbufio->getbuffer();
}


int TAGbufferIO::writedata(void *data, long long dlen)
{
  tagbuflen += fsbufio->write_buf(data, dlen);
  return tagbuflen;
}


int TAGbufferIO::writetag(int tag, void *data, long long dlen)
{
  //int bytewritten = 0;
  tagbuflen += fsbufio->write_int(tag);

  // ??? todo: check if tag is length-less
  if (tag != TAG_OLD_COLORTABLE &&
      tag != TAG_GCAMORPH_GEOM && tag != TAG_GCAMORPH_TYPE && tag != TAG_GCAMORPH_LABELS)
    tagbuflen += fsbufio->write_long(dlen);
  
  tagbuflen += fsbufio->write_buf(data, dlen);

  return tagbuflen;
}


int TAGbufferIO::writematrix(MATRIX *M, int tag)
{
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
  
  return writetag(tag, matbuf, dlen);
}


int TAGbufferIO::write_old_colortable(COLOR_TABLE *ctab)
{
  int ctbufsize = 8 * 1024;

  FSbufferIO ctbufio(ctbufsize);
  
  CTABwriteIntoBinaryBuffer(ctab, &ctbufio);
  
  return writetag(TAG_OLD_COLORTABLE, ctbufio.getbuffer(), ctbufio.getbufferlen());
}


// based on mriio.cpp::znzTAGwriteMRIframes()
int TAGbufferIO::write_mri_frames(MRI *mri)
{
  int mriframebufsize = 4 * 1024;
  FSbufferIO mriframebufio(mriframebufsize);
  
  //long long fstart, fend, here;

  // write some extra space so that we have enough room (can't seek in zz files)
  //long long len = 10 * mri->nframes * sizeof(MRI_FRAME);
  //znzTAGwriteStart(fp, TAG_MRI_FRAME, &fstart, len);
  //here = znztell(fp);
  for (int fno = 0; fno < mri->nframes; fno++) {
    MRI_FRAME *frame = &mri->frames[fno];
    mriframebufio.write_int(frame->type);
    mriframebufio.write_float(frame->TE);
    mriframebufio.write_float(frame->TR);
    mriframebufio.write_float(frame->flip);
    mriframebufio.write_float(frame->TI);
    mriframebufio.write_float(frame->TD);
    mriframebufio.write_float(frame->TM);
    mriframebufio.write_int(frame->sequence_type);
    mriframebufio.write_float(frame->echo_spacing);
    mriframebufio.write_float(frame->echo_train_len);
    for (int i = 0; i < 3; i++) mriframebufio.write_float(frame->read_dir[i]);
    for (int i = 0; i < 3; i++) mriframebufio.write_float(frame->pe_dir[i]);
    for (int i = 0; i < 3; i++) mriframebufio.write_float(frame->slice_dir[i]);
    mriframebufio.write_int(frame->label);
    mriframebufio.write_buf(frame->name, STRLEN);
    mriframebufio.write_int(frame->dof);
    if (frame->m_ras2vox && frame->m_ras2vox->rows > 0) {
      // znzWriteMatrix(fp, frame->m_ras2vox, 0);
      mriframebufio.write_int(0);
      mriframebufio.write_long(MATRIX_STRLEN);
      mriframebufio.write_matrix(frame->m_ras2vox);
    }
    else {
      MATRIX *m = MatrixAlloc(4, 4, MATRIX_REAL);
      // znzWriteMatrix(fp, m, 0);
      mriframebufio.write_int(0);
      mriframebufio.write_long(MATRIX_STRLEN);
      mriframebufio.write_matrix(m);      
      MatrixFree(&m);
    }
    mriframebufio.write_float(frame->thresh);
    mriframebufio.write_int(frame->units);
    if (frame->type == FRAME_TYPE_DIFFUSION_AUGMENTED)  // also store diffusion info
    {
      mriframebufio.write_double(frame->DX);
      mriframebufio.write_double(frame->DY);
      mriframebufio.write_double(frame->DZ);

      mriframebufio.write_double(frame->DR);
      mriframebufio.write_double(frame->DP);
      mriframebufio.write_double(frame->DS);
      mriframebufio.write_double(frame->bvalue);
      mriframebufio.write_double(frame->TM);

      mriframebufio.write_long(frame->D1_ramp);
      mriframebufio.write_long(frame->D1_flat);
      mriframebufio.write_double(frame->D1_amp);

      mriframebufio.write_long(frame->D2_ramp);
      mriframebufio.write_long(frame->D2_flat);
      mriframebufio.write_double(frame->D2_amp);

      mriframebufio.write_long(frame->D3_ramp);
      mriframebufio.write_long(frame->D3_flat);
      mriframebufio.write_double(frame->D3_amp);

      mriframebufio.write_long(frame->D4_ramp);
      mriframebufio.write_long(frame->D4_flat);
      mriframebufio.write_double(frame->D4_amp);
    }
  }
  /*
  fend = znztell(fp);
  len -= (fend - here);  // unused space
  if (len > 0) {
    char *buf = (char *)calloc(len, sizeof(char));
    mriframebufio.write_buf(buf, len);
    free(buf);
  }
  znzTAGwriteEnd(fp, fend);
  */

  return writetag(TAG_MRI_FRAME, mriframebufio.getbuffer(), mriframebufio.getbufferlen());  
}

    
int TAGbufferIO::write_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target)
{
  int bufsize = 2 * 1024;
  FSbufferIO geombufio(bufsize);
  
  geombufio.write_geom(source);
  geombufio.write_geom(target);
  
  return writetag(TAG_GCAMORPH_GEOM, geombufio.getbuffer(), geombufio.getbufferlen());
}


int TAGbufferIO::write_gcamorph_meta(int warpFieldFormat, int gcamorphSpacing, double gcamorphExp_k)
{
  long long buflen = 0;
  int bufsize = 32;
  FSbufferIO metabufio(bufsize);
  
  buflen += metabufio.write_int(warpFieldFormat);

  buflen += metabufio.write_int(gcamorphSpacing);
  
  buflen += metabufio.write_float(gcamorphExp_k);

  return writetag(TAG_GCAMORPH_META, metabufio.getbuffer(), buflen);
}


int TAGbufferIO::write_gcamorph_labels(int x0, int y0, int z0, int   ***gcamorphLabel)
{
  int bufsize = x0 * y0 * z0 * sizeof(int);
  FSbufferIO labelbufio(bufsize);
  
  for (int x = 0; x < x0; x++)
    for (int y = 0; y < y0; y++)
      for (int z = 0; z < z0; z++)
        labelbufio.write_int(gcamorphLabel[x][y][z]);
  
  return writetag(TAG_GCAMORPH_LABELS, labelbufio.getbuffer(), labelbufio.getbufferlen());
}

