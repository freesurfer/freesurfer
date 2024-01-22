#ifndef FSTAGSIO_H
#define FSTAGSIO_H

#include "fio.h"

// forward declarations
struct VOL_GEOM;
class  MRI;
struct MATRIX;
struct COLOR_TABLE;

class FStagsIO
{
public:
  FStagsIO(znzFile fp, bool niftiheaderext=false);
  ~FStagsIO();

  // the following getlen_*() methods return TAG length as one of the two:
  //   TAG w/o  length: tagid + len(tagdata)
  //   TAG w/ a length: tagid + sizeof(long long) + len(tagdata)
  static long long getlen_tag(int tag, long long len, bool niftiheaderext=false);
  static long long getlen_matrix();
  static long long getlen_old_colortable(COLOR_TABLE *ctab, bool niftiheaderext=false);
  static long long getlen_mri_frames(MRI *mri);
  static long long getlen_gcamorph_geom(bool niftiheaderext=false);
  static long long getlen_gcamorph_meta();
  static long long getlen_gcamorph_labels(int x, int y, int z, int len, bool niftiheaderext=false); 

  // methods to write various TAGs including tagid and len(tagdata) if the TAG has a length
  int write_tag(int tag, void *data, long long dlen);
  int write_matrix(MATRIX *M, int tag);
  int write_old_colortable(COLOR_TABLE *ctab);
  int write_mri_frames(MRI *mri);
  
  int write_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target);
  int write_gcamorph_meta(int warpFieldFormat, int gcamorphSpacing, double gcamorphExp_k);
  int write_gcamorph_labels(int x, int y, int z, int ***gcamorphLabel);  

  // retrieve tagid, datalength
  // if the TAG is in 'tagid len data' format, *plen = len(data);
  // otherwise, *plen = 0
  // after the call, the file pointer points to the data
  int read_tagid_len(long long *plen, int tagwithzerolen);

  // the file pointer should pass tagid and datalength, and point to the data
  int read_data(void *databuf, long long len);
  MATRIX* read_matrix();
  COLOR_TABLE* read_old_colortable();
  int read_mri_frames(MRI *mri, long len);
  
  int read_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target);
  int read_gcamorph_meta(int *warpFieldFormat, int *gcamorphSpacing, double *gcamorphExp_k);
  int read_gcamorph_labels(int x0, int y0, int z0, int ***gcamorphLabel);

  // skip tag data (len of bytes)
  int skip_tag(int tag, long long len);
private:
  znzFile fp;

  bool niftiheaderext;
};

#endif  // FSTAGSIO_H

