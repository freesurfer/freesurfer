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
  static long long getlen_tag(int tag, long long len, bool niftiheaderext=false, bool addtaglength=true);
  static long long getlen_matrix(bool niftiheaderext, bool addtaglength=true);
  static long long getlen_old_colortable(COLOR_TABLE *ctab, bool niftiheaderext=false, bool addtaglength=true);
  static long long getlen_mri_frames(MRI *mri, bool addtaglength=true);
  static long long getlen_gcamorph_geom(bool niftiheaderext=false, bool addtaglength=true);
  static long long getlen_gcamorph_meta(bool addtaglength=true);
  static long long getlen_gcamorph_labels(int x, int y, int z, int len, bool niftiheaderext=false, bool addtaglength=true);
  static long long getlen_dof(int dof, bool addtaglength=true);
  static long long getlen_scan_parameters(MRI *mri, bool addtaglength=true);
  static long long getlen_ras_xform(MRI *mri, bool addtaglength=true);

  // methods to write various TAGs including tagid and len(tagdata) if the TAG has a length
  int write_tag(int tag, void *data, long long dlen);
  int write_matrix(MATRIX *M, int tag);
  int write_old_colortable(COLOR_TABLE *ctab);
  int write_mri_frames(MRI *mri);
  
  int write_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target);
  int write_gcamorph_meta(int warpFieldFormat, int gcamorphSpacing, double gcamorphExp_k);
  int write_gcamorph_labels(int x, int y, int z, int ***gcamorphLabel);  

  // these are for nifti header extension only
  int write_dof(int dof);
  int write_scan_parameters(MRI *mri);
  int write_ras_xform(MRI *mri);

  // retrieve tagid, datalength
  // if the TAG is in 'tagid len data' format, *plen = len(data);
  // otherwise, *plen = 0
  // after the call, the file pointer points to the data
  int read_tagid_len(long long *plen, int tagwithzerolen=0);

  // the file pointer should pass tagid and datalength, and point to the data
  int read_data(void *databuf, long long len);
  MATRIX* read_matrix();
  COLOR_TABLE* read_old_colortable();
  int read_mri_frames(MRI *mri, long len);
  
  int read_gcamorph_geom(VOL_GEOM *source, VOL_GEOM *target);
  int read_gcamorph_meta(int *warpFieldFormat, int *gcamorphSpacing, double *gcamorphExp_k);
  int read_gcamorph_labels(int x0, int y0, int z0, int ***gcamorphLabel);

  // for nifti header extension only
  int read_intent_encoded_version(int *version);
  int read_dof(int *dof);
  int read_scan_parameters(MRI *mri, long long dlen);
  int read_ras_xform(MRI *mri);
  
  // skip tag data (len of bytes)
  int skip_tag(int tag, long long len);

private:
  int __write_matrix_niftiheaderext(MATRIX *M, int tag=0);
  MATRIX* __read_matrix_niftiheaderext();
  
private:
  znzFile fp;

  bool niftiheaderext;
};

#endif  // FSTAGSIO_H

