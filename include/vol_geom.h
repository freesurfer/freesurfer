#ifndef VOL_GEOM_H
#define VOL_GEOM_H

//#include "gca.h"
#include "mri.h"

struct VOL_GEOM;

struct VOL_GEOM
{
  int           valid;   /* whether this is a
                                    valid info or not (1 valid, 0 not valid) */
  int           width ;
  int           height ;
  int           depth ;
  float         xsize ;
  float         ysize ;
  float         zsize ;
  float         x_r, x_a, x_s;
  float         y_r, y_a, y_s;
  float         z_r, z_a, z_s;
  float         c_r, c_a, c_s;
  char          fname[STRLEN];  // volume filename

  MATRIX *vox2ras;
  MATRIX *ras2vox;
  MATRIX *tkregvox2ras;
  MATRIX *tkregras2vox;

  VOL_GEOM(const char *srcVol = NULL);  // utils/transform.cpp:void initVolGeom(VOL_GEOM *vg)
  ~VOL_GEOM();

  void init();  // utils/transform.cpp:void initVolGeom(VOL_GEOM *vg)
  void init_average_305();
  void print();  // utils/transform.cpp:void vg_print(const VOL_GEOM *vg)
  void copyFromMRI(const MRI *src); // utils/transform.cpp:void getVolGeom(const MRI *src, VOL_GEOM *dst)
  void copyToMRI(MRI *dst);  // utils/transform.cpp:void useVolGeomToMRI(const VOL_GEOM *src, MRI *dst)
  void write(FILE *fp);  // utils/transform.cpp:void writeVolGeom(FILE *fp, const VOL_GEOM *vg)
  void read(FILE *fp);   // utils/transform.cpp:void readVolGeom(FILE *fp, VOL_GEOM *vg)
  MATRIX *getVox2RAS(int base = 0);     // utils/transform.cpp:MATRIX *vg_i_to_r(const VOL_GEOM *vg)
  MATRIX *getRAS2Vox(int base = 0);     // utils/transform.cpp:MATRIX *vg_r_to_i(const VOL_GEOM *vg)
  MATRIX *getTkregVox2RAS(int base = 0);  // utils/transform.cpp:MATRIX *TkrVox2RASfromVolGeom(const VOL_GEOM *vg)
  MATRIX *getTkregRAS2Vox(int base = 0);  // utils/transform.cpp:MATRIX *TkrRAS2VoxfromVolGeom(const VOL_GEOM *vg)
  int isEqual(const VOL_GEOM* vg2); // utils/transform.cpp:int vg_isEqual(const VOL_GEOM *vg1, const VOL_GEOM *vg2)
  int isNotEqualThresh(const VOL_GEOM *vg2, const double thresh);  // utils/transform.cpp:int vg_isNotEqualThresh(const VOL_GEOM *vg1, const VOL_GEOM *vg2, const double thresh)
  //MATRIX *getVoxelToRasXform(int base = 0);  // utils/transform.cpp:MATRIX *VGgetVoxelToRasXform(VOL_GEOM *vg, MATRIX *m, int base)
  //MATRIX *getRasToVoxelXform(int base = 0);  //utils/transform.cpp:MATRIX *VGgetRasToVoxelXform(VOL_GEOM *vg, MATRIX *m, int base)
  MRI* allocMRI(int type, int nframes, int HeaderOnly);  // utils/transform.cpp:MRI *MRIallocFromVolGeom(VOL_GEOM *vg, int type, int nframes, int HeaderOnly)
  void copy(VOL_GEOM *dst);  // utils/transform.cpp:void copyVolGeom(const VOL_GEOM *src, VOL_GEOM *dst)

  // VOL_GEOM created with srcVol
  // utils/transform.cpp:int mincFindVolume(const char *line, const char *line2, char **srcVol, char **dstVol)
  // utils/transform.cpp:void mincGetVolumeInfo(const char *srcVol, VOL_GEOM *vgSrc)
  // utils/transform.cpp:void mincGetVolInfo(const char *infoline, const char *infoline2, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)

  // utils/transform.cpp:static void LTAgetV2V(MATRIX *mod, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)
  // utils/transform.cpp:static void LTAgetR2R(MATRIX *mod, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)
  // utils/transform.cpp:int TransformGetSrcVolGeom(const TRANSFORM *transform, VOL_GEOM *vg)
  // utils/transform.cpp:int TransformGetDstVolGeom(const TRANSFORM *transform, VOL_GEOM *vg)
  // utils/transform.cpp:MRI *MRIallocFromVolGeom(VOL_GEOM *vg, int type, int nframes, int HeaderOnly)
  // utils/transform.cpp:void copyVolGeom(const VOL_GEOM *src, VOL_GEOM *dst)

  // utils/mri.cpp:int MRIcopyVolGeomToMRI(MRI *mri, const VOL_GEOM *vg)
  // utils/mri.cpp:int MRIcopyVolGeomFromMRI(const MRI *mri, VOL_GEOM *vg)
  //void setfrom(GCA *gca);  // utils/gca.cpp:void GCAsetVolGeom(GCA *gca, VOL_GEOM *vg)

  // mri_modify/mri_modify.cpp:int get_option(int argc, char *argv[], VOL_GEOM &vg)
};

typedef VOL_GEOM VG;
#endif
