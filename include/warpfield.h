#ifndef WARPFIELD_H
#define WARPFIELD_H

// forward declarations
struct GCA_MORPH;
struct VOL_GEOM;
class  MRI;
struct MRIS;
struct MATRIX;

struct WarpfieldDTFMT{
  static const int WARPFIELD_DTFMT_UNKNOWN  = -1;
  // one of these will be saved under TAG_WARPFIELD_DTFMT in mgz    
  static const int WARPFIELD_DTFMT_ABS_CRS  = 0;
  static const int WARPFIELD_DTFMT_DISP_CRS = 1;
  static const int WARPFIELD_DTFMT_ABS_RAS  = 2;
  static const int WARPFIELD_DTFMT_DISP_RAS = 3;
};

class Warpfield
{
public:
  // src = image, dst/trg = atlas
  Warpfield();
  ~Warpfield();

  // convert M3z into 4-frame MRI warp map
  MRI *convert(const char *fname, const int dataformat=WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS, int doGCAMsampleMorph=0);
  MRI *convert(GCA_MORPH *gcam, const int dataformat=WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS, int doGCAMsampleMorph=0);

  // !!!invert functions have not been tested!!!
  // invert M3z into 4-frame MRI warp map
  MRI *invert(const char *fname, const int dataformat=WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS);
  MRI *invert(GCA_MORPH *gcam, const int dataformat=WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS);
  
  // read 4-frame MRI warp map into __warpmap
  GCA_MORPH *read(const char *fname);
  
  // write 4-frame MRI warp map saved in __warpmap to disk
  int write(const char *fname);

  // creat mgz warp with given dimensions, src/dst VOL_GEOM, dataformat
  void create(int width, int height, int depth, const VOL_GEOM& srcGeom, const VOL_GEOM& dstGeom,
	      int spacing, double exp_k, const MATRIX *affine,
	      const int dataformat=WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS);
  // set source coordinates at target [c,r,s] based on dataformat
  void setWarp(int c, int r, int s, float fcs, float frs, float fss, int label);
  
  // apply warpmap to MRI/MRIS
  int applyWarp(const MRI *inmri, MRI *outmri);
  int applyWarp(const MRIS *insurf, MRIS *outsurf);

private:
  int __mgzVersion;               // mgz version

  int __invert;                   // __warpmap is inverted
  
  MATRIX *__srcRAS2Vox;           // source ras2vox
  MATRIX *__srcVox2RAS;           // source vox2ras
  MATRIX *__dstRAS2Vox;           // target ras2vox
  MATRIX *__dstVox2RAS;           // target vox2ras

  MRI *__warpmap;                 // 4-frame MRI warping map (dst => src)
  MRI *__warpmap_inv;             // inverted __warpmap (src => dst)
};

#endif // WARPFIELD_H
