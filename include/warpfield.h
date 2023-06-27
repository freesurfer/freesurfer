#ifndef WARPFIELD_H
#define WARPFIELD_H

// forward declarations
struct GCA_MORPH;
class MRI;
struct MRIS;

class Warpfield
{
public:
  // src = image, dst/trg = atlas
  enum WarpfieldDTFMT {
    WARPFIELD_DTFMT_UNKNOWN = -1,
    // one of these will be saved under TAG_WARPFIELD_DTFMT in mgz    
    WARPFIELD_DTFMT_ABS_CRS,
    WARPFIELD_DTFMT_DISP_CRS,
    WARPFIELD_DTFMT_ABS_RAS,
    WARPFIELD_DTFMT_DISP_RAS
  };

  Warpfield();
  ~Warpfield();

  // convert M3z into 3-frame MRI warp map
  int convert(const char *fname, const WarpfieldDTFMT dataformat=WARPFIELD_DTFMT_ABS_CRS);
  int convert(GCA_MORPH *gcam, const WarpfieldDTFMT dataformat=WARPFIELD_DTFMT_ABS_CRS);

  // invert M3z into 3-fram MRI warp map
  int invert(const char *fname, const WarpfieldDTFMT dataformat=WARPFIELD_DTFMT_ABS_CRS);
  int invert(GCA_MORPH *gcam, const WarpfieldDTFMT dataformat=WARPFIELD_DTFMT_ABS_CRS);
  
  // read 3-frame MRI warp map into __warpmap
  int read(const char *fname);
  
  // write 3-frame MRI warp map saved in __warpmap to disk
  int write(const char *fname);

  // apply warpmap to MRI/MRIS
  int applyWarp(const MRI *inmri, MRI *outmri);
  int applyWarp(const MRIS *insurf, MRIS *outsurf);

private:
  int __mri_version;              // mri version
  WarpfieldDTFMT __dataformat;    // WarpfieldDT
  MRI *__warpmap;                 // 3-frame MRI warping map (dst => src)
  MRI *__warpmap_inv;             // inverted __warpmap (src => dst)
};

#endif // WARPFIELD_H
