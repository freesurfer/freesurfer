#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "warpfield.h"
#include "gcamorph.h"
#include "matrix.h"
#include "mri_circulars.h"
#include "mri_identify.h"

/* This class implements methods
 *   1. reads mgz warp file into GCAM
 *   2. converts GCAM to mgz warp format
 *   3. writes warp in mgz format (version = ((MGZ_INTENT_WARPMAP & 0xffff ) << 8) | MGH_VERSION).
 *
 * The warp file follows mgz format with these tags:
 *   TAG_GCAMORPH_GEOM             followed by gcamorph image (source) geom and gcamorph atlas (target) geom
 *   TAG_GCAMORPH_GEOM_PLUSSHEAR   followed by outputs under TAG_GCAMORPH_GEOM + shear components for image (source) geom and atlas (target) geom
 *   TAG_GCAMORPH_META             followed by data-length, 
 *         WARPFIELD_DTFMT_ABS_CRS|WARPFIELD_DTFMT_DISP_CRS|WARPFIELD_DTFMT_ABS_RAS|WARPFIELD_DTFMT_DISP_RAS
 *         spacing (int)
 *         exp_k   (double)
 *   TAG_GCAMORPH_AFFINE           followed by data-length (1600), matrix data
 * 
 * The data array (width x height x depth x nframes) is indexed by atlas CRS.
 *     frame 0 - image voxel ABS coordinate C, image voxel DISP coordinate C, 
 *               RAS ABS coordinate X, or RAS DISP coordinate X
 *     frame 1 - image voxel ABS coordinate R, image voxel DISP coordinate R,
 *               RAS ABS coordinate Y, or RAS DISP coordinate Y
 *     frame 2 - image voxel ABS coordinate S, image voxel DISP coordinate S,
 *               RAS ABS coordinate Z, or RAS DISP coordinate Z
 *
 * Here are the 4 data formats supported:
 *     WARPFIELD_DTFMT_ABS_CRS   - CRS coordinates in image space
 *     WARPFIELD_DTFMT_DISP_CRS  - displacement CRS, delta = image_CRS - atlas_CRS
 *     WARPFIELD_DTFMT_ABS_RAS   - RAS coordinates in image space
 *     WARPFIELD_DTFMT_DISP_RAS  - displacement RAS, delta = image_RAS - atlas_RAS
 */

// constructor
Warpfield::Warpfield()
{
  __warpmap = NULL;  __warpmap_inv = NULL;
  __freewarpmap = true;
  __invert = 0;
  __mgzVersion = ((MGZ_INTENT_WARPMAP & 0xffff ) << 8) | MGH_VERSION;
  __FSLWarpfield = false;  

  __srcRAS2Vox = NULL;
  __srcVox2RAS = NULL;
  __dstRAS2Vox = NULL;
  __dstVox2RAS = NULL; 
}


// constructor
Warpfield::Warpfield(MRI* mri)
{
  if (mri->intent != MGZ_INTENT_WARPMAP)
  {
    printf("[ERROR] Can't create Warpfield object with non MGZ_INTENT_WARPMAP MRI\n");
    exit(1);
  }
  
  __warpmap = mri;  __warpmap_inv = NULL;
  __freewarpmap = false;
  __invert = 0;
  __mgzVersion = ((MGZ_INTENT_WARPMAP & 0xffff ) << 8) | MGH_VERSION;
  __FSLWarpfield = false;

  // pre-calulated transform matrix taking shears into consideration
  __srcRAS2Vox = __warpmap->gcamorph_image_vg.get_RAS2Vox(0, true);
  __srcVox2RAS = __warpmap->gcamorph_image_vg.get_Vox2RAS(0, true);
  __dstRAS2Vox = __warpmap->gcamorph_atlas_vg.get_RAS2Vox(0, true);
  __dstVox2RAS = __warpmap->gcamorph_atlas_vg.get_Vox2RAS(0, true);  
}


// destructor
Warpfield::~Warpfield()
{
  if (__freewarpmap && __warpmap != NULL)
    MRIfree(&__warpmap);

  if (__freewarpmap && __warpmap_inv != NULL)
    MRIfree(&__warpmap_inv);

  if (__srcRAS2Vox != NULL)
    MatrixFree(&__srcRAS2Vox);
  if (__srcVox2RAS != NULL)
    MatrixFree(&__srcVox2RAS);
  if (__dstRAS2Vox != NULL)
    MatrixFree(&__dstRAS2Vox);
  if (__dstVox2RAS != NULL)
    MatrixFree(&__dstVox2RAS);
}


// convert given MGH_MORPH (.m3z/.m3d) to mgz warp
MRI* Warpfield::convert(const char *fname, const int dataformat, int doGCAMsampleMorph)
{
  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN)
  {
    printf("ERROR: unknown dataformat\n");
    exit(1);
  }
  
  int type = mri_identify(fname);
  if (type != MGH_MORPH)  // .m3z/.m3d
  {
    printf("[ERROR] Warpfield::convert(): %s is not in m3z format\n", fname);
    exit(1);
  }

  GCA_MORPH *gcam = GCAMread(fname);

  return convert(gcam, dataformat);
}

// convert GCAM to mgz warp
//
// similar functionality is also implemented in
//   MRI *GCAMwriteWarpToMRI(const GCA_MORPH *gcam, MRI *mri_warp);         (gcamorph.cpp)
//       2026-02-17 YJH:
//       It looks like that before calling GCAMwriteWarpToMRI(),
//       GCAM needs to be converted to GCAM_VOX. The output MRI is in disp_crs.
//   void write_world(const string& fname, GCAM* gcam, bool is_lps=false);  (mri_warp_convert.cpp)
//   void write_voxel(const string& fname, GCAM* gcam);                     (mri_warp_convert.cpp)
//   MRI *GCAMtoMRI(GCAM *gcam, MRI *mri);                                  (gcamorph.cpp)
//       2026-02-17 YJH:
//       It seems that the output MRI is either abs_crs or abs_ras
//       depending on if gcan->type is GCAM_VOX or GCAM_RAS.
//   MRI *GCAMwriteMRI(GCAM *gcam, MRI *mri, int frame);                    (gcamorph.cpp)
MRI* Warpfield::convert(GCA_MORPH *gcam, const int dataformat, int doGCAMsampleMorph)
{
  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN)
  {
    printf("[ERROR] unknown dataformat\n");
    exit(1);
  }

  // the logic here only work with GCAM_VOX, convert GCAM_RAS to GCAM_VOX first
  if (gcam->type == GCAM_RAS)
  {
    printf("converting GCAM from GCAM_RAS to GCAM_VOX\n");
    GCAMrasToVox(gcam, NULL);
  }
  
  printf("[INFO] Warpfield::convert(): converting GCAM%s ...\n", (doGCAMsampleMorph) ? " (do GCAMsampleMorph)" : "");
  
  printf("[INFO] Warpfield::convert(): gcam       [%d x %d x %d]\n", gcam->width, gcam->height, gcam->depth);
  printf("[INFO] Warpfield::convert(): gcam image [%d x %d x %d]\n", gcam->image.width, gcam->image.height, gcam->image.depth);
  printf("[INFO] Warpfield::convert(): gcam atlas [%d x %d x %d]\n", gcam->atlas.width, gcam->atlas.height, gcam->atlas.depth);  

  // create MRI using gcam dimensions
  // copy geom from gcam->atlas to __warpmap (width, height, deph are not copied)
  // gcam->image vol geom and gcam->atlas vol geom will be saved in mgz under TAG_GCAMORPH_GEOM
  __warpmap = new MRI({gcam->width, gcam->height, gcam->depth, 3}, MRI_FLOAT);
  MRIcopyVolGeomToMRI(__warpmap, &gcam->atlas);
  //__warpmap = new MRI(gcam->atlas, MRI_FLOAT, 3, 0);  //__warpmap = new MRI({gcam->atlas.width, gcam->atlas.height, gcam->atlas.depth, 3}, MRI_FLOAT);

  // TAG_GCAMORPH_META
  __mgzVersion = ((MGZ_INTENT_WARPMAP & 0xffff ) << 8) | MGH_VERSION;
  __warpmap->intent  = MGZ_INTENT_WARPMAP;
  __warpmap->version = __mgzVersion;
  __warpmap->warpFieldFormat = dataformat;
  __warpmap->gcamorphSpacing = gcam->spacing;
  __warpmap->gcamorphExp_k = gcam->exp_k;

  __warpmap->gcamorph_image_vg = gcam->image;
  __warpmap->gcamorph_atlas_vg = gcam->atlas;
  
  if (gcam->m_affine)
  {
    printf("[DEBUG] Warpfield::convert() gcam->m_affine (spacing=%d, exp-k=%.2f, det=%.2f):\n", gcam->spacing, gcam->exp_k, gcam->det);
    MatrixPrint(stdout, gcam->m_affine);
    __warpmap->gcamorphAffine = MatrixCopy(gcam->m_affine, NULL);
    printf("[DEBUG] Warpfield::convert() __warpmap->gcamorphAffine (spacing=%d, exp-k=%.2f)\n", __warpmap->gcamorphSpacing, __warpmap->gcamorphExp_k);
    MatrixPrint(stdout, __warpmap->gcamorphAffine);    
  }
    
  // allocate __warpmap->gcamorphLabel
  __warpmap->initGCAMorphLabel();
  
  // pre-calulated transform matrix
  __srcRAS2Vox = gcam->image.get_RAS2Vox();
  __srcVox2RAS = gcam->image.get_Vox2RAS();
  __dstRAS2Vox = gcam->atlas.get_RAS2Vox();
  __dstVox2RAS = gcam->atlas.get_Vox2RAS();

  // pre-allocated MATRIX
  MATRIX *image_CRS  = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *image_RAS  = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *atlas_CRS0 = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *atlas_RAS0 = MatrixAlloc(4, 1, MATRIX_REAL); 
  
  // ??? what about gcamn->invalid ???
  int out_of_gcam_count = 0;
  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	float fcs = 0, frs = 0, fss = 0;
	if (doGCAMsampleMorph)
	{
          // (c, r, s) is in atlas (target) volume, (fcs, frs, fss) is in image (source) volume
	  // (c, r, s) => (fcs, frs, fss)	
	  int out_of_gcam = GCAMsampleMorph(gcam, (float)c, (float)r, (float)s, &fcs, &frs, &fss);
	  if (out_of_gcam)
	  {
	    out_of_gcam_count++;
	    continue;
	  }
	}
	else
	{
	  // this will work only if gcam and gcam->atlas have the same size
	  fcs = gcam->nodes[c][r][s].x;
	  frs = gcam->nodes[c][r][s].y;
	  fss = gcam->nodes[c][r][s].z;
	}

	setWarp(c, r, s, fcs, frs, fss, gcam->nodes[c][r][s].label);
      }  // s
    }  // r
  }  // c

  printf("[INFO] Warpfield::convert(): total out of range voxel count: %d\n", out_of_gcam_count);

  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS0);
  MatrixFree(&atlas_RAS0);
  
  return __warpmap;
}


// invert M3z into 3-fram MRI warp map
// !!!It has not been tested!!!
MRI* Warpfield::invert(const char *fname, const int dataformat)
{
  printf("Warpfield::invert(const char*, const int) is not implemented\n");
  return NULL;

  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN)
  {
    printf("ERROR: unknown dataformat\n");
    exit(1);
  }
  
  int type = mri_identify(fname);
  if (type != MGH_MORPH)  // .m3z/.m3d
  {
    printf("[ERROR] Warpfield::invert() %s is not in m3z format\n", fname);
    exit(1);
  }
    
  GCA_MORPH *gcam = GCAMread(fname);

  return invert(gcam, dataformat);  
}

// invert GCAM
// !!!It has not been tested!!!
MRI* Warpfield::invert(GCA_MORPH *gcam, const int dataformat)
{
  printf("Warpfield::invert(GCA_MORPH*, const int) is not implemented\n");
  return NULL;
    
  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN)
  {
    printf("ERROR: unknown dataformat\n");
    exit(1);
  }

  // the logic here only work with GCAM_VOX, convert GCAM_RAS to GCAM_VOX first
  if (gcam->type == GCAM_RAS)
  {
    printf("converting GCAM from GCAM_RAS to GCAM_VOX\n");
    GCAMrasToVox(gcam, NULL);
  }
    
  printf("[INFO] Warpfield::invert(): inverting GCAM ...\n");
  __invert = 1;
  
  // create GCAM inverse
  gcam->spacing = 1;

  // purpose of tempMri is just to pass image dimensions to GCAMinvert()
  MRI *tempMri = new MRI(gcam->image, MRI_FLOAT, 3, 0);
  GCAMinvert(gcam, tempMri);
  MRIfree(&tempMri);

  // create MRI using image vol_geom
  __warpmap = new MRI(gcam->image, MRI_FLOAT, 4, 0);

  __mgzVersion = ((MGZ_INTENT_WARPMAP_INV & 0xffff ) << 8) | MGH_VERSION;
  __warpmap->intent  = MGZ_INTENT_WARPMAP_INV;
  __warpmap->version = __mgzVersion;
  __warpmap->warpFieldFormat = dataformat;
  __warpmap->gcamorphSpacing = gcam->spacing;
  __warpmap->gcamorphExp_k   = gcam->exp_k;

  __warpmap->gcamorph_image_vg = gcam->image;
  __warpmap->gcamorph_atlas_vg = gcam->atlas;
  
  if (gcam->m_affine)
    __warpmap->gcamorphAffine = MatrixCopy(gcam->m_affine, NULL);
  
  // pre-calculated transform matrix
  __srcRAS2Vox = gcam->image.get_RAS2Vox();
  __srcVox2RAS = gcam->image.get_Vox2RAS();
  __dstRAS2Vox = gcam->atlas.get_RAS2Vox();
  __dstVox2RAS = gcam->atlas.get_Vox2RAS();  

  // pre-allocated MATRIX
  MATRIX *dst_CRS  = MatrixAlloc(4, 1, MATRIX_REAL);
  MATRIX *dst_RAS  = MatrixAlloc(4, 1, MATRIX_REAL);	  
  MATRIX *src_CRS0 = MatrixAlloc(4, 1, MATRIX_REAL);	    
  MATRIX *src_RAS0 = MatrixAlloc(4, 1, MATRIX_REAL);	    
  
  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	float fct = 0, frt = 0, fst = 0;
        // (c, r, s) is in image (source) volume, (fct, frt, fst) is in atlas (target) volume
	// (c, r, s) => (fct, frt, fst)	
	int out_of_gcam = GCAMsampleInverseMorph(gcam, (float)c, (float)r, (float)s, &fct, &frt, &fst);
	if (out_of_gcam)
	  continue;

	__warpmap->gcamorphLabel[c][r][s] = gcam->nodes[c][r][s].label;
	
        if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
        {
	  // in target (atlas) voxel space
          MRIsetVoxVal(__warpmap, c, r, s, 0, fct);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, frt);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, fst);
	}
	else if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
	{
	  // delta = src_CRS - dst_CRS
	  MRIsetVoxVal(__warpmap, c, r, s, 0, (float)c - fct);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, (float)r - frt);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, (float)s - fst);
	}
	else if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS ||
                 dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
	{
	  // convert (fct, frt, fst) to dst_RAS
	  dst_CRS->rptr[1][1] = fct;
          dst_CRS->rptr[2][1] = frt;
          dst_CRS->rptr[3][1] = fst;
          dst_CRS->rptr[4][1] = 1;

	  MatrixMultiplyD(__dstVox2RAS, dst_CRS, dst_RAS);

	  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
	  {
	    // in target (atlas) RAS space
	    MRIsetVoxVal(__warpmap, c, r, s, 0, dst_RAS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, dst_RAS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, dst_RAS->rptr[3][1]);
	  }
	  else // dataformat == WARPFIELD_DTFMT_DISP_RAS
	  {
	    // convert (c, r, s) to src_RAS
	    src_CRS0->rptr[1][1] = c;
            src_CRS0->rptr[2][1] = r;
            src_CRS0->rptr[3][1] = s;
            src_CRS0->rptr[4][1] = 1;

	    MatrixMultiplyD(__srcVox2RAS, src_CRS0, src_RAS0);

	    // delta = src_RAS - dst_RAS
	    MRIsetVoxVal(__warpmap, c, r, s, 0, src_RAS0->rptr[1][1] - dst_RAS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, src_RAS0->rptr[2][1] - dst_RAS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, src_RAS0->rptr[3][1] - dst_RAS->rptr[3][1]);
	  }
	}  // WARPFIELD_DTFMT_ABS_RAS || WARPFIELD_DTFMT_DISP_RAS
      }  // s
    }  // r
  }  // c

  MatrixFree(&dst_CRS);
  MatrixFree(&dst_RAS);
  MatrixFree(&src_CRS0);
  MatrixFree(&src_RAS0);
  
  return __warpmap;
}


// read 3-frame MRI warp map into __warpmap,
// copy the warp into a GCAM,
// return GCAM created
//
// similar functionality is also implemented in
//   int GCAMreadWarpFromMRI(GCA_MORPH *gcam, const MRI *mri_warp, int DeformationFlag)     (gcamorph.cpp)
//       2026-02-17 YJH:
//       I think GCAMreadWarpFromMRI() requires the input mri_warp to be in abs_crs or disp_crs interpretation
//       because gcamn->origx, origy, origz are initialized to be target crs in GCAMalloc().
//   GCAM* read_voxel(const string& warp_file, const string& src_geom);                     (mri_warp_convert.cpp)
//       [origx, origy, origz], [xn, yn, zn] are set to dst [c, r, s]
//   GCAM* read_world(const string& warp_file, const string& src_geom, bool is_lps=false);  (mri_warp_convert.cpp)
//       [origx, origy, origz], [xn, yn, zn] are set to dst [c, r, s]
GCA_MORPH *Warpfield::read(const char *fname)
{
  int type = mri_identify(fname);
  if (type != MRI_MGH_FILE && type != NII_FILE)
  {
    printf("[ERROR] Warpfield::read(): %s is not in mgz/nii format\n", fname);
    return NULL;
  }
  
  // the function doesn't handle invert warp
  __mgzVersion = ((MGZ_INTENT_WARPMAP & 0xffff ) << 8) | MGH_VERSION;

  __warpmap = MRIread(fname);  //mghRead(fname);
  if (__warpmap == NULL)
  {
    printf("[ERROR] Warpfield::read() failed reading %s\n", fname);
    return NULL;
  }

  if (__warpmap->version != __mgzVersion)
  {
    printf("[ERROR] %s is not mgz warp file\n", fname);
    return NULL;
  }

  GCA_MORPH *gcam = GCAMalloc(__warpmap->width, __warpmap->height, __warpmap->depth);
  if (gcam == NULL)
    return NULL;
  
  gcam->det = 1;
  gcam->spacing = __warpmap->gcamorphSpacing;
  gcam->exp_k   = __warpmap->gcamorphExp_k;
  
  gcam->type = GCAM_VOX;  
  gcam->image = __warpmap->gcamorph_image_vg;
  gcam->atlas = __warpmap->gcamorph_atlas_vg;

  // pre-calulated transform matrix taking shears into consideration
  __srcRAS2Vox = gcam->image.get_RAS2Vox(0, true);
  __srcVox2RAS = gcam->image.get_Vox2RAS(0, true);
  __dstRAS2Vox = gcam->atlas.get_RAS2Vox(0, true);
  __dstVox2RAS = gcam->atlas.get_Vox2RAS(0, true);

  /* YJH 2024-12-19
   * If the function is called from mri_warp_convert, 
   * this will make the source/target geom in the converted warp shearless.
   * Move the removal of shear components from transformed volume geometry to GCAMmorphToAtlas().
   * // remove shear components
   * gcam->image.shearless_components();
   * gcam->atlas.shearless_components();
   */

  if (__warpmap->gcamorphAffine)
  {
    printf("[DEBUG] Warpfield::read() __warpmap->gcamorphAffine (spacing=%d, exp-k=%.2f):\n", __warpmap->gcamorphSpacing, __warpmap->gcamorphExp_k);
    MatrixPrint(stdout, __warpmap->gcamorphAffine);
    gcam->m_affine = MatrixCopy(__warpmap->gcamorphAffine, NULL);
    gcam->det = MatrixDeterminant(gcam->m_affine);
    printf("[DEBUG] Warpfield::read() gcam->m_affine (spacing=%d, exp-k=%.2f, det=%.2f):\n", gcam->spacing, gcam->exp_k, gcam->det);
    MatrixPrint(stdout, gcam->m_affine);
  }

  if (__warpmap->gcamorphLabel)
  {
    printf("[DEBUG] Warpfield::read() gcam->status = GCAM_LABELED\n");
    gcam->status = GCAM_LABELED;
  }
  
  // pre-allocated MATRIX
  MATRIX *image_CRS  = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *image_RAS  = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *atlas_CRS0 = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *atlas_RAS0 = MatrixAlloc(4, 1, MATRIX_REAL);   
  
  // mri_warp_convert::readFSL2() uses GCAMreadWarpFromMRI(gcam, mri_warp, 0)
  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	GCA_MORPH_NODE *gcamn = &gcam->nodes[c][r][s];
        gcamn->origx = (float)c;
        gcamn->origy = (float)r;
        gcamn->origz = (float)s;
	gcamn->xn = c;
        gcamn->yn = r;
        gcamn->zn = s;

	if (__warpmap->gcamorphLabel)
          gcamn->label = __warpmap->gcamorphLabel[c][r][s];
	
	// ??? mark invalid for each node
	// gcamn->invalid = GCAM_POSITION_INVALID, GCAM_AREA_INVALID, GCAM_VALID
	
	if (__warpmap->warpFieldFormat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS ||
	    __warpmap->warpFieldFormat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
	{	  
	  if (__warpmap->warpFieldFormat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
	  {
            gcamn->x = MRIgetVoxVal(__warpmap, c, r, s, 0);
	    gcamn->y = MRIgetVoxVal(__warpmap, c, r, s, 1);
	    gcamn->z = MRIgetVoxVal(__warpmap, c, r, s, 2);
	  }
	  else // WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS
	  {
            gcamn->x = MRIgetVoxVal(__warpmap, c, r, s, 0) + gcamn->origx;
	    gcamn->y = MRIgetVoxVal(__warpmap, c, r, s, 1) + gcamn->origy;
	    gcamn->z = MRIgetVoxVal(__warpmap, c, r, s, 2) + gcamn->origz;	  
	  }	  
	}
	else if (__warpmap->warpFieldFormat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS ||
		 __warpmap->warpFieldFormat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
	{
	  if (__warpmap->warpFieldFormat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
	  {
            image_RAS->rptr[1][1] = MRIgetVoxVal(__warpmap, c, r, s, 0);
            image_RAS->rptr[2][1] = MRIgetVoxVal(__warpmap, c, r, s, 1);
            image_RAS->rptr[3][1] = MRIgetVoxVal(__warpmap, c, r, s, 2);
            image_RAS->rptr[4][1] = 1;
	  }
	  else // WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS
	  {
	    atlas_CRS0->rptr[1][1] = c;
            atlas_CRS0->rptr[2][1] = r;
            atlas_CRS0->rptr[3][1] = s;
            atlas_CRS0->rptr[4][1] = 1;
            MatrixMultiplyD(__dstVox2RAS, atlas_CRS0, atlas_RAS0);
	    
            image_RAS->rptr[1][1] = MRIgetVoxVal(__warpmap, c, r, s, 0) + atlas_RAS0->rptr[1][1];
            image_RAS->rptr[2][1] = MRIgetVoxVal(__warpmap, c, r, s, 1) + atlas_RAS0->rptr[2][1];
            image_RAS->rptr[3][1] = MRIgetVoxVal(__warpmap, c, r, s, 2) + atlas_RAS0->rptr[3][1];
            image_RAS->rptr[4][1] = 1;
	  }

	  // compute image_CRS from image_RAS
	  MatrixMultiplyD(__srcRAS2Vox, image_RAS, image_CRS);
	  gcamn->x = image_CRS->rptr[1][1];
	  gcamn->y = image_CRS->rptr[2][1];
	  gcamn->z = image_CRS->rptr[3][1];
	}
      } // s
    } // r
  } // c

  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS0);
  MatrixFree(&atlas_RAS0);
     
  return gcam;
}


// write out the 3-frame MRI warping map
int Warpfield::write(const char *fname, bool asFSLWarp)
{
  if (asFSLWarp && !__FSLWarpfield)
    toFSLWarpfield();
    
  int type = mri_identify(fname);
  if (type != MRI_MGH_FILE && type != NII_FILE)
  {
    printf("[ERROR] Warpfield::write(): %s is not in mgz/nii format\n", fname);
    exit(1);
  }

  __mgzVersion = ((MGZ_INTENT_WARPMAP & 0xffff ) << 8) | MGH_VERSION;
  if (__invert)
    __mgzVersion = ((MGZ_INTENT_WARPMAP_INV & 0xffff ) << 8) | MGH_VERSION;

  int ret = MRIwrite(__warpmap, fname, NULL, MGZ_INTENT_UNKNOWN, asFSLWarp);  //mghWrite(__warpmap, fname);
  if (ret)
    printf("ERROR: Warpfield::write(%s)\n", fname);
  
  return ret;
}


// set source coordinates at target [c,r,s] based on dataformat
// (fcs, frs, fss) is absolute CRS in source voxel space
void Warpfield::setWarp(int c, int r, int s, float fcs, float frs, float fss, int label)
{
  int dataformat = __warpmap->warpFieldFormat;
  
  // pre-allocated MATRIX
  MATRIX *image_CRS  = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *image_RAS  = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *atlas_CRS0 = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *atlas_RAS0 = MatrixAlloc(4, 1, MATRIX_REAL); 

  __warpmap->gcamorphLabel[c][r][s] = label;

  // ??? what about gcamn->invalid ???
  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
  {
    // in source (unmorphed, image) voxel space
    MRIsetVoxVal(__warpmap, c, r, s, 0, fcs);
    MRIsetVoxVal(__warpmap, c, r, s, 1, frs);
    MRIsetVoxVal(__warpmap, c, r, s, 2, fss);
  }
  else if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
  {
    // set the displacement: delta = image_CRS - atlas_CRS
    MRIsetVoxVal(__warpmap, c, r, s, 0, fcs - (float)c);
    MRIsetVoxVal(__warpmap, c, r, s, 1, frs - (float)r);
    MRIsetVoxVal(__warpmap, c, r, s, 2, fss - (float)s);	     
  }
  else if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS ||
           dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
  {
    // convert (fcs, frs, fss) to image_RAS
    image_CRS->rptr[1][1] = fcs;
    image_CRS->rptr[2][1] = frs;
    image_CRS->rptr[3][1] = fss;
    image_CRS->rptr[4][1] = 1;

    MatrixMultiplyD(__srcVox2RAS, image_CRS, image_RAS);

    if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
    {
      // in source (unmorphed, image) RAS space
      MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1]);
      MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1]);
      MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1]);
    }
    else // dataformat == WARPFIELD_DTFMT_DISP_RAS
    {
      atlas_CRS0->rptr[1][1] = c;
      atlas_CRS0->rptr[2][1] = r;
      atlas_CRS0->rptr[3][1] = s;
      atlas_CRS0->rptr[4][1] = 1;

      MatrixMultiplyD(__dstVox2RAS, atlas_CRS0, atlas_RAS0);
	    
      // set the displacement: delta = image_RAS - atlas_RAS
      MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1] - atlas_RAS0->rptr[1][1]);
      MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1] - atlas_RAS0->rptr[2][1]);
      MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1] - atlas_RAS0->rptr[3][1]);
    }
  }  // WARPFIELD_DTFMT_ABS_RAS || WARPFIELD_DTFMT_DISP_RAS

  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS0);
  MatrixFree(&atlas_RAS0);  
}


// change the warp field format to the given type
void Warpfield::changeFormat(const int newformat)
{
  int format = __warpmap->warpFieldFormat;
  __warpmap->warpFieldFormat = newformat;
  
  if (format == newformat)
    return;

  if (format == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
    __changeFormatFrom_abs_crs(newformat);
  else if (format == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
    __changeFormatFrom_disp_crs(newformat);
  else if (format == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
    __changeFormatFrom_abs_ras(newformat);
  else if (format == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
    __changeFormatFrom_disp_ras(newformat);
}


void Warpfield::__changeFormatFrom_abs_crs(const int newformat)
{
  // pre-allocated MATRIX
  MATRIX *image_CRS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *image_RAS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *atlas_CRS = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *atlas_RAS = MatrixAlloc(4, 1, MATRIX_REAL);
  
  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	// abs_crs (image_crs)
	image_CRS->rptr[1][1] = MRIgetVoxVal(__warpmap, c, r, s, 0);
        image_CRS->rptr[2][1] = MRIgetVoxVal(__warpmap, c, r, s, 1);
        image_CRS->rptr[3][1] = MRIgetVoxVal(__warpmap, c, r, s, 2);
        image_CRS->rptr[4][1] = 1;	

	// atlas_crs
	atlas_CRS->rptr[1][1] = c;
	atlas_CRS->rptr[2][1] = r;
	atlas_CRS->rptr[3][1] = s;
	atlas_CRS->rptr[4][1] = 1;
	if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
	{
	  // abs_crs => disp_crs
	  MRIsetVoxVal(__warpmap, c, r, s, 0, image_CRS->rptr[1][1] - atlas_CRS->rptr[1][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, image_CRS->rptr[2][1] - atlas_CRS->rptr[2][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, image_CRS->rptr[3][1] - atlas_CRS->rptr[3][1]);
	}
	else
	{
	    // abs_crs => abs_ras (image_ras) 
            MatrixMultiplyD(__srcVox2RAS, image_CRS, image_RAS);
	    if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
	    {
	      // abs_ras (image_ras)
	      MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1]);
	      MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1]);
	      MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1]);
	    }
	    else if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
	    {
	      // atlas_ras
	      MatrixMultiplyD(__dstVox2RAS, atlas_CRS, atlas_RAS);

	      // abs_ras => disp_ras	      
	      MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1] - atlas_RAS->rptr[1][1]);
	      MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1] - atlas_RAS->rptr[2][1]);
	      MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1] - atlas_RAS->rptr[3][1]);
	    }
	}
      }
    }
  }
  
  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS);
  MatrixFree(&atlas_RAS);
}


void Warpfield::__changeFormatFrom_disp_crs(const int newformat)
{
  // pre-allocated MATRIX
  MATRIX *delta_CRS = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *image_CRS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *image_RAS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *atlas_CRS = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *atlas_RAS = MatrixAlloc(4, 1, MATRIX_REAL);
  
  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	// disp_crs
	delta_CRS->rptr[1][1] = MRIgetVoxVal(__warpmap, c, r, s, 0);
	delta_CRS->rptr[2][1] = MRIgetVoxVal(__warpmap, c, r, s, 1);
	delta_CRS->rptr[3][1] = MRIgetVoxVal(__warpmap, c, r, s, 2);
	delta_CRS->rptr[4][1] = 1;

	// atlas_crs
	atlas_CRS->rptr[1][1] = c;
	atlas_CRS->rptr[2][1] = r;
	atlas_CRS->rptr[3][1] = s;
	atlas_CRS->rptr[4][1] = 1;
	
	// disp_crs => abs_crs (image_crs)	
	image_CRS->rptr[1][1] = delta_CRS->rptr[1][1] + atlas_CRS->rptr[1][1];
        image_CRS->rptr[2][1] = delta_CRS->rptr[2][1] + atlas_CRS->rptr[2][1];
        image_CRS->rptr[3][1] = delta_CRS->rptr[3][1] + atlas_CRS->rptr[3][1];
        image_CRS->rptr[4][1] = 1;
	if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
	{
	  // abs_crs (image_crs)
	  MRIsetVoxVal(__warpmap, c, r, s, 0, image_CRS->rptr[1][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, image_CRS->rptr[2][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, image_CRS->rptr[3][1]);	  
	}
	else
	{
	  // abs_crs => abs_ras (image_ras)
	  MatrixMultiplyD(__srcVox2RAS, image_CRS, image_RAS);
	  if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
	  {
	    // abs_ras (image_ras)
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1]);
	  }
	  else if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
	  {
	    // atlas_ras
	    MatrixMultiplyD(__dstVox2RAS, atlas_CRS, atlas_RAS);
	    
	    // abs_ras => disp_ras	    
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1] - atlas_RAS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1] - atlas_RAS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1] - atlas_RAS->rptr[3][1]);
	  }
	}
      }
    }
  }

  MatrixFree(&delta_CRS);
  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS);
  MatrixFree(&atlas_RAS);  
}


void Warpfield::__changeFormatFrom_abs_ras(const int newformat)
{
  // pre-allocated MATRIX
  MATRIX *image_CRS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *image_RAS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *atlas_CRS = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *atlas_RAS = MatrixAlloc(4, 1, MATRIX_REAL);
  
  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	// abs_ras (image_ras)
	image_RAS->rptr[1][1] = MRIgetVoxVal(__warpmap, c, r, s, 0);
	image_RAS->rptr[2][1] = MRIgetVoxVal(__warpmap, c, r, s, 1);
	image_RAS->rptr[3][1] = MRIgetVoxVal(__warpmap, c, r, s, 2);
	image_RAS->rptr[4][1] = 1;

	// atlas_crs
	atlas_CRS->rptr[1][1] = c;
	atlas_CRS->rptr[2][1] = r;
	atlas_CRS->rptr[3][1] = s;
	atlas_CRS->rptr[4][1] = 1;	
	if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
	{
	  // atlas_ras
	  MatrixMultiplyD(__dstVox2RAS, atlas_CRS, atlas_RAS);
	  
	  // abs_ras => disp_ras
	  MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1] - atlas_RAS->rptr[1][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1] - atlas_RAS->rptr[2][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1] - atlas_RAS->rptr[3][1]);
	}
	else
	{
	  // abs_ras => abs_crs (image_crs)
	  MatrixMultiplyD(__srcRAS2Vox, image_RAS, image_CRS);
	  if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
	  {
	    // abs_crs (image_crs)
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_CRS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_CRS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_CRS->rptr[3][1]);	    
	  }
	  else if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
	  {
	    // abs_crs => disp_crs
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_CRS->rptr[1][1] - atlas_CRS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_CRS->rptr[2][1] - atlas_CRS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_CRS->rptr[3][1] - atlas_CRS->rptr[3][1]);	    
	  }
	}
      }
    }
  }

  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS);
  MatrixFree(&atlas_RAS);
}


void Warpfield::__changeFormatFrom_disp_ras(const int newformat)
{ 
  // pre-allocated MATRIX
  MATRIX *delta_RAS = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *image_CRS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *image_RAS = MatrixAlloc(4, 1, MATRIX_REAL); 
  MATRIX *atlas_CRS = MatrixAlloc(4, 1, MATRIX_REAL);  
  MATRIX *atlas_RAS = MatrixAlloc(4, 1, MATRIX_REAL);

  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	// disp_ras
	delta_RAS->rptr[1][1] = MRIgetVoxVal(__warpmap, c, r, s, 0);
	delta_RAS->rptr[2][1] = MRIgetVoxVal(__warpmap, c, r, s, 1);
	delta_RAS->rptr[3][1] = MRIgetVoxVal(__warpmap, c, r, s, 2);
	delta_RAS->rptr[4][1] = 1;

	// atlas_crs
	atlas_CRS->rptr[1][1] = c;
	atlas_CRS->rptr[2][1] = r;
	atlas_CRS->rptr[3][1] = s;
	atlas_CRS->rptr[4][1] = 1;

	// atlas_ras
	MatrixMultiplyD(__dstVox2RAS, atlas_CRS, atlas_RAS);

	// image_ras (abs_ras)
	image_RAS->rptr[1][1] = delta_RAS->rptr[1][1] + atlas_RAS->rptr[1][1];
	image_RAS->rptr[2][1] = delta_RAS->rptr[2][1] + atlas_RAS->rptr[2][1];
	image_RAS->rptr[3][1] = delta_RAS->rptr[3][1] + atlas_RAS->rptr[3][1];
	image_RAS->rptr[4][1] = 1;	
	if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
	{
	  // abs_ras
	  MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1]);
	}
	else
	{
	  // abs_ras => abs_crs (image_crs)
	  MatrixMultiplyD(__srcRAS2Vox, image_RAS, image_CRS);
	  if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
	  {
	    // abs_crs (image_crs)
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_CRS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_CRS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_CRS->rptr[3][1]);	    
	  }
	  else if (newformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
	  {
	    // abs_crs => disp_crs
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_CRS->rptr[1][1] - atlas_CRS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_CRS->rptr[2][1] - atlas_CRS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_CRS->rptr[3][1] - atlas_CRS->rptr[3][1]);	    
	  }
	}
      }
    }
  }

  MatrixFree(&delta_RAS);
  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS);
  MatrixFree(&atlas_RAS);
}


// convert Freesurfer 3D morph to FSL warpfield
void Warpfield::toFSLWarpfield()
{
  // FSL warpfields are RAS space
  if (__warpmap->warpFieldFormat != WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
    changeFormat(WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS);
  
  // todo:
  // adjust warpfield values for non-negative determinant
  // ...

  __FSLWarpfield = true;  
}


// apply warpmap to MRI
int Warpfield::applyWarp(const MRI *inmri, MRI *outmri)
{
  printf("Warpfield::applyWarp(const MRI *, MRI*) is not implemented\n");
  return 0;
}


// apply warpmap to surface
// ?? apply the inverted __warpfield from scr to dst??
int Warpfield::applyWarp(const MRIS *insurf, MRIS *outsurf)
{
  printf("Warpfield::applyWarp(const MRIS *, MRIS*) is not implemented\n");  
  return 0;
}
