#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "warpfield.h"
#include "gcamorph.h"
#include "matrix.h"

// constructor
Warpfield::Warpfield()
{
  __warpmap = NULL;  __warpmap_inv = NULL;
  __invert = 0;
  __mgzVersion = ((MGZ_WARPMAP & 0xff ) << 8) | MGH_VERSION;
  __dataformat = WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN;
  
  __srcRAS2Vox = NULL;
  __srcVox2RAS = NULL;
  __dstRAS2Vox = NULL;
  __dstVox2RAS = NULL; 
}


// destructor
Warpfield::~Warpfield()
{
  if (__warpmap != NULL)
    MRIfree(&__warpmap);

  if (__warpmap_inv != NULL)
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


int Warpfield::convert(const char *fname, const int dataformat, int doGCAMsampleMorph)
{
  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN)
  {
    printf("ERROR: unknown dataformat\n");
    exit(1);
  }
  
  int type = TransformFileNameType((char *)fname);
  if (type != MORPH_3D_TYPE)
  {
    printf("ERROR: %s is not in m3z format\n", fname);
    exit(1);
  }
    
  GCA_MORPH *gcam = GCAMread(fname);

  return convert(gcam, __dataformat);
}

// convert GCAM
int Warpfield::convert(GCA_MORPH *gcam, const int dataformat, int doGCAMsampleMorph)
{
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
  
  int (*nintfunc)( double );
  nintfunc = &nint;

  printf("[INFO] Warpfield::convert(): converting GCAM%s ...\n", (doGCAMsampleMorph) ? " (do GCAMsampleMorph)" : "");
  
  printf("[INFO] Warpfield::convert(): gcam       [%d x %d x %d]\n", gcam->width, gcam->height, gcam->depth);
  printf("[INFO] Warpfield::convert(): gcam image [%d x %d x %d]\n", gcam->image.width, gcam->image.height, gcam->image.depth);
  printf("[INFO] Warpfield::convert(): gcam atlas [%d x %d x %d]\n", gcam->atlas.width, gcam->atlas.height, gcam->atlas.depth);  
    
  // create MRI using atlas vol_geom
  __warpmap = new MRI(gcam->atlas, MRI_FLOAT, 3, 0);  //__warpmap = new MRI({gcam->atlas.width, gcam->atlas.height, gcam->atlas.depth, 3}, MRI_FLOAT);
  __dataformat = dataformat;

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
	
        if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
        {
	  // in source (unmorphed, image) voxel space
          MRIsetVoxVal(__warpmap, c, r, s, 0, fcs);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, frs);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, fss);
	}
	else if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
	{
	  // set the displacement: delta = image_CRS - atlas_CRS
	  MRIsetVoxVal(__warpmap, c, r, s, 0, fcs - (float)c);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, frs - (float)r);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, fss - (float)s);	     
#if 0
	  /* the followng logic was copied from mri_warp_convert::write_voxel().
           * This conversion doesn't make sense: image_CRS0 = atlas2image_vox * atlas_CRS0
           */
	  // convert CRS0 in target (atlas) voxel space => source (image) voxel space
	  MATRIX *atlas_CRS0 = MatrixAlloc(4, 1, MATRIX_REAL);
	  atlas_CRS0->rptr[1][1] = c;
          atlas_CRS0->rptr[2][1] = r;
          atlas_CRS0->rptr[3][1] = s;
          atlas_CRS0->rptr[4][1] = 1;

	  MATRIX *atlas2image_vox = MatrixMultiply(gcam->image.get_RAS2Vox(), gcam->atlas.get_Vox2RAS(), NULL);
	  MATRIX *image_CRS0 = MatrixMultiply(atlas2image_vox, atlas_CRS0, NULL);  // CRS0 is now in source (image) voxel space

	  // set the displacement in source (image) voxel space
	  MRIsetVoxVal(__warpmap, c, r, s, 0, fcs - image_CRS0->rptr[1][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, frs - image_CRS0->rptr[2][1]);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, fss - image_CRS0->rptr[3][1]);
#endif
	}
	else if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS ||
                 __dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
	{
	  // convert (fcs, frs, fss) to image_RAS
	  image_CRS->rptr[1][1] = nintfunc(fcs);
          image_CRS->rptr[2][1] = nintfunc(frs);
          image_CRS->rptr[3][1] = nintfunc(fss);
          image_CRS->rptr[4][1] = 1;

	  MatrixMultiply(__srcVox2RAS, image_CRS, image_RAS);

	  if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
	  {
	    // in source (unmorphed, image) RAS space
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1]);
	  }
	  else // __dataformat == WARPFIELD_DTFMT_DISP_RAS
	  {
	    atlas_CRS0->rptr[1][1] = c;
            atlas_CRS0->rptr[2][1] = r;
            atlas_CRS0->rptr[3][1] = s;
            atlas_CRS0->rptr[4][1] = 1;

            MatrixMultiply(__dstVox2RAS, atlas_CRS0, atlas_RAS0);
	    
	    // set the displacement: delta = image_RAS - atlas_RAS
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1] - atlas_RAS0->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1] - atlas_RAS0->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1] - atlas_RAS0->rptr[3][1]);
#if 0
	    /* the logic has the same problem as WARPFIELD_DTFMT_DISP_CRS
             */
	    // convert CRS0 in target (atlas) voxel space => source (image) voxel space
	    MATRIX *atlas_CRS0 = MatrixAlloc(4, 1, MATRIX_REAL);
	    atlas_CRS0->rptr[1][1] = c;
            atlas_CRS0->rptr[2][1] = r;
            atlas_CRS0->rptr[3][1] = s;
            atlas_CRS0->rptr[4][1] = 1;

	    MATRIX *atlas2image_vox = MatrixMultiply(gcam->image.get_RAS2Vox(), gcam->atlas.get_Vox2RAS(), NULL);
	    MATRIX *image_CRS0 = MatrixMultiply(atlas2image_vox, atlas_CRS0, NULL);  // CRS0 is now in source (image) voxel space
            MATRIX *image_RAS0 = MatrixMultiply(gcam->image.get_Vox2RAS(), image_CRS0, NULL); // RAS0 is now in source (image) RAS space
	    
	    // set the displacement in source (image) RAS space
	    MRIsetVoxVal(__warpmap, c, r, s, 0, image_RAS->rptr[1][1] - image_RAS0->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, image_RAS->rptr[2][1] - image_RAS0->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, image_RAS->rptr[3][1] - image_RAS0->rptr[3][1]);
#endif
	  }
	}  // WARPFIELD_DTFMT_ABS_RAS || WARPFIELD_DTFMT_DISP_RAS
      }  // s
    }  // r
  }  // c

  printf("[INFO] Warpfield::convert(): total out of range voxel count: %d\n", out_of_gcam_count);

  MatrixFree(&image_CRS);
  MatrixFree(&image_RAS);
  MatrixFree(&atlas_CRS0);
  MatrixFree(&atlas_RAS0);
  
  return 0;
}


// invert M3z into 3-fram MRI warp map
int Warpfield::invert(const char *fname, const int dataformat)
{
  if (dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN)
  {
    printf("ERROR: unknown dataformat\n");
    exit(1);
  }
  
  int type = TransformFileNameType((char *)fname);
  if (type != MORPH_3D_TYPE)
  {
    printf("ERROR: %s is not in m3z format\n", fname);
    exit(1);
  }
    
  GCA_MORPH *gcam = GCAMread(fname);

  return invert(gcam, __dataformat);  
}

// invert GCAM
int Warpfield::invert(GCA_MORPH *gcam, const int dataformat)
{
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
    
  int (*nintfunc)( double );
  nintfunc = &nint;

  printf("[INFO] Warpfield::invert(): inverting GCAM ...\n");
  __invert = 1;
  
  // create GCAM inverse
  gcam->spacing = 1;
  
  // purpose of tempMri is just to pass image dimensions to GCAMinvert()
  MRI *tempMri = new MRI(gcam->image, MRI_FLOAT, 3, 0);
  GCAMinvert(gcam, tempMri);
  MRIfree(&tempMri);

  // create MRI using image vol_geom
  __warpmap = new MRI(gcam->image, MRI_FLOAT, 3, 0);
  __dataformat = dataformat;

  // pre-calculated transform matrix
  __srcRAS2Vox = gcam->atlas.get_RAS2Vox();
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
	
        if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS)
        {
	  // in target (atlas) voxel space
          MRIsetVoxVal(__warpmap, c, r, s, 0, fct);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, frt);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, fst);
	}
	else if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS)
	{
	  // delta = src_CRS - dst_CRS
	  MRIsetVoxVal(__warpmap, c, r, s, 0, c - fct);
	  MRIsetVoxVal(__warpmap, c, r, s, 1, r - frt);
	  MRIsetVoxVal(__warpmap, c, r, s, 2, s - fst);
	}
	else if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS ||
                 __dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS)
	{
	  // convert (fct, frt, fst) to dst_RAS
	  dst_CRS->rptr[1][1] = nintfunc(fct);
          dst_CRS->rptr[2][1] = nintfunc(frt);
          dst_CRS->rptr[3][1] = nintfunc(fst);
          dst_CRS->rptr[4][1] = 1;

	  MatrixMultiply(__dstVox2RAS, dst_CRS, dst_RAS);

	  if (__dataformat == WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS)
	  {
	    // in target (atlas) RAS space
	    MRIsetVoxVal(__warpmap, c, r, s, 0, dst_RAS->rptr[1][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 1, dst_RAS->rptr[2][1]);
	    MRIsetVoxVal(__warpmap, c, r, s, 2, dst_RAS->rptr[3][1]);
	  }
	  else // __dataformat == WARPFIELD_DTFMT_DISP_RAS
	  {
	    // convert (c, r, s) to src_RAS
	    src_CRS0->rptr[1][1] = c;
            src_CRS0->rptr[2][1] = r;
            src_CRS0->rptr[3][1] = s;
            src_CRS0->rptr[4][1] = 1;

	    src_RAS0 = MatrixMultiply(__srcVox2RAS, src_CRS0, src_RAS0);

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
  
  return 0;
}


// read 3-frame MRI warp map into __warpmap
GCA_MORPH * Warpfield::read(const char *fname)
{  
  __warpmap = MRIread(fname);
  if (__warpmap == NULL)
  {
    printf("ERROR: Warpfield::read(%s)\n", fname);
    return NULL;
  }

  GCA_MORPH *gcam = GCAMalloc(__warpmap->width, __warpmap->height, __warpmap->depth);
  if (gcam == NULL)
    return NULL;
  
  gcam->image = *__warpmap;
  gcam->atlas = *__warpmap;

  for (int c = 0; c < __warpmap->width; c++)
  {
    for (int r = 0; r < __warpmap->height; r++)
    {
      for (int s = 0; s < __warpmap->depth; s++)
      {
	GCA_MORPH_NODE *gcamn = &gcam->nodes[c][r][s];
	
	gcamn->origx = c;
	gcamn->origy = r;
	gcamn->origz = s;
	
        gcamn->x = MRIgetVoxVal(__warpmap, c, r, s, 0);
	gcamn->y = MRIgetVoxVal(__warpmap, c, r, s, 1);
	gcamn->z = MRIgetVoxVal(__warpmap, c, r, s, 2);
      } // s
    } // r
  } // c
  
  return gcam;
}


// write out the 3-frame MRI warping map
int Warpfield::write(const char *fname)
{
  if (__invert)
    __mgzVersion = ((MGZ_WARPMAP_INV & 0xff ) << 8) | MGH_VERSION;

  //printf("[DEBUG] Warpfield::write(): __mgzVersion = %d\n", __mgzVersion);
  __warpmap->setWarpfieldMeta(__warpmap, __mgzVersion, __dataformat, __srcRAS2Vox);

  int ret = MRIwrite(__warpmap, fname);
  if (ret)
    printf("ERROR: Warpfield::write(%s)\n", fname);
  
  return ret;
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
