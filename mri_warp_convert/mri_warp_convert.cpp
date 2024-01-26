/**
 * @brief A program to convert non-linear deformation field file formats
 *
 */

/*
 * Original Author: Oliver Hinds
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <string>
#include <iostream>
#include <fstream>

#include "error.h"
#include "gcamorph.h"
#include "macros.h"
#include "mri.h"
#include "mri_circulars.h"
#include "version.h"
#include "warpfield.h"

using namespace std;

namespace filetypes {
  enum FileType { UNKNOWN, M3Z, FSL, ITK, VOX, RAS, SPM, MGZWARP };
}

struct Parameters
{
  string in_warp;
  string out_warp;
  string in_src_geom;
  string in_warpformat;
  string out_warpformat;
  filetypes::FileType in_type;
  filetypes::FileType out_type;
  bool downsample;
  LTA *lta1;
  LTA *lta2;
};


static struct Parameters P =
  { "", "", "", "abs-crs", "abs-crs", filetypes::UNKNOWN, filetypes::UNKNOWN, false,NULL,NULL};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);

const char *Progname = NULL;

// input is MGH_MORPH
GCAM* readM3Z(const string& warp_file)
// Read an m3z file. Just calls down to GCAMread
{
  GCAM* gcam = GCAMread(warp_file.c_str());
  if (gcam == NULL)
  {
    cerr << "ERROR readM3Z: cannot read " << warp_file << endl;
    exit(1);
  }

  return gcam;
}

// input is MRI_MGH_FILE
GCAM *readMGZWarp(const string& mgzwarpfile)
{
  Warpfield *warpfield = new Warpfield();

  GCAM *gcam = warpfield->read(mgzwarpfile.c_str());
  if (gcam == NULL)
  {
    cerr << "ERROR readMGZWarp: cannot read " << mgzwarpfile << endl;
    exit(1);
  }
  
  return gcam;
}

GCAM* readSPM(const string& warp_file, const string& src_geom)
{
  // This version should properly handle all voxel sizes in the warp and the source
  // See also MRI *MRIapplySpmWarp(MRI *vol, LTA *srclta, MRI *warp, int LRRev, int interp, MRI *out)
  printf("readSPM() as %s\n", (P.in_warpformat.compare("abs-ras") == 0) ? "abs-ras"  : "abs-crs");
  MRI *warp = MRIread(warp_file.c_str()) ;
  if(warp == NULL) ErrorExit(ERROR_NOFILE, "%s: could not read warp volume %s\n", Progname, warp_file.c_str()) ;
  MRI *src = MRIread(src_geom.c_str()) ;
  if(src == NULL) ErrorExit(ERROR_NOFILE, "%s: could not source volume %s\n", Progname, src_geom.c_str()) ;
  
  // it is either abs-ras or abs-crs for --inspm
  if (P.in_warpformat.compare("abs-ras") == 0)
  {
    // this section of codes read SPM as ABS_RAS, and do the RAS2VOX conversion using src_geom
    MATRIX *vox2ras1 = MRIxfmCRS2XYZ(src, 1); // spm crs base=1
    MATRIX *ras2vox1 = MatrixInverse(vox2ras1,NULL);
    MATRIX *spmras = MatrixAlloc(4,1,MATRIX_REAL);
    spmras->rptr[4][1] = 1;
    MATRIX *crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;

    // use 1 to N-1 instead of 1 to N because edge voxels are invalid in the warp
    for(int c=1; c < warp->width-1; c++){
      for(int r=1; r < warp->height-1; r++){
        for(int s=1; s < warp->depth-1; s++){
	  // Get the RAS in the spm input/source space
	  for(int k=0; k<3; k++) spmras->rptr[k+1][1] = MRIgetVoxVal(warp,c,r,s,k);
	  // Get the 1-based CRS in the spm input/source space (usually a conformed space)
	  crs = MatrixMultiplyD(ras2vox1,spmras,crs);
	  // Subtract 1 to make 0-based (could do this in ras2vox1)
	  for(int k=0; k<3; k++) MRIsetVoxVal(warp, c,r,s,k, crs->rptr[k+1][1]-1);
        }
      }
    }
    MatrixFree(&vox2ras1);
    MatrixFree(&ras2vox1);
    MatrixFree(&spmras);
    MatrixFree(&crs);
  }

  //Now copy the warp into a GCAM
  GCA_MORPH* gcam = GCAMalloc(warp->width, warp->height, warp->depth) ;
  GCAMinitVolGeom(gcam, src, warp) ;
  GCAMreadWarpFromMRI(gcam, warp, 0) ; //0 = warp is in absolute source CRS coords

  MRIfree(&warp);
  MRIfree(&src);  

  return(gcam);
}
GCAM* readFSL2(const string& warp_file, const string& src_geom)
{
  // This version should properly handle all voxel sizes in the warp and the source
  printf("readFSL2()\n");
  MRI *warp = MRIread(warp_file.c_str()) ;
  if(warp == NULL) ErrorExit(ERROR_NOFILE, "%s: could not read warp volume %s\n", Progname, warp_file.c_str()) ;
  MRI *src = MRIread(src_geom.c_str()) ;
  if(src == NULL) ErrorExit(ERROR_NOFILE, "%s: could not source volume %s\n", Progname, src_geom.c_str()) ;

  MATRIX *vox2ras = MRIgetVoxelToRasXform(warp) ;
  double det = MatrixDeterminant(vox2ras);
  MatrixFree(&vox2ras);
  printf("det = %g\n",det);
  if(det > 0) printf("non-negative Jacobian determinant -- converting to radiological ordering\n");

  MATRIX *vox2fslras_warp = MRIxfmCRS2XYZfsl(warp);
  MATRIX *vox2fslras_src = MRIxfmCRS2XYZfsl(src);
  MATRIX *fslras2vox_src = MatrixInverse(vox2fslras_src,NULL);
  MatrixFree(&vox2fslras_src);

  // The FSL warp is delta in FSL RAS. Change it to absolute CRS in the source space
  MATRIX *fslras = MatrixAlloc(4,1,MATRIX_REAL);
  fslras->rptr[4][1] = 1;
  MATRIX *crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[4][1] = 1;
  int c=0,r=0,s=0;
  for(c=0; c < warp->width; c++){
    for(r=0; r < warp->height; r++){
      for(s=0; s < warp->depth; s++){
	crs->rptr[1][1] = c;
	crs->rptr[2][1] = r;
	crs->rptr[3][1] = s;
	// Get the FSL RAS coords in the warp volume at this voxel
	fslras = MatrixMultiply(vox2fslras_warp,crs,fslras);
	// Add the warp value at this voxel to get the FSL RAS in the source space
	for(int k=0; k<3; k++) {
	  double w = MRIgetVoxVal(warp, c,r,s,k);
	  if(det > 0) w *= -1; // only flip first frame (by negating relative shifts)
	  fslras->rptr[k+1][1] += w;
	}
	// Now compute vox indices in the source space from the voxel coords
	crs = MatrixMultiply(fslras2vox_src,fslras,crs);
	for(int k=0; k<3; k++) MRIsetVoxVal(warp, c,r,s,k, crs->rptr[k+1][1]);
      }
    }
  }
  MatrixFree(&fslras);
  MatrixFree(&crs);
  MatrixFree(&fslras2vox_src);
  MatrixFree(&vox2fslras_warp);

  //Now copy the warp into a GCAM
  GCA_MORPH* gcam = GCAMalloc(warp->width, warp->height, warp->depth) ;
  GCAMinitVolGeom(gcam, src, warp) ;
  GCAMreadWarpFromMRI(gcam, warp, 0) ; //0 = absolute source CRS coords

  // In the first incarnation of this function, this was run. If it is to be run
  // here, then the function needs to be changed to handle an absolute warp
  // GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri) ;

  MRIfree(&warp);
  MRIfree(&src);

  return(gcam);
}

GCAM* readFSL(const string& warp_file)
// Read in an FSL warp. This is the code that used to reside in
// mri_warp_convert.c.
{
  printf("readFSL()\n");
  MRI* mri = MRIread(warp_file.c_str()) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read warp volume %s\n",
              Progname, warp_file.c_str()) ;

  MATRIX* m = MRIgetVoxelToRasXform(mri) ;

  // NOTE: this assumes a standard siemens image orientation in which
  // case a neurological orientation means that the first frame is
  // flipped

  if ( MatrixDeterminant(m) > 0 )
  {
    fprintf(stdout, "non-negative Jacobian determinant -- converting to radiological ordering\n");
  }
  {
    // 2012/feb/08: tested with anisotropic voxel sizes

    MRI *mri2 = NULL ;
    int c=0,r=0,s=0;
    float v;

    mri2 = MRIcopy(mri,NULL);
    for(c=0; c < mri->width; c++)
    {
      for(r=0; r < mri->height; r++)
      {
        for(s=0; s < mri->depth; s++)
        {
          // only flip first frame (by negating relative shifts)
          v = MRIgetVoxVal(mri, c,r,s,0) / mri->xsize;
          if ( MatrixDeterminant(m) > 0 )
            MRIsetVoxVal(    mri2,c,r,s,0,-v);
          else
            MRIsetVoxVal(    mri2,c,r,s,0, v);

          v = MRIgetVoxVal(mri, c,r,s,1) / mri->ysize;
          MRIsetVoxVal(    mri2,c,r,s,1, v);

          v = MRIgetVoxVal(mri, c,r,s,2) / mri->zsize;
          MRIsetVoxVal(    mri2,c,r,s,2, v);

        }
      }
    }
    MRIfree(&mri);
    mri = mri2;

  }
  MatrixFree(&m) ;


  // this does all the work! (gcamorph.c)
  GCA_MORPH* gcam = GCAMalloc(mri->width, mri->height, mri->depth) ;
  GCAMinitVolGeom(gcam, mri, mri) ;

  // not sure if removing singularities is ever a bad thing. But takes 10min
#if 1
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri) ;
#else
  GCAMreadWarpFromMRI(gcam, mri, 1) ;
#endif

  return gcam;
}

// Read a warp file containing displacements in RAS or LPS space.
// Note src_geom is the geom of the input to the warp space. Most of
// the time, this is just the same as the warp volume.
GCAM* read_world(const string& warp_file, const string& src_geom,
    bool is_lps=false)
{
  // Poor names here: in=warp volume. src defines the geometry of
  // input space of the warp which shares a RAS with the warp. So a
  // given CRS in the input space can be converted to an RAS which can
  // be converted to a CRS in the wapr volume.  Most of the time this
  // is just the same as the warp. It does make it confusing.
  MRI* in = MRIread( warp_file.c_str() );
  if (in == NULL) {
    cerr << "ERROR: couldn't read input warp from " << warp_file << endl;
    return NULL;
  }
  MRI* src = MRIread( src_geom.c_str() );
  if (src == NULL) {
	cerr << "ERROR: couldn't read source/atlas geometry from " << src_geom << endl;
    return NULL;
  }

  GCA_MORPH* out = GCAMalloc(in->width, in->height, in->depth) ;
  //int GCAMinitVolGeom(GCAM *gcam, MRI *mri_image, MRI *mri_atlas)
  GCAMinitVolGeom(out, src, in) ;
  out->type = GCAM_VOX;

  MATRIX* dst_vox2mm = MRIgetVoxelToRasXform(in);
  MATRIX* src_vox2mm = MRIgetVoxelToRasXform(src);
  if (is_lps) {
      MATRIX *ras2lps = MatrixIdentity(4, NULL);
      ras2lps->rptr[1][1] = -1;
      ras2lps->rptr[2][2] = -1;
      dst_vox2mm = MatrixMultiplyD(ras2lps, dst_vox2mm, dst_vox2mm);
      src_vox2mm = MatrixMultiplyD(ras2lps, src_vox2mm, src_vox2mm);
      MatrixFree(&ras2lps);
  }
  MATRIX* src_mm2vox = MatrixInverse(src_vox2mm, NULL);

  VECTOR* src_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR* src_mm = VectorAlloc(4, MATRIX_REAL);
  VECTOR* dst_mm = VectorAlloc(4, MATRIX_REAL);
  VECTOR* wrp_mm = VectorAlloc(4, MATRIX_REAL);

  VECTOR_ELT(wrp_mm, 4) = 0;
  VECTOR_ELT(dst_vox, 4) = 1;
  for(int s=0; s < in->depth; s++) {
    for(int c=0; c < in->width; c++) {
      for(int r=0; r < in->height; r++) {
        GCA_MORPH_NODE* node = &out->nodes[c][r][s];
        node->origx = c;
        node->origy = r;
        node->origz = s;
        node->xn = c;
        node->yn = r;
        node->zn = s;

        VECTOR3_LOAD(dst_vox, c, r, s);
        dst_mm = MatrixMultiplyD(dst_vox2mm, dst_vox, dst_mm);

        VECTOR3_LOAD(wrp_mm, MRIgetVoxVal(in, c, r, s, 0),
            MRIgetVoxVal(in, c, r, s, 1),
            MRIgetVoxVal(in, c, r, s, 2));
        src_mm = VectorAdd(dst_mm, wrp_mm, src_mm);
        src_vox = MatrixMultiplyD(src_mm2vox, src_mm, src_vox);

        node->x = VECTOR_ELT(src_vox, 1);
        node->y = VECTOR_ELT(src_vox, 2);
        node->z = VECTOR_ELT(src_vox, 3);
      }
    }
  }

  MRIfree(&in);
  MRIfree(&src);
  VectorFree(&src_vox);
  VectorFree(&src_mm);
  VectorFree(&wrp_mm);
  VectorFree(&dst_mm);
  VectorFree(&dst_vox);
  MatrixFree(&dst_vox2mm);
  MatrixFree(&src_vox2mm);
  MatrixFree(&src_mm2vox);
  return out;
}

// Read a warp file as displacements in source-voxel space.
GCAM* read_voxel(const string& warp_file, const string& src_geom)
{
  MRI* in = MRIread( warp_file.c_str() );
  if (in == NULL) {
    cerr << "ERROR: couldn't read input warp from " << warp_file << endl;
    return NULL;
  }
  MRI* src = MRIread( src_geom.c_str() );
  if (src == NULL) {
	cerr << "ERROR: couldn't read atlas/source geometry from " << src_geom << endl;
    return NULL;
  }

  GCA_MORPH* out = GCAMalloc(in->width, in->height, in->depth) ;
  GCAMinitVolGeom(out, src, in) ;
  out->type = GCAM_VOX;

  MATRIX* dst_vox2ras = MRIgetVoxelToRasXform(in);
  MATRIX* src_vox2ras = MRIgetVoxelToRasXform(src);
  MATRIX* src_ras2vox = MatrixInverse(src_vox2ras, NULL);
  MATRIX* dst2src_vox = MatrixMultiplyD(src_ras2vox, dst_vox2ras, NULL);
  MATRIX* src_vox = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(dst_vox, 4) = 1;

  for(int s=0; s < in->depth; s++) {
    for(int c=0; c < in->width; c++) {
      for(int r=0; r < in->height; r++) {
        GCA_MORPH_NODE* node = &out->nodes[c][r][s];
        node->origx = c;
        node->origy = r;
        node->origz = s;
        node->xn = c;
        node->yn = r;
        node->zn = s;

        VECTOR3_LOAD(dst_vox, c, r, s);
        src_vox = MatrixMultiplyD(dst2src_vox, dst_vox, src_vox);
        node->x = MRIgetVoxVal(in, c, r, s, 0) + VECTOR_ELT(src_vox, 1);
        node->y = MRIgetVoxVal(in, c, r, s, 1) + VECTOR_ELT(src_vox, 2);
        node->z = MRIgetVoxVal(in, c, r, s, 2) + VECTOR_ELT(src_vox, 3);
      }
    }
  }

  MRIfree(&in);
  MRIfree(&src);
  VectorFree(&src_vox);
  VectorFree(&dst_vox);
  MatrixFree(&dst_vox2ras);
  MatrixFree(&src_vox2ras);
  MatrixFree(&src_ras2vox);
  MatrixFree(&dst2src_vox);
  return out;
}

void writeM3Z(const string& fname, GCAM *gcam, bool downsample=false)
// Write an m3z file. Just calls down to GCAMwrite
{
  GCA_MORPH* out = downsample ? GCAMdownsample2(gcam) : gcam;
  GCAMwrite(out, fname.c_str());
  if (downsample) {
      GCAMfree(&out);
  }
}

void writeMGZWarp(const string& fname, GCAM *gcam)
{
  int dataformat =  WarpfieldDTFMT::WARPFIELD_DTFMT_UNKNOWN;
  if (P.out_warpformat.compare("abs-crs") == 0)
    dataformat = WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_CRS;
  else if (P.out_warpformat.compare("disp-crs") == 0)
    dataformat = WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_CRS;
  else if (P.out_warpformat.compare("abs-ras") == 0)
    dataformat = WarpfieldDTFMT::WARPFIELD_DTFMT_ABS_RAS;
  else if (P.out_warpformat.compare("disp-ras") == 0)
    dataformat = WarpfieldDTFMT::WARPFIELD_DTFMT_DISP_RAS;

  Warpfield *warpfield = new Warpfield();
  warpfield->convert(gcam, dataformat);
  warpfield->write(fname.c_str());
}

void writeFSL(const string& fname, const GCAM *gcam)
// Write an FSL warp file.
// NOT IMPLEMENTED
{
  cerr << "ERROR writeFSL is not implemented, sorry!" << endl;
  exit(1);
}

// Write warp as displacements in RAS or LPS space.
// delta = src_RAS - dst_RAS
// output warp is in ITK_MORPH format
void write_world(const string& fname, GCAM* gcam, bool is_lps=false)
{
  MATRIX* dst_vox2mm = VGgetVoxelToRasXform(&gcam->atlas, NULL, 0);
  MATRIX* src_vox2mm = VGgetVoxelToRasXform(&gcam->image, NULL, 0);

  MRI* out = MRIallocSequence(gcam->atlas.width, gcam->atlas.height,
      gcam->atlas.depth, MRI_FLOAT, 3);
  MRIsetResolution(out, gcam->atlas.xsize, gcam->atlas.ysize,
      gcam->atlas.zsize);
  MRIsetVox2RASFromMatrix(out, dst_vox2mm);
  MRIcopyVolGeomToMRI(out, &gcam->atlas);

  if (is_lps) {
      MATRIX *ras2lps = MatrixIdentity(4, NULL);
      ras2lps->rptr[1][1] = -1;
      ras2lps->rptr[2][2] = -1;
      dst_vox2mm = MatrixMultiplyD(ras2lps, dst_vox2mm, dst_vox2mm);
      src_vox2mm = MatrixMultiplyD(ras2lps, src_vox2mm, src_vox2mm);
      MatrixFree(&ras2lps);
  }

  int x, y, z;
  float xw, yw, zw;
  MATRIX* src_vox = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(src_vox, 4) = 1;
  VECTOR_ELT(dst_vox, 4) = 1;
  MATRIX* src_mm = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_mm = VectorAlloc(4, MATRIX_REAL);
  const bool is_same_size = out->width==gcam->width && out->height==gcam->height
      && out->depth==gcam->depth;
  for (x = 0; x < out->width; x++) {
    for (y = 0; y < out->height; y++) {
      for (z = 0; z < out->depth; z++) {
        if (is_same_size) {
            GCA_MORPH_NODE* node = &gcam->nodes[x][y][z];
            xw = node->x;
            yw = node->y;
            zw = node->z;
        }
        else {
            GCAMsampleMorph(gcam, x, y, z, &xw, &yw, &zw);
        }
        VECTOR3_LOAD(dst_vox, x, y, z);
        MatrixMultiplyD(dst_vox2mm, dst_vox, dst_mm);
        VECTOR3_LOAD(src_vox, xw, yw, zw);
        MatrixMultiplyD(src_vox2mm, src_vox, src_mm);

        MRIsetVoxVal(out, x, y, z, 0, VECTOR_ELT(src_mm,1)-VECTOR_ELT(dst_mm,1));
        MRIsetVoxVal(out, x, y, z, 1, VECTOR_ELT(src_mm,2)-VECTOR_ELT(dst_mm,2));
        MRIsetVoxVal(out, x, y, z, 2, VECTOR_ELT(src_mm,3)-VECTOR_ELT(dst_mm,3));
      }
    }
  }

  if (MRIwriteType(out, fname.c_str(), ITK_MORPH) != 0) {
    cerr << "Error writing warp to " << fname << endl;
  }
  MRIfree(&out);
  MatrixFree(&dst_vox2mm);
  MatrixFree(&src_vox2mm);
  MatrixFree(&src_vox);
  MatrixFree(&dst_vox);
  MatrixFree(&src_mm);
  MatrixFree(&dst_mm);
}

// Write a warp file as displacements in source-voxel space.
// ??? is this right to calculate src_CRS from dst_CRS as following:
//     MatrixMultiplyD(dst2src_vox, dst_vox, src_vox); ???
// ??? GCAMsampleMorph() (node->[x,y,z]) returns src_CRS for given dst_CRS ???
void write_voxel(const string& fname, GCAM* gcam)
{
  MATRIX* dst_vox2ras = VGgetVoxelToRasXform(&gcam->atlas, NULL, 0);
  MATRIX* src_vox2ras = VGgetVoxelToRasXform(&gcam->image, NULL, 0);
  MATRIX* src_ras2vox = MatrixInverse(src_vox2ras, NULL);
  MATRIX* dst2src_vox = MatrixMultiplyD(src_ras2vox, dst_vox2ras, NULL);

  MRI* out = MRIallocSequence(gcam->atlas.width, gcam->atlas.height,
      gcam->atlas.depth, MRI_FLOAT, 3);
  MRIsetResolution(out, gcam->atlas.xsize, gcam->atlas.ysize,
      gcam->atlas.zsize);
  MRIsetVox2RASFromMatrix(out, dst_vox2ras);
  MRIcopyVolGeomToMRI(out, &gcam->atlas);

  MATRIX* src_vox = VectorAlloc(4, MATRIX_REAL);
  MATRIX* dst_vox = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(dst_vox, 4) = 1;
  const bool is_same_size = out->width==gcam->width && out->height==gcam->height
      && out->depth==gcam->depth;
  float x, y, z;
  for (int c = 0; c < out->width; c++) {
    for (int r = 0; r < out->height; r++) {
      for (int s = 0; s < out->depth; s++) {
        if (is_same_size) {
			GCA_MORPH_NODE* node = &gcam->nodes[c][r][s];
			x = node->x;
			y = node->y;
			z = node->z;
		}
		else {
		  GCAMsampleMorph(gcam, c, r, s, &x, &y, &z);
		}
        VECTOR3_LOAD(dst_vox, c, r, s);
        MatrixMultiplyD(dst2src_vox, dst_vox, src_vox);
        MRIsetVoxVal(out, c, r, s, 0, x - VECTOR_ELT(src_vox, 1));
        MRIsetVoxVal(out, c, r, s, 1, y - VECTOR_ELT(src_vox, 2));
        MRIsetVoxVal(out, c, r, s, 2, z - VECTOR_ELT(src_vox, 3));
      }
    }
  }
  if (MRIwrite(out, fname.c_str()) != 0)
  {
    cerr << "Error writing VOX warp to " << fname << endl;
  }
  MRIfree(&out);
  MatrixFree(&dst_vox2ras);
  MatrixFree(&src_vox2ras);
  MatrixFree(&src_ras2vox);
  MatrixFree(&dst2src_vox);
}

int main(int argc, char *argv[])
{
  cout << getVersion() << endl << endl;

  // Default initialization
  int nargs = handleVersionOption(argc, argv, "mri_warp_convert");
  if (nargs && argc - nargs == 1)
  {
    exit(0);
  }
  argc -= nargs;
  Progname = argv[0];
  argc--;
  argv++;
  ErrorInit(NULL, NULL, NULL);

  // Parse command line
  if (!parseCommandLine(argc, argv, P))
  {
    //printUsage();
    exit(1);
  }

  GCA_MORPH* gcam = NULL;
  bool is_lps = false;
  switch (P.in_type) {
    case filetypes::M3Z:
      gcam = readM3Z(P.in_warp.c_str());
      break;
    case filetypes::MGZWARP:
      gcam = readMGZWarp(P.in_warp.c_str());
      break;      
    case filetypes::FSL:
      if(P.in_src_geom.empty()) gcam = readFSL(P.in_warp.c_str());
      else                      gcam = readFSL2(P.in_warp.c_str(),P.in_src_geom);
      break;
    case filetypes::SPM:
      if (P.in_warpformat.compare("abs-ras") != 0 && P.in_warpformat.compare("abs-crs") != 0) {
	printf("ERROR: --inspm only handles abs-ras or abs-crs\n");
	exit(1);
      }
      if(P.in_src_geom.empty()){
	printf("ERROR: --inspm needs source/atlas geometry, use --insrcgeom or -g to specify\n");
	exit(1);
      }
      gcam = readSPM(P.in_warp.c_str(),P.in_src_geom);
      break;
    case filetypes::ITK:
      is_lps = true;
      gcam = read_world(P.in_warp.c_str(), P.in_src_geom, is_lps);
      break;
    case filetypes::VOX:
      gcam = read_voxel(P.in_warp.c_str(), P.in_src_geom);
      break;
    case filetypes::RAS:
      is_lps = false;
      gcam = read_world(P.in_warp.c_str(), P.in_src_geom, is_lps);
      break;
    default:
      ErrorExit(ERROR_BADFILE, "%s: Unknown input type for %s",
                Progname, P.in_warp.c_str());
  }

  if (!gcam)
  {
    ErrorExit(ERROR_BADFILE, "%s: can't read input file %s",
              Progname, P.in_warp.c_str());
  }

  if(P.lta1 || P.lta2){
    // Create composite morph for warping a 
    // source image -> LTA1 -> GCAM -> LTA2 -> atlas/destination image
    // The LTAs can be any type as they will be converted to VOX2VOX inside concat3
    printf("Applying LTAs to the GCAM\n");
    GCA_MORPH *gcam2 = GCAMconcat3(P.lta1, gcam, P.lta2, NULL);
    if(!gcam2) exit(1);
    GCAMfree(&gcam);
    gcam = gcam2;
  }

  switch (P.out_type) {
    case filetypes::M3Z:
      writeM3Z(P.out_warp.c_str(), gcam, P.downsample);
      break;
    case filetypes::MGZWARP:
      writeMGZWarp(P.out_warp.c_str(), gcam);
      break;      
    case filetypes::FSL:
      writeFSL(P.out_warp.c_str(), gcam);
      break;
    case filetypes::ITK:
      is_lps = true;
      write_world(P.out_warp.c_str(), gcam, is_lps);
      break;
    case filetypes::VOX:
      write_voxel(P.out_warp.c_str(), gcam);
      break;
    case filetypes::RAS:
      is_lps = false;
      write_world(P.out_warp.c_str(), gcam, is_lps);
      break;
    default:
      ErrorExit(ERROR_BADFILE, "%s: Unknown output type for %s",
                Progname, P.out_warp.c_str());
  }

  GCAMfree(&gcam);
  printf("%s successful.\n", Progname);
  return (0);
}

#include "mri_warp_convert.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(mri_warp_convert_help_xml, mri_warp_convert_help_xml_len);
}

/*!
 \fn int parseNextCommand(int argc, char **argv)
 \brief Parses the command-line for next command
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       number of used arguments for this command
 */
static int parseNextCommand(int argc, char *argv[], Parameters & P)
{
  bool have_input = false;
  bool have_output = false;

  int nargs = 0;
  char *option;

  option = argv[0] + 1;                     // remove '-'
  if (option[0] == '-')
  {
    option = option + 1;  // remove second '-'
  }
  StrUpper(option);

  if (!strcmp(option, "INM3Z") )
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::M3Z;
    nargs = 1;
    cout << "--inm3z: " << P.in_warp << " input M3Z warp." << endl;
  }
  else if (!strcmp(option, "INMGZWARP"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::MGZWARP;
    nargs = 1;
    cout << "--inmgzwarp: " << P.in_warp << " input MGZWARP." << endl;
  }  
  else if (!strcmp(option, "INFSL"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::FSL;
    nargs = 1;
    cout << "--infsl: " << P.in_warp << " input FSL warp." << endl;
  }
  else if (!strcmp(option, "INSPM"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::SPM;
    nargs = 1;
    cout << "--inspm: " << P.in_warp << " input SPM warp." << endl;
  }
  else if (!strcmp(option, "INITK") || !strcmp(option, "INLPS"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::ITK;
    nargs = 1;
    cout << "--inlps: " << P.in_warp << " input LPS warp." << endl;
  }
  else if (!strcmp(option, "INVOX"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::VOX;
    nargs = 1;
    cout << "--invox: " << P.in_warp << " input VOX warp." << endl;
  }
  else if (!strcmp(option, "INRAS"))
  {
    if (have_input) {
      cerr << endl << endl << "ERROR: Only one input warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_input = true;

    P.in_warp = string(argv[1]);
    P.in_type = filetypes::RAS;
    nargs = 1;
    cout << "--inras: " << P.in_warp << " input RAS warp." << endl;
  }
  else if (!strcmp(option, "OUTMGZWARP"))
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::MGZWARP;
    nargs = 1;
    cout << "--outmgz: " << P.out_warp << " output MGZWARP." << endl;    
  }
  else if (!strcmp(option, "OUTM3Z") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::M3Z;
    nargs = 1;
    cout << "--outm3z: " << P.out_warp << " output M3Z." << endl;
  }
  else if (!strcmp(option, "OUTFSL") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::FSL;
    nargs = 1;
    cout << "--outfsl: " << P.out_warp << " output FSL warp." << endl;
  }
  else if (!strcmp(option, "OUTITK") || !strcmp(option, "OUTLPS"))
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::ITK;
    nargs = 1;
    cout << "--outlps: " << P.out_warp << " output LPS warp." << endl;
  }
  else if (!strcmp(option, "OUTVOX") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::VOX;
    nargs = 1;
    cout << "--outvox: " << P.out_warp << " output VOX warp." << endl;
  }
  else if (!strcmp(option, "OUTRAS") )
  {
    if (have_output) {
      cerr << endl << endl << "ERROR: Only one output warp can be specified"
           << endl << endl;
      printUsage();
      exit(1);
    }
    have_output = true;

    P.out_warp = string(argv[1]);
    P.out_type = filetypes::RAS;
    nargs = 1;
    cout << "--outras: " << P.out_warp << " output RAS warp." << endl;
  }
  else if (!strcmp(option, "INSRCGEOM") || !strcmp(option, "G"))
  {
    // Note src_geom is the geom of the input to the warp space, ie,
    // it shares a RAS with the warp volume. Most of the time, this is
    // just the same as the warp volume.
    P.in_src_geom = string(argv[1]);
    nargs = 1;
    cout << "--insrcgeom: " << P.in_src_geom
         << " atlas/source image (used for geometry)." << endl;
  }
  else if (!strcmp(option, "LTA1")) {
    P.lta1 = LTAread(argv[1]);
    if(!P.lta1) exit(1);
    nargs = 1;
  }
  else if (!strcmp(option, "LTA1-INV")) {
    P.lta1 = LTAread(argv[1]);
    if(!P.lta1) exit(1);
    P.lta1 = LTAinvert(P.lta1,P.lta1);
    nargs = 1;
  }
  else if (!strcmp(option, "LTA2"))  {
    P.lta2 = LTAread(argv[1]);
    if(!P.lta2) exit(1);
    nargs = 1;
  }
  else if (!strcmp(option, "LTA2-INV")) {
    P.lta2 = LTAread(argv[1]);
    if(!P.lta2) exit(1);
    P.lta2 = LTAinvert(P.lta2,P.lta2);
    nargs = 1;
  }
  else if (!strcmp(option, "DOWNSAMPLE") || !strcmp(option, "D"))
  {
    if (!P.downsample)
        cout << "--downsample: save M3Z at half resolution." << endl;
    P.downsample = true;
    nargs = 0;
  }
  else if (!strcmp(option, "INWARPFORMAT"))
  {
    P.in_warpformat = string(argv[1]);
    nargs = 1;
    cout << "--inwarpformat: " << P.in_warpformat
         << " (specify input warp data format: abs-crs (default), disp-crs, abs-ras, or disp-ras)" << endl;
  }
  else if (!strcmp(option, "OUTWARPFORMAT"))
  {
    P.out_warpformat = string(argv[1]);
    nargs = 1;
    cout << "--outwarpformat: " << P.out_warpformat
         << " (specify output warp data format: abs-crs (default), disp-crs, abs-ras, or disp-ras)" << endl;
  } 
  else if (!strcmp(option, "HELP") )
  {
    printUsage();
    exit(1);
  }
  else
  {
    cerr << endl << endl << "ERROR: Option: " << argv[0]
         << " unknown (see --help) !! " << endl << endl;
    exit(1);
  }

  fflush(stdout);

  return (nargs);
}

/*!
 \fn int parseCommandLine(int argc, char **argv)
 \brief Parses the command-line
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       if all necessary parameters were set
 */
static bool parseCommandLine(int argc, char *argv[], Parameters & P)
{
  int nargs;
  int inputargs = argc;
  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv, P);
    argc -= nargs;
    argv += nargs;
  }

  if (inputargs == 0)
  {
    printUsage();
    exit(1);
  }

  bool need_geom;
  switch (P.in_type) {
    case filetypes::ITK:
    case filetypes::RAS:
    case filetypes::VOX:
        need_geom = true;
        break;
    default:
        need_geom = false;
  }
  if (P.in_src_geom.empty() && need_geom) {
    cerr << endl << endl << "ERROR: specified input warp requires atlas --insrcgeom"
         << endl << endl;
    return false;
  }

  if (P.out_type != filetypes::M3Z && P.downsample) {
    cerr << endl << endl;
    cerr << "ERROR: --downsample flag only valid for output type M3Z"
         << endl << endl;
    return false;
  }

  return true;
}
