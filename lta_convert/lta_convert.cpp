/**
 * @brief A programm to convert linear transform file formats
 *
 */

/*
 * Original Author: Martin Reuter
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
#include "macros.h"
#include "mri.h"
#include "transform.h"
#include "resample.h"
#include "registerio.h"
#include "version.h"

using namespace std;

namespace intypes {
enum InputType { UNKNOWN, LTA, REG, FSL, MNI, NIFTYREG, ITK, VOX };
}

struct Parameters
{
  string transin;
  string ltaout;
  string fslout;
  string mniout;
  string regout;
  string niftyregout;
  string itkout;
  string voxout;
  string src;
  string trg;
  bool   invert;
  int    ltaouttype;
  bool   trgconform;
  bool   trgconform_dc;
  bool   trgconform_min;
  float  trgconform_size;
  string subject;
  intypes::InputType intype;
  bool  srcconform;
  bool  srcconform_dc;
  bool  srcconform_min;
  float srcconform_size;
  bool  srcupsample;
  bool  srcdownsample;
  int   srcupsampleFactor;
  int   srcdownsampleFactor;
  bool  trgupsample;
  bool  trgdownsample;
  int   trgupsampleFactor;
  int   trgdownsampleFactor;
};

static struct Parameters P = {
  "",                  // string transin
  "",                  // string ltaout
  "",                  // string fslout
  "",                  // string mniout
  "",                  // string regout
  "",                  // string niftyregout
  "",                  // string itkout
  "",                  // string voxout
  "",                  // string src
  "",                  // string trg
  false,               // bool   invert
  LINEAR_RAS_TO_RAS,   // int    ltaouttype
  false,               // bool   trgconform
  false,               // bool   trgconform_dc
  false,               // bool   trgconform_min
  1.0,                 // float  trgconform_size
  "",                  // string subject
  intypes::UNKNOWN,    // intypes::InputType intype
  false,               // bool  srcconform
  false,               // bool  srcconform_dc
  false,               // bool  srcconform_min
  1.0,                 // float srcconform_size
  false,               // bool  srcupsample
  false,               // bool  srcdownsample
  1,                   // int   srcupsampleFactor
  1,                   // int   srcdownsampleFactor
  false,               // bool  trgupsample
  false,               // bool  trgdownsample  
  1,                   // int   trgupsampleFactor
  1,                   // int   trgdownsampleFactor
};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);

const char *Progname = NULL;

static int do_sqrt = 0 ;

// this function is modified from mri_convert if (conform_flag)
void conformGeom(VOL_GEOM *vg, bool conform_min, float conform_size0, bool confkeepdc)
{
  // make a copy of source VOL_GEOM, input VOL_GEOM will be updated in place
  VOL_GEOM vg_src = *vg;

  int conform_width = 256;
  float conform_size = conform_size0;
  if (conform_min == TRUE)
    conform_size = MRIfindMinSize(vg, &conform_width);
  else
    conform_width = MRIfindRightSize(vg, conform_size);
  
  // the following codes are modified from MRIconformedTemplate()
  vg->width = vg->height = vg->depth = conform_width;
  vg->xsize  = vg->ysize = vg->zsize = conform_size;
  if (confkeepdc)
  {
    char ostr[4];
    int conform_FoV = conform_width * conform_size;
    MRIdircosToOrientationString(&vg_src, ostr);

    int iLR, iIS, iAP;
    for (iLR = 0; iLR < 3; iLR++)
      if (ostr[iLR] == 'L' || ostr[iLR] == 'R') break;
    for (iIS = 0; iIS < 3; iIS++)
      if (ostr[iIS] == 'I' || ostr[iIS] == 'S') break;
    for (iAP = 0; iAP < 3; iAP++)
      if (ostr[iAP] == 'A' || ostr[iAP] == 'P') break;
    
    printf("keeping DC %d %d %d\n", iLR, iIS, iAP);
    printf("ostr %s, width %d, size %g\n", ostr, conform_width, conform_size);

    int Nvox[3], FoV[3];
    double delta[3];
    
    Nvox[0] = vg_src.width;
    Nvox[1] = vg_src.height;
    Nvox[2] = vg_src.depth;
    delta[0] = vg_src.xsize;
    delta[1] = vg_src.ysize;
    delta[2] = vg_src.zsize;
    
    for (int c = 0; c < 3; c++)
      FoV[c] = Nvox[c] * delta[c];

    // K maps voxels in mri to voxels in mri_template
    MATRIX *K = MatrixAlloc(4, 4, MATRIX_REAL);
    K->rptr[4][4] = 1;

    // If the delta=conform_size, then no interpolation will result
    // Otherwise, there will be interpolation that depends on voxel size
    // Using round() forces no interpolation at the edge of the FoV
    // pad is the number of conformed voxels of padding when Nvox != conform_width
    // set pad this way makes the C_RASs be about the same under general conditions

    double step = delta[iLR] / conform_size;
    double pad = round(((conform_FoV - FoV[iLR]) / 2.0) / conform_size);
    if (ostr[iLR] == 'L') {
      K->rptr[1][iLR + 1] = step;
      K->rptr[1][4] = pad;
    }
    else {
      K->rptr[1][iLR + 1] = -step;
      K->rptr[1][4] = conform_width - pad;
    }

    step = delta[iIS] / conform_size;
    pad = round(((conform_FoV - FoV[iIS]) / 2.0) / conform_size);
    if (ostr[iIS] == 'I') {
      K->rptr[2][iIS + 1] = step;
      K->rptr[2][4] = pad;
    }
    else {
      K->rptr[2][iIS + 1] = -step;
      K->rptr[2][4] = conform_width - pad;
    }

    step = delta[iAP] / conform_size;
    pad = round(((conform_FoV - FoV[iAP]) / 2.0) / conform_size);
    if (ostr[iAP] == 'A') {
      K->rptr[3][iAP + 1] = step;
      K->rptr[3][4] = pad;
    }
    else {
      K->rptr[3][iAP + 1] = -step;
      K->rptr[3][4] = conform_width - pad;
    }

    MATRIX *invK = MatrixInverse(K, NULL);
    MATRIX *Smri = MRIxfmCRS2XYZ(&vg_src, 0);
    MATRIX *Stemp = MatrixMultiplyD(Smri, invK, NULL);
    MRIsetVox2RASFromMatrix(vg, Stemp);

    printf("K ---------------\n");
    MatrixPrint(stdout, K);
    printf("Kinv ---------------\n");
    MatrixPrint(stdout, invK);
    printf("Smri ---------------\n");
    MatrixPrint(stdout, Smri);
    printf("Stemp ---------------\n");
    MatrixPrint(stdout, Stemp);
    printf("----------------------\n");

    MatrixFree(&K);
    MatrixFree(&invK);
    MatrixFree(&Smri);
    MatrixFree(&Stemp);
  }
  else
  {
    // replicates old method exactly
    // these are the same as VOL_GEOM initial values
    vg->x_r = -1.0;
    vg->x_a = 0.0;
    vg->x_s = 0.0;
    vg->y_r = 0.0;
    vg->y_a = 0.0;
    vg->y_s = -1.0;
    vg->z_r = 0.0;
    vg->z_a = 1.0;
    vg->z_s = 0.0;
    // ??? what about c_[ras] ???
  }
}

// this function is modified from MRIupsampleN()
void upsampleGeom(VOL_GEOM *vg, int N)
{
  // Computes CRAS based on location of the 1st voxel in upsampled space
  // The new location is 1/(2*Nth) of a voxel from the corner. The RAS of
  // a voxel is at the center of the voxel (unfortunately), so the corner
  // is located at CRS=[-.5 -.5 -.5]
  MATRIX *Vox2RAS = MRIxfmCRS2XYZ(vg, 0);  // scanner vox2ras of source mri
  MATRIX *CRS0 = MatrixZero(4, 1, NULL);
  CRS0->rptr[1][1] = -0.5 + 1.0 / (2 * N);
  CRS0->rptr[2][1] = -0.5 + 1.0 / (2 * N);
  CRS0->rptr[3][1] = -0.5 + 1.0 / (2 * N);
  CRS0->rptr[4][1] = 1.0;
  
  MATRIX *RAS0 = MatrixMultiply(Vox2RAS, CRS0, NULL);

  // Recompute geometry for finer resolution
  // Only the xsize and cras change
  vg->width  *= N;
  vg->height *= N;
  vg->depth  *= N;
  vg->xsize  /= N;
  vg->ysize  /= N;
  vg->zsize  /= N;
  
  MRIp0ToCRAS(vg, RAS0->rptr[1][1], RAS0->rptr[2][1], RAS0->rptr[3][1]);

  MRIreInitCache(vg);

  MatrixFree(&Vox2RAS);
  MatrixFree(&CRS0);
  MatrixFree(&RAS0);
}

// this function is modified from MRIdownsampleN()
void downsampleGeom(VOL_GEOM *vg, int N)
{
  if (vg->width % N != 0) {
    printf("ERROR: MRIdownsampleN: width=%d, N=%d\n", vg->width, N);
    return;
  }
  if (vg->height % N != 0) {
    printf("ERROR: MRIdownsampleN: height=%d, N=%d\n", vg->height, N);
    return;
  }
  if (vg->depth % N != 0) {
    printf("ERROR: MRIdownsampleN: depth=%d, N=%d\n", vg->depth, N);
    return;
  }

  // Compute P0 for dst
  // CRS0 corresponds to the center of the 1st dst vox
  MATRIX *CRS0 = MatrixZero(4, 1, NULL);
  CRS0->rptr[1][1] = -0.5 + N / 2.0;
  CRS0->rptr[2][1] = -0.5 + N / 2.0;
  CRS0->rptr[3][1] = -0.5 + N / 2.0;
  CRS0->rptr[4][1] = 1.0;
  
  MATRIX *Vox2RAS = MRIxfmCRS2XYZ(vg, 0);
  MATRIX *RAS0 = MatrixMultiply(Vox2RAS, CRS0, NULL);

  // Recompute geometry
  // Only the xsize and cras change
  vg->width  /= N;
  vg->height /= N;
  vg->depth  /= N;
  vg->xsize  *= N;
  vg->ysize  *= N;
  vg->zsize  *= N;

  // Compute New CRAS for dst
  MRIp0ToCRAS(vg, RAS0->rptr[1][1], RAS0->rptr[2][1], RAS0->rptr[3][1]);

  MRIreInitCache(vg);

  MatrixFree(&RAS0);
  MatrixFree(&Vox2RAS);
  MatrixFree(&CRS0);
}

LTA * shallowCopyLTA(const LTA * lta)
{
  LTA * ltatmp = LTAalloc(1,NULL);
  ltatmp->xforms[0].m_L=MatrixCopy(lta->xforms[0].m_L,NULL);
  //copyVolGeom(&lta->xforms[0].src,&ltatmp->xforms[0].src);
  ltatmp->xforms[0].src = lta->xforms[0].src;
  //copyVolGeom(&lta->xforms[0].dst,&ltatmp->xforms[0].dst);
  ltatmp->xforms[0].dst = lta->xforms[0].dst;
  ltatmp->type = lta->type;
  ltatmp->fscale = lta->fscale;
  strcpy(ltatmp->subject, lta->subject); 
  return ltatmp;
}

LTA * readLTA(const string& xfname, const string& sname, const string& tname)
// here sname and tname are not necessary
{
  LTA* lta = LTAread(P.transin.c_str());
  if (lta == NULL)
  {
    cerr << "ERROR readLTA: cannot read " << xfname << endl;
    exit(1);
  }
  
  if (lta->type != LINEAR_RAS_TO_RAS)
    LTAchangeType(lta, LINEAR_RAS_TO_RAS);  
  if (lta->type != LINEAR_RAS_TO_RAS)
  {
    cerr << "ERROR readLTA: cannot change type to RAS_TO_RAS." << endl;
    exit(1);  
  }
  
  // if src and trg mri are passed, change geometry
  if (sname != "")
  {
    MRI * src = MRIreadHeader(sname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
    if (src == NULL)
    {
      cerr << "ERROR readLTA: cannot read src MRI" << sname << endl;
      exit(1);
    }
    //getVolGeom(src, &lta->xforms[0].src);
    lta->xforms[0].src = *src;
    // getVolGeom() set valid = 1;
    lta->xforms[0].src.valid = 1;  // ??? valid and ras_good_flag mean the same thing ???
    MRIfree(&src);
  }
  if (tname != "")
  { 
    MRI * trg = MRIreadHeader(tname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
    if (trg == NULL)
    {
      cerr << "ERROR readFSL: cannot read trg MRI" << tname << endl;
      exit(1);
    }
    //getVolGeom(trg, &lta->xforms[0].dst);
    lta->xforms[0].dst = *trg;
    // getVolGeom() set valid = 1;
    lta->xforms[0].dst.valid = 1;  // ??? valid and ras_good_flag mean the same thing ???
    MRIfree(&trg);
  }
  return lta;
}

LTA * readFSL(const string& xfname, const string& sname, const string& tname)
// use lta transform to readin fslreg
// and then an lta change type from FSLREG_TYPE (I implemented that in transform.c)
{

  LTA * lta = LTAreadExType(xfname.c_str(),FSLREG_TYPE);

  // read src and target mri header
  MRI * src = MRIreadHeader(sname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (src == NULL)
  {
    cerr << "ERROR readFSL: cannot read src MRI" << sname << endl;
    exit(1);
  }
  MRI * trg = MRIreadHeader(tname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (trg == NULL)
  {
    cerr << "ERROR readFSL: cannot read trg MRI" << tname << endl;
    exit(1);
  }
  
  getVolGeom(src, &lta->xforms[0].src);
  getVolGeom(trg, &lta->xforms[0].dst);
  //lta->type = FSLREG_TYPE; // necessary before fix type in transform.c
  lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);

  MRIfree(&src);
  MRIfree(&trg);

  return lta;

}


LTA * readMNI(const string& xfname, const string& sname, const string& tname)
// based on regio_read_mincxfm for reading the matrix
// then the matrix should be RAS2RAS
{


  LTA * lta = LTAreadExType(xfname.c_str(),MNI_TRANSFORM_TYPE);
  if (lta == NULL)
    ErrorExit(ERROR_BADFILE, "can't read input file %s", xfname.c_str());

  // read src and target mri header
  MRI * src = MRIreadHeader(sname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (src == NULL)
  {
    cerr << "ERROR readMNI: cannot read src MRI" << sname << endl;
    exit(1);
  }
  MRI * trg = MRIreadHeader(tname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (trg == NULL)
  {
    cerr << "ERROR readMNI: cannot read trg MRI" << tname << endl;
    exit(1);
  }
  
  // NMI XFM matrix should be identical with RAS2RAS?
  getVolGeom(src, &lta->xforms[0].src);
  getVolGeom(trg, &lta->xforms[0].dst);
  //lta->type = LINEAR_RAS_TO_RAS;

  MRIfree(&src);
  MRIfree(&trg);

  return lta;

}

LTA * readREG(const string& xfname, const string& sname, const string& tname)
//based on regio_read_register for reading and then lta change type from REGISTER_DAT
{

  LTA * lta = LTAreadExType(xfname.c_str(),REGISTER_DAT);

  // read src and target mri header
  MRI * src = MRIreadHeader(sname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (src == NULL)
  {
    cerr << "ERROR readREG: cannot read src MRI" << sname << endl;
    exit(1);
  }
  MRI * trg = MRIreadHeader(tname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (trg == NULL)
  {
    cerr << "ERROR readREG: cannot read trg MRI" << tname << endl;
    exit(1);
  }
  
  getVolGeom(src, &lta->xforms[0].src);
  getVolGeom(trg, &lta->xforms[0].dst);
  //lta->type = REGISTER_DAT;  // necessary before fix type in transform.c
  //// uses MRItkReg2Native internally:
  lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);

  MRIfree(&src);
  MRIfree(&trg);

  return lta;
}

LTA * readNIFTYREG(const string& xfname, const string& sname, const string& tname)
// nifty reg writes the inverse RAS2RAS matrix (trg -> src)
// this functionality needs to be moved to LTA (transform) in the future
{
  LTA* lta = LTAalloc(1, NULL) ;
  LINEAR_TRANSFORM * lt = &lta->xforms[0] ;
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0 ;
  lta->type = LINEAR_RAS_TO_RAS;
  
  std::ifstream transfile (xfname.c_str());
  MATRIX* m_L = MatrixAlloc(4,4,MATRIX_REAL);
  if(transfile.is_open())
  {
    int row=1;
    float v1,v2,v3,v4;
    while(!transfile.eof())
    {
      transfile >> v1 >> v2 >> v3 >> v4;
      *MATRIX_RELT(m_L,row,1) = v1;
      *MATRIX_RELT(m_L,row,2) = v2;
      *MATRIX_RELT(m_L,row,3) = v3;
      *MATRIX_RELT(m_L,row,4) = v4;
      row++;
      if(row>4) break;
    }
    transfile.close();
  }
  else
  {
    cerr << "readNIFTYREG Error opening " << xfname << endl;
    exit(1);
  }
  lt->m_L = MatrixInverse( m_L, lt->m_L );
  MatrixFree(&m_L);
  
  // read src and target mri header
  MRI * src = MRIreadHeader(sname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (src == NULL)
  {
    cerr << "ERROR readNIFTYREG: cannot read src MRI" << sname << endl;
    exit(1);
  }
  MRI * trg = MRIreadHeader(tname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (trg == NULL)
  {
    cerr << "ERROR readNIFTYREG: cannot read trg MRI" << tname << endl;
    exit(1);
  }
  
  getVolGeom(src, &lta->xforms[0].src);
  getVolGeom(trg, &lta->xforms[0].dst);
  
  MRIfree(&src);
  MRIfree(&trg);

  return lta;
}

LTA * readITK(const string& xfname, const string& sname, const string& tname)
// ITK (and Ants) uses Left Posterior Superior coordinates
// and stores inverse (map from ref to mov).
// For some Ants output, call ConvertTransformFile (in Ants) to first
// convert to ITK fromat (e.g. from binary .mat). 
// this functionality needs to be moved to LTA (transform) in the future
{

//#Insight Transform File V1.0
//#Transform 0
//Transform: AffineTransform_double_3_3
//Parameters: 1.0899336847275212 -0.039727996262830335 0.03140529159493116 0.03896690477110579 1.0495356962743074 0.24617959701005143 -0.01840643428186128 -0.27758480101278094 1.059354268159089 -0.9666057029459277 -1.6941071744720673 7.8829725769991175
//FixedParameters: 0 0 0

  // read src and target mri header
  MRI * src = MRIreadHeader(sname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (src == NULL)
  {
    cerr << "ERROR readITK: cannot read src MRI " << sname << endl;
    exit(1);
  }
  MRI * trg = MRIreadHeader(tname.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
  if (trg == NULL)
  {
    cerr << "ERROR readITK: cannot read trg MRI " << tname << endl;
    exit(1);
  }
  
  // determine if 2D (and which dim is flat)
  bool src3d = src->width >1 && src->height > 1 && src->depth > 1;
  bool trg3d = trg->width >1 && trg->height > 1 && trg->depth > 1;
  bool mri3d = src3d || trg3d;
  int flatdim =-1;
  if (! mri3d )
  {
    if (src->width == 1 && trg->width ==1) flatdim = 1;
    else if (src->height == 1 && trg->height ==1) flatdim = 2;
    else if (src->depth == 1 && trg->depth ==1) flatdim = 3;
    else
    {
      cerr << "ERROR readITK: 2D images but plane does not agree. " << tname << endl;
      exit(1);
    }
  }
  int o1 = 0; 
  int o2 = 0;
  if (flatdim == 1) o1=1;
  if (flatdim == 2) o2=1;

  LTA* lta = LTAalloc(1, NULL) ;
  LINEAR_TRANSFORM * lt = &lta->xforms[0] ;
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0 ;
  lta->type = LINEAR_RAS_TO_RAS;
  
  std::ifstream transfile (xfname.c_str());
  MATRIX* m_L = MatrixAlloc(4,4,MATRIX_REAL);
  m_L = MatrixIdentity(4,m_L);
  
  MATRIX* m_R = MatrixAlloc(3,3,MATRIX_REAL);
  VECTOR* v_T = VectorAlloc(3,MATRIX_REAL);
  VECTOR* v_F = VectorAlloc(3,MATRIX_REAL);
  
  std::string str;
  bool is3d = true;
  if(transfile.is_open())
  {
    float v1,v2,v3=0.0;
    while(!transfile.eof())
    {
      while (transfile.peek() == '#')
        getline (transfile,str);
      transfile >> str;
      if (str == "Transform:")
      {
        //getline(transfile,str);
        transfile >> str;
        std::cout << "Transform: "<< str <<std::endl;
//        if (str == "AffineTransform_double_2_2" || str == "MatrixOffsetTransformBase_double_2_2")
        if (str == "AffineTransform_double_2_2" )
          is3d = false;
        else if (str != "AffineTransform_double_3_3" && str != "MatrixOffsetTransformBase_double_3_3")
        {
          std::cout << "ERROR readITK: Transform type unknown!"<< std::endl;
          exit(1);    
        }
        if (!is3d && mri3d)
        {
          std::cout << "ERROR readITK: Transform is 2D but images 3D?"<< std::endl;
          exit(1);    
        }
      }
      else if (str == "Parameters:")
      { 
       if (is3d)
       {
        // convert to ras2ras (from lps2lps)
        // read and mult with diag(-1,-1,1,1) from left and right:
        // mR stores unmodified input, mL is the converted matrix
        
        transfile >> v1 >> v2 >> v3;
        *MATRIX_RELT(m_R,1,1) =  v1;
        *MATRIX_RELT(m_R,1,2) =  v2;
        *MATRIX_RELT(m_R,1,3) =  v3;        
        *MATRIX_RELT(m_L,1,1) =  v1;
        *MATRIX_RELT(m_L,1,2) =  v2;
        *MATRIX_RELT(m_L,1,3) = -v3;
        
        transfile >> v1 >> v2 >> v3;
        *MATRIX_RELT(m_R,2,1) =  v1;
        *MATRIX_RELT(m_R,2,2) =  v2;
        *MATRIX_RELT(m_R,2,3) =  v3;
        *MATRIX_RELT(m_L,2,1) =  v1;
        *MATRIX_RELT(m_L,2,2) =  v2;
        *MATRIX_RELT(m_L,2,3) = -v3;
        
        transfile >> v1 >> v2 >> v3;
        *MATRIX_RELT(m_R,3,1) =  v1;
        *MATRIX_RELT(m_R,3,2) =  v2;
        *MATRIX_RELT(m_R,3,3) =  v3;
        *MATRIX_RELT(m_L,3,1) = -v1;
        *MATRIX_RELT(m_L,3,2) = -v2;
        *MATRIX_RELT(m_L,3,3) =  v3;
        
        // mL is only correct, if no FixedParameter is specified below
        transfile >> v1 >> v2 >> v3;
        VECTOR_ELT(v_T,1)     =  v1;
        VECTOR_ELT(v_T,2)     =  v2;
        VECTOR_ELT(v_T,3)     =  v3;
        *MATRIX_RELT(m_L,1,4) = -v1;
        *MATRIX_RELT(m_L,2,4) = -v2;
        *MATRIX_RELT(m_L,3,4) =  v3;
        
        // already set above (identity)
        //*MATRIX_RELT(m_L,4,1) = 0.0;
        //*MATRIX_RELT(m_L,4,2) = 0.0;
        //*MATRIX_RELT(m_L,4,3) = 0.0;
        //*MATRIX_RELT(m_L,4,4) = 1.0;   
       }
       else
       {
         transfile >> v1 >> v2 ;
        *MATRIX_RELT(m_L,1+o1,1+o1) = v1;
        *MATRIX_RELT(m_L,1+o1,2+o1+o2) = v2;
       
         transfile >> v1 >> v2 ;
        *MATRIX_RELT(m_L,2+o1+o2,1+o1) = v1;
        *MATRIX_RELT(m_L,2+o1+o2,2+o1+o2) = v2;
             
         transfile >> v1 >> v2 ;
        *MATRIX_RELT(m_L,1+o1,4) = -v1;
        *MATRIX_RELT(m_L,2+o1+o2,4) = -v2;
        //*MATRIX_RELT(m_L,4,1) = 0.0;
        //*MATRIX_RELT(m_L,4,2) = 0.0;
        //*MATRIX_RELT(m_L,4,3) = 0.0;
        //*MATRIX_RELT(m_L,4,4) = 1.0;   
       }        
      }
      else if (str == "FixedParameters:")
      {
        if (is3d)
        {
          // we need to update the last column in mL
          // see itkMatrixOffsetTransformBase.hxx  MatrixOffsetTransformBase
        
          transfile >> v1>> v2 >> v3; // center of rotation 
          VECTOR_ELT(v_F,1)   =  v1;
          VECTOR_ELT(v_F,2)   =  v2;
          VECTOR_ELT(v_F,3)   =  v3;
          
          // update mT:
          for (int i = 1; i<=3; i++)
          {
            VECTOR_ELT(v_T,i) += VECTOR_ELT(v_F,i);
            for (int j = 1; j<=3; j++)
            {
              VECTOR_ELT(v_T,i) -= (*MATRIX_RELT(m_R,i,j)) * VECTOR_ELT(v_F,j);
            }
          }
          
          // convert to RAS2RAS and replace last column in mL
          *MATRIX_RELT(m_L,1,4) = - VECTOR_ELT(v_T,1);
          *MATRIX_RELT(m_L,2,4) = - VECTOR_ELT(v_T,2);
          *MATRIX_RELT(m_L,3,4) =   VECTOR_ELT(v_T,3);
          
          
        }  
        else
        {
          transfile >> v1 >> v2;
          
          // this can be implemented similarly to the 3d case (just with 2 dimensions)
          // but I don't have a test case


          if (v1!=0 || v2!=0 )
          {
            std::cout << "ERROR readITK: 2D fixedParameters not equal to zero (not implemented): " << v1 << " " << v2 << std::endl;
            exit(1);
          }
        }
      }
      else
      {
         std::cout <<"ERROR readITK: "<< str << " unknown!" <<std::endl;
         exit(1);
      }
    }
    transfile.close();
  }
  else
  {
    cerr << "ERROR readITK:  opening " << xfname << endl;
    exit(1);
  }
  
  lt->m_L = MatrixInverse( m_L, lt->m_L );
  MatrixFree(&m_L);
  
  getVolGeom(src, &lta->xforms[0].src);
  getVolGeom(trg, &lta->xforms[0].dst);
  
  MRIfree(&src);
  MRIfree(&trg);

  return lta;
}

LTA * readVOX(const string& xfname, const string& sname, const string& tname)
// read transform defined in source voxel space going from target to source
// coordinates (i.e. "inverse transform" as NiftyReg, thinking of images)
{
  // src and target mri header
  MRI* src = MRIreadHeader(sname.c_str(), MRI_VOLUME_TYPE_UNKNOWN);
  if (src == NULL)
  {
    cerr << "ERROR readVOX: cannot read src MRI" << sname << endl;
    exit(1);
  }
  MRI* trg = MRIreadHeader(tname.c_str(), MRI_VOLUME_TYPE_UNKNOWN);
  if (trg == NULL)
  {
    cerr << "ERROR readVOX: cannot read trg MRI" << tname << endl;
    exit(1);
  }

  // LTA
  LTA* lta = LTAalloc(1, NULL);
  LINEAR_TRANSFORM * lt = &lta->xforms[0];
  lt->sigma = 1.0f;
  lt->x0 = lt->y0 = lt->z0 = 0;
  lta->type = LINEAR_RAS_TO_RAS;
  getVolGeom(src, &lta->xforms[0].src);
  getVolGeom(trg, &lta->xforms[0].dst);

  // transform in source voxel space
  std::ifstream transfile (xfname.c_str());
  MATRIX* mat = MatrixAlloc(4, 4, MATRIX_REAL);
  if (transfile.is_open())
  {
    int row = 1;
    float v1, v2, v3, v4;
    while (!transfile.eof())
    {
      transfile >> v1 >> v2 >> v3 >> v4;
      *MATRIX_RELT(mat, row, 1) = v1;
      *MATRIX_RELT(mat, row, 2) = v2;
      *MATRIX_RELT(mat, row, 3) = v3;
      *MATRIX_RELT(mat, row, 4) = v4;
      if (++row > 4) break;
    }
    transfile.close();
  }
  else
  {
    cerr << "readVOX error opening " << xfname << endl;
    exit(1);
  }

  // conversion to RAS-to-RAS
  MATRIX* src_to_ras = MRIgetVoxelToRasXform(src);
  MATRIX* ras_to_src = MatrixInverse(src_to_ras /*source*/, NULL /*target*/);
  mat = MatrixMultiplyD(mat, ras_to_src, mat /*target*/);
  mat = MatrixMultiplyD(src_to_ras, mat, mat /*target*/);

  // after the inverse, RAS-to-RAS transforms from src to trg
  MatrixInverse(mat /*source*/, lt->m_L /*target*/);

  // cleanup
  MatrixFree(&src_to_ras);
  MatrixFree(&ras_to_src);
  MatrixFree(&mat);
  MRIfree(&src);
  MRIfree(&trg);
  return lta;
}

void writeFSL(const string& fname, const LTA * lta)
{
  // shallow copy
  LTA * ltatmp = shallowCopyLTA(lta);
     
  // I implemented this in transform.c instead of here
  if (ltatmp->type != FSLREG_TYPE)
    LTAchangeType(ltatmp, FSLREG_TYPE);

  if(LTAwrite(ltatmp, fname.c_str()) != NO_ERROR)
  {
    cerr << "ERROR writeFSL: cannot create file " << fname << endl;
    exit(1);
  }
  LTAfree(&ltatmp);

  return;
}

void writeMNI(const string& fname, const LTA * lta)
// this is the xfm format
{
  if (lta->type != LINEAR_RAS_TO_RAS)
  {
    cerr << "ERROR: lta should be RAS_TO_RAS by now!!!"<< endl;
    exit(1);  
  }
  // shallow copy
  LTA * ltatmp = shallowCopyLTA(lta);
  
  // to force mni output for a RAS2RAS
  ltatmp->type = MNI_TRANSFORM_TYPE;

  if(LTAwrite(ltatmp, fname.c_str()) != NO_ERROR)
  {
    cerr << "ERROR writeFSL: cannot create file " << fname << endl;
    exit(1);
  }
  LTAfree(&ltatmp);

  return;
}

void writeREG(const string& fname, const LTA * lta)
{
  // shallow copy
  LTA * ltatmp = shallowCopyLTA(lta);
     
  // I implemented this in transform.c instead of here
  if (ltatmp->type != REGISTER_DAT)
    LTAchangeType(ltatmp, REGISTER_DAT);

  if(LTAwrite(ltatmp, fname.c_str()) != NO_ERROR)
  {
    cerr << "ERROR writeREG: cannot create file " << fname << endl;
    exit(1);
  }
  LTAfree(&ltatmp);

  return;
}

void writeNIFTYREG(const string& fname, const LTA * lta)
{
  // shallow copy
  LTA * ltatmp = shallowCopyLTA(lta);
  if (ltatmp->type != LINEAR_RAS_TO_RAS)
    LTAchangeType(ltatmp, LINEAR_RAS_TO_RAS);

  // already right-anterior-superior, invert
  MATRIX* m_L = MatrixAlloc(4,4,MATRIX_REAL);
  m_L = MatrixInverse( ltatmp->xforms[0].m_L , m_L);
  LTAfree(&ltatmp);

  std::ofstream transfile (fname.c_str());
  if(transfile.is_open())
  {
    transfile.precision(7);
    for (int row = 1; row <= 4; row++)
    {
      for (int col = 1; col <= 4; col++)
      {
        transfile << *MATRIX_RELT(m_L,row,col);
        if (col < 4) transfile << " ";
      }
      transfile << std::endl;
    }
    transfile.close();
  }
  else
  {
    cerr << "writeNIFTYREG error opening " << fname << endl;
    exit(1);
  }

  MatrixFree(&m_L);
  return;
}

void writeITK(const string& fname, const LTA * lta)
{
  // shallow copy
  LTA * ltatmp = shallowCopyLTA(lta);
  if (ltatmp->type != LINEAR_RAS_TO_RAS)
    LTAchangeType(ltatmp, LINEAR_RAS_TO_RAS);

  // invert
  MATRIX* m_L = MatrixAlloc(4,4,MATRIX_REAL);
  m_L = MatrixInverse( ltatmp->xforms[0].m_L , m_L);
  
  // convert to left-post-superior
  // multiply from left and right with diag(-1,-1,1,1)
  *MATRIX_RELT(m_L,1,3) = - *MATRIX_RELT(m_L,1,3);
  *MATRIX_RELT(m_L,2,3) = - *MATRIX_RELT(m_L,2,3);
  *MATRIX_RELT(m_L,3,1) = - *MATRIX_RELT(m_L,3,1);
  *MATRIX_RELT(m_L,3,2) = - *MATRIX_RELT(m_L,3,2);
  *MATRIX_RELT(m_L,1,4) = - *MATRIX_RELT(m_L,1,4);
  *MATRIX_RELT(m_L,2,4) = - *MATRIX_RELT(m_L,2,4);
  
  LTAfree(&ltatmp);

  std::ofstream transfile (fname.c_str());
  if(transfile.is_open())
  {
    transfile.precision(17);
    transfile << "#Insight Transform File V1.0" << std::endl;
    transfile << "#Transform 0" << std::endl;
    transfile << "Transform: AffineTransform_double_3_3" << std::endl;
    transfile << "Parameters:";
    for (int row = 1; row < 4; row++)
      for (int col = 1; col < 4; col++)
        transfile << " " << *MATRIX_RELT(m_L,row,col);
    for (int row = 1; row < 4 ; row++)
      transfile << " " << *MATRIX_RELT(m_L,row,4);
    transfile << std::endl;
    transfile << "FixedParameters: 0 0 0" << std::endl;
    transfile.close();
  }
  else
  { 
    cerr << "writeITK Error opening " << fname << endl;
    exit(1);
  }
  MatrixFree(&m_L);


  return;

}

void writeVOX(const string& fname, const LTA * lta)
{
  // shallow copy
  LTA* ltatmp = shallowCopyLTA(lta);
  if (ltatmp->type != LINEAR_RAS_TO_RAS)
    LTAchangeType(ltatmp, LINEAR_RAS_TO_RAS);

  // matrix is RAS-to-RAS from source to target coordinates; invert
  MATRIX* mat = MatrixAlloc(4, 4, MATRIX_REAL);
  mat = MatrixInverse(ltatmp->xforms[0].m_L /*source*/, mat /*target*/);
  LTAfree(&ltatmp);

  // conversion to source voxel space
  MATRIX* src_to_ras = VGgetVoxelToRasXform(&lta->xforms[0].src, NULL, 0);
  MATRIX* ras_to_src = VGgetRasToVoxelXform(&lta->xforms[0].src, NULL, 0);
  mat = MatrixMultiplyD(mat, src_to_ras, mat /*target*/);
  mat = MatrixMultiplyD(ras_to_src, mat, mat /*target*/);

  // output
  std::ofstream transfile (fname.c_str());
  if (transfile.is_open())
  {
    transfile.precision(17);
    for (int row = 1; row <= 4; row++)
    {
      for (int col = 1; col <= 4; col++)
      {
        transfile << *MATRIX_RELT(mat, row, col);
        if (col < 4) transfile << " ";
      }
      transfile << std::endl;
    }
    transfile.close();
  }
  else
  {
    cerr << "writeVOX error opening " << fname << endl;
    exit(1);
  }

  // cleanup
  MatrixFree(&mat);
  MatrixFree(&src_to_ras);
  MatrixFree(&ras_to_src);
  return;
}

int main(int argc, char *argv[])
{
  cout << getVersion() << endl << endl;

  // Default initialization
  int nargs = handleVersionOption(argc, argv, "lta_convert");
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
  
  // Read input transform and convert to RAS2RAS:
  LTA * lta = NULL;  
  if (P.intype==intypes::LTA)
    lta = readLTA(P.transin.c_str(),P.src,P.trg);
  else if (P.intype==intypes::FSL)
    lta = readFSL(P.transin.c_str(),P.src,P.trg);
  else if (P.intype==intypes::MNI)
    lta = readMNI(P.transin.c_str(),P.src,P.trg);
  else if (P.intype==intypes::REG)
    lta = readREG(P.transin.c_str(),P.src,P.trg);
  else if (P.intype==intypes::NIFTYREG)
    lta = readNIFTYREG(P.transin.c_str(),P.src,P.trg);
  else if (P.intype==intypes::ITK)
    lta = readITK(P.transin.c_str(),P.src,P.trg);
  else if (P.intype==intypes::VOX)
    lta = readVOX(P.transin.c_str(),P.src,P.trg);
  if (!lta)
  {
    ErrorExit(ERROR_BADFILE, "%s: can't read input file %s",Progname, P.transin.c_str());
  }
  if (lta->type != LINEAR_RAS_TO_RAS)
  {
    cerr << "ERROR: lta should be RAS_TO_RAS by now!!!"<< endl;
    exit(1);  
  }
  cout << " LTA read, type : " << lta->type << endl;
  MatrixPrint(stdout,lta->xforms[0].m_L);
  
  // conform src
  if (P.srcconform)
    conformGeom(&lta->xforms[0].src, P.srcconform_min, P.srcconform_size, P.srcconform_dc);
  if (P.srcupsample)
    upsampleGeom(&lta->xforms[0].src, P.srcupsampleFactor);
  else if (P.srcdownsample)
    downsampleGeom(&lta->xforms[0].src, P.srcdownsampleFactor);
    
  // conform trg
  if (P.trgconform)
    // initVolGeom(&lta->xforms[0].dst);
    conformGeom(&lta->xforms[0].dst, P.trgconform_min, P.trgconform_size, P.trgconform_dc);
  if (P.trgupsample)
    upsampleGeom(&lta->xforms[0].dst, P.trgupsampleFactor);
  else if (P.trgdownsample)
    downsampleGeom(&lta->xforms[0].dst, P.trgdownsampleFactor);
  
  // invert if desired
  if (P.invert)
  {
    VOL_GEOM vgtmp;
    MATRIX *m_tmp = lta->xforms[0].m_L ;
    lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
    MatrixFree(&m_tmp) ;
    LT *lt = &lta->xforms[0];
    if (lt->dst.valid == 0 || lt->src.valid == 0)
    {
      cerr << "WARNING:********************************************************\n";
      cerr << "WARNING: dst or src volume is invalid.  Inverse likely wrong.\n";
      cerr << "WARNING:********************************************************\n";
    }
    //copyVolGeom(&lt->dst, &vgtmp);
    vgtmp = lt->dst;
    //copyVolGeom(&lt->src, &lt->dst);
    lt->dst = lt->src;
    //copyVolGeom(&vgtmp, &lt->src);
    lt->src = vgtmp;
  }

  if(P.subject.size() > 0){
    printf("setting subject to %s\n",P.subject.c_str());
    strcpy(lta->subject,P.subject.c_str());
  }
  lta->fscale = 0.1;
    
  // write final
  if (P.ltaout!="")
  {
    if (lta->type != P.ltaouttype) // can only be ras2ras (default) or vox2vox here, ??? REGISTER_DAT too (--ltatkreg) ???
    {
      LTAchangeType(lta, P.ltaouttype);
    }

    cout << "Writing  LTA to file "<<P.ltaout.c_str()<<"...\n";
    FILE* fo = fopen(P.ltaout.c_str(),"w");
    if (fo==NULL)
      ErrorExit(ERROR_BADFILE,
                "%s: can't create file %s",Progname, P.ltaout.c_str());
    LTAprint(fo, lta);
    fclose(fo);
  }

  if (do_sqrt)
  {
    MATRIX *m_L = lta->xforms[0].m_L, *m_V, *m_U ;
    VECTOR *v_z ;
    
    v_z = VectorAlloc(m_L->cols, MATRIX_REAL);
    m_V = MatrixAlloc(m_L->rows, m_L->rows, MATRIX_REAL);
    m_U = MatrixSVD(m_L, v_z, m_V);
    MatrixPrint(stdout, m_L) ;
    MatrixPrint(stdout, v_z) ;
    MatrixPrint(stdout, m_V) ;
    MatrixPrint(stdout, m_U) ;
      
  }

  if (P.fslout!="")
  {
    writeFSL(P.fslout,lta);
  }
  if (P.mniout!="")
  {
    writeMNI(P.mniout,lta);
  }
  if (P.regout!="")
  {
    writeREG(P.regout,lta);
  }
  if (P.niftyregout!="")
  {
    writeNIFTYREG(P.niftyregout,lta);
  }
  if (P.itkout!="")
  {
    writeITK(P.itkout,lta);
  }
  if (P.voxout!="")
  {
    writeVOX(P.voxout,lta);
  }
  
  LTAfree(&lta);
  printf("%s successful.\n", Progname);
  return (0);
}

#include "lta_convert.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(lta_convert_help_xml, lta_convert_help_xml_len);
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
  int nargs = 0;
  char *option;

  option = argv[0] + 1;                     // remove '-'
  if (option[0] == '-')
  {
    option = option + 1;  // remove second '-'
  }
  StrUpper(option);

  //cout << " option: " << option << endl;

//  if (!strcmp(option, "IN") )
//  {
//    P.transin = string(argv[1]);
//    nargs = 1;
//    cout << "--in: " << P.transin << " input transform." << endl;
//  }
  if (!strcmp(option, "INLTA") )
  {
    P.transin = string(argv[1]);
    P.intype = intypes::LTA;
    nargs = 1;
    cout << "--inlta: " << P.transin << " input LTA transform." << endl;
  }
  else if (!strcmp(option, "INFSL") )
  {
    P.transin = string(argv[1]);
    P.intype = intypes::FSL;
    nargs = 1;
    cout << "--infsl: " << P.transin << " input FSL transform." << endl;
  }
  else if (!strcmp(option, "INMNI") || !strcmp(option, "INXFM"))
  {
    P.transin = string(argv[1]);
    P.intype = intypes::MNI;
    nargs = 1;
    cout << "--inmni: " << P.transin << " input MNI/XFM transform." << endl;
  }
  else if (!strcmp(option, "INREG"))
  {
    P.transin = string(argv[1]);
    P.intype = intypes::REG;
    nargs = 1;
    cout << "--inreg: " << P.transin << " input TK REG transform." << endl;
  }
  else if (!strcmp(option, "SUBJECT"))
  {
    P.subject = string(argv[1]);
    nargs = 1;
    cout << "--s: " << P.subject << " subject name" << endl;
  }
  else if (!strcmp(option, "SQRT"))
  {
    printf(" !!!!!!! WARNING !!!!!!!!!!\n");
    printf(" !!!! MATRIX SQRT not supported yet !!!!!!!!!!\n");
    do_sqrt = 1 ;
    cout << "--sqrt: True " << endl;
  }
  else if (!strcmp(option, "INNIFTYREG"))
  {
    P.transin = string(argv[1]);
    P.intype = intypes::NIFTYREG;
    nargs = 1;
    cout << "--inniftyreg: " << P.transin << " input NiftyReg transform." << endl;
  }
  else if (!strcmp(option, "INITK"))
  {
    P.transin = string(argv[1]);
    P.intype = intypes::ITK;
    nargs = 1;
    cout << "--initk: " << P.transin << " input ITK txt transform." << endl;
  }
  else if (!strcmp(option, "INVOX"))
  {
    P.transin = string(argv[1]);
    P.intype = intypes::VOX;
    nargs = 1;
    cout << "--invox: " << P.transin << " input source voxel space transform." << endl;
  }
  else if (!strcmp(option, "OUTLTA") )
  {
    P.ltaout = string(argv[1]);
    nargs = 1;
    cout << "--outlta: " << P.ltaout << " output LTA." << endl;
  }
  else if (!strcmp(option, "OUTFSL") )
  {
    P.fslout = string(argv[1]);
    nargs = 1;
    cout << "--outfsl: " << P.fslout << " output FSL matrix." << endl;
  }
  else if (!strcmp(option, "OUTMNI") )
  {
    P.mniout = string(argv[1]);
    nargs = 1;
    cout << "--outmni: " << P.mniout << " output MNI/XFM matrix." << endl;
  }
  else if (!strcmp(option, "OUTREG") )
  {
    P.regout = string(argv[1]);
    nargs = 1;
    cout << "--outreg: " << P.regout << " output reg.dat matrix." << endl;
  }
  else if (!strcmp(option, "OUTNIFTYREG"))
  {
    P.niftyregout = string(argv[1]);
    nargs = 1;
    cout << "--outniftyreg: " << P.niftyregout << " output NiftyReg matrix." << endl;
  }
  else if (!strcmp(option, "OUTITK") )
  {
    P.itkout = string(argv[1]);
    nargs = 1;
    cout << "--outitk: " << P.itkout << " output ITK txt matrix." << endl;
  }
  else if (!strcmp(option, "OUTVOX") )
  {
    P.voxout = string(argv[1]);
    nargs = 1;
    cout << "--outvox: " << P.voxout << " output source voxel space matrix." << endl;
  }
  else if (!strcmp(option, "SRC") )
  {
    P.src = string(argv[1]);
    nargs = 1;
    cout << "--src: " << P.src << " src image (geometry)." << endl;
  }
  else   if (!strcmp(option, "SRCCONFORM"))
  {
    P.srcconform = true;
    cout << "--srcconform:: will conform source geometry." << endl;
  }
  else if (!strcmp(option, "SRCCONFORM-DC"))
  {
    P.srcconform = true;
    P.srcconform_dc = true;
    cout << "--srcconform-dc:: will keep source geometry dc." << endl;
  }
  else if (!strcmp(option, "SRCCONFORM-MIN"))
  {
    P.srcconform = true;
    P.srcconform_min = true;
    cout << "--srcconform-min:: will conform source geometry." << endl;
  }
  else if (!strcmp(option, "SRCCONFORM-SIZE"))
  {
    P.srcconform = true;
    P.srcconform_size = atof(argv[1]);
    cout << "--srcconform-size: " << P.srcconform_size << " input src conform size." << endl;
  }
  else if (!strcmp(option, "SRCUPSAMPLE"))
  {
    P.srcupsample = true;
    P.srcupsampleFactor = atoi(argv[1]);
    nargs = 1;
    cout << "--srcupsample: " << P.srcupsampleFactor << " input src upsample factor." << endl;
  }
  else if (!strcmp(option, "SRCDOWNSAMPLE"))
  {
    P.srcdownsample = true;
    P.srcdownsampleFactor = atoi(argv[1]);
    nargs = 1;
    cout << "--srcdownsample: " << P.srcdownsampleFactor << " input src downsample factor." << endl;
  }
  else if (!strcmp(option, "TRG") )
  {
    P.trg = string(argv[1]);
    nargs = 1;
    cout << "--trg: " << P.trg << " trg image (geometry)." << endl;
  }
  else if (!strcmp(option, "TRGCONFORM") )
  {
    P.trgconform = true;
    cout << "--trgconform: will conform target geometry." << endl;
  }
  else if (!strcmp(option, "TRGCONFORM-DC"))
  {
    P.trgconform = true;
    P.trgconform_dc = true;
    cout << "--trgconform-dc:: will keep target geometry dc." << endl;
  }
  else if (!strcmp(option, "TRGCONFORM-MIN"))
  {
    P.trgconform = true;
    P.trgconform_min = true;
    cout << "--trgconform-min:: will conform target geometry." << endl;
  }
  else if (!strcmp(option, "TRGCONFORM-SIZE"))
  {
    P.trgconform = true;
    P.trgconform_size = atof(argv[1]);
    cout << "--trgconform-size: " << P.trgconform_size << " input trg conform size." << endl;
  }
  else if (!strcmp(option, "TRGUPSAMPLE"))
  {
    P.trgupsample = true;
    P.trgupsampleFactor = atoi(argv[1]);
    nargs = 1;
    cout << "--trgupsample: " << P.trgupsampleFactor << " input trg upsample factor." << endl;
  }
  else if (!strcmp(option, "TRGDOWNSAMPLE"))
  {
    P.trgdownsample = true;
    P.trgdownsampleFactor = atoi(argv[1]);
    nargs = 1;
    cout << "--trgdownsample: " << P.trgdownsampleFactor << " input trg downsample factor." << endl;
  }  
  else if (!strcmp(option, "LTAVOX2VOX") )
  {
    P.ltaouttype = LINEAR_VOX_TO_VOX;
    cout << "--ltavox2vox: output LTA as VOX_TO_VOX transform." << endl;
  }
  else if (!strcmp(option, "LTATKREG") )
  {
    P.ltaouttype = REGISTER_DAT;
    cout << "--ltatkreg: output LTA as REGISTER_DAT transform." << endl;
  }
  else if (!strcmp(option, "INVERT") )
  {
    P.invert = true;
    cout << "--invert: will invert transform." << endl;
  }
  else if (!strcmp(option, "HELP") )
  {
    printUsage();
    exit(1);
  }
//  else if (!strcmp(option, "INTYPE") )
//  {
//    char* it = argv[1];
//    StrUpper(it);
//    nargs = 1;
//    string sit = string(it);
//    if (sit == "LTA") P.intype = intypes::LTA;
//    else if (sit == "MNI") P.intype = intypes::MNI;
//    else if (sit == "XFM") P.intype = intypes::MNI;
//    else if (sit == "FSL") P.intype = intypes::FSL;
//    else if (sit == "REG") P.intype = intypes::REG;
//    else
//    {
//      cout << "WARNING: intype " << sit << "unknown, will try to detect input..." << endl;
//      P.intype = intypes::UNKNOWN;
//    }
//    cout << "--intype: " << sit << endl;
//  }
  else
  {
    cerr << endl << endl << "ERROR: Option: " << argv[0] << " unknown (see --help) !! "
        << endl << endl;
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

/*  bool test1 = (P.mov != "" && P.dst != "" && P.lta != "");
  if (!test1)
  {
    printUsage();
    cerr << endl << endl << "ERROR: Please specify --mov --dst and --lta !  "
        << endl << endl;
    exit(1);
  }
  bool test2 = (P.satit || P.sat > 0 || P.cost != Registration::ROB
      || P.leastsquares);
  if (!test2)
  {
    printUsage();
    cerr << endl << endl
        << "ERROR: Please specify either --satit or --sat <float> !  " << endl
        << endl;
    exit(1);
  }
  bool test3 = (P.iscaleout == "" || P.iscale);
  if (!test3)
  {
    printUsage();
    cerr << endl << endl
        << "ERROR: Please specify --iscale together with --iscaleout to compute and output global intensity scaling! "
        << endl << endl;
    exit(1);
  }
  bool test4 = (P.warpout == "" || (P.warpout != P.weightsout));
  if (!test4)
  {
    printUsage();
    cerr << endl << endl
        << "ERROR: Resampled input name (--mapmov) cannot be same as --weights output!"
        << endl << endl;
    exit(1);
  }

  return (test1 && test2 && test3 && test4);*/
  return true;
}


