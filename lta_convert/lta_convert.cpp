/**
 * @brief A programm to convert linear transform file formats
 *
 */

/*
 * Original Author: Martin Reuter
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
enum InputType { UNKNOWN, LTA, REG, FSL, MNI, NIFTYREG, ITK };
}

struct Parameters
{
  string transin;
  string ltaout;
  string fslout;
  string mniout;
  string regout;
  string itkout;
  string src;
  string trg;
  bool   invert;
  int    ltaouttype;
  bool   trgconform;
  string subject;
  intypes::InputType intype;
};

static struct Parameters P =
{ "", "", "", "", "" ,"" ,"" ,"" , false , LINEAR_RAS_TO_RAS, false,"", intypes::UNKNOWN};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);

const char *Progname = NULL;

LTA * shallowCopyLTA(const LTA * lta)
{
  LTA * ltatmp = LTAalloc(1,NULL);
  ltatmp->xforms[0].m_L=MatrixCopy(lta->xforms[0].m_L,NULL);
  copyVolGeom(&lta->xforms[0].src,&ltatmp->xforms[0].src);
  copyVolGeom(&lta->xforms[0].dst,&ltatmp->xforms[0].dst);
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
    getVolGeom(src, &lta->xforms[0].src);
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
    getVolGeom(trg, &lta->xforms[0].dst);
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
  
  // conform trg
  if (P.trgconform)
    initVolGeom(&lta->xforms[0].dst);
  
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
    copyVolGeom(&lt->dst, &vgtmp);
    copyVolGeom(&lt->src, &lt->dst);
    copyVolGeom(&vgtmp, &lt->src);  
  }

  if(P.subject.size() > 0){
    printf("setting subject to %s\n",P.subject.c_str());
    strcpy(lta->subject,P.subject.c_str());
  }
  lta->fscale = 0.1;
    
  // write final
  if (P.ltaout!="")
  {
    if (lta->type != P.ltaouttype) // can only be ras2ras (default) or vox2vox here
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
  if (P.itkout!="")
  {
    writeITK(P.itkout,lta);
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
  else if (!strcmp(option, "INNIFTYREG"))
  {
    P.transin = string(argv[1]);
    P.intype = intypes::NIFTYREG;
    nargs = 1;
    cout << "--inniftyreg: " << P.transin << " input Nifty Reg transform." << endl;
  }
  else if (!strcmp(option, "INITK"))
  {
    P.transin = string(argv[1]);
    P.intype = intypes::ITK;
    nargs = 1;
    cout << "--initk: " << P.transin << " input ITK txt transform." << endl;
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
  else if (!strcmp(option, "OUTITK") )
  {
    P.itkout = string(argv[1]);
    nargs = 1;
    cout << "--outitk: " << P.itkout << " output ITK txt matrix." << endl;
  }
  else if (!strcmp(option, "SRC") )
  {
    P.src = string(argv[1]);
    nargs = 1;
    cout << "--src: " << P.src << " src image (geometry)." << endl;
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


