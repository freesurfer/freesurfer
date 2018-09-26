#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifdef __cplusplus
extern "C"
{
#endif

#include "error.h"
#include "utils.h"
#include "macros.h"
#include "mri.h"
#include "version.h"
#include "transform.h"
#include "gcamorph.h"

#ifdef __cplusplus
}
#endif


struct Parameters
{
  std::string progName;
  std::string outFile;
  std::string srcImage;
  std::string dstImage;
  std::vector<std::string> fileList;
  bool reduce = false;
  bool ras = false;
  bool invert = false;
};


// TODO:
// [X] Test concat GCAM w/ GCAM.
// [X] Test GCAM noop.
// [X] Test LTA noop.
// [X] Test LTA inversion.
// [X] Test GCAM inversion by concat GCAM w/ inverse.
// [X] Test LTA concat with inverse.
// [ ] Change to concat_gcam.
// [ ] Write a test.
// [X] Adapt readme.
// [X] Talk to Andrew and change fun names.


static void parseCommand(int argc, char *argv[], Parameters &par);
static void forward(int &argc, char **&argv);
static void printUsage(void);
TRANSFORM *concat(std::vector<std::string> fileList);


int main(int argc, char *argv[])
{
  Parameters par;
  parseCommand(argc, argv, par);
  TRANSFORM *out = concat(par.fileList);
  
  if (!par.srcImage.empty())
  {
    MRI *mri = MRIreadHeader(par.srcImage.c_str(), MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri)
    {
      ErrorExit(NULL, "ERROR: cannot read image geometry");
    }
    TransformSetSrcVolGeomFromMRI(mri, out); // Before GCAM inversion.
  }
  
  if (!par.dstImage.empty())
  {
    MRI *mri = MRIreadHeader(par.dstImage.c_str(), MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri)
    {
      ErrorExit(NULL, "ERROR: cannot read image geometry");
    }
    TransformSetDstVolGeomFromMRI(mri, out);
  }
  
  if (par.invert)
  {
    TransformInvertReplace(out, NULL);
  }
  
  if (par.reduce && out->type != MORPH_3D_TYPE)
  {
    LTA *tmp = (LTA *)out->xform;
    out->xform = (void *)LTAreduce(tmp);
    LTAfree(&tmp);
  }
  
  switch (out->type)
  {
    case MORPH_3D_TYPE:
      if (par.ras) GCAMvoxToRas((GCAM *)out->xform); // Otherwise GCAM_VOX.
      break;
    default:
      int ltaType = par.ras ? LINEAR_RAS_TO_RAS : LINEAR_VOX_TO_VOX;
      LTAchangeType((LTA *)out->xform, ltaType); // NOOP if same type.
      break;
  }
  
  TransformWrite(out, par.outFile.c_str());
  TransformFree(&out);
  return NO_ERROR;
}


// Try to be remotely memory efficient (should we have several GCAMs).
TRANSFORM *concat(std::vector<std::string> fileList)
{
  TRANSFORM *out = TransformRead(fileList.back().c_str());
  fileList.pop_back();
  if (!out)
  {
    exit(EXIT_FAILURE);
  }
  while (!fileList.empty())
  {
    const int numTrx = 2;
    TRANSFORM *next = TransformRead(fileList.back().c_str());
    fileList.pop_back();
    if (!next)
    {
      exit(EXIT_FAILURE);
    }
    TRANSFORM *trxArray[numTrx];
    trxArray[1] = out;
    trxArray[0] = next;
    out = TransformConcat(trxArray, numTrx);
    TransformFree(&trxArray[1]);
    TransformFree(&trxArray[0]);
  }
  return out;
}


static void parseCommand(int argc, char *argv[], Parameters &par)
{
  par.progName = std::string(argv[0]);
  forward(argc, argv);

  if (argc == 0)
  {
    printUsage();
    exit(EXIT_FAILURE);
  }

  while (argc > 0)
  {
    if (!strcmp(*argv, "--reduce") || !strcmp(*argv, "-r"))
    {
      forward(argc, argv);
      par.reduce = true;
      continue;
    }
    
    if (!strcmp(*argv, "--ras") || !strcmp(*argv, "-w"))
    {
      forward(argc, argv);
      par.ras = true;
      continue;
    }
    
    if (!strcmp(*argv, "--invert") || !strcmp(*argv, "-i"))
    {
      forward(argc, argv);
      par.invert = true;
      continue;
    }
    
    if (!strcmp(*argv, "--source") || !strcmp(*argv, "-s"))
    {
      forward(argc, argv);
      if (argc==0 || ISOPTION(*argv[0]))
      {
        ErrorExit(ERROR_BADPARM, "ERROR: no volume for source image geometry.");
      }
      par.srcImage = std::string(*argv);
      forward(argc, argv);
      continue;
    }
    
    if (!strcmp(*argv, "--dest") || !strcmp(*argv, "-d"))
    {
      forward(argc, argv);
      if (argc==0 || ISOPTION(*argv[0]))
      {
        ErrorExit(ERROR_BADPARM, "ERROR: no volume for destination geometry.");
      }
      par.dstImage = std::string(*argv);
      forward(argc, argv);
      continue;
    }
    
    if (ISOPTION(*argv[0]))
    {
      ErrorExit(ERROR_BADPARM, "ERROR: unknown option %s", *argv);
    }
    
    par.fileList.push_back(std::string(*argv));
    forward(argc, argv);
  }
  
  if (par.fileList.size() < 2)
  {
    ErrorExit(ERROR_BADPARM, "ERROR: specify output and at least one input");
  }
  par.outFile = par.fileList.back();
  par.fileList.pop_back(); // Remove last element.
}


#include "mri_concatenate_gcam.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(mri_concatenate_gcam_help_xml,
                mri_concatenate_gcam_help_xml_len);
}


static void forward(int &argc, char **&argv)
{
  argc--;
  argv++;
}

