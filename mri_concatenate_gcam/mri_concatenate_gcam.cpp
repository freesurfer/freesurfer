#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "error.h"
#include "utils.h"
#include "macros.h"
#include "mri.h"
#include "version.h"
#include "transform.h"
#include "gcamorph.h"

struct Parameters
{
  std::string progName;
  std::string outFile;
  std::string srcImage;
  std::string dstImage;
  std::vector<std::string> fileList;
  bool reduce = false;
  bool invert = false;
  bool downsample = false;
};


static void parseCommand(int argc, char *argv[], Parameters &par);
static void forward(int &argc, char **&argv);
static void printUsage(void);
static TRANSFORM *concat(std::vector<std::string> fileList);


int main(int argc, char *argv[])
{
  vg_isEqual_Threshold = 10e-4; // Override, include/transform.h.
  Parameters par;
  parseCommand(argc, argv, par);
  TRANSFORM *out = concat(par.fileList);
  MRI *mri_src = NULL;
  MRI *mri_dst = NULL;
  
  if (!par.srcImage.empty()) {
    mri_src = MRIreadHeader(par.srcImage.c_str(), MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_src) {
      exit(EXIT_FAILURE);
    }
  }
  
  if (!par.dstImage.empty()) {
    mri_dst = MRIreadHeader(par.dstImage.c_str(), MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_dst) {
      exit(EXIT_FAILURE);
    }
  }
  
  // Before inversion, in case the required source image was moved.
  if ((mri_src || mri_dst) && out->type==MORPH_3D_TYPE) {
    GCAM *gcam = (GCAM *)out->xform;
    out->xform = (void *)GCAMchangeVolGeom(gcam, mri_src, mri_dst);
    GCAMfree(&gcam);
  }
  
  if (par.invert) {
    TransformInvertReplace(out, NULL);
  }
  
  if (par.reduce && out->type!=MORPH_3D_TYPE) {
    LTA *tmp = (LTA *)out->xform;
    out->xform = (void *)LTAreduce(tmp);
    LTAfree(&tmp);
  }
  
  if (par.downsample && out->type==MORPH_3D_TYPE) {
    GCAM *gcam = (GCAM *)out->xform;
    if (gcam->spacing == 1) {
      out->xform = (void *)GCAMdownsample2(gcam);
      GCAMfree(&gcam);
    }
    else
      printf("INFO: spacing is %d, no downsampling needed\n", gcam->spacing);
  }
  
  TransformWrite(out, par.outFile.c_str());
  TransformFree(&out);
  return (NO_ERROR);
}


// Try to be remotely memory efficient (should we have several GCAMs).
TRANSFORM *concat(std::vector<std::string> fileList)
{
  TRANSFORM *out = TransformRead(fileList.back().c_str());
  fileList.pop_back();
  if (!out) {
    exit(EXIT_FAILURE);
  }
  while (!fileList.empty()) {
    const int numTrx = 2;
    TRANSFORM *next = TransformRead(fileList.back().c_str());
    fileList.pop_back();
    if (!next) {
      exit(EXIT_FAILURE);
    }
    TRANSFORM *trxArray[numTrx];
    trxArray[1] = out;
    trxArray[0] = next;
    out = TransformConcat(trxArray, numTrx);
    TransformFree(&trxArray[1]);
    TransformFree(&trxArray[0]);
  }
  return (out);
}


static void parseCommand(int argc, char *argv[], Parameters &par)
{
  par.progName = std::string(argv[0]);
  forward(argc, argv);

  if (argc == 0) {
    printUsage();
    exit(EXIT_FAILURE);
  }

  while (argc > 0) {
    if (!strcmp(*argv, "--reduce") || !strcmp(*argv, "-r")) {
      forward(argc, argv);
      par.reduce = true;
      continue;
    }
    
    if (!strcmp(*argv, "--invert") || !strcmp(*argv, "-i")) {
      forward(argc, argv);
      par.invert = true;
      continue;
    }
    
    if (!strcmp(*argv, "--downsample") || !strcmp(*argv, "-d")) {
      forward(argc, argv);
      par.downsample = true;
      continue;
    }
    
    if (!strcmp(*argv, "--change-source") || !strcmp(*argv, "-s")) {
      forward(argc, argv);
      if (argc==0 || ISOPTION(*argv[0])) {
        ErrorExit(ERROR_BADPARM, "ERROR: no volume for source image geometry");
      }
      par.srcImage = std::string(*argv);
      forward(argc, argv);
      continue;
    }
    
    if (!strcmp(*argv, "--change-target") || !strcmp(*argv, "-t")) {
      forward(argc, argv);
      if (argc==0 || ISOPTION(*argv[0])) {
        ErrorExit(ERROR_BADPARM, "ERROR: no volume for target image geometry");
      }
      par.dstImage = std::string(*argv);
      forward(argc, argv);
      continue;
    }
    
    if (ISOPTION(*argv[0])) {
      ErrorExit(ERROR_BADPARM, "ERROR: unknown option %s", *argv);
    }
    
    par.fileList.push_back(std::string(*argv));
    forward(argc, argv);
  }
  
  if (par.fileList.size() < 2) {
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

