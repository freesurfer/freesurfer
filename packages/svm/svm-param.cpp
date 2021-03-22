/*
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


#include <string.h>
#include "svm-io-format.h"
#include "svm-param.h"

const double SvmParam::dC = 100;
const Kernel SvmParam::dKernel = Kernel("1");
const int    SvmParam::dVerbose = 1;
const double SvmParam::dAlphaEpsilon = -1;
const double SvmParam::dClassEpsilon = 1e-6;
const int    SvmParam::dMaxIterations = 500;
const double SvmParam::dSigDig = 12;
const double SvmParam::dOptEpsilon = 1e-4;


const char* SvmParam::nameC = "-C";
const char* SvmParam::nameKernel = "-k";
const char* SvmParam::nameVerbose = "-v";
const char* SvmParam::nameAlphaEpsilon = "-alpha_epsilon";
const char* SvmParam::nameClassEpsilon = "-class_epsilon";
const char* SvmParam::nameMaxIterations = "-max_iterations";
const char* SvmParam::nameSigDig = "-sig_digits";
const char* SvmParam::nameOptEpsilon = "-opt_epsilon";

using namespace std;

SvmParam::SvmParam() {
  C = dC;
  kernel = dKernel;
  verbose = dVerbose;
  alphaEpsilon = dAlphaEpsilon;
  classEpsilon = dClassEpsilon;
  maxIterations = dMaxIterations;
  sigDig = dSigDig;
  optEpsilon = dOptEpsilon;
}


int find_param(const char* param, const char* const* argv, int argc) {
  for (int i = 0; i < argc; i++ )
    if (!strcmp(param,argv[i]) )
      return i;
  return -1;
}

template <class T>
bool scan_param (T& param, const char* valStr, const char *format = "%d " ) {
  return ( sscanf(valStr,format,&param) == 1 );

}

template <class T>
bool scan_single_param (T& param, const char *paramString,
                        const char* const* argv, int argc) {
  int index;
  if ( (index = find_param(paramString,argv,argc)) >= 0 )
    return scan_param(param,argv[index+1],readFormat(param));
  return false;
}




void SvmParam::parse(const char* const* argv, int argc) {
  // Kernel is a special case, it has two parameters.
  int index = find_param(nameKernel,argv,argc);
  if ( index >= 0 ) {
    char *kernelStr = new char[strlen(argv[index+1])+strlen(argv[index+2])+3];
    sprintf(kernelStr,"%s %s", argv[index+1],argv[index+2]);
    kernel.parse(kernelStr);
    delete [] kernelStr;
  }



  // Parse all the single valued parameters
  scan_single_param(C,nameC,argv,argc);

  scan_single_param(verbose,nameVerbose,argv,argc);
  scan_single_param(alphaEpsilon,nameAlphaEpsilon,argv,argc);
  scan_single_param(classEpsilon,nameClassEpsilon,argv,argc);

  scan_single_param(maxIterations,nameMaxIterations,argv,argc);
  scan_single_param(sigDig,nameSigDig,argv,argc);
  scan_single_param(optEpsilon,nameOptEpsilon,argv,argc);
}




ostream& SvmParam::printUsage(ostream& s) {
  s << "\nSVM options:\n"
  << " " << nameKernel << " kernel_type [kernel_param]"
  << " : kernel used in the training step\n"
  << "\t kernel_type: 1-linear, 2 - rbf, 3 - polynomial (default = "
  << dKernel.getString() << ")\n"
  << "\t kernel_parameter: width for rbf, degree for polynomial\n"
  << " " << nameC << " double (default = " << dC << ")"
  << " : soft margin parameter in the training step\n"
  << "\nAdditional options:\n"
  << " " << nameVerbose << " int (0-3, default = " << dVerbose << ")"
  << " : verbose mode, controls the level of log detail\n"
  << " " << nameAlphaEpsilon << " double (default = " << dAlphaEpsilon << ")"
  << " : a (positive) threshold below which\n"
  << "\t an alpha coefficient is considered zero; if -1, the threshold\n"
  << "\t is set adaptively to 10^-5 of the largest alpha\n"
  << " " << nameClassEpsilon << " double (default = " << dClassEpsilon << ")"
  << " : threshold for the classification\n"
  << "\t result; has to exceed the threshold to count as correct\n"
  << "\nInterior poit optimization options:\n"
  << " " << nameMaxIterations << " int (default = " << dMaxIterations << ")"
  << " : maximal numer of optimization iterations\n"
  << " " << nameSigDig << " int (default = " << dSigDig << ") "
  << " : number of significant digits to which\n"
  << "\t the primal and the dual cost functions have to agree\n"
  << " " << nameOptEpsilon << " double (default = " << dOptEpsilon << ")"
  << " : epsilon region around\n"
  << "\t the constraints that is enforced during the optimization\n"
  << "\n";

  return s;
}






