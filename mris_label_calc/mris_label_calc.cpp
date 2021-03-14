/**
 * @brief A programm to calculate stuff on surface labels
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

#include "mrimorph.h"
#include "version.h"

const char *Progname = NULL;

void printUsage() {
  std::cout << std::endl;
  std::cout << "> mris_label_calc command input1 input2 output" << std::endl;
  std::cout << std::endl;
  std::cout << "   To calculate stuff on surface labels..." << std::endl;
  std::cout << std::endl;
  std::cout << "   commands: " << std::endl;
  std::cout << "      union         union (OR) of both input labels"
            << std::endl;
  std::cout << "      intersect     intersection (AND) of both input labels"
            << std::endl;
  std::cout << "      invert        inverse (NOT) of label on surface (input2)"
            << std::endl;
  std::cout << "      erode <n>     erode  label <n> times on surface (input2)"
            << std::endl;
  std::cout << "      dilate <n>    dilate label <n> times on surface (input2)"
            << std::endl;
  std::cout << std::endl;
}

int main(int argc, char *argv[]) {

  Progname = argv[0];

  if (argc < 5) {
    printUsage();
    exit(1);
  }
  std::string comm = argv[1];

  if (comm == "union") {
    if (argc != 5) {
      std::cerr << "Command 'union' needs 4 arguments: union inlabel1 inlabel2 "
                   "outlabel"
                << std::endl;
      exit(1);
    }
    std::string if1 = argv[2];
    std::string if2 = argv[3];
    std::string of  = argv[4];
    LABEL *     l1  = LabelRead(nullptr, if1.c_str());
    LABEL *     l2  = LabelRead(nullptr, if2.c_str());
    assert(l1 != nullptr);
    assert(l2 != nullptr);
    LABEL *ret = LabelCombine(l1, l2);
    LabelRemoveDuplicates(ret);
    ret->subject_name[0] = '\0';
    LabelWrite(ret, of.c_str());
  } else if (comm == "intersect") {
    if (argc != 5) {
      std::cerr << std::endl
                << "  Command 'intersect' needs 4 arguments: " << std::endl
                << std::endl;
      ;
      std::cerr << "> mris_label_calc intersect inlabel1 inlabel2 outlabel"
                << std::endl;
      std::cerr << std::endl;
      exit(1);
    }
    std::string if1 = argv[2];
    std::string if2 = argv[3];
    std::string of  = argv[4];
    LABEL *     l1  = LabelRead(nullptr, if1.c_str());
    LABEL *     l2  = LabelRead(nullptr, if2.c_str());
    assert(l1 != nullptr);
    assert(l2 != nullptr);
    LabelIntersect(l1, l2);
    l1->subject_name[0] = '\0';
    LabelWrite(l1, of.c_str());

  } else if (comm == "invert") {
    if (argc != 5) {
      std::cerr << std::endl
                << "  Command 'invert' needs 4 arguments:" << std::endl
                << std::endl;
      std::cerr << "> mris_label_calc invert inlabel insurface outlabel"
                << std::endl
                << std::endl;
      exit(1);
    }
    std::string if1       = argv[2];
    std::string if2       = argv[3];
    std::string of        = argv[4];
    LABEL *     l1        = LabelRead(nullptr, if1.c_str());
    MRIS *      surf      = MRISread(if2.c_str());
    LABEL *     linv      = MRISlabelInvert(surf, l1);
    linv->subject_name[0] = '\0';
    LabelWrite(linv, of.c_str());
  } else if (comm == "erode") {
    if (argc != 6) {
      std::cerr << std::endl
                << "  Command 'erode' needs 5 arguments:" << std::endl
                << std::endl;
      std::cerr
          << "> mris_label_calc erode iterations inlabel insurface outlabel"
          << std::endl
          << std::endl;
      exit(1);
    }
    int         it   = atoi(argv[2]);
    std::string if1  = argv[3];
    std::string if2  = argv[4];
    std::string of   = argv[5];
    LABEL *     l1   = LabelRead(nullptr, if1.c_str());
    MRIS *      surf = MRISread(if2.c_str());
    if (LabelErode(l1, surf, it) == NO_ERROR) {
      l1->subject_name[0] = '\0';
      LabelWrite(l1, of.c_str());
    }
  } else if (comm == "dilate") {
    if (argc != 6) {
      std::cerr << std::endl
                << "  Command 'dilate' needs 5 arguments:" << std::endl
                << std::endl;
      std::cerr
          << "> mris_label_calc dilate iterations inlabel insurface outlabel"
          << std::endl
          << std::endl;
      exit(1);
    }
    int         it   = atoi(argv[2]);
    std::string if1  = argv[3];
    std::string if2  = argv[4];
    std::string of   = argv[5];
    LABEL *     l1   = LabelRead(nullptr, if1.c_str());
    MRIS *      surf = MRISread(if2.c_str());
    if (LabelDilate(l1, surf, it, CURRENT_VERTICES) == NO_ERROR) {
      l1->subject_name[0] = '\0';
      LabelWrite(l1, of.c_str());
    }
  } else {
    std::cerr << std::endl
              << "  Command: " << comm << " unknown !" << std::endl;
    printUsage();
    exit(1);
  }
}
