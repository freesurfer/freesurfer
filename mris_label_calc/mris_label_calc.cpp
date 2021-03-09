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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"
#include "label.h"

using namespace std;
const char *Progname = NULL;


void printUsage()
{
    cout << endl;
    cout << "> mris_label_calc command input1 input2 output" << endl;
    cout << endl;
		cout << "   To calculate stuff on surface labels..." << endl;
		cout << endl;
    cout << "   commands: " << endl;
    cout << "      union         union (OR) of both input labels" << endl;
    cout << "      intersect     intersection (AND) of both input labels" << endl;
    cout << "      invert        inverse (NOT) of label on surface (input2)" << endl;
    cout << "      erode <n>     erode  label <n> times on surface (input2)" << endl;
    cout << "      dilate <n>    dilate label <n> times on surface (input2)" << endl;
    cout << endl;
}

int main(int argc, char *argv[])
{

  Progname = argv[0] ;

  if (argc < 5)
  {
    printUsage();
    exit(1);
  }
  string comm  = argv[1];
	
	if (comm == "union")
	{
	  if (argc != 5)
		{
		   cerr << "Command 'union' needs 4 arguments: union inlabel1 inlabel2 outlabel" << endl;
			 exit(1);
		}
    string if1 = argv[2];
    string if2 = argv[3];
	  string of  = argv[4];
    LABEL *l1  = LabelRead(NULL,if1.c_str());
	  LABEL *l2  = LabelRead(NULL,if2.c_str());
	  assert (l1 != NULL);
	  assert (l2 != NULL);
	  LABEL * ret = LabelCombine(l1,l2);
		LabelRemoveDuplicates(ret);
		ret->subject_name[0]='\0';
		LabelWrite(ret,of.c_str());
	}
	else if (comm == "intersect")
	{
	  if (argc != 5)
		{
		   cerr << endl << "  Command 'intersect' needs 4 arguments: " << endl << endl;;
			 cerr << "> mris_label_calc intersect inlabel1 inlabel2 outlabel" << endl;
			 cerr << endl;
			 exit(1);
		}
    string if1 = argv[2];
    string if2 = argv[3];
	  string of  = argv[4];
    LABEL *l1  = LabelRead(NULL,if1.c_str());
	  LABEL *l2  = LabelRead(NULL,if2.c_str());
	  assert (l1 != NULL);
	  assert (l2 != NULL);
		LabelIntersect(l1,l2);
		l1->subject_name[0]='\0';
		LabelWrite(l1,of.c_str());
	
	}
	else if (comm == "invert")
	{
	  if (argc != 5)
		{
		   cerr << endl << "  Command 'invert' needs 4 arguments:" << endl << endl;
			 cerr << "> mris_label_calc invert inlabel insurface outlabel" << endl<< endl;
			 exit(1);
		}
    string if1  = argv[2];
    string if2  = argv[3];
	  string of   = argv[4];
	  LABEL *l1   = LabelRead(NULL,if1.c_str());
		MRIS *surf  = MRISread(if2.c_str());
	  LABEL *linv = MRISlabelInvert(surf,l1);
		linv->subject_name[0]='\0';
		LabelWrite(linv,of.c_str());
	}
	else if (comm == "erode")
	{
	  if (argc != 6)
		{
		   cerr << endl << "  Command 'erode' needs 5 arguments:"<< endl << endl;
			 cerr << "> mris_label_calc erode iterations inlabel insurface outlabel" << endl<< endl;
			 exit(1);
		}
		int it      = atoi(argv[2]);
    string if1  = argv[3];
    string if2  = argv[4];
	  string of   = argv[5];
	  LABEL *l1   = LabelRead(NULL,if1.c_str());
		MRIS *surf  = MRISread(if2.c_str());
		if (LabelErode(l1,surf,it) == NO_ERROR)
		{
		  l1->subject_name[0]='\0';
		  LabelWrite(l1,of.c_str());		  
		}
  }
	else if (comm == "dilate")
	{
	  if (argc != 6)
		{
		   cerr << endl << "  Command 'dilate' needs 5 arguments:" << endl << endl;
			 cerr << "> mris_label_calc dilate iterations inlabel insurface outlabel" << endl << endl;
			 exit(1);
		}
		int it      = atoi(argv[2]);
    string if1  = argv[3];
    string if2  = argv[4];
	  string of   = argv[5];
	  LABEL *l1   = LabelRead(NULL,if1.c_str());
		MRIS *surf  = MRISread(if2.c_str());
		if (LabelDilate(l1,surf,it, CURRENT_VERTICES) == NO_ERROR)
		{
		  l1->subject_name[0]='\0';
		  LabelWrite(l1,of.c_str());		  
		}
  }
	else
	{
	  cerr << endl<< "  Command: " << comm << " unknown !" <<  endl;
		printUsage();
		exit(1);
  }

}
