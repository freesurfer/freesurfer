#include "topology/globals.h"

#include <iostream>
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif
	//#include <stdio.h>
	//#include <stdlib.h>
	//#include <string.h>
	//#include <math.h>
	//#include <ctype.h>
	//#include "macros.h"
	#include "error.h"
	//#include "tags.h"
	//#include "diag.h"
	//#include "proto.h"
	//#include "timer.h"
	//#include "mrisurf.h"
	//#include "mri.h"
	#include "macros.h"
	//#include "icosahedron.h"
	//#include "mrishash.h"
	//#include "version.h"
#ifdef __cplusplus
}
#endif


void check(bool exp){
	if(exp==false)	cout << "e";   
}

void ErrorExit(string s){
	cout << endl << "ERROR: " << s << endl;
	exit(-1);
}

int Random(int nmax){
	//  return rand()*nmax/RAND_MAX;
	return nint(randomNumber(0.0, (double)nmax-1));
}

