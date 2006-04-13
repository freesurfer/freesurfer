#include "topology/globals.h"

#include <iostream>
using namespace std;

void check(bool exp){
	if(exp==false)	cout << "e";   
}

void ErrorExit(string s){
	cout << endl << "ERROR: " << s << endl;
	exit(-1);
}
