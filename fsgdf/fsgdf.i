// To build fsgdf_wrap.i, do something like:
// /home/kteich/local/bin/swig -I/home/kteich/local/swig -I/home/kteich/local/swig/tcl -tcl fsgdf.i

%module fsgdf
%{
  #include "fsgdf.h"
%}

%include typemaps.i

%typemap(in) char *OUTSTRING {
     $1 = (char*)calloc(256,sizeof(char));
}
%typemap(argout) char *OUTSTRING {
     int len = strlen($1);
     Tcl_Obj *o = Tcl_NewStringObj($1,len);
     Tcl_ListObjAppendElement(interp,$result,o);
}
%typemap(freearg) char *OUTSTRING {
     free($1);
}

extern FSGD *gdfRead(char *gdfname, int loaddata);
extern int gdfReadRegistration(FSGD *gd, int type, char *regname,
			MATRIX* tkregmat,
			mriTransformRef client_transform,
			mriVolumeRef client_volume);
extern int gdfPrintStdout(FSGD *gd);
extern int gdfGetTitle(FSGD *gd, char *OUTSTRING);
extern int gdfGetMeasurementName(FSGD *gd, char *OUTSTRING);
extern int gdfGetSubjectName(FSGD *gd, char *OUSTRING);
extern int gdfGetDataFileName(FSGD *gd, char *OUTSTRING);
extern int gdfGetNumClasses(FSGD *gd, int *OUTPUT);
extern int gdfGetNthClassLabel(FSGD *gd, int nclass, char *OUTSTRING);
extern int gdfGetNthClassMarker(FSGD *gd, int nclass, char *OUTSTRING);
extern int gdfGetNthClassColor(FSGD *gd, int nclass, char *OUTSTRING);
extern int gdfGetNumVariables(FSGD *gd, int *OUTPUT);
extern int gdfGetNthVariableLabel(FSGD *gd, int nvariable, char *OUTSTRING);
extern int gdfGetDefaultVariable(FSGD *gd, char *OUTSTRING);
extern int gdfGetDefaultVariableIndex(FSGD *gd, int *OUTPUT);
extern int gdfGetNumSubjects(FSGD *gd, int *OUTPUT);
extern int gdfGetNthSubjectID(FSGD *gd, int nsubject, char *OUTSTRING);
extern int gdfGetNthSubjectClass(FSGD *gd, int nsubject, int *OUTPUT);
extern int gdfGetNthSubjectNthValue(FSGD *gd, int nsubject, 
				     int nvariable, float *OUTPUT);
extern int gdfGetNthSubjectMeasurement(FSGD *gd, int nsubject, 
					int x, int y, int z, float *OUTPUT);
extern int gdfOffsetSlope(FSGD *gd, int nclass, int nvar,
			  int x, int y, int z, float *OUTPUT, float *OUTPUT);
