%module fsgdf
%{
  #include "fsgdf.h"
%}

%include typemaps.i

%typemap(in) char *OUTSTRING {
     $1 = (char*)calloc(256,sizeof(char));
}
%typemap(argout) char *OUTSTRING {
     Tcl_Obj *o = Tcl_NewStringObj($1,256);
     Tcl_ListObjAppendElement(interp,$result,o);
}
%typemap(freearg) char *OUTSTRING {
     free($1);
}


extern FSGD *gdfRead(char *gdfname);
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

