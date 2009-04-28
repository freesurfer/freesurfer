//
// RegPowell is a class to compute a robust registration using Powell
//
// written by Martin Reuter
// Apr. 22th ,2009
//

#ifndef RegPowell_H
#define RegPowell_H

#ifdef __cplusplus
extern "C" {
#endif
#include "matrix.h"
#include "mri.h"
#ifdef __cplusplus
}
#endif

#include <utility>
#include <string>
#include <vector>
#include "Registration.h"


class RegPowell : public Registration
{
  public:
    RegPowell():Registration(){};
    RegPowell(MRI * s, MRI *t):Registration(s,t){};
    virtual ~RegPowell(){};
    
    virtual std::pair <MATRIX*, double> computeIterativeRegistration( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);

  protected:
  
    static float costFunction(float p[] );
    static RegPowell* tocurrent;
    static MRI * scf;
    static MRI * tcf;
    static int pcount;
    static MATRIX * mh1;
    static MATRIX * mh2;
    static int icount;
    
};


#endif

