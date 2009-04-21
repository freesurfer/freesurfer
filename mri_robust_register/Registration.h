//
// Registration is a class to compute a robust registration
//
// written by Martin Reuter
// Nov. 4th ,2008
//

#ifndef Registration_H
#define Registration_H

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
//#include <iostream>

class Registration
{
  public:
    Registration():transonly(false),rigid(true),robust(true), sat(-1),
                   iscale(false),rtype(1),subsamplesize(-1),debug(0),
                   mri_source(NULL),mri_target(NULL), Minit(NULL),Mfinal(NULL),
		   mri_weights(NULL), mov2weights(NULL),dst2weights(NULL),
		   lastp(NULL), mri_indexing(NULL) {};
    Registration(MRI * s, MRI *t):transonly(false),rigid(true),robust(true),
                   sat(-1),iscale(false),rtype(1),subsamplesize(-1),debug(0),
                   mri_source(MRIcopy(s,NULL)),mri_target(MRIcopy(t,NULL)),
		   Minit(NULL),Mfinal(NULL),mri_weights(NULL),
		   mov2weights(NULL),dst2weights(NULL),lastp(NULL),
		   mri_indexing(NULL) {};
  
    ~Registration()
    { // we cleanup our private variables
       if (mri_source) MRIfree(&mri_source);
       if (mri_target) MRIfree(&mri_target);
       if (Minit) MatrixFree(&Minit);
       if (Mfinal) MatrixFree(&Minit);
       if (lastp) MatrixFree(&lastp);
       if (mri_indexing) MRIfree(&mri_indexing);
       if (mri_weights) MRIfree(&mri_weights);
       if (gpS.size() > 0) freeGaussianPyramid(gpS);
       if (gpT.size() > 0) freeGaussianPyramid(gpT);
       if (mov2weights) MatrixFree(&mov2weights);
       if (dst2weights) MatrixFree(&dst2weights);
    }
    
    void clear(); // initialize registration (keep source and target and gauss pyramid)
  
    // Set parameters:
    void setTransonly(bool r) {transonly = r;};
    void setRigid(bool r)  {rigid = r;};
    void setRobust(bool r) {robust = r;};
    void setSaturation(double d) {sat = d;};
    void setDebug(int d) {debug = d;};
    void setIscale(bool i) {iscale = i;};
    void setRtype(int r)   {rtype = r;};
    void setMinit(MATRIX* m){Minit = MatrixCopy(m,Minit);};
    void setSource (MRI * s, bool fixvoxel = false, bool fixtype = false);
    void setTarget (MRI * t, bool fixvoxel = false, bool fixtype = false);
    void setSubsamplesize (int sss){subsamplesize = sss;};
    void setName(const std::string &n);

    bool isIscale()        {return iscale;};
    std::string  getName() {return name;};
    MRI * getWeights() {return mri_weights;};
    std::pair <MATRIX*, MATRIX*> getHalfWayMaps() {std::pair <MATRIX*, MATRIX*> md2w(mov2weights,dst2weights); return md2w;};
  
    // compute registration
    std::pair <MATRIX*, double> computeIterativeRegistration( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);
    std::pair <MATRIX*, double> computeIterativeRegSat( int n,double epsit,MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);
    std::pair <MATRIX*, double> computeMultiresRegistration (int stopres, int n,double epsit, MRI * mriS= NULL, MRI* mriT= NULL, MATRIX* Minit = NULL, double iscaleinit = 1.0);
 
    bool warpSource(const std::string & fname, MATRIX* M = NULL, double is = -1);
    bool warpSource(MRI* orig, MRI* target, const std::string &fname, MATRIX* M = NULL, double is = -1);
 
    // testing
    void testRobust(const std::string & fname, int testno);

   double computeSatEstimate (int reslevel, int n,double epsit, MRI * mriS=NULL, MRI* mriT=NULL, MATRIX* mi=NULL, double scaleinit=1.0 );

   MRI *  MRIvalscale(MRI *mri_src, MRI *mri_dst, double s);
   double RigidTransDistSq(MATRIX *a, MATRIX *b = NULL);
   double AffineTransDistSq(MATRIX *a, MATRIX *b = NULL, double r=100);
   MATRIX * MatrixSqrt(MATRIX * m, MATRIX * sqrtm=NULL);
   MRI* makeConform(MRI *mri, MRI *out, bool fixvoxel = true, bool fixtype = true);

  protected:

    //   returns weights:
    std::pair < MATRIX*, MRI* > computeRegistrationStepW(MRI * mriS = NULL, MRI* mriT=NULL);
    //   returns param vector:
    MATRIX* computeRegistrationStepP(MRI * mriS = NULL, MRI* mriT = NULL);
    //   returns 4x4 matrix and iscale:
    std::pair <MATRIX*, double> computeRegistrationStep(MRI * mriS = NULL, MRI* mriT = NULL);
  
  
  private:
  
  // construct Ab and R:
   MATRIX* constructR(MATRIX* p);
   std::pair < MATRIX*, VECTOR* > constructAb(MRI *mriS, MRI *mriT);
   std::pair < MATRIX*, VECTOR* > constructAb2(MRI *mriS, MRI *mriT);

   // my MRI operations (helpers for construct Ab):
   MRI * convolute(MRI * mri, MRI * filter, int dir);
   MRI * getPrefilter();
   MRI * getDerfilter();
   MRI * subSample(MRI * mri);
   MRI * getBlur(MRI* mriS);
   MRI * getPartial(MRI* mriS, int dir);
   bool  getPartials(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur);
   MRI * getBlur2(MRI* mri);
   bool  getPartials2(MRI* mri, MRI* & outfx, MRI* & outfy, MRI* &outfz, MRI* &outblur);
   
   int findRightSize(MRI *mri, float conform_size);
   
   MATRIX* MRIgetZslice(MRI * mri, int slice);
   
   // conversions
   MATRIX* getMatrix(std::vector < double > d, int r, int c=-1, MATRIX* m=NULL);
   MATRIX * rt2mat(MATRIX * r, MATRIX * t, MATRIX *outM); // uses global rtype flag
   MATRIX * p2mat(MATRIX * p6, MATRIX *outM);
   MATRIX * aff2mat(MATRIX * aff, MATRIX *outM);
   std::pair < MATRIX*, double > convertP2Md(MATRIX* p);
   MATRIX * getHalfRT (MATRIX * m, MATRIX * mhalf=NULL);
   double RotMatrixLogNorm(MATRIX * m);
   double RotMatrixGeoDist(MATRIX * a, MATRIX *b = NULL);
   
   // gaussian pyramid:
   std::vector < MRI* > buildGaussianPyramid (MRI * mri_in, int n);
   void freeGaussianPyramid(std::vector< MRI* >& p);

   // tools
   double getFrobeniusDiff(MATRIX *m1, MATRIX *m2);

   // initial registration using pca DOES NOT WORK!
   MATRIX * initialize_transform(MRI *mri_in, MRI *mri_ref);
   int init_scaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L);
   int order_eigenvectors(MATRIX *m_src_evectors, MATRIX *m_dst_evectors);
   MATRIX * compute_pca(MRI *mri_in, MRI *mri_ref);
   MATRIX * pca_matrix(MATRIX *m_in_evectors, double in_means[3], MATRIX *m_ref_evectors, double ref_means[3]);
   
   
    bool transonly;
    bool rigid;
    bool robust;
    double sat;
    bool iscale;
    int rtype;
    int subsamplesize;
    std::string name;
    std::string nbase;
    int debug;
    //bool outweights;
    //std::string weightsname;
    
    MRI * mri_source;
    std::vector < MRI* > gpS;
    MRI * mri_target;
    std::vector < MRI* > gpT;
    MATRIX * Minit;
    MATRIX * Mfinal;
    double iscalefinal;
    MRI * mri_weights;
    MATRIX * mov2weights;
    MATRIX * dst2weights;
    
    // help vars
    MATRIX* lastp;
    double zeroweights;
    
    MRI * mri_indexing;
};


#endif
