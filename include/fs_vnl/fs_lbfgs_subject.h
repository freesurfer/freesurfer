#ifndef FS_LBFGS_SUBJECT_H_
#define FS_LBFGS_SUBJECT_H_

#include "fs_vnl/fs_lbfgs_observer.h"

class fs_lbfgs_subject
{
public:
  fs_lbfgs_subject();
  ~fs_lbfgs_subject();

  void setObserver( fs_lbfgs_observer* observer );
  void notify( double bestF, vnl_vector< double >* bestX );

private:
  fs_lbfgs_observer* mObserver;

};

#endif /*FS_LBFGS_SUBJECT_H_*/
