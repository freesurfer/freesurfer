#include "fs_vnl/fs_lbfgs_subject.h"

fs_lbfgs_subject::fs_lbfgs_subject() {
  mObserver = NULL;
}

fs_lbfgs_subject::~fs_lbfgs_subject() {
}

void fs_lbfgs_subject::setObserver( fs_lbfgs_observer* observer ) {
  mObserver = observer;
}

void fs_lbfgs_subject::notify( double bestF, vnl_vector< double >* bestX ) {
  // we might want a number of observers, but right now we only have one
  if( mObserver != NULL ) {
    mObserver->update( bestF, bestX );
  }
}
