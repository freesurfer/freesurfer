#define export
#include <vnl/vnl_nonlinear_minimizer.h>
#undef export

#include <stdio.h>
#include <stdlib.h>

#if 0
    // THESE ARE DEFINED IN vnl_nonlinear_minimizer.cxx
    // centos6-x86_64-packages/vxl/current/lib/libvnl.a

  //: Return the name of the class.
  //  Used by polymorphic IO
vcl_string vnl_nonlinear_minimizer::is_a() const
{
    fprintf(stderr, "%s:%d called vnl_nonlinear_minimizer::is_a\n",
    	__FILE__,__LINE__);
    exit(1);
    return NULL;
}

  //: Return true if the name of the class matches the argument.
  //  Used by polymorphic IO
bool vnl_nonlinear_minimizer::is_class(vcl_string const& s) const
{
    fprintf(stderr, "%s:%d called vnl_nonlinear_minimizer::is_class\n",
    	__FILE__,__LINE__);
    exit(1);
    return false;
}

#endif
