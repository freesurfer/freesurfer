#ifndef _PVL_SPARSE_SYSTEM_H_
#define _PVL_SPARSE_SYSTEM_H_
#include <sbl/math/Vector.h>
using namespace sbl;
namespace pvl {


// fixed-length equation structure 
struct FixedEquation9 {
    int i1;
    int i2;
    int i3;
    int i4;
    int i5;
    int i6;
    int i7;
    int i8;
    int i9;
    float a1;
    float a2;
    float a3;
    float a4;
    float a5;
    float a6;
    float a7;
    float a8;
    float a9;
    float a1Inv;
};


// fixed-length equation structure 
struct FixedEquation6 {
    int i1;
    int i2;
    int i3;
    int i4;
    int i5;
    int i6;
    float a1;
    float a2;
    float a3;
    float a4;
    float a5;
    float a6;
    float a1Inv;
};


// fixed-length equation structure 
struct FixedEquation5 {
    int i1;
    int i2;
    int i3;
    int i4;
    int i5;
    float a1;
    float a2;
    float a3;
    float a4;
    float a5;
    float a1Inv;
};


// maximum number of terms in variable-length equation structure
#define MAX_TERM_COUNT 40


// variable-length equation structure
struct VarEquation {
    int termCount;
    int index[ MAX_TERM_COUNT ];
    float a[ MAX_TERM_COUNT ];
    float a1Inv;
};


/// The SparseSystem class represents a system of linear equations with sparse coefficients.
class SparseSystem {
public:

    // constructor / destructor
    SparseSystem( int equationCount );
    ~SparseSystem();

    /// optionally set maximum number of iterations
    inline void setMaxIter( int maxIter ) { m_maxIter = maxIter; }

    /// set SOR factor (0 to 2)
    inline void setSORFactor( float sorFactor ) { m_sorFactor = sorFactor; }

    /// enable adaptive SOR
    inline void enableAdaptive( int adaptIterMask, float adaptUpdateThresh ) { m_adaptIterMask = adaptIterMask; m_adaptUpdateThresh = adaptUpdateThresh; }

    /// return L2 norm of Ax - b
    double residual( const VectorF &x );

    /// display all equations
    void dispAll();

    /// add eqn of form: a_i x_i + a_j x_j = b;
    /// for SOR method, first variable should be on diagonal and have non-zero coef
    // fix(clean): unify fixed-length cases
    void addEquation( int i1, int i2, int i3, int i4, int i5, float a1, float a2, float a3, float a4, float a5, float b );
    void addEquation( int i1, int i2, int i3, int i4, int i5, int i6, float a1, float a2, float a3, float a4, float a5, float a6, float b );
    void addEquation( int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8, float a9, float b );

    // variable length equations
    inline int moreTermsAllowed( int eqnIndex ) { return MAX_TERM_COUNT - m_varEqn[ eqnIndex ].termCount; }
    void addTerm( int eqnIndex, int varIndex, float a );
    void setB( int eqnIndex, float b ) { m_b[ eqnIndex ] = b; }
    void addToB( int eqnIndex, float b ) { m_b[ eqnIndex ] += b; }

    /// set initial value of x (copies arg)
    void setInit( const VectorF &init ) { m_x = init; }

    /// solve and return resulting x vector
    VectorF solve();

    /// set indentation of displayed info (or 0 for no diplayed info, the default)
    void setDisplayIndent( int indent ) { m_indent = indent; } 

private:

    /// solve and return resulting x vector using SOR algorithm
    VectorF solveSOR();

    /// check whether equations are valid
    void checkEquations();

    // equations 
    VectorF m_b;
    VectorF m_x;
    FixedEquation5 *m_fixedEqn5;
    FixedEquation6 *m_fixedEqn6;
    FixedEquation9 *m_fixedEqn9;
    VarEquation *m_varEqn;
    int m_eqnCount;
    int m_nextEqnIndex;

    // optimization parameters
    int m_maxIter;
    float m_sorFactor;
    int m_indent;

    // parameters for adaptive SOR; cuts short optimization of parameters that aren't changing
    int m_adaptIterMask;
    float m_adaptUpdateThresh;

    // disable copy constructor and assignment operator
    SparseSystem( const SparseSystem &x );
    SparseSystem &operator=( const SparseSystem &x );
};


} // end namespace pvl
#endif // _PVL_SPARSE_SYSTEM_H_

