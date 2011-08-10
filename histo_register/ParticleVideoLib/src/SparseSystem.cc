#include <pvl/SparseSystem.h>
#include <math.h>
namespace pvl {


//-------------------------------------------
// CONSTRUCTOR / DESTRUCTOR
//-------------------------------------------


// basic constructor
SparseSystem::SparseSystem( int equationCount ) {

    // equations
    m_fixedEqn5 = NULL;
    m_fixedEqn6 = NULL;
    m_fixedEqn9 = NULL;
    m_varEqn = NULL;
    m_b.setLength( equationCount );
    m_eqnCount = equationCount;
    m_nextEqnIndex = 0;

    // solver parameters
    m_maxIter = 0;
    m_sorFactor = 1.5f;
    m_indent = 0;
    m_adaptIterMask = 0;
    m_adaptUpdateThresh = 0;
}


// deallocate the equations
SparseSystem::~SparseSystem() {
    if (m_fixedEqn5) delete [] m_fixedEqn5;
    if (m_fixedEqn6) delete [] m_fixedEqn6;
    if (m_fixedEqn9) delete [] m_fixedEqn9;
    if (m_varEqn) delete [] m_varEqn;
}


//-------------------------------------------
// BUILD EQUATIONS
//-------------------------------------------


/// adds term to variable length equation 
void SparseSystem::addTerm( int eqnIndex, int varIndex, float a ) {

    // if first time, allocate and initialize equations
    if (m_varEqn == NULL) {
        m_varEqn = new VarEquation[ m_eqnCount ];
        if (m_varEqn == NULL)
            fatalError( "SparseSystem::AddTerm: allocate failed" );
        for (int i = 0; i < m_eqnCount; i++) {
            m_varEqn[ i ].termCount = 0;
        }
    }

    // get current equation
    VarEquation *eqn = m_varEqn + eqnIndex;

    // sanity checks
    if (eqnIndex >= m_eqnCount)
        fatalError( "SparseSystem::AddTerm: eqnIndex too large" );
    if (eqnIndex < 0)
        fatalError( "SparseSystem::AddTerm: eqnIndex less than 0" );
    if (varIndex >= m_eqnCount)
        fatalError( "SparseSystem::AddTerm: varIndex too large" );
    if (varIndex < 0)
        fatalError( "SparseSystem::AddTerm: varIndex less than 0" );
    if (eqn->termCount >= MAX_TERM_COUNT)
        fatalError( "SparseSystem::AddTerm: too many terms" );
    if (eqn->termCount < 0)
        fatalError( "SparseSystem::AddTerm: invalid eqn" );

    // check for existing term
    for (int i = 0; i < eqn->termCount; i++) {
        if (eqn->index[ i ] == varIndex) {
            eqn->a[ i ] += a;
            if (i == 0) // update a1Inv if necessary
                eqn->a1Inv = 1.0f / eqn->a[ 0 ];
            return;
        }
    }

    // add new term
    int termIndex = eqn->termCount;
    eqn->a[ termIndex ] = a;
    eqn->index[ termIndex ] = varIndex;
    if (termIndex == 0) { // update a1Inv if necessary
        if (varIndex != eqnIndex)
            fatalError( "must add diagonal terms first" );
        eqn->a1Inv = 1.0f / a;
    }
    eqn->termCount++;
}


/// add eqn of form: a_i x_i + a_j x_j = b
void SparseSystem::addEquation( int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8, float a9, float b ) {
    if (m_fixedEqn9 == NULL)
        m_fixedEqn9 = new FixedEquation9[ m_eqnCount ]; 
    if (m_fixedEqn9 == NULL)
        fatalError( "SparseSystem::AddEquation: out of memory" );
    FixedEquation9 *eqn = m_fixedEqn9 + m_nextEqnIndex;
    eqn->i1 = i1;
    eqn->i2 = i2;
    eqn->i3 = i3;
    eqn->i4 = i4;
    eqn->i5 = i5;
    eqn->i6 = i6;
    eqn->i7 = i7;
    eqn->i8 = i8;
    eqn->i9 = i9;
    eqn->a1 = a1;
    eqn->a1Inv = 1.0f / a1;
    eqn->a2 = a2;
    eqn->a3 = a3;
    eqn->a4 = a4;
    eqn->a5 = a5;
    eqn->a6 = a6;
    eqn->a7 = a7;
    eqn->a8 = a8;
    eqn->a9 = a9;
    m_b[ m_nextEqnIndex ] = b;
    m_nextEqnIndex++;
}


/// add eqn of form: a_i x_i + a_j x_j = b
void SparseSystem::addEquation( int i1, int i2, int i3, int i4, int i5, int i6, float a1, float a2, float a3, float a4, float a5, float a6, float b ) {
    if (m_fixedEqn6 == NULL)
        m_fixedEqn6 = new FixedEquation6[ m_eqnCount ]; 
    if (m_fixedEqn6 == NULL)
        fatalError( "SparseSystem::AddEquation: out of memory" );
    FixedEquation6 *eqn = m_fixedEqn6 + m_nextEqnIndex;
    eqn->i1 = i1;
    eqn->i2 = i2;
    eqn->i3 = i3;
    eqn->i4 = i4;
    eqn->i5 = i5;
    eqn->i6 = i6;
    eqn->a1 = a1;
    eqn->a1Inv = 1.0f / a1;
    eqn->a2 = a2;
    eqn->a3 = a3;
    eqn->a4 = a4;
    eqn->a5 = a5;
    eqn->a6 = a6;
    m_b[ m_nextEqnIndex ] = b;
    m_nextEqnIndex++;
}


/// add eqn of form: a_i x_i + a_j x_j = b
void SparseSystem::addEquation( int i1, int i2, int i3, int i4, int i5, float a1, float a2, float a3, float a4, float a5, float b ) {
    if (m_fixedEqn5 == NULL)
        m_fixedEqn5 = new FixedEquation5[ m_eqnCount ]; 
    if (m_fixedEqn5 == NULL)
        fatalError( "SparseSystem::AddEquation: out of memory" );
    FixedEquation5 *eqn = m_fixedEqn5 + m_nextEqnIndex;
    eqn->i1 = i1;
    eqn->i2 = i2;
    eqn->i3 = i3;
    eqn->i4 = i4;
    eqn->i5 = i5;
    eqn->a1 = a1;
    eqn->a1Inv = 1.0f / a1;
    eqn->a2 = a2;
    eqn->a3 = a3;
    eqn->a4 = a4;
    eqn->a5 = a5;
    m_b[ m_nextEqnIndex ] = b;
    m_nextEqnIndex++;
}



//-------------------------------------------
// DIAGNOSTICS
//-------------------------------------------


/// check whether equations are valid
void SparseSystem::checkEquations() {

    // if using length 5 equations
    if (m_fixedEqn5) {
        if (m_eqnCount != m_x.length())
            fatalError( "system is not square" );
        for (int i = 0; i < m_eqnCount; i++) {
            if (m_fixedEqn5[ i ].i1 != i) {
                fatalError( "first var is not on diagonal" );
                break;
            }
            if (m_fixedEqn5[ i ].a1 < 1e-8 && m_fixedEqn5[ i ].a1 > -1e-8) {
                disp( 1, "diagonal: %f", m_fixedEqn5[ i ].a1 );
                fatalError( "diagonal too small" );
                break;
            }
        }
    }

    // if using length 6 equations
    if (m_fixedEqn6) {
        if (m_eqnCount != m_x.length())
            fatalError( "system is not square" );
        for (int i = 0; i < m_eqnCount; i++) {
            if (m_fixedEqn6[ i ].i1 != i) {
                fatalError( "first var is not on diagonal (%d, %d)", m_fixedEqn6[ i ].i1, i );
                break;
            }
            if (m_fixedEqn6[ i ].a1 < 1e-8 && m_fixedEqn6[ i ].a1 > -1e-8) {
                disp( 1, "diagonal: %f", m_fixedEqn6[ i ].a1 );
                fatalError( "diagonal too small" );
                break;
            }
        }
    }

    // if using length 9 equations
    if (m_fixedEqn9) {
        if (m_eqnCount != m_x.length())
            fatalError( "system is not square" );
        for (int i = 0; i < m_eqnCount; i++) {
            if (m_fixedEqn9[ i ].i1 != i) {
                fatalError( "first var is not on diagonal (%d, %d)", m_fixedEqn9[ i ].i1, i );
                break;
            }
            if (m_fixedEqn9[ i ].a1 < 1e-8 && m_fixedEqn9[ i ].a1 > -1e-8) {
                disp( 1, "diagonal: %f", m_fixedEqn9[ i ].a1 );
                fatalError( "diagonal too small" );
                break;
            }
        }
    }

    // if using var-length equations
    if (m_varEqn) {
        if (m_eqnCount != m_x.length())
            fatalError( "system is not square" );
        int minTermCount = 1000, maxTermCount = 0;
        float minCoef = 1e8, maxCoef = -1e8;
        float minInv = 1e8, maxInv = -1e8;
        for (int i = 0; i < m_eqnCount; i++) {
            VarEquation *eqn = m_varEqn + i;
            int termCount = eqn->termCount;
            if (termCount < minTermCount)
                minTermCount = termCount;
            if (termCount > maxTermCount)
                maxTermCount = termCount;
            for (int j = 1; j < eqn->termCount; j++) {
                float coef = eqn->a[ j ];
                if (coef < minCoef) minCoef = coef;
                if (coef > maxCoef) maxCoef = coef;
            }
            if (eqn->a[ 0 ] < 1e-4 && eqn->a[ 0 ] > -1e-4)
                fatalError( "diagonal too small" );
            float inv = eqn->a1Inv;
            if (inv < minInv) minInv = inv;
            if (inv > maxInv) maxInv = inv;
            if (i != eqn->index[ 0 ])
                fatalError( "first var is not on diagonal" );
        }
        if (m_indent) 
            disp( m_indent, "var-length system: eqns: %d, min terms: %d, max terms: %d, min coef: %f, max coef: %f, min inv: %f, max inv: %f", 
                  m_eqnCount, minTermCount, maxTermCount, minCoef, maxCoef, minInv, maxInv );
    }
}


/// return L2 norm of Ax - b
double SparseSystem::residual( const VectorF &xVect ) {
    double sumSqDiff = 0;
    const float *x = xVect.dataPtr();

    // if using fixed-length equations
    if (m_fixedEqn5) {
        for (int i = 0; i < m_eqnCount; i++) {
            FixedEquation5 *eqn = m_fixedEqn5 + i;
            double axSum = eqn->a1 * x[ eqn->i1 ] + eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ];
            double diff = axSum - m_b[ i ];
            sumSqDiff += diff * diff;
        }
    }

    // if using fixed-length equations
    if (m_fixedEqn6) {
        for (int i = 0; i < m_eqnCount; i++) {
            FixedEquation6 *eqn = m_fixedEqn6 + i;
            double axSum = eqn->a1 * x[ eqn->i1 ] + eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ] + eqn->a6 * x[ eqn->i6 ];
            double diff = axSum - m_b[ i ];
            sumSqDiff += diff * diff;
        }
    }

    // if using fixed-length equations
    if (m_fixedEqn9) {
        for (int i = 0; i < m_eqnCount; i++) {
            FixedEquation9 *eqn = m_fixedEqn9 + i;
            double axSum = eqn->a1 * x[ eqn->i1 ] + eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ] + eqn->a6 * x[ eqn->i6 ] + eqn->a7 * x[ eqn->i7 ] + eqn->a8 * x[ eqn->i8 ] + eqn->a9 * x[ eqn->i9 ];
            double diff = axSum - m_b[ i ];
            sumSqDiff += diff * diff;
        }
    }

    // if using var-length equations
    if (m_varEqn) {
        for (int i = 0; i < m_eqnCount; i++) {
            VarEquation *eqn = m_varEqn + i;
            double axSum = 0;
            for (int j = 0; j < eqn->termCount; j++) {
                axSum += eqn->a[ j ] * x[ eqn->index[ j ] ];
            }
            double diff = axSum - m_b[ i ];
            sumSqDiff += diff * diff;
        }
    }

    // return L2 norm of residual
    return sqrtf( (float) (sumSqDiff / m_eqnCount) );
}


/// display all equations 
void SparseSystem::dispAll() {
    if (m_varEqn) {
        for (int eqnIndex = 0; eqnIndex < m_eqnCount; eqnIndex++) {
            VarEquation *eqn = m_varEqn + eqnIndex;
            disp( 1, "equation %d:", eqnIndex );
            for (int i = 0; i < eqn->termCount; i++) {
                disp( 2, "%f * x%d", eqn->a[ i ], eqn->index[ i ] );
            }
            disp( 3, "= %f", m_b[ eqnIndex ] );
        }
    }
    // fix(later): add other cases
}


//-------------------------------------------
// SOLVERS
//-------------------------------------------


/// solve and return resulting x vector using SOR algorithm
VectorF SparseSystem::solveSOR() {
    int cols = m_x.length(); 
    checkEquations();

    // copy the initial guess into the solution vector
    VectorF x( m_x );

    // if using length-5 fixed SOR equations
    if (m_fixedEqn5) {

        // iterate 
        for (int iter = 0; iter < m_maxIter; iter++) {

            // loop over elements of x
            for (int i = 0; i < cols; i++) {
                FixedEquation5 *eqn = m_fixedEqn5 + i;

                // compute Ax with off diagonal elements
                float ax = eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ];
                ax = (m_b[ i ] - ax) * eqn->a1Inv;

                // update x
                x[ i ] += m_sorFactor * (ax - x[ i ]);
            }
        }

    // if using length-6 fixed SOR equations
    } else if (m_fixedEqn6) {

        // iterate
        VectorF xUpdate( x.length() );
        VectorI xActive( x.length() );
        xActive.clear( 1 );
        for (int iter = 0; iter < m_maxIter; iter++) {

            // if adaptive
            if (m_adaptIterMask) {
            
                // loop over elements of x
                for (int i = 0; i < cols; i++) {
                    if (xActive[ i ]) {
                        FixedEquation6 *eqn = m_fixedEqn6 + i;

                        // compute Ax with off diagonal elements
                        float ax = eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ] + eqn->a6 * x[ eqn->i6 ];
                        ax = (m_b[ i ] - ax) * eqn->a1Inv;

                        // update x
                        xUpdate[ i ] = m_sorFactor * (ax - x[ i ]);
                        x[ i ] += xUpdate[ i ];
                    }
                }

                // every N iterations, check updates; disable vars with small updates
                // N = m_adaptIterThresh + 1; assumes N is power of 2
                if ((iter & m_adaptIterMask) == m_adaptIterMask) {
                    for (int i = 0; i < cols; i++) {
                        xActive[ i ] = xUpdate[ i ] > m_adaptUpdateThresh || xUpdate[ i ] < -m_adaptUpdateThresh;
                    }
                }

            // if not adaptive
            } else {

                // loop over elements of x
                for (int i = 0; i < cols; i++) {
                    FixedEquation6 *eqn = m_fixedEqn6 + i;

                    // compute Ax with off diagonal elements
                    float ax = eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ] + eqn->a6 * x[ eqn->i6 ];
                    ax = (m_b[ i ] - ax) * eqn->a1Inv;

                    // update x
                    x[ i ] += m_sorFactor * (ax - x[ i ]);
                }    
            }
        }

    // if using length-9 fixed SOR equations
    } else if (m_fixedEqn9) {

        // iterate
        VectorF xUpdate( x.length() );
        VectorI xActive( x.length() );
        xActive.clear( 1 );
        for (int iter = 0; iter < m_maxIter; iter++) {

            // if adaptive
            if (m_adaptIterMask) {
            
                // loop over elements of x
                for (int i = 0; i < cols; i++) {
                    if (xActive[ i ]) {
                        FixedEquation9 *eqn = m_fixedEqn9 + i;

                        // compute Ax with off diagonal elements
                        float ax = eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ] + eqn->a6 * x[ eqn->i6 ] + eqn->a7 * x[ eqn->i7 ] + eqn->a8 * x[ eqn->i8 ] + eqn->a9 * x[ eqn->i9 ];
                        ax = (m_b[ i ] - ax) * eqn->a1Inv;

                        // update x
                        xUpdate[ i ] = m_sorFactor * (ax - x[ i ]);
                        x[ i ] += xUpdate[ i ];
                    }
                }

                // every N iterations, check updates; disable vars with small updates
                // N = m_adaptIterThresh + 1; assumes N is power of 2
                if ((iter & m_adaptIterMask) == m_adaptIterMask) {
                    for (int i = 0; i < cols; i++) {
                        xActive[ i ] = xUpdate[ i ] > m_adaptUpdateThresh || xUpdate[ i ] < -m_adaptUpdateThresh;
                    }
                }

            // if not adaptive
            } else {

                // loop over elements of x
                for (int i = 0; i < cols; i++) {
                    FixedEquation9 *eqn = m_fixedEqn9 + i;

                    // compute Ax with off diagonal elements
                    float ax = eqn->a2 * x[ eqn->i2 ] + eqn->a3 * x[ eqn->i3 ] + eqn->a4 * x[ eqn->i4 ] + eqn->a5 * x[ eqn->i5 ] + eqn->a6 * x[ eqn->i6 ] + eqn->a7 * x[ eqn->i7 ] + eqn->a8 * x[ eqn->i8 ] + eqn->a9 * x[ eqn->i9 ];
                    ax = (m_b[ i ] - ax) * eqn->a1Inv;

                    // update x
                    x[ i ] += m_sorFactor * (ax - x[ i ]);
                }    
            }
        }

    // if using variable-length equations
    } else {

        // iterate 
        for (int iter = 0; iter < m_maxIter; iter++) {

            // loop over elements of x
            for (int i = 0; i < cols; i++) {
                VarEquation *eqn = m_varEqn + i;

                // compute Ax with off diagonal elements
                float ax = 0;
                for (int j = 1; j < eqn->termCount; j++) 
                    ax += eqn->a[ j ] * x[ eqn->index[ j ] ];
                ax = (m_b[ i ] - ax) * eqn->a1Inv;

                // update x
                x[ i ] += m_sorFactor * (ax - x[ i ]);
            }
        }

        // display info
        if (m_indent)
            disp( m_indent, "solved SOR system: eqns: %d, iters: %d, factor: %f, resid: %lf", 
                  m_eqnCount, m_maxIter, m_sorFactor, residual( x ) );
    }

    // return solution vector
    return x; 
}


/// solve and return resulting x vector
VectorF SparseSystem::solve() {
    return solveSOR();
}


} // end namespace pvl 

