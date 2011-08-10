#ifndef _SBL_OPTIMIZER_H_
#define _SBL_OPTIMIZER_H_
#include <sbl/core/String.h>
#include <sbl/math/Vector.h>
#include <sbl/math/Matrix.h>
namespace sbl {


//----------------------------------
// OBJECTIVE CLASS
//----------------------------------


/// The Objective class provides an interface to an objective function to be optimized.
class Objective {
public:

    // basic constructor/destructor
    Objective() {}
    virtual ~Objective() {}

    /// evaluate objective function at given point
    virtual double eval( const VectorD &point ) = 0;

private:

    // disable copy constructor and assignment operator
    Objective( const Objective &x );
    Objective &operator=( const Objective &x );
};

    
//----------------------------------
// OPTIMIZER CLASS
//----------------------------------


/// The Optimizer class defines an interface to numerical optimization methods.
class Optimizer {
public:

    // basic constructor
    Optimizer( Objective &objective );
    virtual ~Optimizer() {}

    /// set the optimization starting point
    inline void setStart( const VectorD &startPoint ) { m_startPoint = startPoint; }

    /// set optimization bounds; points outside the bounds receive a penality
    inline void setBounds( const VectorD &lBound, const VectorD &uBound ) { m_lBound = lBound; m_uBound = uBound; }

    /// set scale factor for penalty applied to points outside bounds
    inline void setPenaltyFactor( double penaltyFactor ) { m_penaltyFactor = penaltyFactor; }

    /// set termination tolerance (relative difference between objective values)
    inline void setFinalTolerance( double finalTolerance ) { m_finalTolerance = finalTolerance; }

    /// run the optimization until termination (or user cancel); returns best point found
    virtual VectorD run( double *finalObjValue = 0 ) = 0;

    /// repeatedly run the optimization until termination (or user cancel); returns best point found
    VectorD repeatRun( double *finalObjValue = 0 );

    /// enable recording used by plotting functions below
    inline void storeHistory() { m_storeHistory = true; }

    /// plot sequence of evaluated objective values 
    void plotEvalHistory( const String &fileName ) const;

    /// plot sequence of best objective values (values that were the best at the time they were evaluated)
    void plotBestHistory( const String &fileName ) const;

    /// plot all 1D slices of the objective function, along with evaluation history
    void plotObjectiveSlices( const VectorD &center, int pointsPerDim, const String &filePrefix ) const;

    /// plot a 1D slice of the objective function, along with evaluation history
    void plotObjectiveSlice( const VectorD &center, int pointsPerDim, int dim, const String &fileName ) const;

    /// plot a 2D map of the objective function, along with evaluation history
    void plotObjectiveMap( const VectorD &center, int dim1, int dim2, int points1, int points2, int scaleFactor, const String &fileName ) const;

    /// display statistics about optimization run
    void dispHistoryStats( int indent ) const;

protected:

    /// evaluate objective function at given point, applying a penality to out-of-bounds points
    double evalWithPenalty( const VectorD &point );

    // represents objective function being optimized
    Objective &m_objective;

    // the optimization starting point
    VectorD m_startPoint;

    // optimization bounds; points outside the bounds receive a penality
    VectorD m_lBound;
    VectorD m_uBound;

    // set scale factor for penalty applied to points outside bounds
    double m_penaltyFactor;

    // termination tolerance (relative difference between objective values)
    double m_finalTolerance;

    // history of best and eval points
    bool m_storeHistory;
    Array<VectorD> m_bestPoints;
    VectorD m_bestValues;
    Array<VectorD> m_evalPoints;
    VectorD m_evalValues;

private:

    // disable copy constructor and assignment operator
    Optimizer( const Optimizer &x );
    Optimizer &operator=( const Optimizer &x );
};


//----------------------------------
// SIMPLEX OPTIMIZER CLASS
//----------------------------------


/// The SimplexOptimizer class performs a downhill simplex optimization (Nelder and Mead, 1965)
class SimplexOptimizer : public Optimizer {
public:

    // basic constructor
    SimplexOptimizer( Objective &objective ) : Optimizer( objective ) { m_startFrac = 0.1; m_reductionCount = 0; }
    virtual ~SimplexOptimizer() {}

    /// set the fracion of the search space (defined by bounds) to use for the initial simplex
    void setStartFrac( double startFrac ) { m_startFrac = startFrac; }

    /// run the optimization until termination (or user cancel); returns best point found
    VectorD run( double *finalObjValue  = 0 );

private:

    /// perform a single step of the optimization
    void step();

    // the fracion of the search space (defined by bounds) to use for the initial simplex
    double m_startFrac;

    // the current simplex
    Array<VectorD> m_points;

    // the objective funciton value associated with each point in the current simplex
    VectorD m_values;

    // number of consecutive reduction steps
    int m_reductionCount;

    // disable copy constructor and assignment operator
    SimplexOptimizer( const SimplexOptimizer &x );
    SimplexOptimizer &operator=( const SimplexOptimizer &x );
};


} // end namespace hb
#endif // _SBL_OPTIMIZER_H_

