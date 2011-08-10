// Licensed under MIT license; see license.txt.

#include <sbl/math/Optimizer.h>
#include <sbl/core/Command.h> // for checkCommandEvents
#include <sbl/core/StringUtil.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/image/ImageUtil.h> // for diagnostic plotting
#include <sbl/image/ImageTransform.h> // for diagnostic plotting
#include <sbl/other/Plot.h> // for diagnostic plotting
namespace sbl {


//-------------------------------------------
// OPTIMIZER CLASS
//-------------------------------------------


// basic constructor
Optimizer::Optimizer( Objective &objective ) : m_objective( objective ) { 
    m_penaltyFactor = 1.0; 
    m_finalTolerance = 1e-6; 
    m_storeHistory = false;
}


/// repeatedly run the optimization until termination (or user cancel); returns best point found
VectorD Optimizer::repeatRun( double *finalObjValue ) {
    double bestValue = 0;
    VectorD bestPoint;

    // iterate until convergence or cancel
    for (int iter = 0; iter < 100; iter++) {

        // run the optimization
        double runValue = 0;
        VectorD runPoint = run( &runValue );
        disp( 1, "iter: %d, result value: %f", iter, runValue );

        // if new result is better, continue from that point
        if (runValue < bestValue || bestPoint.length() == 0) {
            bestValue = runValue;
            bestPoint = runPoint;
            setStart( runPoint );
        } else {
            break;
        }

        // if new result is within final tolerance of previous best, stop
        if ((runValue - bestValue) / ((runValue + bestValue) * 0.5f + 1e-8) < m_finalTolerance) {
            break;
        }

        // check for user cancel
        if (checkCommandEvents())
            break;
    }

    // return best point found
    if (finalObjValue)
        *finalObjValue = bestValue;
    return bestPoint;
}


/// plot sequence of evaluated objective values 
void Optimizer::plotEvalHistory( const String &fileName ) const {
    aptr<Plot> plot( new Plot );
    plot->add( m_evalValues );
    plot->save( fileName );
    Plot::disp( plot );
}


/// plot sequence of best objective values (values that were the best at the time they were evaluated)
void Optimizer::plotBestHistory( const String &fileName ) const {
    aptr<Plot> plot( new Plot );
    plot->add( m_bestValues );
    plot->save( fileName );
    Plot::disp( plot );
}


/// plot all 1D slices of the objective function, along with evaluation history
void Optimizer::plotObjectiveSlices( const VectorD &center, int pointsPerDim, const String &filePrefix ) const {
    int dimCount = m_startPoint.length();
    for (int i = 0; i < dimCount; i++) {
        String fileName = filePrefix + sprintF( "-%d.svg", i );
        plotObjectiveSlice( center, pointsPerDim, i, fileName );
    }
}


/// plot a 1D slice of the objective function, along with evaluation history
void Optimizer::plotObjectiveSlice( const VectorD &center, int pointsPerDim, int dim, const String &fileName ) const {
    assertAlways( pointsPerDim > 1 );

    // compute objective values
    VectorD query( center );
    VectorD objPoints;
    VectorD objValues;
    double min = m_lBound[ dim ];
    double max = m_uBound[ dim ];
    for (int i = 0; i < pointsPerDim; i++) {
        double v = min + (double) i / (double) (pointsPerDim - 1) * (max - min);
        query[ dim ] = v;
        objPoints.append( v );
        objValues.append( m_objective.eval( query ) );
    }

    // extract history coords along desired dim
    VectorD bestPoints( m_bestPoints.count() );
    for (int i = 0; i < m_bestPoints.count(); i++) 
        bestPoints[ i ] = m_bestPoints[ i ][ dim ];
    VectorD evalPoints( m_evalPoints.count() );
    for (int i = 0; i < m_evalPoints.count(); i++) 
        evalPoints[ i ] = m_evalPoints[ i ][ dim ];

    // create plot
    Plot plot;
    plot.add( objPoints, objValues );
    plot.setStyle( PLOT_DOTS );
    plot.setColor( 255, 200, 0 );
    plot.add( evalPoints, m_evalValues );
    plot.setColor( 255, 0, 0 );
    plot.add( bestPoints, m_bestValues );

    // save plot
    plot.save( fileName );
}


/// plot a 2D map of the objective function, along with evaluation history
void Optimizer::plotObjectiveMap( const VectorD &center, int dim1, int dim2, int points1, int points2, int scaleFactor, const String &fileName ) const {
    assertAlways( points1 > 4 );
    assertAlways( points2 > 4 );

    // get bounds for objective map
    double min1 = m_lBound[ dim1 ];
    double min2 = m_lBound[ dim2 ];
    double max1 = m_uBound[ dim1 ];
    double max2 = m_uBound[ dim2 ];

    // create image
    VectorD query( center );
    ImageGrayF objValues( points1, points2 );
    for (int j = 0; j < points2; j++) {
        for (int i = 0; i < points1; i++) {
            query[ dim1 ] = min1 + (double) i / (double) (points1 - 1) * (max1 - min1);
            query[ dim2 ] = min2 + (double) j / (double) (points2 - 1) * (max2 - min2);
            objValues.data( i, j ) = (float) m_objective.eval( query );
        }
    }

    // resize image
    int size1 = points1 * scaleFactor;
    int size2 = points2 * scaleFactor;
    aptr<ImageGrayF> objLarge = resize( objValues, size1, size2, true );
    aptr<ImageGrayU> objGray = toUChar( *objLarge );
    aptr<ImageColorU> objColor = toColor( *objGray );
    
    // add history points
    for (int i = 0; i < m_evalPoints.count(); i++) {
        double p1 = m_evalPoints[ i ][ dim1 ];
        double p2 = m_evalPoints[ i ][ dim2 ];
        // fix(later): this is probably off by scaleFactor/2 or so
        int x = round( (p1 - min1) / (max1 - min1) * (double) size1 );
        int y = round( (p2 - min2) / (max2 - min2) * (double) size2 );
        if (objColor->inBounds( x, y )) 
            objColor->setRGB( x, y, 255, 200, 0 );
    }
    for (int i = 0; i < m_bestPoints.count(); i++) {
        double p1 = m_bestPoints[ i ][ dim1 ];
        double p2 = m_bestPoints[ i ][ dim2 ];
        // fix(later): this is probably off by scaleFactor/2 or so
        int x = round( (p1 - min1) / (max1 - min1) * (double) size1 );
        int y = round( (p2 - min2) / (max2 - min2) * (double) size2 );
        if (objColor->inBounds( x, y )) 
            objColor->setRGB( x, y, 255, 0, 0 );
    }

    // save image
    saveImage( *objColor, fileName );
}


/// display statistics about optimization run
void Optimizer::dispHistoryStats( int indent ) const {
    disp( indent, "eval count: %d, eval min: %f, eval max: %f", m_evalValues.length(), m_evalValues.min(), m_evalValues.max() );
    disp( indent, "best count: %d, best min: %f, best max: %f", m_bestValues.length(), m_bestValues.min(), m_bestValues.max() );
}


/// evaluate objective function at given point, applying a penality to out-of-bounds points
double Optimizer::evalWithPenalty( const VectorD &point ) {

    // compute objective function value
    double obj = m_objective.eval( point );

    // add penalty
    int dim = point.length();
    for (int j = 0; j < dim; j++) {
        double v = point[ j ];
        if (v < m_lBound[ j ]) 
            obj += m_penaltyFactor * (m_lBound[ j ] - v) / (m_uBound[ j ] - m_lBound[ j ]);
        if (v > m_uBound[ j ])
            obj += m_penaltyFactor * (v - m_uBound[ j ]) / (m_uBound[ j ] - m_lBound[ j ]);
    }

    // if requested, add to history
    if (m_storeHistory) {
        m_evalPoints.appendCopy( point );
        m_evalValues.append( obj );
        if (m_bestValues.length() == 0 || obj < m_bestValues.endValue()) {
            m_bestPoints.appendCopy( point );
            m_bestValues.append( obj );
        }
    }
    return obj;
}


//-------------------------------------------
// SIMPLEX OPTIMIZER CLASS
//-------------------------------------------


/// run the optimization until termination (or user cancel); returns best point found
VectorD SimplexOptimizer::run( double *finalObjValue ) {

    // create initial simplex
    int dim = m_startPoint.length();
    m_points.appendCopy( m_startPoint );
    for (int i = 0; i < dim; i++) {
        VectorD *point = new VectorD( m_startPoint );
        (*point)[ i ] += (m_uBound[ i ] - m_lBound[ i ]) * m_startFrac;
        m_points.append( point );
    }

    // evaluate each point in initial simplex
    for (int i = 0; i < m_points.count(); i++) {
        m_values.append( evalWithPenalty( m_points[ i ] ) );
    }

    // iterate until convergence or user cancel
    bool done = false;
    while (done == false) {

        // run one optimization step
        step();

        // check for convergence
        double minValue = m_values.min();
        double maxValue = m_values.max();
        if ((maxValue - minValue) / ((minValue + maxValue) * 0.5f + 1e-8) < m_finalTolerance) 
            break;;

        // check for getting stuck (results in a series of reductions)
        // fix(later): fix this problem
        if (m_reductionCount > 10) {
            warning( "SimplexOptimizer: repeated reductions; halting optimization" );
            break;
        }

        // check for user cancel
        if (checkCommandEvents())
            break;
    }

    // return best point found
    int bestIndex = argMin( m_values );
    if (finalObjValue)
        *finalObjValue = m_values[ bestIndex ];
    return m_points[ bestIndex ];
}


/// perform a single step of the optimization (implementation based on Wikipedia description)
void SimplexOptimizer::step() {

    // data dimensions
    int count = m_points.count();
    assertDebug( count );
    int dim = m_points[ 0 ].length();
    assertDebug( count == dim + 1 );

    // sort the points by value
    VectorI sortInd = sortIndex( m_values );
    Array<VectorD> newPoints;
    VectorD newValues( count );
    for (int i = 0; i < count; i++) {
        int index = sortInd[ i ];
        newValues[ i ] = m_values[ index ];
        newPoints.appendCopy( m_points[ index ] );
    }
    m_points = newPoints;
    m_values = newValues;
    const VectorD &worstPoint = m_points[ count - 1 ];

    // compute centroid of first N - 1 points
    VectorD center( count );
    center.clear( 0 );
    for (int i = 0; i < count - 1; i++) {
        for (int j = 0; j < dim; j++)
            center[ j ] += m_points[ i ][ j ];
    }
    double factor = 1.0 / (double) (count - 1);
    for (int j = 0; j < dim; j++) {
        center[ j ] *= factor;
    }

    // compute point by reflecting across centroid
    VectorD reflectedPoint( dim );
    for (int j = 0; j < dim; j++) 
        reflectedPoint[ j ] = center[ j ] + (center[ j ] - worstPoint[ j ]);
    double reflectedValue = evalWithPenalty( reflectedPoint );
    if (m_values[ 0 ] <= reflectedValue && reflectedValue < m_values[ count - 2 ]) {
        m_values[ count - 1 ] = reflectedValue;
        m_points[ count - 1 ] = reflectedPoint;
        m_reductionCount = 0;
        return;
    }

    // if reflected is best so far, compute expanded point
    if (reflectedValue < m_values[ 0 ]) {
        VectorD expandedPoint( dim );
        for (int j = 0; j < dim; j++) {
            expandedPoint[ j ] = center[ j ] + 2.0 * (center[ j ] - worstPoint[ j ]);
        }
        double expandedValue = evalWithPenalty( expandedPoint );
        if (expandedValue < reflectedValue) {
            m_values[ count - 1 ] = expandedValue;
            m_points[ count - 1 ] = expandedPoint;
        } else {
            m_values[ count - 1 ] = reflectedValue;
            m_points[ count - 1 ] = reflectedPoint;
        }
        m_reductionCount = 0;
        return;
    }

    // next try the contracted point
    VectorD contractedPoint( dim );
    for (int j = 0; j < dim; j++) {
        contractedPoint[ j ] = worstPoint[ j ] + 0.5 * (center[ j ] - worstPoint[ j ]);
    }
    double contractedValue = evalWithPenalty( contractedPoint );
    if (contractedValue < m_values[ count - 1 ]) {
        m_values[ count - 1 ] = contractedValue;
        m_points[ count - 1 ] = contractedPoint;
        m_reductionCount = 0;
        return;        
    }

    // if none of these produced a better point, perform reduction
    const VectorD &bestPoint = m_points[ 0 ];
    for (int i = 1; i < count; i++) {
        for (int j = 0; j < dim; j++) {
            m_points[ i ][ j ] = bestPoint[ j ] + 0.5f * (m_points[ i ][ j ] - bestPoint[ j ]);
        }
    }
    m_reductionCount++;
}


} // end namespace sbl

