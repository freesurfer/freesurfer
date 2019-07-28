import numpy as np


def backprojectKroneckerProductBasisFunctions(kroneckerProductBasisFunctions, coefficients):
    numberOfDimensions = len(kroneckerProductBasisFunctions)
    Ms = np.zeros(numberOfDimensions, dtype=np.uint32)  # Number of basis functions in each dimension
    Ns = np.zeros(numberOfDimensions, dtype=np.uint32)  # Number of basis functions in each dimension
    transposedKroneckerProductBasisFunctions = []
    for dimensionNumber in range(numberOfDimensions):
        Ms[dimensionNumber] = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
        Ns[dimensionNumber] = kroneckerProductBasisFunctions[dimensionNumber].shape[0]
        transposedKroneckerProductBasisFunctions.append(kroneckerProductBasisFunctions[dimensionNumber].T)
    y = projectKroneckerProductBasisFunctions(transposedKroneckerProductBasisFunctions, coefficients.reshape(Ms, order='F') )
    Y = y.reshape(Ns, order='F')
    return Y


def projectKroneckerProductBasisFunctions(kroneckerProductBasisFunctions, T):
    #
    # Compute
    #   c = W' * t
    # where
    #   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
    # and
    #   t = T( : )
    numberOfDimensions = len(kroneckerProductBasisFunctions)
    currentSizeOfT = list(T.shape)
    for dimensionNumber in range(numberOfDimensions):
        # Reshape into 2-D, do the work in the first dimension, and shape into N-D
        T = T.reshape((currentSizeOfT[0], -1), order='F')
        T = ( kroneckerProductBasisFunctions[dimensionNumber] ).T @ T
        currentSizeOfT[0] = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
        T = T.reshape(currentSizeOfT, order='F')
        # Shift dimension
        currentSizeOfT = currentSizeOfT[1:] + [currentSizeOfT[0]]
        T = np.rollaxis(T, 0, 3)
    # Return result as vector
    coefficients = T.flatten(order='F')
    return coefficients


def computePrecisionOfKroneckerProductBasisFunctions(kroneckerProductBasisFunctions, B):
    #
    # Compute
    #   H = W' * diag( B ) * W
    # where
    #   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
    # and B is a weight matrix
    numberOfDimensions = len( kroneckerProductBasisFunctions )

    # Compute a new set of basis functions (point-wise product of each combination of pairs) so that we can
    # easily compute a mangled version of the result
    Ms = np.zeros( numberOfDimensions , dtype=np.uint32) # Number of basis functions in each dimension
    hessianKroneckerProductBasisFunctions = {}
    for dimensionNumber in range(numberOfDimensions):
        M = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
        A = kroneckerProductBasisFunctions[dimensionNumber]
        hessianKroneckerProductBasisFunctions[dimensionNumber] = np.kron( np.ones( (1, M )), A ) * np.kron( A, np.ones( (1, M) ) )
        Ms[dimensionNumber] = M
    result = projectKroneckerProductBasisFunctions( hessianKroneckerProductBasisFunctions, B )
    new_shape = list(np.kron( Ms, [ 1, 1 ] ))
    new_shape.reverse()
    result = result.reshape(new_shape)
    permutationIndices = np.hstack((2 * np.r_[: numberOfDimensions ], 2 * np.r_[: numberOfDimensions ] +1))
    result = np.transpose(result, permutationIndices)
    precisionMatrix = result.reshape( ( np.prod( Ms ), np.prod( Ms ) ) )
    return precisionMatrix

