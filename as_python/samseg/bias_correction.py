import numpy as np


def backprojectKroneckerProductBasisFunctions(kroneckerProductBasisFunctions, coefficients):
    #
    # numberOfDimensions = length( kroneckerProductBasisFunctions );
    numberOfDimensions = len(kroneckerProductBasisFunctions)
    # Ms = zeros( 1, numberOfDimensions ); % Number of basis functions in each dimension
    Ms = np.zeros(numberOfDimensions)  # Number of basis functions in each dimension
    # Ns = zeros( 1, numberOfDimensions ); % Number of data points in each dimension
    Ns = np.zeros(numberOfDimensions)  # Number of basis functions in each dimension
    # transposedKroneckerProductBasisFunctions = cell( 0, 0 );
    transposedKroneckerProductBasisFunctions = []
    # for dimensionNumber = 1 : numberOfDimensions
    for dimensionNumber in range(numberOfDimensions):
        #   Ms( dimensionNumber ) = size( kroneckerProductBasisFunctions{ dimensionNumber }, 2 );
        Ms[dimensionNumber] = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
        #   Ns( dimensionNumber ) = size( kroneckerProductBasisFunctions{ dimensionNumber }, 1 );
        Ns[dimensionNumber] = kroneckerProductBasisFunctions[dimensionNumber].shape[0]
        #   transposedKroneckerProductBasisFunctions{ dimensionNumber } = kroneckerProductBasisFunctions{ dimensionNumber }';
        transposedKroneckerProductBasisFunctions.append(kroneckerProductBasisFunctions[dimensionNumber].T)
        # end
        #
        #
    # y = projectKroneckerProductBasisFunctions( transposedKroneckerProductBasisFunctions, reshape( coefficients, Ms ) );
    y = projectKroneckerProductBasisFunctions(transposedKroneckerProductBasisFunctions, coefficients.reshape(Ms) )
    # Y = reshape( y, Ns );
    Y = y.reshape(Ns)
    return Y


def projectKroneckerProductBasisFunctions(kroneckerProductBasisFunctions, T):
    # %
    # % Compute
    # %   c = W' * t
    # % where
    # %   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
    # % and
    # %   t = T( : )
    #
    numberOfDimensions = len(kroneckerProductBasisFunctions)
    # currentSizeOfT = size( T );
    currentSizeOfT = list(T.shape)
    # for dimensionNumber = 1 : numberOfDimensions
    for dimensionNumber in range(numberOfDimensions):
        #   % Reshape into 2-D, do the work in the first dimension, and shape into N-D
        #   T = reshape( T, currentSizeOfT( 1 ), [] );
        T = T.reshape((currentSizeOfT[0], -1))
        #   T = ( kroneckerProductBasisFunctions{ dimensionNumber } )' * T;
        T = ( kroneckerProductBasisFunctions[dimensionNumber] ).T @ T
        #   currentSizeOfT( 1 ) = size( kroneckerProductBasisFunctions{ dimensionNumber }, 2 );
        currentSizeOfT[0] = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
        #   T = reshape( T, currentSizeOfT );
        T = T.reshape(currentSizeOfT)
        #
        #   % Shift dimension
        #   currentSizeOfT = [ currentSizeOfT( 2 : end ) currentSizeOfT( 1 ) ];
        currentSizeOfT = currentSizeOfT[1:] + [currentSizeOfT[0]]
        #   T = shiftdim( T, 1 );
        T = np.rollaxis(T, 0, 3)
        # end
    #
    # % Return result as vector
    # coefficients = T(:);
    coefficients = T.flatten()
    return coefficients


def computePrecisionOfKroneckerProductBasisFunctions(kroneckerProductBasisFunctions, B):
    # %
    # % Compute
    # %   H = W' * diag( B ) * W
    # % where
    # %   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
    # % and B is a weight matrix
    #
    #
    # numberOfDimensions = length( kroneckerProductBasisFunctions );
    numberOfDimensions = len( kroneckerProductBasisFunctions )
    #
    # % Compute a new set of basis functions (point-wise product of each combination of pairs) so that we can
    # % easily compute a mangled version of the result
    # Ms = zeros( 1, numberOfDimensions ); % Number of basis functions in each dimension
    Ms = np.zeros( numberOfDimensions ) # % Number of basis functions in each dimension
    # hessianKroneckerProductBasisFunctions = cell( 0, 0 );
    hessianKroneckerProductBasisFunctions = {}
    # for dimensionNumber = 1 : numberOfDimensions
    for dimensionNumber in range(numberOfDimensions):
        #   M = size( kroneckerProductBasisFunctions{ dimensionNumber }, 2 );
        M = kroneckerProductBasisFunctions[dimensionNumber].shape[1]
        #
        #   if 0
        #     N = size( kroneckerProductBasisFunctions{ dimensionNumber }, 1 );
        #
        #     A = kroneckerProductBasisFunctions{ dimensionNumber };
        #     A2 = zeros( N, M^2 );
        #     for i = 1 : M
        #       shapeI = A( :, i );
        #       for j = 1 : M
        #         shapeJ = A( :, j );
        #         A2( :, i + (j-1) * M ) = shapeI .* shapeJ;
        #       end
        #     end
        #
        #     hessianKroneckerProductBasisFunctions{ dimensionNumber } = A2;
        #   else
        #     A = kroneckerProductBasisFunctions{ dimensionNumber };
        A = kroneckerProductBasisFunctions[dimensionNumber]
        #     hessianKroneckerProductBasisFunctions{ dimensionNumber } = kron( ones( 1, M ), A ) .* kron( A, ones( 1, M ) );
        hessianKroneckerProductBasisFunctions[dimensionNumber] = np.kron( np.ones( (1, M )), A ) * np.kron( A, np.ones( (1, M) ) )
        #   end
        #
        #   Ms( dimensionNumber ) = M;
        Ms[dimensionNumber] = M
    #
    # end
    # result = projectKroneckerProductBasisFunctions( hessianKroneckerProductBasisFunctions, B );
    result = projectKroneckerProductBasisFunctions( hessianKroneckerProductBasisFunctions, B );
    # %result = reshape( result, [ Ms.^2 ] );
    result = result.reshape(np.kron( Ms, [ 1, 1 ] ) )
    # permutationIndices = [ 2 * [ 1 : numberOfDimensions ]-1 2 * [ 1 : numberOfDimensions ] ];
    permutationIndices = np.hstack((2 * np.r_[: numberOfDimensions ], 2 * np.r_[: numberOfDimensions ] +1))
    # result = permute( result, permutationIndices );
    result = np.transpose(result, permutationIndices)
    # precisionMatrix = reshape( result, [ prod( Ms ) prod( Ms ) ] );
    precisionMatrix = result.reshape( ( np.prod( Ms ), np.prod( Ms ) ) )
    return precisionMatrix
