function precisionMatrix = computePrecisionOfKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, B )
%
% Compute 
%   H = W' * diag( B ) * W 
% where 
%   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
% and B is a weight matrix


numberOfDimensions = length( kroneckerProductBasisFunctions );

% Compute a new set of basis functions (point-wise product of each combination of pairs) so that we can
% easily compute a mangled version of the result
Ms = zeros( 1, numberOfDimensions ); % Number of basis functions in each dimension
hessianKroneckerProductBasisFunctions = cell( 0, 0 );
for dimensionNumber = 1 : numberOfDimensions
  M = size( kroneckerProductBasisFunctions{ dimensionNumber }, 2 );
  
  if 0
    N = size( kroneckerProductBasisFunctions{ dimensionNumber }, 1 );
    
    A = kroneckerProductBasisFunctions{ dimensionNumber };
    A2 = zeros( N, M^2 );
    for i = 1 : M
      shapeI = A( :, i );
      for j = 1 : M
        shapeJ = A( :, j );
        A2( :, i + (j-1) * M ) = shapeI .* shapeJ;
      end
    end

    hessianKroneckerProductBasisFunctions{ dimensionNumber } = A2;
  else
    A = kroneckerProductBasisFunctions{ dimensionNumber };
    hessianKroneckerProductBasisFunctions{ dimensionNumber } = kron( ones( 1, M ), A ) .* kron( A, ones( 1, M ) );
  end
  
  Ms( dimensionNumber ) = M;
  
end
result = projectKroneckerProductBasisFunctions( hessianKroneckerProductBasisFunctions, B );
%result = reshape( result, [ Ms.^2 ] );
result = reshape( result, kron( Ms, [ 1 1 ] ) );
permutationIndices = [ 2 * [ 1 : numberOfDimensions ]-1 2 * [ 1 : numberOfDimensions ] ];
result = permute( result, permutationIndices );
precisionMatrix = reshape( result, [ prod( Ms ) prod( Ms ) ] );
