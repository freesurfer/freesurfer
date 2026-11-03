function coefficients = projectKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, T )
%
% Compute 
%   c = W' * t 
% where 
%   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }
% and 
%   t = T( : )

numberOfDimensions = length( kroneckerProductBasisFunctions );
currentSizeOfT = size( T );
for dimensionNumber = 1 : numberOfDimensions
  % Reshape into 2-D, do the work in the first dimension, and shape into N-D
  T = reshape( T, currentSizeOfT( 1 ), [] );
  T = ( kroneckerProductBasisFunctions{ dimensionNumber } )' * T;
  currentSizeOfT( 1 ) = size( kroneckerProductBasisFunctions{ dimensionNumber }, 2 );
  T = reshape( T, currentSizeOfT ); 
  
  % Shift dimension 
  currentSizeOfT = [ currentSizeOfT( 2 : end ) currentSizeOfT( 1 ) ];
  T = shiftdim( T, 1 );
end

% Return result as vector
coefficients = T(:);
