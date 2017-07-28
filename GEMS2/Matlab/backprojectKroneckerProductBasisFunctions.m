function Y = backprojectKroneckerProductBasisFunctions( kroneckerProductBasisFunctions, coefficients )
%
% Compute 
%   y = W * c 
% where 
%   W = W{ numberOfDimensions } \kron W{ numberOfDimensions-1 } \kron ... W{ 1 }.
%
% The results is returned as a numberOfDimensions-dimensional matrix

numberOfDimensions = length( kroneckerProductBasisFunctions );
Ms = zeros( 1, numberOfDimensions ); % Number of basis functions in each dimension
Ns = zeros( 1, numberOfDimensions ); % Number of data points in each dimension
transposedKroneckerProductBasisFunctions = cell( 0, 0 );
for dimensionNumber = 1 : numberOfDimensions
  Ms( dimensionNumber ) = size( kroneckerProductBasisFunctions{ dimensionNumber }, 2 );
  Ns( dimensionNumber ) = size( kroneckerProductBasisFunctions{ dimensionNumber }, 1 );
  transposedKroneckerProductBasisFunctions{ dimensionNumber } = kroneckerProductBasisFunctions{ dimensionNumber }';
end

y = projectKroneckerProductBasisFunctions( transposedKroneckerProductBasisFunctions, reshape( coefficients, Ms ) );
Y = reshape( y, Ns );

