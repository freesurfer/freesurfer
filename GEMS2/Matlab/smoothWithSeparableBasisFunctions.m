function [ smoothedData, C ] = smoothWithSeparableBasisFunctions( data, numberOfBasisFunctions, mask, maximumNumberOfBasisFunctions )
%
%

persistent cache
if ( nargin == 0 )
  disp( 'Usage: data = smoothWithSeparableBasisFunctions( data, numberOfBasisFunctions, mask, maximumNumberOfBasisFunctions )' )
  disp( '' )
  disp( '' )
  disp( 'Let''s test drive this thing!' )
  testDriveOurselves();
  return

end
if ( nargin > 2 )
  % Calculate everything related to the setup
  [ Nx Ny Nz Ax Ay Az hessian ] = calculateEverythingForSetup( mask, maximumNumberOfBasisFunctions );

  % Cholesky decomposition
  R = chol( hessian );

  % Put it in the cache
  if ( exist( 'cache' ) ~= 1 )
    cache = struct();
  end
  cache.mask = mask;
  cache.Nx = Nx;
  cache.Ny = Ny;
  cache.Nz = Nz;
  cache.maximumNumberOfBasisFunctions = maximumNumberOfBasisFunctions;
  cache.hessian = hessian;
  cache.resolution{maximumNumberOfBasisFunctions} = ...
       struct( 'Ax', Ax, 'Ay', Ay, 'Az', Az, 'R', R, 'Rtranspose', R' );



end
if ( nargin == 1 )
  % We assume mask has already been cached. Test if that's
  % indeed the case
  if ( exist( 'cache' ) ~= 1 )
    error( 'Can''t use cache because nothing there yet' )
  end
  if ( ~isfield( cache, 'Nx' ) )
    error( 'Can''t use cache because nothing there yet' )
  end

  % Use maximum number of basis functions available
  numberOfBasisFunctions = cache.maximumNumberOfBasisFunctions;

end

% Ready for the real work. Make sure we having things in cache because we're gonna use it
if ( exist( 'cache' ) ~= 1 )
  error( 'Can''t use cache because nothing there yet' )
end

disp( [ '   Smoothing with ' num2str( numberOfBasisFunctions ) ' basis functions ' ] );
tic

if ~isstruct( cache.resolution{ numberOfBasisFunctions } )
  % We haven't done this number of basis functions before. Pre-compute
  % stuff we will need and put it in the cache
  disp( [ '      Caching stuff for ' num2str( numberOfBasisFunctions ) ' basis functions' ] )

  % 1-D basis functions
  maximumNumberOfBasisFunctions = cache.maximumNumberOfBasisFunctions;
  Ax = cache.resolution{maximumNumberOfBasisFunctions}.Ax( :, 1 : numberOfBasisFunctions );
  Ay = cache.resolution{maximumNumberOfBasisFunctions}.Ay( :, 1 : numberOfBasisFunctions );
  Az = cache.resolution{maximumNumberOfBasisFunctions}.Az( :, 1 : numberOfBasisFunctions );

  % Cholesky decomposition
  C = zeros( maximumNumberOfBasisFunctions, maximumNumberOfBasisFunctions, maximumNumberOfBasisFunctions );
  C( 1 : numberOfBasisFunctions, 1 : numberOfBasisFunctions, 1 : numberOfBasisFunctions ) = 1;
  playingIndices = find( C );
  hessian = cache.hessian;
  hessian = hessian( playingIndices, : );
  hessian = hessian( :, playingIndices );
  R = chol( hessian );

  % Store away things
  cache.resolution{numberOfBasisFunctions} = ...
       struct( 'Ax', Ax, 'Ay', Ay, 'Az', Az, 'R', R, 'Rtranspose', R' );

end



% Calculate coefficients using seperability of basis functions to speed things up
Nx = cache.Nx;
Ny = cache.Ny;
Nz = cache.Nz;
%disp( '      Pushing the masked data through the basis functions' )
%tic
tmp = cache.resolution{numberOfBasisFunctions}.Ax' * reshape( data, [ Nx Ny*Nz ] );
tmp = permute( reshape( tmp, [ numberOfBasisFunctions Ny Nz ] ), [ 2 1 3 ] );
tmp = cache.resolution{numberOfBasisFunctions}.Ay' * reshape( tmp, [ Ny numberOfBasisFunctions*Nz ] );
tmp = permute( reshape( tmp, [ numberOfBasisFunctions numberOfBasisFunctions Nz ] ), [ 2 1 3 ] );
tmp = reshape( tmp, [ numberOfBasisFunctions^2 Nz ] );
C = reshape( tmp * cache.resolution{numberOfBasisFunctions}.Az, [ numberOfBasisFunctions numberOfBasisFunctions numberOfBasisFunctions ] );
%toc

% Need to still correct for non-orthonormality of the basis functions within the mask by inverting the hessian. If the mask was one everywhere, the hessian would be unity matrix (since wer're using DCT basis functions) so this step wouldn't be needed.
%disp( '      Inverting the Hessian' )
%tic
C(:) = cache.resolution{numberOfBasisFunctions}.R \ ...
       ( cache.resolution{numberOfBasisFunctions}.Rtranspose \ C(:) );
%toc

% Expand back into image domain
% (this code is actually the same as before, but consistenlty swap A' with A, and numberOfBasisFunctions
% with N)
%disp( '      Expanding back into image domain' )
%tic
tmp = cache.resolution{numberOfBasisFunctions}.Ax * reshape( C, [ numberOfBasisFunctions numberOfBasisFunctions^2 ] );
tmp = permute( reshape( tmp, [ Nx numberOfBasisFunctions numberOfBasisFunctions ] ), [ 2 1 3 ] );
tmp = cache.resolution{numberOfBasisFunctions}.Ay * reshape( tmp, [ numberOfBasisFunctions Nx*numberOfBasisFunctions ] );
tmp = permute( reshape( tmp, [ Ny Nx numberOfBasisFunctions ] ), [ 2 1 3 ] );
tmp = reshape( tmp, [ Nx*Ny numberOfBasisFunctions ] );
smoothedData = reshape( tmp * cache.resolution{numberOfBasisFunctions}.Az', [ Nx Ny Nz ] );

%
smoothedData = smoothedData .* cache.mask;
elapsedTimeSmoothing = toc;
disp( [ '       ' num2str( elapsedTimeSmoothing ) ] )

end


%
%
%
function  [ Nx Ny Nz Ax Ay Az hessian ] = calculateEverythingForSetup( mask, maximumNumberOfBasisFunctions )
%
%

  % Get dimensions from mask
  [ Nx Ny Nz ] = size( mask );

  % Construct DCT basis functions Let's use the definition of Matlab: 
  % http://www.mathworks.com/access/helpdesk/help/toolbox/signal/dct.html
  disp( '      Building DCT basis functions' )
  numberOfBasisFunctions = maximumNumberOfBasisFunctions;
  
  tmpx = ( 2 * [ 1 : Nx ]' - 1 ) * pi / 2 / Nx;
  Ax = cos( tmpx * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / Nx );
  Ax( :, 1 ) = Ax( :, 1 ) / sqrt( 2 );
  
  tmpy = ( 2 * [ 1 : Ny ]' - 1 ) * pi / 2 / Ny;
  Ay = cos( tmpy * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / Ny );
  Ay( :, 1 ) = Ay( :, 1 ) / sqrt( 2 );
  
  tmpz = ( 2 * [ 1 : Nz ]' - 1 ) * pi / 2 / Nz;
  Az = cos( tmpz * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / Nz );
  Az( :, 1 ) = Az( :, 1 ) / sqrt( 2 );
  
  
  
  % Build Hessian as fast as possible by exploiting the separability of the basis functions
  disp( '      Computing Hessian' )
  tic
  A2x = zeros( Nx, numberOfBasisFunctions^2 );
  for i = 1 : numberOfBasisFunctions
    shapeI = Ax( :, i );
    for j = 1 : numberOfBasisFunctions
      shapeJ = Ax( :, j );
      A2x( :, i + (j-1) * numberOfBasisFunctions ) = shapeI .* shapeJ;
    end
  end
  
  A2y = zeros( Ny, numberOfBasisFunctions^2 );
  for i = 1 : numberOfBasisFunctions
    shapeI = Ay( :, i );
    for j = 1 : numberOfBasisFunctions
      shapeJ = Ay( :, j );
      A2y( :, i + (j-1) * numberOfBasisFunctions ) = shapeI .* shapeJ;
    end
  end
  
  A2z = zeros( Nz, numberOfBasisFunctions^2 );
  for i = 1 : numberOfBasisFunctions
    shapeI = Az( :, i );
    for j = 1 : numberOfBasisFunctions
      shapeJ = Az( :, j );
      A2z( :, i + (j-1) * numberOfBasisFunctions ) = shapeI .* shapeJ;
    end
  end
  
  % Do the operations direction-wise
  tmp = A2x' * reshape( mask, [ Nx Ny*Nz ] );
  tmp = permute( reshape( tmp, [ numberOfBasisFunctions^2 Ny Nz ] ), [ 2 1 3 ] );
  tmp = A2y' * reshape( tmp, [ Ny numberOfBasisFunctions^2*Nz ] );
  tmp = permute( reshape( tmp, [ numberOfBasisFunctions^2 numberOfBasisFunctions^2 Nz ] ), [ 2 1 3 ] );
  tmp = reshape( tmp, [ numberOfBasisFunctions^4 Nz ] );
  Gamma = reshape( tmp * A2z, [ numberOfBasisFunctions^2 numberOfBasisFunctions^2 numberOfBasisFunctions^2 ] );
  
  % The matrix Gamma above contains the right entries, but in the incorrect place. Let's construct the Hessian by reordering the elements. The correction of the indices is given by:
  %
  %   ( n - m' ) * M + ( p - n ) * M^2 + ( m' - n' ) * M^3 + ( n' - p ) * M^4
  %
  % where the needed (n,m',p,n') values are calculated w.r.t. Gamma
  %
  n = repmat( [ 1 : numberOfBasisFunctions ], ...
              [ numberOfBasisFunctions^2 numberOfBasisFunctions ] );
  n = repmat( n, [ 1 1 numberOfBasisFunctions^2 ] );
  
  mAccent = repmat( [ 1 : numberOfBasisFunctions ], [ numberOfBasisFunctions 1 ] );
  mAccent = repmat( mAccent(:), [ 1 numberOfBasisFunctions^2 ] );
  mAccent = repmat( mAccent, [ 1 1 numberOfBasisFunctions^2 ] );
  
  nAccent = repmat( [ 1 : numberOfBasisFunctions ], [ numberOfBasisFunctions 1 ] );
  nAccent = repmat( ( nAccent(:) )', [ numberOfBasisFunctions^2 1 ] );
  nAccent = repmat( nAccent, [ 1 1 numberOfBasisFunctions^2 ] );
  
  p = repmat( [ 1 : numberOfBasisFunctions ], ...
              [ numberOfBasisFunctions^2 numberOfBasisFunctions ] );
  p = reshape( p, [ numberOfBasisFunctions^2 1 numberOfBasisFunctions^2 ] );
  p = repmat( p, [ 1 numberOfBasisFunctions^2 1 ] );
  
  correctionNeeded = ( n(:) - mAccent(:) ) * numberOfBasisFunctions + ...
                    ( p(:) - n(:) ) * numberOfBasisFunctions^2 + ...
                    ( mAccent(:) - nAccent(:) ) * numberOfBasisFunctions^3 + ...
                    ( nAccent(:) - p(:) ) * numberOfBasisFunctions^4;
  
  hessian = zeros( numberOfBasisFunctions^3, numberOfBasisFunctions^3 );
  hessian( [ 1 : numberOfBasisFunctions^6 ]' + correctionNeeded ) = Gamma(:);
  
  elapsedTimeSmoothing = toc;
  disp( [ '       ' num2str( elapsedTimeSmoothing ) ] )

  
  testValidityOfHessian = false; % For debugging purposes, explicitly build 3-D basis functions and compute Hessian from those as a reference. Don't do this for more than 6 or more basis functions because it's just too slow (e.g., 10 never returns)
  if testValidityOfHessian
    disp( 'Computing Hessian using explicit 3-D basis functions...' )
    tic
  
    basisFunctions = zeros( Nx*Ny*Nz, numberOfBasisFunctions^3 );
    counter = 1;
    for i = 1 : numberOfBasisFunctions
    shapeI = Ax( :, i );
    for j = 1 : numberOfBasisFunctions
    shapeJ = Ay( :, j );
    shapeIJ = shapeI * shapeJ';
    for k = 1 : numberOfBasisFunctions
    shapeK = Az( :, k );
    shapeIJK = reshape( shapeIJ( : ) * shapeK', [ Nx Ny Nz ] );
    basisFunctions( :, ...
                    i + ...
                    ( j-1 ) * numberOfBasisFunctions + ...
                    ( k - 1 ) * numberOfBasisFunctions^2 ) = shapeIJK( : );
    end
    end
    end
    basisFunctions = basisFunctions( find( mask ), : );
    directHessian = basisFunctions' * basisFunctions;
    toc
  
    maxError = max( abs( hessian(:) - directHessian(:) ) )
  
  end
  
  
end


%
%
%
function  testDriveOurselves()
%
%
  % 
  close all
  clear all

  % Let's generate some data
  Nx = 62
  Ny = 64
  Nz = 61
  maximumNumberOfBasisFunctions = 10
  data = randn( Nx, Ny, Nz );

  % Create some mask
  center = [ Nx/2 Ny/2 Nz/2 ]';
  radius = [ Nx*.5 Ny*.5 Nz*.5 ]';
  [ x1, x2, x3 ] = ndgrid( 1:Nx, 1:Ny, 1:Nz );
  distance = ( ( x1 - center( 1 ) ) / radius( 1 ) ).^2 + ...
             ( ( x2 - center( 2 ) ) / radius( 2 ) ).^2 + ...
             ( ( x3 - center( 3 ) ) / radius( 3 ) ).^2;
  mask = distance <= 1;

  % Show what we have
  figure
  data = data .* mask;
  showThreeDImage( data );


  numberOfBasisFunctions = 1;
  smoothedData = smoothWithSeparableBasisFunctions( data, ...
                 numberOfBasisFunctions, mask, maximumNumberOfBasisFunctions );
  for numberOfBasisFunctions = 1 : maximumNumberOfBasisFunctions
    disp( '========================================' )
    disp( [ num2str( numberOfBasisFunctions ) ' basis functions' ] )
    disp( '========================================' )
    for i = 1 : 4
      disp( '  +++++' )
      disp( [ '  run ' num2str( i ) ] )
      disp( '  +++++' )
      smoothedData = smoothWithSeparableBasisFunctions( data, numberOfBasisFunctions );
    end

    % Show result
    figure
    showThreeDImage( smoothedData );
    title( [ num2str( numberOfBasisFunctions ) ' basis functions' ] )
  end


end

