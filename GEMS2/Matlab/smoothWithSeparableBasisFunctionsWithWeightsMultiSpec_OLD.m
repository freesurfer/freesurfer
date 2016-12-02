function [ smoothedData, C ] = smoothWithSeparableBasisFunctionsWithWeightsMultiSpec( data, weights, numberOfBasisFunctions, maximumNumberOfBasisFunctions )
%
%

persistent cache
if ( nargin == 0 )
  disp( 'Usage: data = smoothWithSeparableBasisFunctionsWithWeights( data, weights, numberOfBasisFunctions, maximumNumberOfBasisFunctions )' )
  disp( '' )
  disp( '' )
  disp( 'Let''s test drive this thing!' )
  testDriveOurselves();
  return

end
if ( nargin == 4 )
  % Pre-calculate everything there is to pre-calculate at full resolution
  [ Nx Ny Nz numberOfImages ] = size( data );
  numberOfImages = sqrt(numberOfImages);
  [ Ax Ay Az A2x A2y A2z correctionNeeded ] = ...
             calculateEverythingForThisResolution( Nx, Ny, Nz, ...
                                                   maximumNumberOfBasisFunctions );
  % Put it in the cache
  if ( exist( 'cache' ) ~= 1 )
    cache = struct();
  end
  cache.Nx = Nx;
  cache.Ny = Ny;
  cache.Nz = Nz;
  cache.numberOfImages = numberOfImages;
  cache.maximumNumberOfBasisFunctions = maximumNumberOfBasisFunctions;
  cache.resolution{maximumNumberOfBasisFunctions} = ...
       struct( 'Ax', Ax, 'Ay', Ay, 'Az', Az, ...
               'A2x', A2x, 'A2y', A2y, 'A2z', A2z, ...
               'correctionNeeded', correctionNeeded );

end
if ( nargin == 1 )
  weights = ones( size( data ) );
end
if ( nargin == 2 )
  % We assume the numberOfBasisFunctions has already been cached. Test if that's
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
  %disp( [ '      Caching stuff for ' num2str( numberOfBasisFunctions ) ' basis functions' ] )
  [ Ax Ay Az A2x A2y A2z correctionNeeded ] = ...
             calculateEverythingForThisResolution( cache.Nx, cache.Ny, cache.Nz, ...
                                                   numberOfBasisFunctions );
  cache.resolution{numberOfBasisFunctions} = ...
       struct( 'Ax', Ax, 'Ay', Ay, 'Az', Az, ...
               'A2x', A2x, 'A2y', A2y, 'A2z', A2z, ...
               'correctionNeeded', correctionNeeded );

end



% Calculate coefficients using seperability of basis functions to speed things up
Nx = cache.Nx;
Ny = cache.Ny;
Nz = cache.Nz;
numberOfImages = cache.numberOfImages;
%disp( '      Pushing the weighted data through the basis functions' )
%tic
%numberOfBasisFunctions
%numberOfImages
C = zeros(numberOfBasisFunctions, numberOfBasisFunctions,numberOfBasisFunctions,numberOfImages);
%size(C)
for channel = 1:numberOfImages
    tmp = cache.resolution{numberOfBasisFunctions}.Ax' * reshape( sum(weights(:,:,:,2*(channel-1)+1:2*channel) .* data(:,:,:,2*(channel-1)+1:2*channel),4), [ Nx Ny*Nz ] );
    tmp = permute( reshape( tmp, [ numberOfBasisFunctions Ny Nz ] ), [ 2 1 3 ] );
    tmp = cache.resolution{numberOfBasisFunctions}.Ay' * reshape( tmp, [ Ny numberOfBasisFunctions*Nz ] );
    tmp = permute( reshape( tmp, [ numberOfBasisFunctions numberOfBasisFunctions Nz ] ), [ 2 1 3 ] );
    tmp = reshape( tmp, [ numberOfBasisFunctions^2 Nz ] );
    
    C(:,:,:,channel) =...
        reshape( tmp * cache.resolution{numberOfBasisFunctions}.Az, [ numberOfBasisFunctions numberOfBasisFunctions numberOfBasisFunctions ] );
end
%toc
%C1 = C(:,:,:,1);
%C2 = C(:,:,:,2);
%C = cat(1,C1(:),C2(:));
% Need to still correct for non-orthonormality of the basis functions by inverting the hessian. If the mask was one everywhere, the hessian would be unity matrix (since wer're using DCT basis functions) so this step wouldn't be needed.
%disp( '      Computing and inverting the Hessian' )
%tic
hessian = getHessian( weights, ...
                      cache.resolution{numberOfBasisFunctions}.A2x, ...
                      cache.resolution{numberOfBasisFunctions}.A2y, ...
                      cache.resolution{numberOfBasisFunctions}.A2z, ...
                      cache.resolution{numberOfBasisFunctions}.correctionNeeded );

C(:) = hessian \ C(:);
%toc
%hessian

%coeffs = C;
%C = zeros(numberOfBasisFunctions, numberOfBasisFunctions, numberOfBasisFunctions, numberOfImages);
%for channel = 1:numberOfImages
%    C(:,:,:,channel) = coeffs((channel-1)*numberOfBasisFunctions + [1:numberOfBasisFunctions],:,:);
%end
%clear coeffs;
% Expand back into image domain
% (this code is actually the same as before, but consistenlty swap A' with A, and numberOfBasisFunctions
% with N)
%disp( '      Expanding back into image domain' )
%tic
smoothedData = zeros(Nx,Ny,Nz,numberOfImages);
for channel = 1:numberOfImages
    tmp = cache.resolution{numberOfBasisFunctions}.Ax * reshape( C(:,:,:,channel), [ numberOfBasisFunctions numberOfBasisFunctions^2 ] );
    tmp = permute( reshape( tmp, [ Nx numberOfBasisFunctions numberOfBasisFunctions ] ), [ 2 1 3 ] );
    tmp = cache.resolution{numberOfBasisFunctions}.Ay * reshape( tmp, [ numberOfBasisFunctions Nx*numberOfBasisFunctions ] );
    tmp = permute( reshape( tmp, [ Ny Nx numberOfBasisFunctions ] ), [ 2 1 3 ] );
    tmp = reshape( tmp, [ Nx*Ny numberOfBasisFunctions ] );
    smoothedData(:,:,:,channel) = reshape( tmp * cache.resolution{numberOfBasisFunctions}.Az', [ Nx Ny Nz ] );
end
%
elapsedTimeSmoothing = toc;
disp( [ '       ' num2str( elapsedTimeSmoothing ) ] )

end


%
%
%
function  [ Ax Ay Az A2x A2y A2z correctionNeeded ] = ...
             calculateEverythingForThisResolution( Nx, Ny, Nz, ...
                                                   numberOfBasisFunctions )
%
%
  %disp( [ '      Precalculating stuff for ' num2str( numberOfBasisFunctions ) ' basis functions' ] )
  %tic

  % Construct DCT basis functions Let's use the definition of Matlab: 
  % http://www.mathworks.com/access/helpdesk/help/toolbox/signal/dct.html

  tmpx = ( 2 * [ 1 : Nx ]' - 1 ) * pi / 2 / Nx;
  Ax = cos( tmpx * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / Nx );
  Ax( :, 1 ) = Ax( :, 1 ) / sqrt( 2 );
  
  tmpy = ( 2 * [ 1 : Ny ]' - 1 ) * pi / 2 / Ny;
  Ay = cos( tmpy * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / Ny );
  Ay( :, 1 ) = Ay( :, 1 ) / sqrt( 2 );
  
  tmpz = ( 2 * [ 1 : Nz ]' - 1 ) * pi / 2 / Nz;
  Az = cos( tmpz * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / Nz );
  Az( :, 1 ) = Az( :, 1 ) / sqrt( 2 );
  
  
  
  % Build other set of basis functions needed to quickly compute the Hessian
  %disp( '      Building Hessian basis functions' )
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
  
  % We are going to use these Hessian basis functions A2x, A2y, A2z as follows:
  %  % Do the operations direction-wise
  %  tmp = A2x' * reshape( weights, [ Nx Ny*Nz ] );
  %  tmp = permute( reshape( tmp, [ numberOfBasisFunctions^2 Ny Nz ] ), [ 2 1 3 ] );
  %  tmp = A2y' * reshape( tmp, [ Ny numberOfBasisFunctions^2*Nz ] );
  %  tmp = permute( reshape( tmp, [ numberOfBasisFunctions^2 numberOfBasisFunctions^2 Nz ] ), [ 2 1 3 ] );
  %  tmp = reshape( tmp, [ numberOfBasisFunctions^4 Nz ] );
  %  Gamma = reshape( tmp * A2z, [ numberOfBasisFunctions^2 numberOfBasisFunctions^2 numberOfBasisFunctions^2 ] );
  %  
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
  
  %  hessian = zeros( numberOfBasisFunctions^3, numberOfBasisFunctions^3 );
  %  hessian( [ 1 : numberOfBasisFunctions^6 ]' + correctionNeeded ) = Gamma(:);
  
  %elapsedTime = toc;
  %disp( [ '       ' num2str( elapsedTime ) ] )

  
  
end



%
%
%
function hessian = getHessian( weights, A2x, A2y, A2z, correctionNeeded )
%
%
  
  [ Nx Ny Nz channels] = size( weights );
  channels = sqrt(channels);
  numberOfBasisFunctions = sqrt( size( A2x, 2 ) );
  
  hessian = zeros( channels*numberOfBasisFunctions^3, channels*numberOfBasisFunctions^3 );
  blockHessian = zeros(numberOfBasisFunctions^3, numberOfBasisFunctions^3);
  counter = 1;
  for channel1 = 1:channels
      for channel2 = 1:channels
        % Do the operations direction-wise
        tmp = A2x' * reshape( weights(:,:,:,counter), [ Nx Ny*Nz ] );
        tmp = permute( reshape( tmp, [ numberOfBasisFunctions^2 Ny Nz ] ), [ 2 1 3 ] );
        tmp = A2y' * reshape( tmp, [ Ny numberOfBasisFunctions^2*Nz ] );
        tmp = permute( reshape( tmp, [ numberOfBasisFunctions^2 numberOfBasisFunctions^2 Nz ] ), [ 2 1 3 ] );
        tmp = reshape( tmp, [ numberOfBasisFunctions^4 Nz ] );
        Gamma = reshape( tmp * A2z, [ numberOfBasisFunctions^2 numberOfBasisFunctions^2 numberOfBasisFunctions^2 ] );
  
        % The matrix Gamma above contains the right entries, but in the incorrect place. 
        % Let's construct the Hessian by reordering the elements. 
        
        %save('biasCorrectionStuff.mat','Gamma','numberOfBasisFunctions','correctionNeeded','blockHessian');
       
        blockHessian( [ 1 : numberOfBasisFunctions^6 ]' + correctionNeeded ) = Gamma(:);
        hessian((channel1-1)*numberOfBasisFunctions^3 + [1:numberOfBasisFunctions^3],...
                (channel2-1)*numberOfBasisFunctions^3 + [1:numberOfBasisFunctions^3]) = blockHessian;
        counter = counter + 1;
      end
  end
  % Do the operations direction-wise
  


  %
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
    directHessian = basisFunctions' * ...
                    (basisFunctions .* repmat( weights( : ), [ 1 numberOfBasisFunctions^3 ] ) );
    
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
  data = randn( Nx, Ny, Nz, 4 );
  
  size(data)

  % Create some weights
  center = [ Nx/2 Ny/2 Nz/2 ]';
  radius = [ Nx*.5 Ny*.5 Nz*.5 ]';
  [ x1, x2, x3 ] = ndgrid( 1:Nx, 1:Ny, 1:Nz );
  distance = ( ( x1 - center( 1 ) ) / radius( 1 ) ).^2 + ...
             ( ( x2 - center( 2 ) ) / radius( 2 ) ).^2 + ...
             ( ( x3 - center( 3 ) ) / radius( 3 ) ).^2;
  weights = exp( - 7.5 *( distance - 1 ) );
  weights = weights ./ ( 1 + weights );
  weightsTMP(:,:,:,1) = weights;
  weightsTMP(:,:,:,2) = weights;
  weightsTMP(:,:,:,3) = weights;
  weightsTMP(:,:,:,4) = weights;
  weights = weightsTMP;
  clear weightsTMP

  % Show what we have
  for n = 1:4
    figure
    %showThreeDImage( data );
    showImage(data(:,:,:,n));
  end
  
  figure
  %showThreeDImage( weights, false, [ 0 1 ] );
  showImage(weights(:,:,:,1))


  numberOfBasisFunctions = 1;
  smoothedData = smoothWithSeparableBasisFunctionsWithWeightsMultiSpec( data, weights, ...
                                                               numberOfBasisFunctions, ...
                                                               maximumNumberOfBasisFunctions );
  for numberOfBasisFunctions = 1 : maximumNumberOfBasisFunctions
    disp( '========================================' )
    disp( [ num2str( numberOfBasisFunctions ) ' basis functions' ] )
    disp( '========================================' )
    for i = 1 : 4
      disp( '  +++++' )
      disp( [ '  run ' num2str( i ) ] )
      disp( '  +++++' )
      smoothedData = smoothWithSeparableBasisFunctionsWithWeightsMultiSpec( data, weights, numberOfBasisFunctions );
    end

    % Show result
    for n = 1:2
        figure
        showImage( smoothedData(:,:,:,n) );
        title( [ num2str( numberOfBasisFunctions ) ' basis functions' ] )
    end
  end


end

