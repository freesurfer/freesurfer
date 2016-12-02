function smoothField = expandBasisFunctions( coefficients, Nx, Ny, Nz, downSamplingFactor )
%
%
%


  % Construct DCT basis functions Let's use the definition of Matlab: 
  % http://www.mathworks.com/access/helpdesk/help/toolbox/signal/dct.html
  disp( '      Building DCT basis functions' )
  numberOfBasisFunctions = size( coefficients, 1 ); % Assuming here that all dimensions have equal amount of basis funcs
  
  tmpx = ( 2 * ( ( [ 1 : Nx ]' - 1 ) / downSamplingFactor + 1 ) - 1 ) * pi / 2 /  length( [ 1 : downSamplingFactor : Nx ] );
  Ax = cos( tmpx * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / length( [ 1 : downSamplingFactor : Nx ] ) );
  Ax( :, 1 ) = Ax( :, 1 ) / sqrt( 2 );
  
  tmpy = ( 2 * ( ( [ 1 : Ny ]' - 1 ) / downSamplingFactor + 1 ) - 1 ) * pi / 2 / length( [ 1 : downSamplingFactor : Ny ] );
  Ay = cos( tmpy * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / length( [ 1 : downSamplingFactor : Ny ] ) );
  Ay( :, 1 ) = Ay( :, 1 ) / sqrt( 2 );
  
  tmpz = ( 2 * ( ( [ 1 : Nz ]' - 1 ) / downSamplingFactor + 1 ) - 1 ) * pi / 2 / length( [ 1 : downSamplingFactor : Nz ] );
  Az = cos( tmpz * [ 0 : numberOfBasisFunctions-1 ] ) * sqrt( 2 / length( [ 1 : downSamplingFactor : Nz ] ) );
  Az( :, 1 ) = Az( :, 1 ) / sqrt( 2 );

  %  figure
  %  subplot( 3, 1, 1 );
  %  plot( Ax( :, 1 : 4 ) )
  %  subplot( 3, 1, 2 );
  %  plot( Ay( :, 1 : 4 ) )
  %  subplot( 3, 1, 3 );
  %  plot( Az( :, 1 : 4 ) )
 

  % Expand the basis into the image domain
  tmp = Ax * reshape( coefficients, [ numberOfBasisFunctions numberOfBasisFunctions^2 ] );
  tmp = permute( reshape( tmp, [ Nx numberOfBasisFunctions numberOfBasisFunctions ] ), [ 2 1 3 ] );
  tmp = Ay * reshape( tmp, [ numberOfBasisFunctions Nx*numberOfBasisFunctions ] );
  tmp = permute( reshape( tmp, [ Ny Nx numberOfBasisFunctions ] ), [ 2 1 3 ] );
  tmp = reshape( tmp, [ Nx*Ny numberOfBasisFunctions ] );
  smoothField = reshape( tmp * Az', [ Nx Ny Nz ] );


