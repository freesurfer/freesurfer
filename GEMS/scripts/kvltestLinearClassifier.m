%
close all
clear all


% Generate some data
x1s = randn( 30, 2 ) + repmat( [ 2.5 0 ], [ 30, 1 ] );
x2s = randn( 30, 2 ) + repmat( [ 0.1 0 ], [ 30, 1 ] );
angle = 45
mixingMatrix = [ cos( angle / 180 * pi )  -sin( angle  / 180 * pi ); ...
                 sin( angle / 180 * pi )   cos( angle / 180 * pi ) ];
mixingMatrix = diag( [ 4.5 0.2 ] ) * mixingMatrix;
mixingOffset = [ 20.0 -5.0 ]';
x1s = mixingMatrix * x1s' + repmat( mixingOffset, [ 1 30 ] );
x2s = mixingMatrix * x2s' + repmat( mixingOffset, [ 1 30 ] );



%
figure
l1 = line( x1s( 1, : ), x1s( 2, : ), 'marker', 'o', 'linestyle', 'none', 'color', 'b' );
hold on
l2 = line( x2s( 1, : ), x2s( 2, : ), 'marker', 'o', 'linestyle', 'none', 'color', 'r' );



% Normalize the data
numberOfFeatures = size( x1s, 1 );
for featureNumber = 1 : numberOfFeatures
  data = [ x1s( featureNumber, : ) x2s( featureNumber, : ) ];
  mean = sum( data ) / length( data );
  variance = sum( ( data - mean ).^2 ) / length( data );
  x1s( featureNumber, : ) = ( x1s( featureNumber, : ) - mean ) / sqrt( variance );
  x2s( featureNumber, : ) = ( x2s( featureNumber, : ) - mean ) / sqrt( variance );
end

%
figure
l1 = line( x1s( 1, : ), x1s( 2, : ), 'marker', 'o', 'linestyle', 'none', 'color', 'b' );
hold on
l2 = line( x2s( 1, : ), x2s( 2, : ), 'marker', 'o', 'linestyle', 'none', 'color', 'r' );



% Bishop's book p. 186-189
N1 = size( x1s, 2 );
N2 = size( x2s, 2 );
m1 = sum( x1s, 2 ) / N1;
m2 = sum( x2s, 2 ) / N2;
Sw = ( x1s - repmat( m1, [ 1 N1 ] ) ) * ( x1s - repmat( m1, [ 1 N1 ] ) )' + ...
     ( x2s - repmat( m2, [ 1 N2 ] ) ) * ( x2s - repmat( m2, [ 1 N2 ] ) )';
w = Sw \ ( m2 - m1 );
w = w / sqrt( sum( w.^2 ) );

estimatedAngle = acos( w(1) ) / pi * 180;
if ( w(2) < 0 )
  estimatedAngle = 360 - estimatedAngle;
end
if ( estimatedAngle > 180 )
  estimatedAngle = estimatedAngle - 180;
end
estimatedAngle







%  % The following doesn't say much, as the feature values are not normalized -> not
%  % comparable in magnitude...
%  for i = 1 : numberOfFeatures
%    disp( [ 'optimal weight for ' deblank( labels( i, : ) ) ': ' num2str( w( i ) ) ] )
%  end


% Compare to simply looking at one feature at a time
ws = [ eye( numberOfFeatures, numberOfFeatures ) w ];
titles = strvcat( 'dimension 1', 'dimension 2', 'best projection' ); 

for i = 1 : size( ws, 2 )
  % 
  w = ws( :, i );

  % Projections
  p1s = w' * x1s;
  p2s = w' * x2s;
  pm1 = w' * m1; % projection of m1
  pm2 = w' * m2; % projection of m2
  middle = ( pm1 + pm2 ) / 2;
  p1s = p1s - middle;
  p2s = p2s - middle;

  if 0
    % Instead of just setting threshold in the middle (i.e., at 0 in the current
    % normalized data), try to find the best value
  end

  figure
  if 0
    [ n1 x1 ] = hist( p1s', 10 );
    bar( x1, n1 );
    hold on
    [ n2 x2 ] = hist( p2s', 10 );
    bar( x2, n2, 'r' );
  else
    line( p1s', zeros( N1, 1 ) - .3, 'marker', 'o', 'linestyle', 'none', 'color', 'b' );
    hold on
    line( p2s', zeros( N2, 1 ) + .3, 'marker', 'o', 'linestyle', 'none', 'color', 'r' );
    set( gca, 'ylim', [ -1 1 ] );
    trainingCorrectClassificationRate = ( sum( p1s < 0 ) + sum( p2s > 0 ) ) / ( N1 + N2 );
    trainingCorrectClassificationRate = max( trainingCorrectClassificationRate, 1 - trainingCorrectClassificationRate );
    title( [ deblank( titles( i, : ) ) ' (' num2str( trainingCorrectClassificationRate * 100 ) '%)' ] )
  end
end



% Let's do a proper leave-one-out validation
useSimpleThreshold = true;

xs = [ x1s x2s ];
classifications = [ zeros( 1, size( x1s, 2 ) ) ones( 1, size( x2s, 2 ) ) ];
numberOfSubjects = length( classifications );

numberOfClassificationErrors = 0;
for subjectNumber = 1 : numberOfSubjects
  useSubjects = ones( 1, numberOfSubjects );
  useSubjects( subjectNumber ) = 0;

  trainingXs = xs( :, find( useSubjects ) );
  trainingClassifications = classifications( :, find( useSubjects ) );

  x1s = trainingXs( :, find( trainingClassifications == 0 ) ); % MCI
  x2s = trainingXs( :, find( trainingClassifications == 1 ) ); % AD

  % Bishop's book p. 186-189
  N1 = size( x1s, 2 );
  N2 = size( x2s, 2 );
  m1 = sum( x1s, 2 ) / N1;
  m2 = sum( x2s, 2 ) / N2;
  Sw = ( x1s - repmat( m1, [ 1 N1 ] ) ) * ( x1s - repmat( m1, [ 1 N1 ] ) )' + ...
      ( x2s - repmat( m2, [ 1 N2 ] ) ) * ( x2s - repmat( m2, [ 1 N2 ] ) )';
  w = Sw \ ( m2 - m1 );
  w = w / sqrt( sum( w.^2 ) );


  % Check if we would classify it correctly
  testingX = xs( :, subjectNumber );
  testingClassification = classifications( :, subjectNumber );
  if useSimpleThreshold
    % Use middle between means as the threshold to classify
    pm1 = w' * m1; % projection of m1
    pm2 = w' * m2; % projection of m2
    middle = ( pm1 + pm2 ) / 2;
    predictedClassification = ( ( w' * testingX ) > middle );
  else
    % Use proper Gaussian mixture model to determine the threshold
    px1s = w' * x1s;
    px2s = w' * x2s;
    pm1 = sum( px1s ) / N1;
    pm2 = sum( px2s ) / N2;
    pv1 = sum( ( px1s - pm1 ).^2 ) / N1;
    pv2 = sum( ( px2s - pm2 ).^2 ) / N2;
    propToPosterior1 = 1 / sqrt( 2 * pi * pv1 ) * exp( -( w' * testingX - pm1 )^2 / 2 / pv1 ) * N1;
    propToPosterior2 = 1 / sqrt( 2 * pi * pv2 ) * exp( -( w' * testingX - pm2 )^2 / 2 / pv2 ) * N2;
    predictedClassification = ( propToPosterior2 > propToPosterior1 );
  end


  if ( testingClassification ~= predictedClassification )
    numberOfClassificationErrors = numberOfClassificationErrors + 1;
  end

end % End loop over all subjects

testingCorrectClassificationRate = ( numberOfSubjects - numberOfClassificationErrors ) ...
                                  / numberOfSubjects;
disp( [ 'testingCorrectClassificationRate: ' num2str( testingCorrectClassificationRate * 100 ) ] )


% Compare with just summing
xs = sum( xs );
numberOfClassificationErrors = 0;
for subjectNumber = 1 : numberOfSubjects
  useSubjects = ones( 1, numberOfSubjects );
  useSubjects( subjectNumber ) = 0;

  trainingXs = xs( :, find( useSubjects ) );
  trainingClassifications = classifications( :, find( useSubjects ) );

  x1s = trainingXs( :, find( trainingClassifications == 0 ) ); % MCI
  x2s = trainingXs( :, find( trainingClassifications == 1 ) ); % AD

  % Bishop's book p. 186-189
  N1 = size( x1s, 2 );
  N2 = size( x2s, 2 );
  m1 = sum( x1s, 2 ) / N1;
  m2 = sum( x2s, 2 ) / N2;
  Sw = ( x1s - repmat( m1, [ 1 N1 ] ) ) * ( x1s - repmat( m1, [ 1 N1 ] ) )' + ...
      ( x2s - repmat( m2, [ 1 N2 ] ) ) * ( x2s - repmat( m2, [ 1 N2 ] ) )';
  w = Sw \ ( m2 - m1 );
  w = w / sqrt( sum( w.^2 ) );


  % Check if we would classify it correctly
  testingX = xs( :, subjectNumber );
  testingClassification = classifications( :, subjectNumber );
  if useSimpleThreshold
    % Use middle between means as the threshold to classify
    pm1 = w' * m1; % projection of m1
    pm2 = w' * m2; % projection of m2
    middle = ( pm1 + pm2 ) / 2;
    predictedClassification = ( ( w' * testingX ) > middle );
  else
    % Use proper Gaussian mixture model to determine the threshold
    px1s = w' * x1s;
    px2s = w' * x2s;
    pm1 = sum( px1s ) / N1;
    pm2 = sum( px2s ) / N2;
    pv1 = sum( ( px1s - pm1 ).^2 ) / N1;
    pv2 = sum( ( px2s - pm2 ).^2 ) / N2;
    propToPosterior1 = 1 / sqrt( 2 * pi * pv1 ) * exp( -( w' * testingX - pm1 )^2 / 2 / pv1 ) * N1;
    propToPosterior2 = 1 / sqrt( 2 * pi * pv2 ) * exp( -( w' * testingX - pm2 )^2 / 2 / pv2 ) * N2;
    predictedClassification = ( propToPosterior2 > propToPosterior1 );
  end


  if ( testingClassification ~= predictedClassification )
    numberOfClassificationErrors = numberOfClassificationErrors + 1;
  end

end % End loop over all subjects

testingCorrectClassificationRate = ( numberOfSubjects - numberOfClassificationErrors ) ...
                                  / numberOfSubjects;
disp( [ 'testingCorrectClassificationRate [1 1]: ' num2str( testingCorrectClassificationRate * 100 ) ] )



