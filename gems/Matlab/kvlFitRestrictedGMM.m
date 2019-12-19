function [ means, sigma, cost ] = kvlFitRestrictedGMM( d, W, n, boundaryEnforcement, runSelfTest )
%
% USAGE: [ means, sigma, cost ] = kvlFitRestrictedGMM( d, W, n )
%
% Given a Ix1 vector d of intensities and an associated IxK weight matrix W with weights associating each 
% intensity to each of K classes, optimize the cost function 
%
%  \sum_{i=1}^I \sum_{k=1}^K -\log[ N( d_i ; means_k, sigma^2 ) ]
%
% with respect to means = (means_1, ..., means_k )' and sigma, under the condition that for any pair (k,k') we have that
%
%    abs( means_k - means_{k'} ) <= n * sigma  
%
% for some user-specified hyperparameter n.
%
% For multicontrast data, d is a matrix of size IxJ with J the number of contrasts, and the problem is solved for 
% each contrast separately, and returned as means with dimensions KxJ and sigma with dimensions 1xJ. The returned
% cost is simply the sum of the J individual costs.
%
%  
% The solution implemented here is as follows. First the unconstrainted problem is solved -- if this satisifies the 
% constraint, then we're lucky and nothing else needs to be done. Otherwise we look for the solution on *some* boundary
% of the "allowable" domain defined by the constraints (obtained by replacing "<=" with "=="). Since these boundaries 
% are not perpendicular to each other w.r.t. the various components of means, the effect of imposing one on a specific
% pair of means components may have non-trivial consequences on the other pairs. For this reason (for better or worse)
% my implementation here simply tests all the possible boundary combinations, and retains the one with minimal cost that
% lies within the "allowable" domain.
% 
%
% For a given set of boundary combinations, the solution proceeds as follows. Defining the variables
%  
%  N_k = \sum_i w_{ik}
%  N = \sum_k N_k
%  mu_k = \sum_i w_{ik} d_i / N_k
%  mu = ( mu_1, ..., mu_K )^T
%  Q = diag( N_1, ... N_K )
%  
%  we follow an optimization scheme where x is optimized out for any given value of sigma: 
%  min_{x,sigma} cost( x, sigma ) = min_sigma [ min_x cost( x, sigma ) ], where the inner optimization (to x, for a 
%  a given value of sigma) is given in closed form. Indeed, the inner optimization involves optimizing 
%  
%    ( x - mu )' Q ( x - mu )   subject to E * x = sigma * f
%    
%  where E is of the form [ -1 1 0 0; 0 -1 1 0; ... ] and f = [ \pm n, \pm n, ... ]. Now by parameterizing x as a 
%  function of a lower-dimensional vector y (every row in E removes one degree of freedom) as follows:
%  
%    x = Z * y + sigma * v, where Z and v are chosen so that E*Z=0 and E*v=f, 
%  
%  we get that for any value of y
%  
%    E * x = ( E * Z ) * y + sigma * E * v 
%          = 0 * y + sigma * f 
%          = sigma * f
%    
%  i.e. the condition "E * x = sigma * f" is automatically satisified, so that we now simply have to solve an 
% *unconstrainted* problem in y!
%  
%  Plugging the form "x = Z * y + sigma * v" into "( x - mu ) Q ( x - mu )", taking the derivative and equating to 
%  zero yields as solution of y
%  
%    yHat = inv( Z' * Q * Z ) * Z' * Q * ( mu - sigma * v ), 
%    
%  which is nothing but a weighted least-squares fit of regression model with basis functions Z to the unconstrainted
%  solution mu.
%  
%  Plugging to resulting solution of x 
%    
%    xHat = Z * yHat + sigma * v
%  
%  and therefore, defining the smoother matrix 
%  
%    S = Z * inv( Z' * Q * Z ) * Z' * Q
%  
%  we have that
%  
%  xHat = [ S - eye( K ) ] * ( mu - sigma * v ) + mu;
%       = S * mu - sigma * [ S - eye( K ) ] * v;
%       = p - sigma * q
%    
% with 
%  
%  p = S * mu 
%  q = [ S - eye( K ) ] * v
%    
%  (note that for no constraints S=I so that p=mu and q=0, giving us the expected xHat=mu in that case). 
% 
% With this solution of xHat as a function of sigma, we can now turn our attention to the (outer) optimization 
% of sigma. Plugging "xHat = p - sigma * q" into the cost function, we eneed to minimize
%  
%    \sum_k \sum_i w_{ik} [ ( d_i - p_k + sigma * q_k )^2 / sigma^2 + log( sigma^2 ) ]  
%  
% The derivative w.r.t. sigma is
%  
%    \sum_k \sum_i w_{ik} [ 2 * q_k * ( d_i - p_k + sigma * q_k ) / sigma^2
%                           +
%                           -2 * ( d_i - p_k + sigma * q_k )^2 / sigma^3 
%                           + 
%                           2 / sigma ]  
%    
% Equating to zero yields 
%  
%    \sum_k \sum_i w_{ik} [
%      q_k * ( d_i - p_k ) * sigma    + q_k^2 * sigma^2
%      + 
%      - q_k^2 * sigma^2     - 2 * ( d_i - p_k ) * q_k * sigma   - (d_i - p_k)^2
%      +
%      sigma^2
%    ] = 0
%    
%  Reshuffling yields the quadratic equation 
%  
%    a * sigma^2 + b * sigma  + c = 0 
%    
%    a = N;
%    b = -\sum_k \sum_i w_{ik} ( d_i - p_k ) * q_k;
%    c = -[ \sum_k \sum_i w_{ik} (d_i - p_k)^2 ] 
%  
%  with solution 
%  
%    sigmaHat = ( -b + sqrt( b^2 - 4*a*c ) ) / ( 2 * a )
%    
%  Note that for no constraints p=mu and q=0, so that for that situtation we readily 
%  get the expected "[ \sum_k \sum_i w_{ik} (d_i - mu_k)^2 ] / N" as the solution.
%  
% 
% One more note: in the code, the variable "boundaryEnforcement" is a 1x(K-2) variable that encodes 
% whether and how each of the K-2 boundaries is enforced. Each boundary has three possibilities:
%   0: boundary is not enforced
%   1: boundary is enforced with a positive value
%   2: boundary is enforced with a negative value
%




if ( nargin == 0 )
  % 
  % Run a self-test
  % 
                  
  % Generate some data
  I = 1000;
  means = [ 10 30 60 ]';
  variance = 10^2;
  weights = [ 0.5 0.3 0.2 ];

  K = length( means );
  d = zeros( I, 1 );
  W = zeros( I, K );
  IcreatedSoFar = 0;
  for k=1:K
    %
    Ik = round( I * weights( k )  );
    if ( Ik > ( I - IcreatedSoFar ) )
      Ik = ( I - IcreatedSoFar );
    end
    
    d( IcreatedSoFar + [ 1 : Ik ] ) = sqrt(variance) * randn( Ik, 1 ) + means( k );
    W( IcreatedSoFar + [ 1 : Ik ], k ) = 1;
    
    %
    IcreatedSoFar = IcreatedSoFar + Ik;
  end
  W = W .* [ rand( I, 1 ) * .3 + .7 ];  % Test case where sum( W, 2 ) smaller than 1 for some rows


  numberOfBins = 100;
  [ ~, binEdges, binNumbers ] = histcounts( d, numberOfBins );
  h = zeros( numberOfBins, 1 );
  sumW = sum( W, 2 );
  for binNumber = 1 : numberOfBins
    h( binNumber ) = sum( sumW( binNumbers == binNumber ) );
  end
  binWidth = binEdges( 2 ) - binEdges( 1 );
  binCenters = binEdges( 1 : end-1 ) + binWidth/2;
  figure
  bar( binCenters, h/sum(h) )
  hold on
  model = zeros( size( binCenters ) );
  for k = 1 : K
    gauss = 1/sqrt( 2 * pi * variance ) * exp( -( binCenters - means( k ) ).^2 / 2 / variance );
    weight = sum( W( :, k ) ) / sum( W(:) );
    model = model + gauss * weight;
    
    plot( binCenters, gauss * weight * binWidth )
  end  
  plot( binCenters, model * binWidth )

  %
  n = 2;
  for testNumber = 1 : 3
    switch testNumber
      case 1
        % Test one specific boundary enforcement
        boundaryEnforcement = [ 0 2 ]
        runSelfTest = true;
        [ xHat, sigmaHat, costHat ] = kvlFitRestrictedGMM( d, W, n, boundaryEnforcement, runSelfTest );
      case 2
        %  Normal user-mode: give me the best solution among all possible boundary enforcements
        [ xHat, sigmaHat, costHat ] = kvlFitRestrictedGMM( d, W, n );
      % case 3
      %   % Normal user-mode with only one class
      %  [ xHat, sigmaHat, costHat ] = kvlFitRestrictedGMM( d, W( :, 1 ), n )
      otherwise
        % Normal user-mode with multicontrast data
        d = [ d 100-3*d ];
        [ xHat, sigmaHat, costHat ] = kvlFitRestrictedGMM( d, W, n );
    end
    
    % Analyze 
    E_all = zeros( K-1, K );
    template = [ -1 1 zeros( 1, K-2 ) ];
    for k=1:K-1
      E_all( k, : ) = circshift( template, k-1 );
    end
    xHat
    sigmaHat
    costHat
    ( E_all * xHat ) ./ sigmaHat

    
    
    %
    figure
    bar( binCenters, h/sum(h) )
    hold on
    model = zeros( size( binCenters ) );
    xHat = xHat( :, 1 ); sigmaHat = sigmaHat( 1 ); % Only show for 1D 
    for k = 1 : K
      gauss = 1/sqrt( 2 * pi * sigmaHat^2 ) * exp( -( binCenters - xHat( k ) ).^2 / 2 / sigmaHat^2 );
      weight = sum( W( :, k ) ) / sum( W(:) );
      model = model + gauss * weight;
      
      plot( binCenters, gauss * weight * binWidth )
    end  
    plot( binCenters, model * binWidth )

    %
    %disp( '------ Press ENTER to continue' )
    %pause
    
  end % End loop over tests
  
  
  return
end % End self-test 
  

if ( nargin < 5 )
  runSelfTest = false;
end
  
  
if ( nargin == 3 )
  %
  % Normal end-user usage
  %

  %
  K = size( W, 2 );
  J = size( d, 2 );
  
  % Catch multicontrast input 
  if ( J > 1 )
    %
    means = zeros( K, J );
    sigma = zeros( 1, J )
    cost = 0;
    for j = 1 : J 
      [ currentMeans, currentSigma, currentCost ] = kvlFitRestrictedGMM( d( :, j ), W, n );
      
      means( :, j ) = currentMeans;
      sigma( j ) = currentSigma;
      cost = cost + currentCost;
      
    end
  
    return
  end
  
  % Catch trivial case of only a single class
  if ( K == 1 )
    means = W' * d / sum( W );
    variance = W' * ( d - means ).^2 / sum( W ); 
    sigma = sqrt( variance );
    cost = W' * ( ( d - means ).^2 / sigma^2 + log( sigma^2 ) );
    return
  end
  
  
  % Set up test variables to determine whether solutions fall inside the allowable parameter domain later
  E_all = zeros( K-1, K );
  template = [ -1 1 zeros( 1, K-2 ) ];
  for k=1:K-1
    E_all( k, : ) = circshift( template, k-1 );
  end
  
  % Now loop over all possible boundaryEnforcements, each time computing the corresponding solution and 
  % retaining the one that is in the allowable domain with the minimum cost. 
  numberOfBoundaries = K-1;
  numberOfBoundaryCombinations = 3^(K-1); 
  bestCostSoFar = Inf;
  for boundaryCombinationNumber = 0 : numberOfBoundaryCombinations-1
    % Construct the boundaryEnforcement code
    power = 3.^[ numberOfBoundaries-1 : -1 : 0 ];
    boundaryEnforcement = floor( rem( boundaryCombinationNumber * ones( 1, numberOfBoundaries ), 3 * power ) ./ power );
    
    % Find the corresponding solution
    runSelfTest = false;
    [ currentMeans, currentSigma, currentCost ] = kvlFitRestrictedGMM( d, W, n, boundaryEnforcement, runSelfTest );
    
    % Check if solution is in allowable space
    numericalZeroThreshold = 1e-5; % Theoretically zero
    % isAllowable = ( max( abs( E_all * currentMeans ) - currentSigma * n * ones( K-1, 1 ) ) < numericalZeroThreshold )
    normalizedAbsoluteDistanceBetweenMeans = abs( E_all * currentMeans ) / currentSigma;
    isAllowable = ( max( normalizedAbsoluteDistanceBetweenMeans - n * ones( K-1, 1 ) ) < numericalZeroThreshold );
  
  
    % 
    if isAllowable
      if ( currentCost < bestCostSoFar )
        % Rememer this solution
        means = currentMeans;
        sigma = currentSigma;
        cost = currentCost;
        bestCostSoFar = cost;
    
        % First code (no boundaries enforced whatsoever) is special
        if ( boundaryCombinationNumber == 0 )
          disp( 'immediately hit the jack pot' )
          break;
        end

      end % End check if solution is best so far
     
    end % End check if allowable solution

    %
    % disp( '-------- Press Enter to continue-------' )
    % pause
    
  end % End loop over all boundary combinations  

  
  return
end

    
  
%
K = size( W, 2 );
constraintsToConsider = find( boundaryEnforcement );
E_all = zeros( K-1, K );
template = [ -1 1 zeros( 1, K-2 ) ];
for k=1:K-1
  E_all( k, : ) = circshift( template, k-1 );
end
f_all = zeros( K-1, 1 );  % Just one example - doesn't matter
f_all( constraintsToConsider ) = n * ( -1 ).^( boundaryEnforcement( constraintsToConsider ) == 2 );
v_all = E_all \ f_all;


% 
E = E_all( constraintsToConsider, : );
Z = null( E );
f = f_all( constraintsToConsider );
v = v_all;


% Make sure we've defined Z and v correctly
if runSelfTest
  assert( max( max( abs( E*Z ) ) ) < 1e-12 )
  assert( max( abs( E*v - f ) ) < 1e-12 )
end


% Solve for sigma
Ns = sum( W );
N = sum( Ns );
mu = ( sum( d .* W ) ./ Ns )';
Q = diag( Ns );
S = Z * inv( Z' * Q * Z ) * Z' * Q;
p = S * mu;
q = [ S - eye( K ) ] * v;

a = N;
b = -sum( ( d - p' ) .* W ) * q;
c = -sum( sum( ( d - p' ).^2 .* W ) ); 
sigmaHat = ( -b + sqrt( b^2 - 4*a*c ) ) / ( 2 * a );


% Solve for x
xHat = p - sigmaHat * q;


% Make sure the conditions are satisfied
if runSelfTest
  assert( max( abs( E * xHat - sigmaHat * f ) ) < 1e-12 )
end

%
cost = sum( sum( ( ( d - xHat' ).^2 / sigmaHat^2 + log( sigmaHat^2 ) ) .* W ) );


if runSelfTest
  %
  % Make sure we couldn't have found any better solution. For x, do it by finding y, and changing it a bit 
  % and verifying that the cost is optimal. For sigma, do it by varying it a bit, computing corresponding x,
  % and verifying effect on cost.
  %
  yHat = ( Z' * Z ) \ ( Z' * ( xHat - sigmaHat * v ) );
  %  y = yHat;
  %  x = Z * y + sigmaHat * v;
  %  cost = sum( sum( ( ( d - x' ).^2 / sigmaHat^2 + log( sigmaHat^2 ) ) .* W ) );
  for dimensionNumber = 1 : length( yHat )
    for signNumber = 1 : 2
      y = yHat; y( dimensionNumber ) = y( dimensionNumber ) + .001 * (-1)^signNumber;
      x = Z * y + sigmaHat * v;
      costNew = sum( sum( ( ( d - x' ).^2 / sigmaHat^2 + log( sigmaHat^2 ) ) .* W ) );
    
      assert( costNew > cost )
    end
  end
  for signNumber = 1 : 2
    sigma = sigmaHat + .001 * (-1)^signNumber;
    x = p - sigma * q;
    costNew = sum( sum( ( ( d - x' ).^2 / sigma^2 + log( sigma^2 ) ) .* W ) );
    assert( costNew > cost )
  end
  
end


% Give variables the correct names to return
means = xHat;
sigma = sigmaHat;

