function significance = permutationTest( data1, data2, numberOfPerturbations );
%
%  data1 = randn( 50, 1 ) + 40;
%  data2 = randn( 20, 1 ) + 41;
%  numberOfPerturbations = 10000;
%

if ( nargin < 3 )
  numberOfPerturbations = 10000;
end



totalData = [ data1; data2 ];
numberOfDataPoints = length( totalData );
numberOfData1Points = length( data1 );

features = zeros( numberOfPerturbations, 1 );
for i = 1 : numberOfPerturbations
  permutedIndices = randperm( numberOfDataPoints );
  permutatedData1 = totalData( permutedIndices( 1 : numberOfData1Points ) );
  permutatedData2 = totalData( permutedIndices( numberOfData1Points+1 : end ) );

  meanOfPermutatedData1 = sum( permutatedData1 ) / length( permutatedData1 );
  meanOfPermutatedData2 = sum( permutatedData2 ) / length( permutatedData2 );
  feature = meanOfPermutatedData2 - meanOfPermutatedData1;

  features( i ) = feature;
end

%  figure
%  hist( features, 30 )


meanOfData1 = sum( data1 ) / length( data1 );
meanOfData2 = sum( data2 ) / length( data2 );
observedFeature = meanOfData2 - meanOfData1;


significanceOfBeingSoBig = 1 - sum( observedFeature > features ) / numberOfPerturbations;
significanceOfBeingSoSmall = 1 - sum( observedFeature < features ) / numberOfPerturbations;

significance = min( significanceOfBeingSoBig, significanceOfBeingSoSmall );




