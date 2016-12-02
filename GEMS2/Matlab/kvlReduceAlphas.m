function [ reducedAlphas reducingLookupTable ] = kvlReduceAlphas( alphas, compressionLookupTableFileName, sameGaussianParameters )
% 
%
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );

numberOfNodes = size( alphas, 1 );
numberOfLabels = size( alphas, 2 );
numberOfReducedLabels = length( sameGaussianParameters );
reducedAlphas = zeros( numberOfNodes, numberOfReducedLabels, 'single' );
reducingLookupTable = zeros( numberOfLabels, 1 );

for reducedLabel = 1 : numberOfReducedLabels
  disp( '--------------' )
  sameGaussians = sameGaussianParameters{ reducedLabel };

  for FreeSurferLabel = sameGaussians
    compressedLabel = find( FreeSurferLabels == FreeSurferLabel ) - 1;
    name = deblank( names( find( FreeSurferLabels == FreeSurferLabel ), : ) );
    disp( [ 'Making ' name ' map to reduced label ' num2str( reducedLabel ) ] )

    reducedAlphas( :, reducedLabel ) = reducedAlphas( :, reducedLabel ) + alphas( :, compressedLabel+1 );
    reducingLookupTable( compressedLabel+1 ) = reducedLabel;
  end

end



