function sumOfAlphas = sumAlphas( alphas, FreeSurferLabels, compressedLabels, names, sameGaussians )
% 
%

numberOfNodes = size( alphas, 1 );
sumOfAlphas = zeros( numberOfNodes, 1 );
for FreeSurferLabel = sameGaussians
  compressedLabel = compressedLabels( find( FreeSurferLabels == FreeSurferLabel ) );
  name = deblank( names( find( FreeSurferLabels == FreeSurferLabel ), : ) );
  disp( [ 'Adding ' name ' to global class' ] )

  sumOfAlphas = sumOfAlphas + alphas( :, compressedLabel+1 );
end

