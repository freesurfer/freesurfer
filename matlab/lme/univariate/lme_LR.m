function lrstats = lme_LR(lrmlfull,lrmlred,q)
% lrstats = lme_LR(lrmlfull,lrmlred,q)
%
% Likelihood ratio test for the random effects. It can be used to test if a
% model with q+1 random effects is significantly better than a model with q
% random effects.
%
% Input
% lrmlfull: Maximum value for the restricted log-likelihood of the full 
% model (the one with q+1 random effects).
% lrmlred: Maximum value for the restricted log-likelihood of the reduced 
% model (the one with q random effects).
% q: Number of random effects in the reduced model.
%
% Output
% lrstats.G: Likelihood ratio test statistic.
% lrstats.pval: P-value of the test (based on a 50:50 mixture of 
% chi-squared distributions with q and q+1 degrees of freedom). 
% lrstats.df: Degrees of freedom.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:11 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:11 $
%    $Revision: 1.1.2.2 $
% References: Fitzmaurice GM, Laird NM, Ware JH (2004): Applied longitudinal
% analysis. Hoboken, N.J.: Wiley-Interscience. xix, 506 p. p.
%

if nargin < 3 
    error('Too few inputs');   
end;

G = 2*(lrmlfull-lrmlred);
lrstats.G = G;
lrstats.pval = lrpval(q,G);
lrstats.df = [q q+1];
end





function pval = lrpval(q,G)
%Critical points for a 50:50 mixture of chi-squared distributions with q
%and q+1 degrees of freedom; right hand tail probabilities.
if q > 10 
   error('The maximum supported number of random effects is 10');   
end;
sigs = zeros(1,10);
sigs(1) = 0.3; sigs(2) = 0.2; sigs(3) = 0.1; sigs(4) = 0.05;
sigs(5) = 0.025; sigs(6) = 0.01; sigs(7) = 0.005; sigs(8) = 0.0025;
sigs(9) = 0.001; sigs(10) = 0.0005;

tb = zeros(11,10);
tb(1,1) = 0.28; tb(1,2) = 0.71; tb(1,3) = 1.64; tb(1,4) = 2.71;
tb(1,5) = 3.84; tb(1,6) = 5.41; tb(1,7) = 6.63; tb(1,8) = 7.88;
tb(1,9) = 9.55; tb(1,10) = 10.83; 
tb(2,1) = 1.76; tb(2,2) = 2.50; tb(2,3) = 3.81; tb(2,4) = 5.14;
tb(2,5) = 6.48; tb(2,6) = 8.27; tb(2,7) = 9.63; tb(2,8) = 11.00;
tb(2,9) = 12.81; tb(2,10) = 14.18; 
tb(3,1) = 3.06; tb(3,2) = 3.98; tb(3,3) = 5.53; tb(3,4) = 7.05;
tb(3,5) = 8.54; tb(3,6) = 10.50; tb(3,7) = 11.97; tb(3,8) = 13.43;
tb(3,9) = 15.36; tb(3,10) = 16.80; 
tb(4,1) = 4.29; tb(4,2) = 5.36; tb(4,3) = 7.09; tb(4,4) = 8.76;
tb(4,5) = 10.38; tb(4,6) = 12.48; tb(4,7) = 14.04; tb(4,8) = 15.59;
tb(4,9) = 17.61; tb(4,10) = 19.13; 
tb(5,1) = 5.49; tb(5,2) = 6.68; tb(5,3) = 8.57; tb(5,4) = 10.37;
tb(5,5) = 12.10; tb(5,6) = 14.32; tb(5,7) = 15.97; tb(5,8) = 17.59;
tb(5,9) = 19.69; tb(5,10) = 21.27; 
tb(6,1) = 6.66; tb(6,2) = 7.96; tb(6,3) = 10.00; tb(6,4) = 11.91;
tb(6,5) = 13.7; tb(6,6) = 16.07; tb(6,7) = 17.79; tb(6,8) = 19.47;
tb(6,9) = 21.66; tb(6,10) = 23.29; 
tb(7,1) = 7.82; tb(7,2) = 9.21; tb(7,3) = 11.38; tb(7,4) = 13.40;
tb(7,5) = 15.32; tb(7,6) = 17.76; tb(7,7) = 19.54; tb(7,8) = 21.29;
tb(7,9) = 23.55; tb(7,10) = 25.23; 
tb(8,1) = 8.97; tb(8,2) = 10.44; tb(8,3) = 12.74; tb(8,4) = 14.85;
tb(8,5) = 16.86; tb(8,6) = 19.38; tb(8,7) = 21.23; tb(8,8) = 23.04;
tb(8,9) = 25.37; tb(8,10) = 27.10; 
tb(9,1) = 10.10; tb(9,2) = 11.66; tb(9,3) = 14.07; tb(9,4) = 16.27;
tb(9,5) = 18.35; tb(9,6) = 20.97; tb(9,7) = 22.88; tb(9,8) = 24.74;
tb(9,9) = 27.13; tb(9,10) = 28.91; 
tb(10,1) = 11.23; tb(10,2) = 12.87; tb(10,3) = 15.38; tb(10,4) = 17.67;
tb(10,5) = 19.82; tb(10,6) = 22.52; tb(10,7) = 24.49; tb(10,8) = 26.40;
tb(10,9) = 28.86; tb(10,10) = 30.68; 
tb(11,1) = 12.35; tb(11,2) = 14.06; tb(11,3) = 16.67; tb(11,4) = 19.04;
tb(11,5) = 21.27; tb(11,6) = 24.05; tb(11,7) = 26.07; tb(11,8) = 28.02;
tb(11,9) = 30.54; tb(11,10) = 32.40; 

pos = find(G>tb(q+1,:),1,'last');
if isempty(pos)
   pval = sigs(1);
else
   pval = sigs(pos);
end
end
