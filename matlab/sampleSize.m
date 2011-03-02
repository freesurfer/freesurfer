function N = sampleSize(stdS, effectsize, alpha, targetPOWER, alternative);
% function N = sampleSize(stdS, effectsize, alpha, power, alternative);
% 
% compute the required sample-size for significance level alpha and power
% assuming parallel study with two groups of equal-size
%
% stdS is the standard deviation of measurement (equal for both groups)
% effectsize is the effect to look for (difference in means of two groups)
% alpha is the significance level, default = 0.01
% power is the desired power, default = 0.8
% alternative is a string for either 'one-sided' or 'two-sided'; 
% default is 'two-sided'
% tested against http://www.stat.uiowa.edu/~rlenth/Power/, 
% and the results match perfectly

%
% sampleSize.m
%
% Original Author: Xaio Han
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.5 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if nargin < 2,
    error('Insufficient Inputs');
    return;
end;

if nargin < 3,
   alpha = 0.01;
end;

if nargin < 4
    targetPOWER = 0.8;
end;

if nargin < 5,
    alternative = 'two-sided';
end;
    
if strcmp(alternative, 'one-sided'),
    disp('assuming one-sided t-test');
    alpha = alpha;
else,
    disp('assuming two-sided t-test');
    alpha = 0.5*alpha;
end;
    
N1 = 1;
N2 = 1000;

%first test whether requiring N > 2*1000
%i.e., see whether N2 satisfies the required power and significance level
t_alpha = tinv(1-alpha, N2 + N2 - 2);
noncentrality = effectsize*sqrt(0.5*N2)/stdS;

%if the following is true, the significant level won't be achieved 
%if t_alpha > noncentrality,
%    disp(sprintf('Requires more than %d subjects for each group.', N2));
%    N = N2;
%    return;
%end;

%if the following test is true, then the power level won't be achieved
power_N2 = 1 - nctcdf(t_alpha, N2 + N2 -2, noncentrality);
if power_N2 < targetPOWER,
    disp(sprintf('Requires more than %d subjects for each group.', N2));
    N = N2;
    return; 
end;

N = ceil((N1+N2)*0.5);

%find N by an iterative search
while N2 - N1 > 1,
        
    t_alpha = tinv(1-alpha, N + N -2);
    noncentrality = effectsize*sqrt(0.5*N)/stdS;

    %evaluate power
    power_N = 1 - nctcdf(t_alpha, N + N -2, noncentrality);
    
    if power_N >= targetPOWER,
        N2 = N;
    else
        N1 = N;
    end;

    N = ceil((N1+N2)*0.5);
end;

%requires N samples for each group
disp(sprintf('The required sample size for each group is %d.',N));
return;
