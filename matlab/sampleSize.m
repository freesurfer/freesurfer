function N = sampleSize(stdS, effectsize, alpha, targetPOWER);
% function N = sampleSize(stdS, effectsize, alpha, power);
% 
%                author: Xiao Han
%
% compute the sample-size for significance level alpha (1-sided) and parallel study
%
% stdS is the standard deviation of measurement
% effectsize is the effect to look for (difference in means of two groups)
% alpha is the significance level, default = 0.01
% power is the desired power, default = 0.8

%for two-sided, use
%alpha = 0.5*alpha ?

if nargin < 2,
    error('Insufficient Inputs');
    return;
end;

if nargin < 3,
   alpha = 0.01;
   targetPOWER = 0.8;
elseif nargin < 4
    targetPOWER = 0.8;
end;


N1 = 1;
N2 = 1000;

%first test whether requiring N > 2*1000
%i.e., see whether 2*N2 satisfies the required power and significance level
t_alpha = tinv(1-alpha, N2 + N2 - 2);
t_delta = effectsize*sqrt(0.5*N2)/stdS;

%if the following is true, the significant level won't be achieved 
if t_alpha > t_delta,
    disp(sprintf('Requires more than %d subjects.', 2*N2));
    N = 2*N2;
    return;
end;

%if the following test is true, then the power level won't be achieved
power_N2 = 1 - tcdf(t_alpha - t_delta, N2 + N2 -2);
if power_N2 < targetPOWER,
    disp(sprintf('Requires more than %d subjects.', 2*N2));
    N = 2*N2;
    return; 
end;

N = ceil((N1+N2)*0.5);

%find N by an iterative search
while N2 - N1 > 1,
        
    t_alpha = tinv(1-alpha, N + N -2);
    t_delta = effectsize*sqrt(0.5*N)/stdS;

    %evaluate power
    power_N = 1 - tcdf(t_alpha - t_delta, N + N -2);
    
    if power_N >= targetPOWER,
        N2 = N;
    else
        N1 = N;
    end;

    N = ceil((N1+N2)*0.5);
end;

%reuires N samples for each group, so total size = 2N
N = N + N;

disp(sprintf('The required total sample size (two groups) is %d.',N));
return;