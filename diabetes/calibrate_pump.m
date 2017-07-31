function [Ibest, Sbest, ebest, Cbest] = calibrate_pump(A, M, F, dt)
% [Ibest, Sbest, ebest, Cbest] = calibrate_pump(A, M, F, dt, diag)
% A is active insulin at times 1 and 2
% M is pairs of measured blood glucose level
% F is carbs (in grams) given at time 1 that are metabolized by time 2
% dt is time between measures
%
% outputs:
% Ibest  - BG sensitivity to carbs (amount BG is raised/gram carb)
% Sbest - sensitivity factor (amount BG drops per/unit insulin)
% ebest - error in basal rate
% Cbest - carb ratio (grams of carbohdyrates/unit insulin


min_rms = 10000;

G = M(:,1)-M(:,2);
X = [(A(:,1)-A(:,2))';dt ; -F];
p = inv(X*X')*X*G;
Ibest = p(3) ;
Sbest = p(1) ;
ebest = p(2)/Sbest ;
Cbest = Sbest / Ibest ;


