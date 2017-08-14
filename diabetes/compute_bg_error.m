function [rms_error, BG_predicted] = compute_bg_error(M0, M,A,F,dt,S,basal_absorbed,I)
% rms_error = compute_bg_error(M,A,D,F,dt,S,e,I,c)
% M - BG measures
% A - insulin absorbed since carbs were given
% F - Carbs given
% dt- time since previous measure
% S - sensitivity
% basal_absorbed - the cumulative amount of basal absorbed at each time
% I - carbs/BG ratio

% C - carbs/insulin ratio   not needed - ratio of S and I

ntps = size(M,2);

rms_error = 0 ;
for i=1:ntps
    M1 = predict_bg(M0, S, A(i), basal_absorbed(i), dt(i)-dt(1), F(i), I);
    BG_predicted(i) = M1 ;
%    error = abs(M(i)-M1) ;
    error = (M(i)-M1).^2 ;
%    error = (M(i)-M1).^6 ;
    rms_error = rms_error + error ;
end

%rms_error = (rms_error/ntps) ;
rms_error = sqrt(rms_error/ntps) ;
