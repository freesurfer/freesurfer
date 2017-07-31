function rms_error = compute_bg_error(M,A,D,F,dt,S,e,I,C)
% rms_error = compute_bg_error(M,A,D,F,dt,S,e,I,c)
% M - BG measures
% A - active insulin measures
% D - doses given
% F - Carbs given
% dt- time since previous measure
% S - sensitivity
% e - basal error
% I - carbs/BG ratio
% C - carbs/insulin ratio

ntps = size(M,2);

rms_error = 0 ;
for i=2:ntps
    insulin = A(i-1) + D(i-1) - A(i) ;
    M1 = predict_bg(M(i-1), S, insulin, e, dt(i), F(i-1), C, I);
    error = abs(M(i)-M1) ;
    rms_error = rms_error + error ;
end

rms_error = (rms_error/(ntps-1)) ;
%rms_error = sqrt(rms_error/(ntps-1)) ;
