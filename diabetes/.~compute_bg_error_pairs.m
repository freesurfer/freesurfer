function L1_error = compute_bg_error(M,D,F,dt,S,e,I,C)
% rms_error = compute_bg_error(M,D,F,dt,S,e,I,c)

% M is pairs of BG meauurements, (:,1) when carbs are given
%                   and (:,2) when they are fully absorbed
% D is vector of insulin absorbed for pair i 
% F is carbs given
% dt is time interval

% A used to be active, not used anymore
npairs = size(M,1);

L1_error = 0 ;
for i=1:npairs
%    M1 = predict_bg(M(i,1), S, A(i,1)-A(i,2), e, dt, F(i), C, I);
    M1 = predict_bg(M(i,1), S, D(i), e, dt(i), F(i), C, I);
    error = abs(M(i,2)-M1) ;
    %disp(sprintf('M %2.0f --> %2.0f, predicted = %2.0f', M(i,1), M(i,2), M1));
    L1_error = L1_error + error ;
end

L1_error = (L1_error/npairs) ;
