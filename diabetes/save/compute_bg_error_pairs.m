function L1_error = compute_bg_error(M,D,F,dt,S,e,I,debug)
% rms_error = compute_bg_error(M,D,F,dt,S,e,I,c,debug)

% M is pairs of BG meauurements, (:,1) when carbs are given
%                   and (:,2) when they are fully absorbed
% D is vector of insulin absorbed for pair i 
% F is carbs given
% dt is time interval

% A used to be active, not used anymore
npairs = size(M,1);

if (nargin < 8)
   debug = 0 ;
end

min_error = 10 ;
max_error = 40 ;
num_correct = 0 ;
L1_error = 0 ;
for i=1:npairs
%    M1 = predict_bg(M(i,1), S, A(i,1)-A(i,2), e, dt, F(i), I);
    M1 = predict_bg(M(i,1), S, D(i), e, dt(i), F(i), I);
    error = abs(M(i,2)-M1) ;
    if (error > max_error)
       error = max_error ;
    end
    if (error <= min_error)
       num_correct = num_correct+1 ;
    end

    if (debug)
        disp(sprintf('M %2.0f --> %2.0f, predicted = %2.0f, error = %2.0f', M(i,1), M(i,2), M1,error));
    end
    L1_error = L1_error + error ;
end

%L1_error = 100*(npairs-num_correct)/npairs ;
L1_error = (L1_error/npairs) ;
