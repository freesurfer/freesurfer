function [G, Ip] = compute_BG_from_insulin(I, parms)
% function [G, Ip] = compute_BG_from_insulin(I, parms)

T1  = parms.T1/parms.dt;
T2  = parms.T2/parms.dt;

%% Define operators
t = parms.dt*(0:(T2-1))';
allOnes = ones(T2,1);

OP_Id = eye(T2);
OP_d = (diag(ones(T2-1,1),1) - diag(ones(T2-1,1),-1))/(2*parms.dt);
OP_d2 = (diag(ones(T2-1,1),1) + diag(ones(T2-1,1),-1) - 2*eye(T2))/parms.dt^2;
OP_int = (tril(ones(T2), 0) - .5*eye(T2) - .5*[ones(T2,1) zeros(T2,T2-1)])*parms.dt;

OP_Id = OP_Id(2:end-1, :); % Ignoring the boundary
OP_d = OP_d(2:end-1, :);
OP_d2 = OP_d2(2:end-1, :);
OP_I = (1/parms.k3)*OP_d2 + (1 + (parms.k2+parms.k1)/parms.k3)*OP_d + parms.k1*OP_Id; % Pump insulin

exe = 2*(abs(t-200)<50);
exe = 0*(abs(t-200)<50);

%% Optimize
% Assuming the blood glucose is G = l - K*Ip, we minimize ||G-parms.G_target||_2
% by quadratic programming, or ||G-parms.G_target||_1 by linear programming, or
% simply G by linear programming/
K = parms.insulin_sensitivity*parms.k1*OP_int;
l = parms.G0*allOnes + parms.k6*t + exp(-parms.k7*t).*(OP_int*(exe.*exp(parms.k7*t))) - OP_int*exe;
for nCarb = 1:length(parms.carb_grams)
    l = l + parms.carb_sensitivity*parms.carb_grams(nCarb)*(1-exp(-parms.k4*(t-parms.carb_delay(nCarb)))).*(t>=parms.carb_delay(nCarb));
end

Aeq1 = [1 zeros(1,T2-1); 0 1 zeros(1,T2-2); zeros(1,T2-1) 1];
beq1 = [1; 1; 1]*parms.k6/(parms.k1*parms.insulin_sensitivity);

Ip = [Aeq1(1:2,:); OP_I] \ [beq1(1:2); I(2:end-1)];
G = l - K*Ip;
