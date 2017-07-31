function [BG, insulin] = compute_BG_and_insulin(parms, plotit)
% function [BG, insulin] = compute_BG_and_insulin(parms, plotit)

%% Set variables

% Optimization cost function, minimizing:
% 1: blood glucose,  2: L1 error,  3: L2 error
optimCF = 2;
T1  = parms.T1/parms.dt;
T2  = parms.T2/parms.dt;

%% Define operators
t = parms.dt*(0:(T2-1))';
allOnes = ones(T2,1);
totalInsulin = (parms.k6*parms.T2 + sum(parms.carb_grams)*parms.carb_sensitivity + parms.G0-parms.G_target)/parms.insulin_sensitivity;

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
% Subject to positive insulin
A1 = -OP_I(1:(T1-1),:);
b1 = zeros(T1-1,1);
% Subject to the limits of G
A2 = K; % This is just the lower limit. For the higher limit, use the following two commented lines instead.
b2 = l - parms.low_limit*allOnes;
% A2 = [-K; K];
% b2 = [high_limit*allOnes - l; l - parms.low_limit*allOnes];
% Subject to the total amount of insulin (currently disabled)
A3 = [];%ones(1,T2-2)*OP_I;
b3 = [];%totalInsulin;
% Subject to the positive Ip
Ip_min = zeros(T2,1);
% Steady state at the beginning and end
Aeq1 = [1 zeros(1,T2-1); zeros(1,T2-1) 1];
beq1 = [1; 1]*parms.k6/(parms.k1*parms.insulin_sensitivity);

Aeq1 = [1 zeros(1,T2-1); 0 1 zeros(1,T2-2); zeros(1,T2-1) 1];
beq1 = [1; 1; 1]*parms.k6/(parms.k1*parms.insulin_sensitivity);

% Glucose target at the end
Aeq2 = -K(end,:);
beq2 = parms.G_target - l(end);
% Insulin schedule to be just basal from parms.T1 to parms.T2
Aeq3 = OP_I(T1:end,:);
beq3 = ones(T2-T1-1,1)*parms.k6/parms.insulin_sensitivity;

switch optimCF
    case 1
        % Linear programming, without L1 error
        A = [A1; A2; A3];
        Aeq = [Aeq1; Aeq2; Aeq3];
        Ip = linprog(-K'*allOnes, A, [b1; b2; b3], Aeq, [beq1; beq2; beq3], Ip_min, [], [], optimoptions('linprog', 'Algorithm', 'interior-point'));
    case 2
        % Linear programming
        A = [A1; A2; A3];
        Aeq = [Aeq1; Aeq2; Aeq3];
        Ip_absIp = linprog([0*allOnes; allOnes], [A zeros(size(A,1),T2); K, -eye(T2); -K, -eye(T2)], [b1; b2; b3; l-parms.G_target*allOnes; -l+parms.G_target*allOnes], [Aeq zeros(size(Aeq,1),T2)], [beq1; beq2; beq3], [Ip_min; 0*allOnes], [], [], optimoptions('linprog', 'Algorithm', 'interior-point'));
        Ip = Ip_absIp(1:T2);
    case 3
        % Quadratic programming
        Ip = quadprog(K'*K, -K'*(l-parms.G_target*allOnes), [A1; A2; A3], [b1; b2; b3], [Aeq1; Aeq2; Aeq3], [beq1; beq2; beq3], Ip_min, [], [], optimoptions('quadprog', 'Algorithm', 'interior-point-convex'));
end

% Other quantities
G = l - K*Ip;
BG = G ;
I = [parms.k6/parms.insulin_sensitivity; OP_I*Ip; parms.k6/parms.insulin_sensitivity];
If = ((OP_d + (parms.k2+parms.k1)*OP_Id)/parms.k3) * Ip;
ind = find(abs(I) < 1e-6) ;
I(ind) = 0 ;  % makes looking at insulin easier
insulin = I ;




% only plotting stuff after here
if (nargin >= 2 && plotit)
% baseline.Gp_t doesn't exit, and the lines containing it have been commented.
%disp(['Delay: ' num2str(parms.carb_delay) ', mean(G): ' num2str(mean(baseline.Gp_t(1:T1))) ', mean(opt): ' num2str(mean(G(1:T1))) ', improvement: ', num2str(mean(baseline.Gp_t(1:T1))-mean(G(1:T1)))])
%disp(sprintf('delay %2.1f, mean(G)=%2.1f, mean(opt)=%2.1f, improvement = %2.1f', parms.carb_delay,mean(baseline.Gp_t(1:T1)),mean(G(1:T1)), mean(baseline.Gp_t(1:T1))-mean(G(1:T1))));

%% Plot the results
%figure('name', sprintf(['carb delay=' num2str(parms.carb_delay) ', G0=%2.0f, mean improvement=%2.0f-->%2.0f, %2.1f'], parms.G0, mean(baseline.Gp_t(1:T1)),mean(G(1:T1)),mean(baseline.Gp_t(1:T1))-mean(G(1:T1))))
hold on ;
%plot(baseline.Gp_t, 'r') ;
hx = plotyy(t, G, t, I);
set(get(gca, 'children'), 'linewidth', 6)
%plot(hx(1), t, Ip, 'r');
%legend('Insulin schedule', 'Plasma insulin', 'Blood glucose')
xlabel('Time (mins)', 'fontsize', 16, 'fontweight', 'bold');
ylabel(hx(1), 'Blood glucose (mg/dL)', 'fontsize', 16, 'fontweight', 'bold')
ylabel(hx(2), 'Insulin (U)', 'fontsize', 16, 'fontweight', 'bold')
set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [0 parms.T1],'ytick', [0:50:400],'ylim', [0 425]) ;
ln1 = line([0 parms.T2], [ parms.G_target parms.G_target]);
ln2 = line([0 parms.T2], [parms.low_limit parms.low_limit]);
set(ln1, 'linewidth', 6, 'linestyle', '-.', 'color', 'g') ;
set(get(gca, 'children'), 'linewidth', 6)
set(ln2, 'linewidth', 6, 'linestyle', '-.', 'color', 'r') ;
%l = legend('normal BG response', 'optimal BG response','optimal insulin schedule');
%set(l, 'fontsize', 16, 'fontweight', 'bold');
axes(hx(2))
hold on
set(get(gca, 'children'), 'linewidth', 6)
set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [0 parms.T1],'ytick', [0:.5:4.0],'ylim', [0 4.25]) ;
hold off ;

end



if 0

clf
plot(Ip, 'm') ;
hold on ;
plot(If, 'c') ;
plot(I, 'g') ;
ch = get(gca, 'children')  ;
set(ch, 'linewidth', 6) ;
set(ch(1), 'color', [0 .6 0]) ;
set(gca, 'fontsize', 16, 'fontweight', 'bold', 'xlim', [0 parms.T1]) ;
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('Insulin (U)', 'fontsize', 16, 'fontweight', 'bold') ;

legend('plasma insulin', 'interstitial insulin', 'insulin delivery');
hold off;

end
