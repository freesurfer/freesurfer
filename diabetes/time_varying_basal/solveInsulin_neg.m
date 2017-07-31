low_limit = 70;
high_limit = 180;
carb_delay = 30;
carb_grams = 50;
G_target = 100;   % desired blood glucose (Gp) level

% compute response without optimal schedule
parms.insulin_delay = 15 ;
parms.k1 = .021 ;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
parms.k2 = .001 ;   % rate at which insulin moves from plasma to interstitial fluid (units/min)
parms.k3 = .021 ;   % rate at which insulin moves from fluid to plasma (units/min)
parms.k4 = .05 ;    % rate at which carbs are metabolized from stomach to blood (grams/min)
parms.k5 = parms.insulin_sensitivity ;   % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)
parms.k6 = .5 ;    % rate at which liver drips glucose into blood plasma (glucose/min)
parms.low_limit = 70 ;
parms.high_limit = 180 ;
parms.basal = parms.k6/parms.insulin_sensitivity;   % basal insulin
parms.insulin = parms.basal*ones(size(parms.time_series)) ;
parms.G_target = G_target ;   % desired blood glucose (Gp) level
parms.carb_sensitivity = 7 ;
parms.insulin_sensitivity = 180 ;  % glucose/unit insulin
parms.carb_delay = carb_delay ;
parms.carb_grams = carb_grams ;
parms.time_series = 1:T2;
parms.G0 = G_target ;
carb_ratio = parms.insulin_sensitivity/parms.carb_sensitivity ;
bolus = parms.carb_grams/carb_ratio;                    % insulin to counteract carbs
bolus = bolus + (parms.G0-parms.G_target)/parms.insulin_sensitivity ; % insulin to correct for current BG level

parms.insulin(parms.insulin_delay) = parms.insulin(parms.insulin_delay) + bolus ;

baseline = simulate_timecourse(parms) ;





%% Set variables
T1 = 300;
T2 = 600;
k4 = .021;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
k3 = .001;   % rate at which insulin moves from plasma to interstitial fluid (units/min)
k5 = .021;   % rate at which insulin moves from fluid to plasma (units/min)
k1 = .05;    % rate at which carbs are metabolized from stomach to blood (grams/min)
SI = 180;    % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)
SC = 7;      % Carb sensitivity
k2 = .5 ;    % rate at which liver drips glucose into blood plasma (glucose/min)
G0 = G_target;

%% Define operators
t = (0:(T2-1))';
allOnes = ones(T2,1);
totalInsulin = (k2*T2 + carb_grams*SC + G0-G_target)/SI;

OP_Id = eye(T2);
OP_d = (diag(ones(T2-1,1),1) - diag(ones(T2-1,1),-1))/2;
OP_d2 = diag(ones(T2-1,1),1) + diag(ones(T2-1,1),-1) - 2*eye(T2);
OP_int = tril(ones(T2), 0) - .5*eye(T2) - .5*[ones(T2,1) zeros(T2,T2-1)];

OP_Id = OP_Id(2:end-1, :); % Ignoring the boundary
OP_d = OP_d(2:end-1, :);
OP_d2 = OP_d2(2:end-1, :);
OP_I = (1/k5)*OP_d2 + (1 + (k3+k4)/k5)*OP_d + k4*OP_Id; % Pump insulin

%% Optimize
% Assuming the blodd glucose is G = l - K*Ip, we minimize ||G-G-target||^2
K = SI*k4*OP_int;
l = G0*allOnes + k2*t + SC*carb_grams*(1-exp(-k1*(t-carb_delay))).*(t>=carb_delay);
% Subject to positive insulin
A1 = -OP_I(1:(T1-1),:);
b1 = zeros(T1-1,1);
A1=[];
b1=[];
%b1 = -1*ones(T1-1,1);
% Subject to the limits of G
A2 = K; % This is just the lower limit. For the higher limit, use the following two commented lines instead.
b2 = l - low_limit*allOnes;
% A2 = [-K; K];
% b2 = [high_limit*allOnes - l; l - low_limit*allOnes];
% Subject to the total amount of insulin (currently disabled)
A3 = [];%ones(1,T2-2)*OP_I;
b3 = [];%totalInsulin;
% Subject to the positive Ip
Ip_min = zeros(T2,1);
% Steady state at the beginning and end
Aeq1 = [1 zeros(1,T2-1); zeros(1,T2-1) 1];
beq1 = [1; 1]*k2/(k4*SI);

Aeq1 = [1 zeros(1,T2-1); 0 1 zeros(1,T2-2); zeros(1,T2-1) 1];
beq1 = [1; 1; 1]*k2/(k4*SI);

% Glucose target at the end
Aeq2 = -K(end,:);
beq2 = G_target - l(end);
% Insulin schedule to be just basal from T1 to T2
Aeq3 = OP_I(T1:end,:);
beq3 = ones(T2-T1-1,1)*k2/SI;

% Quadratic programming
Ip = quadprog(K'*K, -K'*(l-G_target*allOnes), [A1; A2; A3], [b1; b2; b3], [Aeq1; Aeq2; Aeq3], [beq1; beq2; beq3], Ip_min, [], [], optimoptions('quadprog', 'Algorithm', 'interior-point-convex'));
G = l - K*Ip;
I = [k2/SI; OP_I*Ip; k2/SI];
If = ((OP_d + (k3+k4)*OP_Id)/k5) * Ip;

%% Plot the results
clf
hold on ;
plot(baseline.Gp_t, 'r') ;
hx = plotyy(t, G, t, I);
set(get(gca, 'children'), 'linewidth', 6)
%plot(hx(1), t, Ip, 'r');
%legend('Insulin schedule', 'Plasma insulin', 'Blood glucose')
l = legend('normal BG response', 'optimal BG response','optimal insulin schedule');
set(l, 'fontsize', 16, 'fontweight', 'bold');
xlabel('Time (mins)', 'fontsize', 16, 'fontweight', 'bold');
ylabel(hx(1), 'Insulin (U)', 'fontsize', 16, 'fontweight', 'bold')
ylabel(hx(2), 'Blood glucose (mg/dL)', 'fontsize', 16, 'fontweight', 'bold')
set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [0 300],'ytick', [0:50:350],'ylim', [0 350]) ;
ln1 = line([0 T2], [ G_target G_target]);
ln2 = line([0 T2], [low_limit low_limit]);
set(ln1, 'linewidth', 6, 'linestyle', '-.', 'color', 'g') ;
set(get(gca, 'children'), 'linewidth', 6)
set(ln2, 'linewidth', 6, 'linestyle', '-.', 'color', 'r') ;
axes(hx(2))
hold on
set(get(gca, 'children'), 'linewidth', 6)
set(gca, 'fontsize', 16, 'fontweight', 'bold','xlim', [0 300]) ;
hold off ;


disp(sprintf('mean(G)=%2.1f, mean(opt)=%2.1f, improvement = %2.1f', mean(baseline.Gp_t(1:300)),mean(G(1:300)), mean(baseline.Gp_t(1:300))-mean(G(1:300))));
