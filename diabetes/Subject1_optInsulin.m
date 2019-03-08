%% Set variables

% Optimization cost function, minimizing:
% 1: blood glucose,  2: L1 error,  3: L2 error
optimCF = 1;

a = 1.001;
w = .009;
k1 = .012;    % rate at which carbs are metabolized from stomach to blood (grams/min)
SI = 89;     % insulin sensitivity (fit from medronic IOB table- amount BG is lowered by insulin)
SC = 5.2;    % Carb sensitivity
k2 = .5;%1e-6;    % rate at which liver drips glucose into blood plasma (glucose/min)

T1 = 120;
T2 = 600;
k4 = sqrt(a)*w;   % rate at which insulin moves from plasma to cells (fit from medronic IOB table, unit/min)
k5 = (a*w^2)/k4;   % rate at which insulin moves from fluid to plasma (units/min)
k3 = (a+1)*w - k4 - k5;   % rate at which insulin moves from plasma to interstitial fluid (units/min)
k6 = .05;    % relative rate at which muscle glucose is refilled (1/min)
G_target = 100;   % desired blood glucose (Gp) level
low_limit = 70;
high_limit = 180;

carb_delay = 10;
carb_grams = 50.1;
allG0 = 70:5:180;
clear dispStr
for index = 1:length(allG0)
    G0 = allG0(index);
    
    %% Define operators
    t = (0:(T2-1))';
    allOnes = ones(T2,1);
    %totalInsulin = (k2*T2 + sum(carb_grams)*SC + G0-G_target)/SI;
    
    OP_Id = eye(T2);
    OP_d = (diag(ones(T2-1,1),1) - diag(ones(T2-1,1),-1))/2;
    OP_d2 = diag(ones(T2-1,1),1) + diag(ones(T2-1,1),-1) - 2*eye(T2);
    OP_int = tril(ones(T2), 0) - .5*eye(T2) - .5*[ones(T2,1) zeros(T2,T2-1)];
    
    OP_Id = OP_Id(2:end-1, :); % Ignoring the boundary
    OP_d = OP_d(2:end-1, :);
    OP_d2 = OP_d2(2:end-1, :);
    OP_I = (1/k5)*OP_d2 + (1 + (k3+k4)/k5)*OP_d + k4*OP_Id; % Pump insulin
    
    %exe = 2*(abs(t-200)<50);
    exe = zeros(T2,1);
    
    %% Optimize
    % Assuming the blood glucose is G = l - K*Ip, we minimize ||G-G_target||_2
    % by quadratic programming, or ||G-G_target||_1 by linear programming, or
    % simply G by linear programming/
    K = SI*k4*OP_int;
    l = G0*allOnes + k2*t + exp(-k6*t).*(OP_int*(exe.*exp(k6*t))) - OP_int*exe;
    for nCarb = 1:length(carb_grams)
        l = l + SC*carb_grams(nCarb)*(1-exp(-k1*(t-carb_delay(nCarb)))).*(t>=carb_delay(nCarb));
    end
    % Subject to positive insulin
    A1 = -OP_I(1:(T1-1),:);
    b1 = zeros(T1-1,1);
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
    %     Aeq1 = [1 zeros(1,T2-1); zeros(1,T2-1) 1];
    %     beq1 = [1; 1]*k2/(k4*SI);
%     Aeq1 = [1 zeros(1,T2-1); 0 1 zeros(1,T2-2); zeros(1,T2-1) 1];
%     beq1 = [1; 1; 1]*k2/(k4*SI);
    Aeq1 = [1 zeros(1,T2-1); 0 1 zeros(1,T2-2)];
    beq1 = [1; 1]*k2/(k4*SI);
    % Glucose target at the end
    Aeq2 = -K(end,:);
    beq2 = G_target - l(end);
    % Insulin schedule to be just basal from T1 to T2
    Aeq3 = OP_I(T1:end,:);
    beq3 = ones(T2-T1-1,1)*k2/SI;
    
    switch optimCF
        case 1
            % Linear programming, without L1 error
            A = [A1; A2; A3];
            Aeq = [Aeq1; Aeq2; Aeq3];
            Ip = linprog(-K'*allOnes, A, [b1; b2; b3], Aeq, [beq1; beq2; beq3], Ip_min, [], [], optimoptions('linprog', 'Algorithm', 'interior-point', 'MaxIterations', 1000, 'OptimalityTolerance', 1e-7));
        case 2
            % Linear programming
            A = [A1; A2; A3];
            Aeq = [Aeq1; Aeq2; Aeq3];
            Ip_absIp = linprog([0*allOnes; allOnes], [A zeros(size(A,1),T2); K, -eye(T2); -K, -eye(T2)], [b1; b2; b3; l-G_target*allOnes; -l+G_target*allOnes], [Aeq zeros(size(Aeq,1),T2)], [beq1; beq2; beq3], [Ip_min; 0*allOnes], [], [], optimoptions('linprog', 'Algorithm', 'interior-point'));
            Ip = Ip_absIp(1:T2);
        case 3
            % Quadratic programming
            Ip = quadprog(K'*K, -K'*(l-G_target*allOnes), [A1; A2; A3], [b1; b2; b3], [Aeq1; Aeq2; Aeq3], [beq1; beq2; beq3], Ip_min, [], [], optimoptions('quadprog', 'Algorithm', 'interior-point-convex'));
    end
    
    % Other quantities
    G = l - K*Ip;
    I = [k2/SI; OP_I*Ip; k2/SI];
    If = ((OP_d + (k3+k4)*OP_Id)/k5) * Ip;
    I(abs(I) < 1e-6) = 0 ;  % makes looking at insulin easier
    
    %     delay_BG(index,:) = G ;
    %     delay_I(index,:) = I ;
    %     disp(['Delay: ' num2str(carb_delay) ', mean(G): ' num2str(mean(baseline.Gp_t(1:T1))) ', mean(opt): ' num2str(mean(G(1:T1))) ', improvement: ', num2str(mean(baseline.Gp_t(1:T1))-mean(G(1:T1)))])
    %disp(sprintf('delay %2.1f, mean(G)=%2.1f, mean(opt)=%2.1f, improvement = %2.1f', carb_delay,mean(baseline.Gp_t(1:T1)),mean(G(1:T1)), mean(baseline.Gp_t(1:T1))-mean(G(1:T1))));
    if 1
        figure('name', sprintf(['carb delay=' num2str(carb_delay) ', G0=%2.0f, mean G = %2.0f'], G0, mean(G(1:T1))))
        hold on ;
        %plot(baseline.Gp_t, 'r') ;
        hx = plotyy(t, G, t, I);
        set(get(gca, 'children'), 'linewidth', 6)
        %plot(hx(1), t, Ip, 'r');
        %legend('Insulin schedule', 'Plasma insulin', 'Blood glucose')
        xlabel('Time (mins)', 'fontsize', 16, 'fontweight', 'bold');
        ylabel(hx(1), 'Blood glucose (mg/dL)', 'fontsize', 16, 'fontweight', 'bold')
        ylabel(hx(2), 'Insulin (U)', 'fontsize', 16, 'fontweight', 'bold')
        set(gca, 'fontsize', 16, 'fontweight', 'bold','ytick', 0:50:400,'ylim', [0 425]) ;
        ln1 = line([0 T2], [ G_target G_target]);
        ln2 = line([0 T2], [low_limit low_limit]);
        set(ln1, 'linewidth', 6, 'linestyle', '-.', 'color', 'g') ;
        set(get(gca, 'children'), 'linewidth', 6)
        set(ln2, 'linewidth', 6, 'linestyle', '-.', 'color', 'r') ;
        %l = legend('normal BG response', 'optimal BG response','optimal insulin schedule');
        %set(l, 'fontsize', 16, 'fontweight', 'bold');
        axes(hx(2))
        hold on
        set(get(gca, 'children'), 'linewidth', 6)
        set(gca, 'fontsize', 16, 'fontweight', 'bold','ytick', 0:.5:4.0,'ylim', [0 4.25]) ;
        %xlim([0 T1])
        hold off ;
    end
    [~, ind] = max(I);
    bolus = sum(I(max(1, ind-10):min(ind+10, length(I))));
    dispStr{index} = ['For G0 = ' num2str(G0) ', use a bolus of ' num2str(bolus) ' units at t = ' num2str(t(ind)) ' minutes.'];
    disp(dispStr{index})
end
for index = 1:length(allG0)
    disp(dispStr{index})
end
