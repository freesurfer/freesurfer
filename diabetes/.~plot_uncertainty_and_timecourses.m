BG_opt = computeUpsampled(insulin_opt3, parms) ;
BG_init = computeUpsampled(insulin_init, parms) ;

%%%%% plot the risk minimizing uncertainty
figure('name', 'risk minimizing', 'position', [786 604 1379  689]) ;
h2 = imagesc(BGopt_prob);
set(gca, 'ylim', [0 350], 'xlim', [0 parms.T1], 'fontsize', 16, 'fontweight', 'bold');
axis xy
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold') ;
colorbar;
title('uncertainty curves for risk-minimizing schedule', 'fontsize', 16, 'fontweight', 'bold') ;
print -dtiff risk_minimizing_uncertainty.tif

%%%%% plot the deteriministic uncertainty
chigh = max(BGopt_prob(:)) ;
figure('name', 'L1 minimizing', 'position', [2174   601 1248  691]) ;
h1 = imagesc(BG_prob, [0 chigh]);
set(gca, 'ylim', [0 350], 'xlim', [0 parms.T1], 'fontsize', 16, 'fontweight', 'bold');
axis xy
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold') ;
colorbar;
title('uncertainty curves for deterministic schedule', 'fontsize', 16, 'fontweight', 'bold') ;
print -dtiff deterministic_L1_uncertainty.tif

%%%% plot the timecourses
figure('name', 'BG for risk and L1-minimizing schedules') ;
plot(BG_opt, 'g')
hold on;
plot(BG_init, 'r')
hold off;
legend('risk minimizing', 'deterministic');
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold') ;
ylabel('BG (mg/dL)', 'fontsize', 16, 'fontweight', 'bold') ;
set(gca, 'xlim', [0 parms.T2], 'fontsize', 16, 'fontweight','bold') ;
ch = get(gca, 'children') ;
set(ch, 'linewidth', 4);
%print -dtiff BG_timecourses.tif

