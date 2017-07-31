function loss = compute_loss(BG_t, G_target, risk_matrix)
% function loss = compute_loss(BG_t, G_target)


loss_matrix = risk_matrix .* BG_t ;
loss_matrix(:, 300:end) = 0 ;
loss = mean(mean(loss_matrix));

if 0
plot(G,loss);
set(gca,'xlim', [50 400], 'fontsize', fs, 'fontweight', 'bold')
xlabel('Blood Glucose Level (mg/dL)', 'fontsize', fs, 'fontweight', 'bold')
ylabel('Loss (A/U)', 'fontsize', fs, 'fontweight', 'bold')
ch = get(gca, 'children');
set(ch, 'linewidth', 10)
end
