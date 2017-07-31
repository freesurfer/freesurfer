fs = 16;
G = 20:500;
Gt=100;
erate = 10 ;
Glow = 50 ;
sigma = 20;
loss = (G-Gt).^2./sigma.^2;
escale = 1000;
eg = escale*exp(-(G-Glow)/erate) ;
loss = loss+eg ;
plot(G,loss);
set(gca,'xlim', [50 400], 'fontsize', fs, 'fontweight', 'bold')
xlabel('Blood Glucose Level (mg/dL)', 'fontsize', fs, 'fontweight', 'bold')
ylabel('Loss (A/U)', 'fontsize', fs, 'fontweight', 'bold')
ch = get(gca, 'children');
set(ch, 'linewidth', 10)
