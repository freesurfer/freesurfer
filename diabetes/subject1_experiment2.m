BG1 = [173 167.5 167.5 178.5 164 217.5 222.5 243 265.5 241.5 237.5 238.5 231.5 233.5 215.5 200.5 189 184.5 179.5 174.5 162.5 159 153.5 147.5 140 135];

BG0 = [143.5 147.5 155.5 166 216 222.5 240 248.5 264.5 294.5 309.5 327.5 332 330.5 326 325 326 326.5 324 310.5 312.5 306.5 296 293.5 280.5 273.5];
time = [-10 0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240];

plot(time, BG0, 'r');
hold on ;
plot(time, BG1, 'g');
ylabel('Blood glucose level (mg/dL)',  'fontsize', 16, 'fontweight', 'bold');
xlabel('time (min)', 'fontsize', 16, 'fontweight', 'bold');
set(gca,  'fontsize', 16, 'fontweight', 'bold');
set(get(gca, 'children'), 'linewidth', 6);
