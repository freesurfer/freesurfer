% **** global optimum 10.571 found at B1=0.300, B2=0.575 k_carbs=0.0095, k_insulin=0.0100, SC=9.6, SI=96.0, C=10.00, BCI=3.5 hours, carb delay 5 min *****
% **** new optimum 10.604 found at B1=0.300, B2=0.575 k_carbs=0.0095, k_insulin=0.0100, SC=9.6, SI=96.0, C=10.00, BCI=3.0 *****
searching BCI=3.0 (79), basal1 = 0.650

define_insulin_constants

protein = 22 ;
carbs = 50 + protein * protein_to_carb_ratio ;

BG0 = [179 183 192 218 241 281 300.5 313 332 330 347 357 361 352 337 ...
       339.5 334.5 317.5 298.5 289 278.5 270.5 259 249.5 244 236 ];

% from experiment 2
BG1 = [155 159.5 160.5 176.5 242 245.5 253 261.5 268.5 283 268.5 ...
       283.5 295.5 291.5 286 280.5 268.5 260 249 233 233.5 227.5 219 ...
       207.5 202 191.5];


carb_delay = 5 / min_per_hour ;
% subject specific stuff
basal_time = -10:dt:time(end) ;  % in hours
[foo, zero_ind] = min(abs(basal_time*60-0));  % time 0 in spreadsheet

% time0 = 9:23
% basals: 
% 1: 12:00AM-8:00AM  0.6
% 2: 8:00 AM- 1:00PM 0.825

% no first bolus
first_bolus = 0 ;

zero_time = 9+(23/60) ; % 9:23 or so
basal = zeros(size(basal_time));

[foo, basal1_index] = min(abs(zero_time+basal_time-0));
[foo, basal2_index] = min(abs(zero_time+basal_time-8));
% this one doesn't matter as there was no bolus prior to experiment
b1_start_ind = zero_ind ;  % basal one will go from  b1_start_end:bci
[foo, first_bolus_index] = min(abs(zero_time+basal_time-5.75));  
basal1 = (.6) ;
basal2 = (.825) ;
basal(1:basal2_index-1) = basal1 ;
basal(basal2_index:end) = basal2 ;
basal1 = basal1 * dt ;
basal2 = basal2 * dt ;
basal = basal * dt ;  % divide it up into 10min intervals
basal(first_bolus_index) = basal(first_bolus_index) + first_bolus ;
% bolus times
% 9:  5.7
pct_insulin_absorbed_at_end = 0.6 ;
pct_carbs_absorbed_at_end = 0.6;

% these are hourly and will be converted to per-minute later
basal1_min = basal1/dt-.3 ;
basal1_max = basal1/dt+.1;
basal2_min = basal2/dt-.25 ;
basal2_max = basal2/dt+.25 ;


bolus = 5.0;
bolus_schedule = zeros(size(time)) ;
[foo, bolus1_index] = min(abs(time-0.0));
bolus_schedule(bolus1_index) = bolus_schedule(bolus1_index)+bolus ;
[foo,time_start_index] = min(abs(basal_time-time(1)));
[foo,time_end_index] = min(abs(basal_time-time(end)));

SI_best = 50 ; % insulin sensitivity, BG/insulin
b1_best = basal1 ;
b2_best = basal2 ;
C_best = 10 ; % carb ratio, carbs/insulin
SC_best = SI_best / C_best ;
k_carbs_best = .005 ;
bci_best = 1 ;
k_insulin_best = .025 ;


bci_range = (zero_ind+(round(0.0/dt)):round(0.5/dt):(zero_ind+(round(3.5/dt))));
C_range = 5:1:12 ;
SI_range = 60:2:150 ;
SC_range = 5:1:15;

fit_insulin_parameters
