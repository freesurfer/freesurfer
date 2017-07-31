read_medtronic_data

estimateIOBTable;
sind = 1 ;
iob_table = activeInsulin_3Hr ;
time_series = extract_time(S,sind) ;
carb_vec = extract_carbs(S,sind,time_series) ;
bolus_vec = extract_insulin_bolus(S,sind,time_series) ;
basal_vec = extract_basals(S,sind,time_series) ;
active_insulin = compute_active_insulin(S,sind,time_series,bolus_vec,basal_vec,iob_table3) ;
bg_readins = S(sind).rows(:,bg_reading_ind) ;

