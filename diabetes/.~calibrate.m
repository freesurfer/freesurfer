read_medtronic_data

estimateIOBTable;
iob_table = activeInsulin_3Hr ;
time_series = extract_time(S,1) ;
carb_vec = extract_carbs(S,1,time_series) ;
bolus_vec = extract_insulin_bolus(S,1,time_series) ;
basal_vec = extract_basals(S,1,time_series) ;
active_insulin = extract_active_insulin(S,1,time_series,bolus_vec,basal_vec,iob_table3) ;
