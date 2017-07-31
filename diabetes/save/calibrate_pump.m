
% A is active insulin
% M is measured blood glucose level
% F is carbs (in grams)
% D is insulin dose given
% I  - BG sensitivity to carbs (amount BG is raised/gram carb)
% S - sensitivity factor (amount BG drops per/unit insulin)
% e - error in basal rate
% C - carb ratio (grams of carbohdyrates/unit insulin

init_medtronic

for subject_index=1:nsubjects
%for subject_index=10:10

extract_subject_data ;

% find pairs of measurements with carbs given then at least two hours
% to let them be fully absorbed
find_carb_pairs ;

min_rms = 10000;

for S=20:300
  if (mod(S,50)==0)
%    disp(sprintf('processing S = %d', S)) 
  end
  for e=-.5:.025:.5
     for I=5:30
       	   rms = compute_bg_error_pairs(M,D,F,dt,S,e,I);
	   if (rms < min_rms)
	      min_rms = rms ;
	      Sbest(subject_index) = S ;
	      ebest(subject_index) = e ;
	      Ibest(subject_index) = I ;
	      Cbest(subject_index) = S/I ;
%	      disp(sprintf('new optimum %2.0f found at S = %d, e=%2.2f, I=%d, C=%d', rms,S,e,I,C)) 
%	      drawnow;
	end
     end
  end
end

ind = find(sensitivity_vec > 0) ;
Spump = sensitivity_vec(ind(1)) ;
ind = find(carb_ratio_vec > 0) ;
Cpump = carb_ratio_vec(ind(1)) ;
Ipump = Spump / Cpump ;

disp(sprintf('subject %d: global optimum %2.0f found at S = %2.0f (%2.0f), e=%2.3f, I=%2.0f (%2.0f), C=%2.1f (%2.1f)', subject_index,min_rms,Sbest(subject_index),Spump,ebest(subject_index),Ibest(subject_index),Ipump,Cbest(subject_index),Cpump)) ;

if (0)
   compute_bg_error_pairs(M,D,F,dt,Sbest(subject_index),ebest(subject_index),Ibest(subject_index),1)
   compute_bg_error_pairs(M,D,F,dt,Spump,0,Ipump,1)
   compute_bg_error_pairs(M,D,F,dt,80,0,7,1)

end

end  % subject loop
