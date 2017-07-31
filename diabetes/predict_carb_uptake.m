function carbs_timecourse = predict_carb_uptake(carbs_t0, k_blood_to_stomach, k_stomach_to_blood, time,carb_time)
% carb_timecourse = predict_carb_uptake(carbs_t0, k_blood_to_stomach, k_stomach_to_blood, time,carb_time)

if (nargin < 5)
    carb_time = 0 ;  
end

carbs_timecourse = zeros(size(time)) ;

carb_time = carb_time * 60;
time = time * 60 ;   % change it to minutes
t0 = min(time) ;
t1 = max(time) ;
dt = .05 ;
i = 1 ;
ic = 2 ;
stomach = 0 ;
blood = 0 ;
carbs_timecourse(ic) = 0 ;
for t=t0:dt:t1
   carbs_highres(i) = blood ;
   if (((t-dt) < carb_time) & (t >= carb_time))
     stomach = carbs_t0 ;
   end

   dblood =   dt *  k_blood_to_stomach * blood ;
   dstomach = dt *k_stomach_to_blood * stomach ;
   stomach = stomach + dblood - dstomach ;
   blood = blood + dstomach - dblood ;
    
   if ( ( (t-dt) < time(ic)) & (t >= time(ic)))
     carbs_timecourse(ic) = blood ;
     ic = ic + 1 ;
     if (ic > length(time))
       break ;
   end
  end
  i = i+1;
end
