function time_series = extract_time(S, sind)
% time_series = extract_time(S, subject_index)
% convert time strings into a vector that is the # of minutes since
% the start of the trial
time_ind = 4 ;

time_strings = S(sind).rows(:,time_ind) ;

month = 1 ;
day = 2 ;
year = 3 ;
hour = 4 ;
min = 5 ;
sec = 6 ;
AM_PM = 7 ;

ntps = length(time_strings) ;
min_index = -1 ;

days_to_minutes = 1.0 / (datenum(0,0,0,0,1,0) - datenum(0,0,0,0,0,0)) ;
time_series = zeros(size(time_strings)) ;

for t=1:ntps
    elts = sscanf(char(time_strings(t,:)), '%d/%d/%d %d:%d:%d %cM') ;
    if (t == 15)
       disp('break') ;
    end
    if (length(elts) < sec)
       time_series(t) = NaN ;
    else
       if (elts(AM_PM) == 'A')    % change 12AM times to be 0 hour
       	  if (elts(hour) == 12)
	      elts(hour) =0 ;
          end
       end
       if (elts(AM_PM) == 'P')  % change >12PM times to be plus 12
       	  if (elts(hour) < 12)
	    elts(hour) = elts(hour) + 12 ;
          end
       end
       days_since_beginning = datenum(elts(year), elts(month), elts(day),elts(hour), elts(min), elts(sec));
       if (min_index < 0)
          min_index = t ;
       	  minutes_at_start = days_to_minutes * days_since_beginning ;
       end
       time_series(t) = (days_since_beginning*days_to_minutes - minutes_at_start)  ; ;
    end
end 

disp(sprintf('min index %d, first time %d (%s)\n', min_index,time_series(min_index), char(time_strings(min_index,:))));
