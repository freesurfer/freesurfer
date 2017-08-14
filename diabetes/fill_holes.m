function [filled, orig_edited replace] = fill_holes(orig_data, filled_data)

% filled = fille_holes(orig_data, filled_data)

[nvec, ntps] = size(orig_data) ;


filled = orig_data ;
if (nargin == 1)
   filled_data = orig_data ;
end

% first remove rows that have less than 3 actual measurements
for row=1:nvec
    ngood = 0 ; 
    replace(row) = 0;
    for t=1:ntps
        if (isnan(orig_data(row,t)) == 0)
	   ngood = ngood + 1 ;
	end
    end
    if (ngood < 3)
       replace(row) = 1 ;
    else
	if (isnan(orig_data(row,1)) && isnan(orig_data(row,2)))
	   replace(row) = 1;
	end
	if (isnan(orig_data(row,end)) && isnan(orig_data(row,end-1)) && isnan(orig_data(row,end-2)))
	   replace(row) = 1;
	end
    end
end
dind = find(replace > 0) ;
total_deleted = length(dind) ;

filled_compacted = zeros(nvec-total_deleted, ntps) ;
orig_compacted = zeros(nvec-total_deleted, ntps) ;
skipped = 0 ;
for row=1:nvec
    if (replace(row) == 0)  % copy this row
       orig_compacted(row-skipped,:) = orig_data(row,:) ;
       filled_compacted(row-skipped,:) = filled_data(row,:) ;
    else
	skipped = skipped+1 ;
    end
end
orig_data = orig_compacted ; filled_data = filled_compacted ;
filled = filled_data ;
[nvec, ntps] = size(orig_data) ;

    
if (nvec > 1)
   filled = filled_data ;
   for n=1:nvec
      filled(n,:) = fill_holes(orig_data(n,:), filled_data(n,:));
   end
else
   ctrl = zeros(size(orig_data)) ;
   for t=1:ntps
      if(isnan(orig_data(t)) == 0)
        ctrl(t) = 1 ;
      end
   end

   filled = soap_bubble(orig_data, ctrl, 100) ;

   % fill in values at left end if NaNs there (soap bubble doesn't extrapolate well
   if (isnan(orig_data(1)))
     for t=1:ntps
        if ((isnan(orig_data(t)) == 0) && t >= 2)
       	  deriv = filled(t) - filled(t+1) ;
	  for t1=t-1:-1:1
	      filled(t1) = filled(t1+1) + deriv ;
	  end
       	  break ;
       end
     end
   end

   % fill in values at right end if NaNs there (soap bubble doesn't extrapolate well
   if (isnan(orig_data(end)))
     for t=ntps:-1:1
       if ((isnan(orig_data(t)) == 0) && t < ntps)
       	  deriv = filled(t) - filled(t-1) ;
	  for t1=t+1:ntps
	      filled(t1) = filled(t1-1) + deriv ;
	  end
       	  break ;
       end
     end
   end
end


orig_edited = orig_data ;
