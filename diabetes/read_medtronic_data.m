function medtronic_subject_data = read_medtronic_data(fname, debug)

if (nargin == 1)
   debug = 0 ;
end

tmp = load(fname) ;
medtronic_subject_data = tmp.S ;

if (debug)
   for f=1:size(medtronic_subject_data(1).fields,2)
       disp(sprintf('col %d: %s', f, cell2mat(medtronic_subject_data(1).fields(:,f)))) ;
   end
end

