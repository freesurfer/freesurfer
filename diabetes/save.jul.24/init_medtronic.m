if (exist('medtronic_subject_data') == 0)
   medtronic_subject_data = read_medtronic_data('allCSV.mat')  ;
end

nsubjects = size(medtronic_subject_data,2) ;
define_medtronic_columns
estimateIOBTable;
