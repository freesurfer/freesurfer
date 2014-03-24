function [bins1, bins2, counts] = plot_2dhisto(fname)
% [bins1, bins2] =  plot_2dhisto(fname)

fid = fopen(fname, 'r') ;

line = fgetl(fid);
bins1 = sscanf(line, '%f') ;
line = fgetl(fid);
bins2 = sscanf(line, '%f') ;

counts = zeros(length(bins1),length(bins2));
for bin1=1:length(bins1)
    line = fgetl(fid);
    counts(bin1,:) = sscanf(line, '%f');
end

fclose(fid) ;
