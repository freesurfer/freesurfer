function lme_timePlot(time,y,ni,subjs)
% lme_timePlot(time,y,ni,subjs)
%
% Spaghetti time plot of longitudinal data.
%
% Input
% time: Time covariate.
% y: Ordered data vector (according to time for each subject).
% ni: Vector whose entries are the number of repeated measures for each
% subject in the study (ordered according to y).
% subjs: A vector with the indices of a particular subset of subjects. If
% not specify, time series data for all subjects will be ploted.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:11 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:11 $
%    $Revision: 1.1.2.2 $
%  
if nargin < 3 
    error('Too few inputs');   
elseif nargin < 4
    t = time;
    dat = y;
else
    n = sum(ni(subjs));
    t = zeros(n,1);
    dat = zeros(n,1);
    pos1 = 0; 
    for i=1:length(subjs)
        pos2 = sum(ni(1:subjs(i)-1));
        t(pos1+1:pos1+ni(subjs(i))) = time(pos2+1:pos2+ni(subjs(i)));
        dat(pos1+1:pos1+ni(subjs(i))) = y(pos2+1:pos2+ni(subjs(i)));
        pos1 = pos1+ni(subjs(i));
    end;
    ni = ni(subjs);
end;
figure('Name','Time plot');
hold on;
pos1 = 0;
for i=1:length(ni)
   plot(t(pos1+1:pos1+ni(i)),dat(pos1+1:pos1+ni(i)),'-ok');
   pos1 = pos1+ni(i);
end; 
hold off;   
   