function lme_lowessPlot(time,y,bw,group)
% lme_lowessPlot(time,y,bw,group)
%
% Time plot of the lowess regression curve estimating the trend in the
% mean response over time (Depends on the curve fitting toolbox).
%
% Input
% time: Time covariate (ordered according to time for each subject) .
% y: Ordered data vector (according to time for each subject).
% bw: Fraction of data points to be included in the sliding window of the
% lowess smoothing method (band width, default 0.6).
% group: A covariate indicating the corresponding group for each subject.
% Default ones(length(time),1) (in this case will also plot the sample data
% points).
%
%
% Original Author: Jorge Luis Bernal Rusiel 
% References:  W.S.Cleveland, "Robust Locally Weighted Regression and 
% Smoothing Scatterplots", _J. of the American Statistical Ass._, Vol 74, 
% No. 368 (Dec.,1979), pp. 829-836.
% http://www.math.tau.ac.il/~yekutiel/MA%20seminar/Cleveland%201979.pdf
%  
if nargin < 2
    error('Too few inputs');
else
    if nargin < 4
        group = ones(length(time),1);
        if nargin < 3
            bw = 0.6;
        end;
    end;
end;
un = unique(group);
figure('Name','Lowess plot');
hold on;
colors = {'b','k','g','c','m','r','y'};
if length(un)==1
    plot(time,y,'ok');    
end;
for i=1:length(un)
    t = time(group == un(i));
    dat = y(group == un(i));
    [st,ix] = sort(t);
    ybw = smooth(st,dat(ix),bw,'rlowess');
    plot(st,ybw,['-' colors{mod(i,7)+1}],'LineWidth',4);
end;
hold off;   

