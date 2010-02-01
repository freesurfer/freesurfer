function [s] = syslog_prefix(varargin)
%
% NAME
%
%  function [s] = syslog_prefix()
%
% ARGUMENTS
% INPUT
%
% OPTIONAL
%	
% OUTPUT
%	s	string		syslog-style prefix string
%
% DESCRIPTION
%
%	This function generates a syslog-style prefix string,
%	suitable for time stamping log-type messages
%
% PRECONDITIONS
%
%	o un*x runtime that understands the 'host' command
%
% POSTCONDITIONS
%
%	o s = "YYYY mmm dd hh:mm:ss <host>"
%
% NOTE:
%
% HISTORY
% 02 August 2007
% o Initial design and coding.
%

[str_ret str_hostname]	= system('hostname');
[str_host str_rem]	= strtok(str_hostname, char(10));

clk			= fix(clock);
year			= clk(1);
month			= clk(2);
day			= clk(3);
hour			= clk(4);
minute			= clk(5);
seconds			= clk(6);

switch(month)
    case 1
	str_month	= 'Jan';
    case 2
	str_month	= 'Feb';
    case 3
	str_month	= 'Mar';
    case 4
	str_month	= 'Apr';
    case 5
	str_month	= 'May';
    case 6
	str_month	= 'Jun';
    case 7
	str_month	= 'Jul';
    case 8
	str_month	= 'Aug';
    case 9
	str_month	= 'Sep';
    case 10
	str_month	= 'Oct';
    case 11
	str_month	= 'Nov';
    case 12
	str_month	= 'Dec';
end

s = sprintf('%d %s %02d %02d:%02d:%02d %s', year, str_month, 	...
	day, hour, minute, seconds, str_host);