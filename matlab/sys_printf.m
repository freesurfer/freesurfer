function [s] = sys_printf(varargin)
%
% NAME
%
%  function [s] = sys_printf(format, ...)
%
% ARGUMENTS
% INPUT
%	format, ...	string		C-style format string to print
% OPTIONAL
%	
% OUTPUT
%	s               string		syslog-style prefixed string
%
% DESCRIPTION
%
%	This function prepends <format> with syslog_prefix().
%
% PRECONDITIONS
%
% POSTCONDITIONS
%	o string is returned and printed.
%
% NOTE:
%
% HISTORY
% 12 December 2008
% o Initial design and coding.
%
sfrmt   = sprintf(varargin{:});
s       = sprintf('%s %s', syslog_prefix(), sfrmt);
fprintf(1, '%s', s);
