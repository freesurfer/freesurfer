function [] = colprintf(aC, astr_LC, varargin)
%
% NAME
%
%  function [] = colprintf(aC, astr_LC, format, ...)
%
% ARGUMENTS
% INPUT
%       aC              class/col width         if class - queried for col
%                                               widths, else parsed for 
%                                               width spec
%       astr_LC         string                  Left-column 'intro' text
%       format, ...     string                  C-style format string to print
%                                               in right column.
%
% OPTIONAL
%
% DESCRIPTION
%
%       Prints two-tuple text inputs in two-columns, with column widths
%       defined in the hosting class <aC> or alternatively if aC is a 
%       string of form "<LC>.<RC>"
%
% NOTE:
%
% HISTORY
% 18 September 2009
% o Initial design and coding.
%

	LC              = 40;
	RC              = 40;
        verbosity       = 1;

        if isobject(aC)
            LC          = aC.m_LC;
            RC          = aC.m_RC;
            verbosity   = aC.m_verbosity;
        else
            [str_LC, str_RC]    = strtok(aC, ';');
            LC                  = str2num(str_LC);
            RC                  = str2num(str_RC);
        end

        sfrmt   = sprintf(varargin{:});
        
        if length(astr_LC) & verbosity
            fprintf(1, '%s', sprintf('%*s', LC, astr_LC));
        end
        if length(sfrmt)   & verbosity
            fprintf(1, '%s', sprintf('%*s', RC, sfrmt));
        end
end
