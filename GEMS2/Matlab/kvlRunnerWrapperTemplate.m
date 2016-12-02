function varargout = kvl@runnerName@( varargin )
%
%

numberOfOutputs = @runnerNumberOfOutputs@;
if nargout > numberOfOutputs
  error( 'Too many outputs requested' )
end

output = cell( 1, numberOfOutputs );
[ output{:} ] = kvlGEMSMatlab( '@runnerName@', varargin{:} );

for i = 1 : nargout
  varargout{ i } = output{ i };
end
