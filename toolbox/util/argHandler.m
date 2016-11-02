%% Matlab Library
%  Bryan Urquhart
%
%  This function is designed to parse and chop up input arguments that are
%  delivered via varargin. This function is useful only when the arguments
%  come in string-value pairs - meanging the length must be even. This
%  function will throw an error if the length is not even.
%
%  This function operates in two modes, depending on whether defaults are specified.
%
%	When no defaults are specified, argHandler will split the args into a
%	two column cell array from a one column cell array, and it will convert
%	all incoming strings to lower case to remove any case sensitivity
%	issues.
%
%	When defaults are specified (struct of defaults expected), argHandler
%	will return a struct with final values for all of the parameters in the
%	defaults struct.  Matches to field names are performed
%	case-insensitively.
%
function [ argv ] = argHandler( args , defs )
%% Verify input length is even
%
% Get length of args
N = length( args );
%
% Check to make sure its even
if( mod(N,2) ~= 0 )
  error( 'Not an even number of input arguments' );
end
%
%% Reshape the arguments if there are an even number
%

% reshape and store in argv
args = reshape( args , 2 , N/2 )';

if( nargin < 2)
	% Lower case conversion
	argv = args;
	argv(:,1) = strtrim(lower(args(:,1)));
else
	argv = defs;
	fieldNs = fieldnames(defs);
	for i = 1:size(args,1); fn = args{i,1};
		idx = strcmpi(fn,fieldNs);
		if(sum(idx) > 1)
			warning('argHandler:cases','argHandler can''t be used to process options that are identical except for case');
		elseif(sum(idx) == 1)
			fn = fieldNs{idx};
			argv.(fn) = args{i,2};
		else
			warning('argHandler:unusedArg','Parameter %s isn''t valid for the calling function',fn);
		end
	end
end
