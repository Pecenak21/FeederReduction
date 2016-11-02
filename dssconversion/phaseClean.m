function [ numPhase ] = phaseClean( phase )
%PHASECLEAN Clean up phase input and return number of phases for OpenDSS
%conversion. Default numPhase: 3.

% process input
if isempty(phase)
	warning('DSSConversion:phaseClean','Need to specify input. Set to 3.');
	phase = 3;
end

if ischar(phase)
	% remove spaces, neutral and ground notation
	phase = regexprep(lower(phase),'[ng\s]','');
	
	% count number of phases excluding neutral line
	phase = length(phase);
elseif ~isnumeric(phase)
	warning('DSSConversion:phaseClean','Invalid input for phases. Set to 3.');
	phase = 3;
end

numPhase = floor(phase);

if numPhase < 0
	warning('DSSConversion:phaseClean','Invalid numeric input for phases. Set to 3.');
	numPhase = 3;
end
if numPhase > 3
	warning('DSSConversion:phaseClean',['NumPhase = ' num2str(numPhase) ' > 3: More than 3 phases.']);
end

end