function s = parseConf( conf )
%PARSECONF parses output from readConf()
%   parseConf looks at each config string and tries to convert to numbers, and if it failes it leaves them as strings.
%
% Usage:
%	conf = parseConf( readConf( conf_path ) )
%	s will be a struct with fields of varying class
%	conf_path should be the path to the file
fn = fieldnames(conf);

% Limit the possible characters for a number to a smaller set. This will
% exclude the time format string (include ":" character) from being converted incorrectly.
% If you try str2num("2012-10-11 08:00:00") you will get 1991 !!!
numericCharacters = '0123456789ij.- ';

for idx = 1:length(fn)
    convertFlag = 1;
    
    if all(ismember(conf.(fn{idx}),numericCharacters))
        % extra step to carefully check the string including '-' character
        % (e.g. '2012-01-02') that could be a date format instead of a
        % negative number (e.g. '-2' or '1-2j'). If datenum can convert it,
        % then it's a date format. Otherwise will try to convert it as a
        % number.
        if find(conf.(fn{idx})=='-',1)
            try
                datenum(conf.(fn{idx}));
                convertFlag = 0;
            catch
            end
        end
        
        if convertFlag
            try
                [ s.(fn{idx}) success ] = str2num(conf.(fn{idx})); %#ok<ST2NM>
                if(success), continue; end
            catch %#ok<*CTCH>
            end
        end
    end
    
    s.(fn{idx}) = conf.(fn{idx});
end