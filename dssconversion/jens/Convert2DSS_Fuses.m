function f = Convert2DSS_Fuses(fs)
% f = Convert2DSS_Fuses(fs)
%
% PURPOSE : Converts fuse from object to OpenDSS
%
%
% INPUT :   fs: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS fuse object
%           Not documented in OpenDSS manual V7.6(November 2012)
%           Information in help file of OpenDSS V7.6.2.1:

f(length(fs)) = dssfuse;

% Name without space
            f.Name = fs.DeviceNumber;

            % 
%             if(fs(i).FuseIsOpen), f(i).Action = 'Open';
%             else f(i).Action = 'Close';
%             end
            
            if ~isempty(strfind(fs.ClosedPhase,'A')) || ~isempty(strfind(fs.ClosedPhase,'B')) || ~isempty(strfind(fs.ClosedPhase,'C'))
                f.Action = 'Close';
            else
                f.Action = 'Open';
            end


    % 		% TODO: check if it can actually work when RatedCurrent is actual phase amps
            f.RatedCurrent = str2double(fs.Amps);
    % 		f(i).RatedCurrent = str2double(fs(i).AmpRating);
    % 
            % Handle spaces in sectionId
            f.MonitoredObj = ['Line.' fs.SectionID];
%             if(~fs(i).NearFromNode)
%                 f(i).SwitchedTerm = 2;
%             end
            if ~strcmp(fs.Location,'L')
                f.SwitchedTerm = 2;
            end