function [r tr] = Convert2DSS_Regulators(rs,trs,phasemap)
% r = Convert2DSS_Regulators(rs)
%
% PURPOSE : Converts regulators from object to OpenDSS, integrated in
% OpenDSSS as transformers
%
%
% INPUT :   rs: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Regulator Controler object (part of transformer)
%           see page 111 in OpenDSS manual V7.4.3(March 2012)

% If the database has a table describing regulator devices, issue a
% warning to let the user know that our script doesn't work with those
% at this time.
% js: looks like there is no db at this time, DevRegulators needs to be
% populated if we get our hands on a db
if(isfield(rs,'DevRegulators'))
    [haveRdevs haveRdevs] = ismember({rs.RegulatorType},{rs.DevRegulators.RegulatorName});
else
    haveRdevs = false(size(rs));
end

r(length(rs)) = dssregcontrol;

if isempty(trs)
    clear trs;
    clear tr;
end

% [busid, idx_bid] = ismember({rs.SectionID},{d.Section.SectionID});

if ~exist('trs','var')
    tr(length(rs)) = dsstransformer;
    tid = [];
    % initialize transformer index;
    tidex = 1;
else
    [busid, tid] = ismember(trs.SectionID,{rs.SectionID});
    % number of existing transformer that match 
    ntr = sum(nonzeros(busid));

    tr( (length(tr)+1) : (length(tr)+ length(rs) -ntr) ) = dsstransformer;

    % initialize transformer index;
    tidex = length(tr) + 1;
end

for i = 1:length(rs)
    ri = r(i);
    rsi = rs(i);

    % create transformer if not exists
    if ~ismember(i,tid)
        tri = tr(tidex);
        % Assign basic properties for a two-winding transformer
        tri.Name = rsi.UniqueDeviceId;
        tri.Buses = {rs.FromNodeID rs.ToNodeID};
        tri.Phases = rsi.ConnectedPhases;
        % don't forget to specify connection style for both windings
        tri.Conns = {rsi.RegulatorConfiguration, rsi.RegulatorConfiguration};
        % select the winding that we're going to put some taps on
        if(rsi.TapsNearFromNode == 1)
            tri.wdg = 1;
        else
            tri.wdg = 2;
        end
        % Fill in some additional transformer properties if they're
        % present in the model
        if(haveRdevs(i))
            RegDev = d.DevRegulators(haveRdevs(i));
            tri.MaxTap = 1+RegDev.RaiseAndLowerMaxPercentage/100;
            tri.MinTap = 2-tri.MaxTap;
            tri.NumTaps = RegDev.NumberOfTaps;
            tri.kVAs = RegDev.RegulatorRatedKva*[1 1];
            tri.kVs = RegDev.RegulatorRatedVoltage*[1 1];
            % I don't think the following two are quite right
            warning('dssconversion:attention','Using some regulator settings that haven''t been tested');
            tri.Rs = RegDev.PercentZOnRegulatorBase*[1 1];
            tri.imag = tri.R *RegDev.RegulatorXRRatio;
        else
            % Mostly the defaults are at least as good as what I can do
            % because I'd be making it up, but there are a couple:
            % if the tap limiter is being used, its range is a
            % reasonable number of taps to use
            % This also doesn't let us guess maxtaps or mintaps because
            % they're pu, not integer.
            if(rsi.TapLimiterActive)
                tri.NumTaps = rsi.TapLimiterHighSetting - rsi.TapLimiterLowSetting;
            end
            % we can try to extract the current/kv ratings from the
            % regulator type name.  In the examples I've seen, these
            % are often supplied as, e.g. '200A 12' or '12.47KV 350'.
            % But this will just end up being a guess
            kva = regexp(rsi.RegulatorType,'([\d.]+)A ([\d.]+)','tokens','once');
            if(isempty(kva))
                kva = regexpi(rsi.RegulatorType,'([\d.]+)KV ([\d.]+)','tokens','once');
                kva = kva(end:-1:1);
            end
            if(~isempty(kva))
                tri.kVAs = prod(str2double(kva))*[1 1];
                tri.kVs = str2double(kva{2})*[1 1];
            else
                % and finally, if all that's failed, try to use the
                % feeder voltage
                try
                    tri.kVs = d.InstFeeders.NominalKvll*[1 1];
                    if(strcmp(tri.Conn,'wye'))
                        tri.kVs = tri.kVs/sqrt(3);
                    end
                catch
                end
            end
        end
        % grab the phase info; 'phase' var is used further down, so
        % don't mess with it.
        if(rsi.TapsAreGangOperated)
            phase = rsi.GangOperatingPhase;
            % let the transformer object convert the integer tap
            % position into a pu tap position; it does this when we use
            % integer data types
% 				tri.tap = int16(rsi.(['TapPositionPhase' phasemap.(phase)(2)]));
        else
            phase = rsi.ConnectedPhases;
            if(sum(phase~=' ' & phase~='N')>1), warning('dssconversion:ThreePhase','Converting a multiphase regulator that is not gang operated may lose info as implemented'); end
            % adopt the most common tap position as the one to use
            tri.tap = int16(mode(nonzeros([rsi.TapPositionPhase1 rsi.TapPositionPhase2 rsi.TapPositionPhase3])));
        end

        tr(tidex) = tri;
        tidex = tidex + 1;
    else
        % load old tranformer
        tri = tr(ismember(rsi.SectionId,{trs.SectionId}));
        % This really shouldn't ever happen in my understanding of the
        % SynerGee model format.  I'm therefore opting not to spend
        % time writing the code correctly to fill out the transformer
        % in this case.
        warning('dssconversion:duplicateRegulator','Found an existing transformer with the same name ''%s'' used by a regulator.',tri.Name);
    end

    % Start setting properties on the regcontrol object
    ri.Name = rsi.UniqueDeviceId;
    ri.transformer = tri.Name;

    % Get bus ID from SectionID ('RG 05201317_05201317A') and
    % NearFromNode properties
    if rsi.NearFromNode == 1
        bus = rs.FromNodeID;
    else
        bus = rs.ToNodeID;
    end
    % select the Phases to monitor on that bus
    % doing it this way actually implies something like Synergee's
    % TapsAreGangOperated, so we'd want to look back at that if they're
    % not.
    [x y phase] = dataclean(phase,'monitoredphase');
    ri.PTphase = phase;
    %bus = [bus '.' phase];
    %ri.bus = bus;

    % select the appropriate winding for the regulator controller
    if rsi.TapsNearFromNode == 1
        ri.winding = 1;
    else
        ri.winding = 2;
    end

    % Set (or guess) the pt ratio
    if(haveRdevs(i))
        ri.ptratio = RegDev.PTRatio;
        % I think these two CT numbers are the same thing, but OpenDSS
        % doesn't document very well what it's doing.  Based on reading
        % the code, I think it uses this as the ratio between line
        % current and control current, which is exactly what synergee
        % does
        ri.CTprim = RegDev.CTRating;

    else
        % PT ratio is supposed to be assigned so that the nominal
        % line-line voltage corresponds to a control Line-Neutral
        % voltage of 120 afaict
        ri.ptratio = tri.kVs(1)/sqrt(3)/120;
        % skip CT rating I guess?  I'm not sure how I'd pick one, and
        % it shouldn't matter anyway as long as we don't use the LDC
        % features in OpenDSS.
    end

    % Use either the specified control values (for gang operated taps)
    % or the most common control values for the forward direction
    if(rsi.TapsAreGangOperated)
        phase = phasemap.(rsi.GangOperatingPhase)(2);
        ri.band = rsi.(['ForwardBWDialPhase' phase]);
        ri.vreg = rsi.(['ForwardVoltageSettingPhase' phase]);
    else
        ri.band = mode(nonzeros([rsi.ForwardBWDialPhase1 rsi.ForwardBWDialPhase2 rsi.ForwardBWDialPhase3]));
        ri.vreg = mode(nonzeros([rsi.ForwardVoltageSettingPhase1 rsi.ForwardVoltageSettingPhase2 rsi.ForwardVoltageSettingPhase3]));
    end
    if(any([rsi.ForwardRDialPhase1 rsi.ForwardRDialPhase2 rsi.ForwardRDialPhase3]))
        warning('dssconversion:unimplemented','Regulator %s uses LDC, which is not implemented',ri.Name);
    end

    % convert the reverse control mode
    % see chart on p167 of synergee tech manual
    switch(rsi.ReverseSensingMode)
        case 'LF'%Lock Forward
        case 'IN'%Neutral Idle
            ri.revNeutral = 'Yes';
        case {'NR','CG'} %no reverse or cogeneration?
        case {'IR','BD','LR'}% idle reverse
            warning('dssconversion:unimplemented','I don''t do anything with the reverse sensing mode %s yet',rsi.ReverseSensingMode)
        otherwise % means I guessed the acronyms above incorrectly
            error('dssconversion:implfailure','Bad Reverse Sensing mode: %s',rsi.ReverseSensingMode);
    end
    % setting for reverse direction operation
    if any([rsi.ReverseVoltageSettingPhase1 rsi.ReverseVoltageSettingPhase2 rsi.ReverseVoltageSettingPhase3])
        ri.reversible = 'Yes';
        % grab the voltage band
        if(rsi.TapsAreGangOperated)
            ri.revband = rsi.(['ReverseBWDialPhase' phase]);
            ri.revvreg = rsi.(['ReverseVoltageSettingPhase' phase]);
        else
            ri.revband = mode(nonzeros([rsi.ReverseBWDialPhase1 rsi.ReverseBWDialPhase2 rsi.ReverseBWDialPhase3]));
            ri.revvreg = mode(nonzeros([rsi.ReverseVoltageSettingPhase1 rsi.ReverseVoltageSettingPhase2 rsi.ReverseVoltageSettingPhase3]));
        end
    end

    if(~(rsi.RegulatorIsOn) || rsi.ManualTapOperation)
        ri.maxtapchange = 0; %disable tap changes
    end
    if(rsi.FirstHouseActive)
        ri.vlimit = rsi.FirstHouseHighVoltage;
        if(rsi.FirstHouseLowVoltage~=0)
            warning('dssconversion:fhlv','First House Low Voltage protection does not exist in opendss');
        end
    end

    if(rsi.NomKvMult~=1)
        warning('dssconversion:unimplemented','Regulator uses a nonunity kv multiplier');
    end

    r(i) = ri;		
end

% %TODO: finally, remove lines that represent these regulator(s)
% [val l_idx] = ismember(dataclean({d.InstRegulators.SectionId},'name'),{l.Name});
% l(nonzeros(l_idx)) = [];



    