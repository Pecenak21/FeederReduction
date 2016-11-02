function tr = Convert2DSS_Transformers(trs)
% tr = Convert2DSS_Transformers(trs)
%
% PURPOSE : Converts transformer from object to OpenDSS
%
%
% INPUT :   trs: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Transformer object
%           see page 111 in OpenDSS manual V7.4.3(March 2012)

% js: probably does not work with SynerGEE, yet. SynerGEE systems converted
% so far did not have Primary Transformers - only regulators, which are
% modeled as transformers

tr(length(trs)) = dsstransformer;

if ~isempty(trs(1).SectionID)
    % make phase mapping easier:
    phasemap.A = '.1';
    phasemap.B = '.2';
    phasemap.C = '.3';
    phasemap.N = '.0';


    % 	[busid, idx_bid] = ismember({trs.SectionID},{d.Section.SectionID});
    %      % js: At this point, I include only Section IDs of lines. Consequently, 
    %      % only the tranformers that have a Section ID that also exists as line 
    %      % are found. Probably need to fix this later.
    %     idx_bid_Transformers=find(idx_bid);
    %     trs=d.Transformers(idx_bid_Transformers);
    %     tr(length(trs)) = dsstransformer;
    %     idx_bid=idx_bid(idx_bid~=0);
    for i = 1:length(trs)
        tr(i).Name = trs(i).UniqueDeviceID;
    %     tr(i).Phases = trs(i).Phases;


        % we only know the properties of the low winding, so we skip directly to set those
        tr(i).wdg = 2;
    %     tr(i).Kv = trs(i).SpecNomKv;
        kV_prim=trs(i).KVLLprim;
        kV_sec=trs(i).KVLLsec;    
        if length(trs(i).Phase)==1
            kV_prim=num2str(str2double(kV_prim)/sqrt(3)); % based on OpenDSS manual, use LL voltage as rating for 2- and 3-phase xfmr and winding voltage otherwise. 
    %         kV_sec=num2str(str2double(kV_sec)/sqrt(3)); %looks like secondary side is on phase to ground voltage already, but not clear
        end
        tr(i).kVs = {kV_prim kV_sec};

    % 	tr(i).kVAs = {'100000' '100000'}; % ignore CYME value and always use 100 MVA base
    %                                       % needed to match OpenDSS fault currents with
    %                                       % customer provided fault currents
        kva=num2str(str2num(trs(i).KVA));
        tr(i).kVAs = {kva kva}; % ignore CYME value and always use 100 MVA base
                                          % needed to match OpenDSS fault currents with
                                          % customer provided fault currents

        tr(i).XHL = trs(i).Z1; 
    %     tr(i).XHL=15;
        secphases = trs(i).Phase;
        secphases(secphases==' ') = [];
    % 	if(~s(i).NeutIsGrounded)
    % 		secphases(secphases=='N') = [];
    % 	end
        tr(i).Phases = length(secphases);
        secphstr = '';
        for phase_idx = 1:length(secphases)
            secphstr = [secphstr phasemap.(secphases(phase_idx))];
        end
        bus1 = [trs(i).FromNodeID secphstr];
        bus2 = [trs(i).ToNodeID secphstr];
        tr(i).Buses = {bus1 bus2};
            if(~trs(i).HighSideNearFromNode)
            tr(i).Buses = tr(i).Buses(end:-1:1);
        end
        % Set connections
        tr(i).Conns = {trs(i).HighSideConnectionCode trs(i).LowSideConnectionCode};
    end
end

end
    
% obj.defaults = struct('Name',{''}, ...
%         'Phases', 3, ...
% 		'Windings',2, ...
% 		'Buses',{{}}, ...
% 		'Conns',{{'wye', 'wye'}}, ...
% 		'kVs',[12.47 12.47], ...
% 		'kVAs',[1000 1000], ...
% 		'Taps',[1, 1], ...
% 		'Rs',[0.2 0.2], ...
% 		'XHL',[], ...
%         'XLT',[], ...
%         'XHT',[], ...
%         'XscArray',[], ...
%         'Thermal',2, ...
%         'n',0.8, ...
%         'm',0.8, ...
%         'flrise',65, ...
%         'hsrise',15, ...
%         'Loadloss',[], ...
%         'Noloadloss',0, ...
%         'imag',0, ...
%         'Ppm_Antifloat',1, ...
%         'NormHKVA',[], ...
%         'EmergHKVA',[], ...
%         'Faultrate',0.007, ...
%         'Like','', ...
%         'sub','No', ... %act as substation
% 		'Basefreq',60, ...
% 		'Wdg',1, ...
% 		'Bus',{{'',''}}, ...
% 		'Conn',{{'wye','wye'}}, ...
%         'kV',{{12.47,12.47}}, ...
%         'kVA',{{1000,1000}}, ...
% 		'tap',{{1,1}}, ...
% 		'MaxTap',{{1.1,1.1}}, ...
% 		'MinTap',{{0.9,0.9}}, ...
% 		'NumTaps',{{32,32}}, ...
% 		'R',{{0.2,0.2}}, ...
% 		'rneut',{{-1,-1}}, ...
% 		'xneut',{{0,0}} );

% The Transfomer is implemented as a multi?terminal (two or more) power delivery element.
% A transfomer consists of two or more Windings, connected in somewhat arbitrary fashion (with
% a default Wye?Delta connection). You can specify the parameters one winding at a time or use
% arrays to set all the winding values at once. Use the "wdg=…" parameter to select a winding for
% editing.
% Transformers have one or more phases. The number of conductors per terminal is always one
% more than the number of phases. For wye? or star?connected windings, the extra conductor is
% the neutral point. For delta?connected windings, the extra terminal is open internally (you
% normally leave this connected to node 0).
% Properties, in order, are:
% Phases= Number of phases. Default is 3.
% Windings= Number of windings. Default is 2.
% For defining the winding values one winding at a time, use the following parameters. Always
% start the winding definition with "wdg = …" when using this method of defining transformer
% parameters. The remainder of the tags are optional as usual if you keep them in order.
% Wdg= Integer representing the winding which will become the active winding for subsequent
% data.
% Bus= Definition for the connection of this winding (each winding is connected to one terminal of
% the transformer and, hence, to one bus).
% Conn= Connection of this winding. One of {wye | ln} for wye connected banks or {delta | ll} for
% delta (line?line) connected banks. Default is wye.
% Kv= Rated voltage of this winding, kV. For transformers designated 2? or 3?phase, enter phaseto?
% phase kV. For all other designations, enter actual winding kV rating. Two?phase
% transfomers are assumed to be employed in a 3?phase system. Default is 12.47 kV.
% Kva= Base kVA rating (OA rating) of this winding.
% Tap = Per unit tap on which this winding is set.
% %R = Percent resistance of this winding on the rated kVA base. (Reactance is between two
% windings and is specified separately ?? see below.)
% rneut = Neutral resistance to ground in ohms for this winding. Ignored if delta winding. For
% open ungrounded neutral, set to a negative number. Default is –1 (capable of being
% ungrounded). The DSS defaults to connecting the neutral to node 0 at a bus, so it
% will still be ground when the system Y is built. To make the neutral floating,
% explicitly connect it to an unused node at the bus, e.g., Bus=Busname.1.2.3.4, when
% node 4 will be the explicit neutral node.
% xneut = Neutral reactance in ohms for this winding. Ignored if delta winding. Assumed to be in
% series with neutral resistance. Default is 0.
% Use the following properties to set the winding values using arrays (setting of wdg= … is
% ignored). The names of these properties are simply the plural form of the property name above.
% Buses = Array of bus definitions for windings [1, 2. …].
% Conns = Array of winding connections for windings [1, 2. …].
% KVs = Array of kV ratings following rules stated above for the kV field for windings [1,2,…].
% KVAs = Array of base kVA ratings for windings [1,2,…].
% Taps = Array of per unit taps for windings [1,2,…].
% %Rs = Array of percent resistances for windings [1, 2. …]
% Use the following propertis to define the reactances of the transformer. For 2? and 3?winding
% transformers, you may use the conventional XHL, XLT, and XHT parameters. You may also put
% the values in an array (xscarray), which is required for higher phase order transformers. There
% are always n*(n?1)/2 different short circuit reactances, where n is the number of windings.
% Always use the kVA base of the first winding for entering impedances. Impedance values are
% entered in percent.
% XHL = Percent reactance high?to?low (winding 1 to winding 2).
% XLT = Percent reactance low?to?tertiary (winding 2 to winding 3).
% XHT = Percent reactance high?to?tertiary (winding 1 to winding 3).
% XscArray = Array of n*(n?1)/2 short circuit reactances in percent on the first winding's kVA base.
% "n" is number of windings. Order is (12, 13, 14, …1n, 23, 24, … 34, …)
% General transformer rating data:
% Thermal = Thermal time constant, hrs. Default is 2.
% n = Thermal exponent, n, from IEEE/ANSI C57. Default is 0.8.
% m = Thermal exponent, m, from IEEE/ANSI C57. Default is 0.8.
% flrise = Full?load temperature rise, degrees centigrade. Default is 65.
% hsrise = Hot?spot temperatire rise, degrees centigrade. Default is 15.
% %Loadloss = Percent Losses at rated load.. Causes the %r values to be set for windings 1 and 2.
% %Noloadloss = Percent No load losses at nominal voltage. Default is 0. Causes a resistive branch
% to be added in parallel with the magnetizing inductance.
% %imag = Percent magnetizing current. Default is 0. An inductance is used to represent the
% magnetizing current. This is embedded within the transformer model as the
% primitive Y matrix is being computed.
% Ppm_Antifloat = Parts per million for anti floating reactance to be connected from each terminal
% to ground. Default is 1. That is, the diagonal of the primitive Y matrix is increased by
% a factor of 1.000001. Prevents singular matrix if delta winding left floating. Set this
% to zero if you don’t need it and the resulting impedance to ground is affecting the
% results. Is inconsequential for most cases.
% NormHKVA = Normal maximum kVA rating for H winding (1). Usually 100 ? 110% of maximum
% nameplate rating.
% EmergHKVA = Emergency maximum kVA rating for H winding (1). Usually 140 ? 150% of
% maximum nameplate rating. This is the amount of loading that will cause 1% loss of
% life in one day.
% Faultrate = Failure rate for transformer. Defaults to 0.007 per year. All are considered
% permanent.
% Basefreq= Base frequency, Hz. Default is 60.0
% Like= Name of another Transformer object on which to base this one.
% Sub = Yes/No. Designates whether this transformer is to be treated as a substation. Default is
% No.
