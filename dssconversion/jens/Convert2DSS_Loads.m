function [ld ge] = Convert2DSS_Loads(ldi)
% l = Convert2DSS_Loads(s)
%
% PURPOSE : Converts loads from object to OpenDSS
%
%
% INPUT :   s: struct that has all entries (e.g., o.R1)
%
% OUTPUT :  OpenDSS Load object
%           see page 122 in OpenDSS manual V7.4.3(March 2012)



if isempty(ldi(1).SectionID)==0
    
    
    
    % process load
    if(length(ldi)==1 && iscell(ldi.SectionID)), ldi = structconv(ldi); end
    ld = dssload(); ld(length(ldi)).Name = '';
    idx = 1;
    phases = [[str2nump(ldi.Phase1Kw)]' [str2nump(ldi.Phase2Kw)]' [str2nump(ldi.Phase3Kw)]'];
%     [tmpvar sectionind] = ismember({ldi.SectionID},{d.Section.SectionID}); % make sure that the Section ID of each load bus is included in the Section structure

    for i = 1:length( ldi )

%         busid = sectionind(i);
        % we're currently opting to model all loads as spot loads at the 'to'
        % node unless expicitly specified as 'F' (from node), which doesn't
        % actually happen in any of our models.
        if(ldi(i).LocationToModelSpotLoads == 'F')
            busid = ldi(i).FromNodeID;
        else
            busid = ldi(i).ToNodeID;
        end
        % From the OpenDSS help file in V 7.6.1.0
        % ZIPV
        %  Array of 7 coefficients:
        %  First 3 are ZIP weighting factors for real power (should sum to 1)
        %  Next 3 are ZIP weighting factors for reactive power (should sum to 1)
        %  Last 1 is cut-off voltage in p.u. of base kV; load is 0 below this cut-off
        %  No defaults; all coefficients must be specified if using model=8.
        zipv_P = [str2nump(ldi(i).ConstantImpedance) str2nump(ldi(i).ConstantCurrent) str2nump(ldi(i).ConstantPower)]./100;
        zipv_Q = zipv_P;
        zipv = [zipv_P zipv_Q .5]; % assume real and reactive parts same, set cutoff to 0.5pu (may be way too low?

        % for now we ignore the Kwh field for each phase, since it's not clear
        % what that would give to the model, and they're all zero anyway.
        if(ldi(i).IsSpotLoad || true)
            ld(idx).bus1 = busid;
            % get a list of the power on each of the phases
            kw = phases(i,phases(i,:)>0);
            % Right now I'm assuming only LN connected loads.  Some loads might
            % be LL, in which case that would be weird, but our current model
            % specifies all the sections as having a neutral connection, so we
            % don't need to worry about that here.
            if(length(kw)<3 || length(unique(kw))>1) %if we have unequal power on the phases or less than three phases, we model power use with one load per phase
                for ph=find(phases(i,:)>0)
                    ld(idx).Name = ['l' busid '_' ldi(i).(sprintf('Phase%i',ph))]; % not sure why this results in XXX for Phase 1, YYY for Phase 2, etc., I think set this way by means of load class
                    ld(idx).bus1 = sprintf('%s.%i',busid,ph);
                    ld(idx).Phases = 1;
                    ld(idx).NumCust = ldi(i).(sprintf('Phase%iCustomers',ph));
                    ld(idx).Kw = phases(i,ph);
                    ld(idx).Kvar = ldi(i).(sprintf('Phase%iKvar',ph));
                    ld(idx).Model = 8; %ZIP load
                    ld(idx).ZIPV = zipv;
                    idx = idx+1;
                end
                if isempty(find(phases(i,:)))
                    ld(idx).Phases = 0;
                    idx = idx+1;
                end
            else % if we have three phases with equal power, we can use a single 3-phase load
                ld(idx).Name = ['l' busid];
                ld(idx).Phases = 3;
                ld(idx).NumCust = str2nump(ldi(i).Phase1Customers) + str2nump(ldi(i).Phase2Customers) + str2nump(ldi(i).Phase3Customers);
                ld(idx).Kw = sum(kw);
                ld(idx).Kvar = str2nump(ldi(i).Phase1Kvar)+str2nump(ldi(i).Phase2Kvar)+str2nump(ldi(i).Phase3Kvar);
                ld(idx).Model = 8; %ZIP load
                ld(idx).ZIPV = zipv;
                idx = idx + 1;
            end
        else
            % The treatment above is actually correct only for spot loads.
            % Non spot loads are indicated to be modeled at 'center' which is
            % not what we're doing right now.
%             ld(idx).bus1 = busid;
        end
        
        ge=[];
        % process generation (if load has generation), applies only to SynerGEE
        if isfield(ldi(i),'GenTotalKw') && ~isempty(ldi(i).GenTotalKw)

            

            % add generator if customer is generating electricity
            if str2num(ldi(i).GenTotalKw)>0
                                
                % intialize generator classes
                genClasses = {'unknown'};
                %new class
                nc = unique(lower({ldi(i).GenCustClass}));
                genClasses = [genClasses setdiff(nc,genClasses)];
                
                % initialize generators
                ge = dssgenerator;
                
                % add data
                GenTotalKw_num=str2num(ldi(i).GenTotalKw);
                GenTotalKvar_num=str2num(ldi(i).GenTotalKvar);
                GenPct_num=str2num(ldi(i).GenPct);
                ge.Name = ldi(i).UniqueDeviceId; 
                ge.Kw = GenTotalKw_num * GenPct_num/100;
                ge.Kvar = GenTotalKvar_num;
                ge.bus1 = busid;

                % PV generator
                if ~isempty(strfind('pv',lower(ldi(i).GenCustClass)))
                    ge.Model = 3;
                end

                % set class. Class 0 is 'unknown' (our convention).
                [val c] = ismember(lower(ldi(i).GenCustClass),genClasses);
                ge.Class = c - 1;

                if strcmp(ldi(i).GenStatus,'P') %powered
                    ge.Enabled = 'yes';
                elseif strcmp(ldi(i).GenStatus,'O') %offline
                    ge.Enabled = 'no';
                else
                    warning('dssconversion:generatorSetup','check generator status again. You might need to change the codes according to your data on how they specify generator status. Set to YES by default.');
                    ge.Enabled = 'yes';
                end
            end
        end
        
    end
    
end

% 	obj.defaults = struct( ...
%         'Name',{''},...
% 		'Phases',3,...
% 		'bus1',{''},...
%         'Kv',12.47,...
%         'Kw',10,...
%         'Pf',0.88,...
%         'Model',1,...
%         'Yearly',{[]},...
%         'Daily',{[]},...
%         'Duty',{[]},...
%         'Growth',{[]},...
%         'Conn','wye',...
%         'Kvar',5.4,...
%         'Rneut',-1,...
%         'Xneut',0,...
%         'Status','variable',...
%         'Class',1,...
%         'Vminpu',0.95,...
%         'Vmaxpu',1.05,...
%         'VminNorm',0,...
%         'VminEmerg',0,...
%         'XfkVA',0.0,...
%         'AllocationFactor',.5,...
%         'kVA',hypot(12.47,5.4),...
%         'mean',50,... % percents
%         'stddev',10,...
%         'CVRwatts',1,...
%         'CVRvars',2,...
%         'kWh',0,...
%         'kWhDays',30,...
%         'CFactor',4.0,...
%         'CVRCurve',{[]},...
%         'NumCust',1,...
% 		'ZIPV',[],... % 7 values: ZIP fractions for real power, ZIP fractions for reactive power, cutoff voltage in terms of pu of basekv (load is 0 below cutoff)
%         'spectrum','defaultload',...
%         'BaseFreq',60,...
%         'Like'      ,{''}...
%     );

% A Load is a complicated Power Conversion element that is at the heart of many analyses. It is
% basically defined by its nominal kW and PF or its kW and kvar. Then it may be modified by a
% number of multipliers, including the global circuit load multiplier, yearly load shape, daily load
% shape, and a dutycycle load shape.
% The default is for the load to be a current injection source. Thus, its primitive Y matrix contains
% only the impedance that might exist from the neutral of a wye?connected load to ground.
% However, if the load model is switched to Admittance from PowerFlow (see Set LoadModel
% command), the load is converted to an admittance and included in the system Y matrix. This
% would be the model used for fault studies where convergence might not be achieved because of
% low voltages.
% Loads are assumed balanced for the number of phases specified. If you would like unbalanced
% loads, enter separate single?phase loads.
% There are three legal ways to specify the base load:
% 1. kW, PF
% 2. kw, kvar
% 3. kVA, PF
% If you sent these properties in the order shown, the definition should work. If you deviate from
% these procedures, the result may or may not be what you want. (To determine if it has
% accomplished the desired effect, execute the Dump command for the desired load(s) and
% observe the settings.)
% The properties, in order, are:
% bus1= Name of bus to which the load is connected. Include node definitions if the
% terminal conductors are connected abnormally. 3?phase Wye?connected loads
% have 4 conductors; Delta?connected have 3. Wye?connected loads, in general, have
% one more conductor than phases. 1?phase Delta has 2 conductors; 2?phase has 3.
% The remaining Delta, or line?line, connections have the same number of conductors
% as phases.
% Phases= No. of phases this load.
% Kv= Base voltage for load. For 2? or 3?phase loads, specified in phase?to?phase kV. For
% all other loads, the actual kV across the load branch. If wye (star) connected, then
% specify phase?to?neutral (L?N). If delta or phase?to?phase connected, specify the
% phase?to?phase (L?L) kV.
% Kw= nominal kW for load. Total of all phases. See kVA.
% Pf= nominal Power Factor for load. Negative PF is leading. Specify either PF or kvar
% (see below). If both are specified, the last one specified takes precedence.
% Model= Integer defining how the load will vary with voltage. Presently defined models are:
% 1: Normal load?flow type load: constant P and Q.
% 2: Constant impedance load
% 3: Constant P, Quadratic Q (somewhat like a motor)
% 4: Linear P, Quadratic Q (Mixed resistive, motor)
% 5: Rectifier load (Constant P, constant current)
% 6: Constant P; Q is fixed at nominal value
% 7: Constant P; Q is fixed impedance at nominal value
% 8: ZIP model
% "Constant" values may be modified by loadshape multipliers. "Fixed" values are
% always the same ?? at nominal, or base, value.
% Yearly= Name of Yearly load shape.
% Daily= Name of Daily load shape.
% Duty= name of Duty cycleload shape. Defaults to Daily load shape if not defined.
% Growth= Name of Growth Shape Growth factor defaults to the circuit's default growth rate
% if not defined. (see Set %Growth command)
% Conn= {wye | y | LN} for Wye (Line?Neutral) connection; {delta | LL} for Delta (Line?Line)
% connection. Default = wye.
% Kvar= Base kvar. If this is specified, supercedes PF. (see PF)
% Rneut= Neutral resistance, ohms. If entered as negative, non?zero number, neutral is
% assumed open, or ungrounded. Ignored for delta or line?line connected loads.
% Xneut= Neutral reactance, ohms. Ignored for delta or line?line connected loads. Assumed
% to be in series with Rneut value.
% Status= {fixed| variable}. Default is variable. If fixed, then the load is not modified by
% multipliers; it is fixed at its defined base value.
% Class= Integer number segregating the load according to a particular class.
% Vminpu = Default = 0.95. Minimum per unit voltage for which the MODEL is assumed to
% apply. Below this value, the load model reverts to a constant impedance model.
% Vmaxpu = Default = 1.05. Maximum per unit voltage for which the MODEL is assumed to
% apply. Above this value, the load model reverts to a constant impedance model.
% VminNorm = Minimum per unit voltage for load EEN evaluations, Normal limit. Default = 0,
% which defaults to system "vminnorm" property (see Set Command under
% Executive). If this property is specified, it ALWAYS overrides the system
% specification. This allows you to have different criteria for different loads. Set to
% zero to revert to the default system value.
% VminEmerg = Minimum per unit voltage for load UE evaluations, Emergency limit. Default = 0,
% which defaults to system "vminemerg" property (see Set Command under
% Executive). If this property is specified, it ALWAYS overrides the system
% specification. This allows you to have different criteria for different loads. Set to
% zero to revert to the default system value.
% XfkVA = Default = 0.0. Rated kVA of service transformer for allocating loads based on
% connected kVA at a bus. Side effect: kW, PF, and kvar are modified. See
% PeakCurrent property of EnergyMeter. See also AllocateLoads Command. See kVA
% property below.
% AllocationFactor = Default = 0.5. Allocation factor for allocating loads based on connected kVA
% at a bus. Side effect: kW, PF, and kvar are modified by multiplying this factor times
% the XFKVA (if > 0). See also AllocateLoads Command.
% kVA = Definition of the Base load in kVA, total all phases. This is intended to be used in
% combination with the power factor (PF) to determine the actual load. Legal ways to
% define base load (kW and kvar):
% kW, PF
% kW, kvar
% kVA, PF
% XFKVA * Allocationfactor, PF
% kWh/(kWhdays*24) * Cfactor, PF
% %mean = Percent mean value for load to use for monte carlo studies if no loadshape is
% assigned to this load. Default is 50.
% %stddev = Percent Std deviation value for load to use for monte carlo studies if no loadshape
% is assigned to this load. Default is 10.
% CVRwatts = Percent reduction in active power (watts) per 1% reduction in voltage from 100%
% rated. Default=1. Typical values range from 0.4 to 0.8. Applies to Model=4 only.
% Intended to represent conservation voltage reduction or voltage optimization
% measures.
% CVRvars = Percent reduction in reactive power (vars) per 1% reduction in voltage from 100%
% rated. Default=2. Typical values range from 2 to 3. Applies to Model=4 only.
% Intended to represent conservation voltage reduction or voltage optimization
% measures.
% kWh = kWh billed for this period. Default is 0. See help on kVA and Cfactor and kWhDays.
% kWhDays = Length of kWh billing period in days (24 hr days). Default is 30. Average demand is
% computed using this value.
% CFactor = Factor relating average kW to peak kW. Default is 4.0. See kWh and kWhdays. See
% kVA.
% CVRCurve = Default is NONE. Curve describing both watt and var factors as a function of time.
% Refers to a LoadShape object with both Mult and Qmult defined. Define a
% Loadshape to agree with yearly or daily curve according to the type of analysis
% being done. If NONE, the CVRwatts and CVRvars factors are used and assumed
% constant for the simulation period.
% NumCust = Number of customers, this load. Default is 1.
% spectrum = Name of harmonic current spectrum for this load. Default is "defaultload", which is
% defined when the DSS starts.
% Basefreq = Base frequency for which this load is defined. Default is 60.0.
% Like = Name of another Load object on which to base this one.

