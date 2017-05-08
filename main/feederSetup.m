function [ c, opt, d, glc] = feederSetup( feederName, configId, overwrite)
% feederSetup: This function generates any specified feeder configuration for simulation purpose. 
% This function currently supports 5 SDG&E feeders with a few configuration each.
%
% input:
%           feederName  : Name/ Id of the feeder. E.g. 'fallbrook' or '520'
%           configId    : id of the configuration to use. E.g. 'original' or 'validated'
%           overwrite   : (default: 0) overwrite saved circuit data if needed
%
% output:
%           c           : core circuit data
%           opt         : customized feeder option with fields 
%                           'fcProfileId' if specified forecast profile(s) should be used
%                           'loadProfId' if specified load profile(s) should be used
%                           'excludeScalingUpIds' if some pv systems should be excluded from the scaling up process for higher PV penetration levels.
%           d           : original feeder data
%           glc         : original linecode data
%
% examples of use:
%           [ c, opt, d, glc] = feederSetup( 'fallbrook', 'validated');
%           [ c, opt, d, glc] = feederSetup( 'alpine', 'wpv_existing_balanced3phase', 1);
%
% supported feeders and configs (TO BE UPDATED, look into the function for more details):
%           Fallbrook feeder: {'fallbrook','520','avocado','a'}
%               configs:    'original' : original converted OpenDSS-like data format from Synergy data.
%                           'validated': validated circuit using short circuit and power flow data given from utility
%                           'wpv_existing': validated circuit with existing pv systems added
%                           'wpv_virtual': validated circuit with both existing and virtual pv systems created for higher pv penetration levels.
%                           {'s1b','balanced_virtualPV_vr_begin'}: balanced-3phase circuit with all existing and virtual pv systems. Voltage regulators' set points are at the beginning of the branches.

% initialization
def_addpath;
global conf; if isempty(conf), conf = getConf; end
global kk;
outDir = conf.outputDir; if ~exist(outDir,'dir'); mkdir(outDir); end
fn = [feederName '_' configId '.mat'];
fp = [outDir '/' fn];

% initialization of customized feeder option parameters
opt.fcProfileId = [];
opt.loadProfId = [];

% check if saved circuit file exists and load to use if available
if exist(fp,'file') && (~exist('overwrite','var') || ~overwrite)
    c = load(fp); d = c.d; glc = c.glc;
    if isfield(c,'opt'), opt = c.opt; end
    c = c.c;
    return;
end
close all;
validatePQ = 0; % run power and voltage validation after constructing the circuit by defaults
excludeScalingUpIds = []; % exclude some generator ids during the scaling up process if needed
switch lower(feederName)
	  case {'ieee8500','8500'}
        %% IEEE8500 feeder
        excludeScalingUpIds = [];
        switch lower(configId)
            case 'original'
				c=dssparse('c:\users\zactus\gridIntegration\masterCircuit_PV.dss');
				glc=[]; ValidatePQ = 0; opt=[]; d=[];
				o = actxserver('OpendssEngine.dss');
				dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
				p='C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\8500-Node\Master.dss';
				dssText.Command = ['Compile "' p '"']; 
				dssCircuit = o.ActiveCircuit;
				c.buslist_orig=c.buslist;
				c.buslist.id=dssCircuit.AllBUSNames;
				c.buslist.coord=zeros(length(c.buslist.id),2);
				
			case 'reduced'
				[c,d,glc] = feederSetup( feederName, 'original');
								
				critical_nodes=c.buslist.id(round((length(c.buslist.id)-2)*rand(kk,1)+1));
				critical_nodes=unique(critical_nodes);
				while length(critical_nodes)<kk
					cnode_tmp=c.buslist.id(round((length(c.buslist.id)-2)*rand(20,1)+1));
					critical_nodes=[critical_nodes;cnode_tmp];
					critical_nodes=unique(critical_nodes);
				end	% 			   [c,~] = FeederReduction(critical_nodes,c);
				
				for ii=1:length(critical_nodes)
					critical_nodes(ii)=regexprep(critical_nodes(ii),'-','_');
				end
				[c] = FeederReduction_SE(critical_nodes,c,[]);
		end
		case {'ieee13','13'}
		excludeScalingUpIds = [];
        switch lower(configId)
        %% IEEE13 feeder
            case 'original'
				c = dssparse('C:\Users\Zactus\Documents\OpenDSS\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt_zack.dss')
				glc=[]; ValidatePQ = 0; opt=[]; d=[];
				o = actxserver('OpendssEngine.dss');
				dssText = o.Text; dssText.Command = 'Clear'; cDir = pwd;
								
			case 'reduced'
				[c,d,glc] = feederSetup( feederName, 'original');
								
% 				critical_nodes=c.buslist.id(round((length(c.buslist.id)-2)*rand(kk,1)+1));
% 				critical_nodes=unique(critical_nodes);
% 				while length(critical_nodes)<kk
% 					cnode_tmp=c.buslist.id(round((length(c.buslist.id)-2)*rand(20,1)+1));
% 					critical_nodes=[critical_nodes;cnode_tmp];
% 					critical_nodes=unique(critical_nodes);
% 				end	% 			   [c,~] = FeederReduction(critical_nodes,c);
% 				
% 				for ii=1:length(critical_nodes)
% 					critical_nodes(ii)=regexprep(critical_nodes(ii),'-','_');
% 				end
				critical_nodes=c.buslist.id([13,14]);

% 				critical_nodes={'611','652','675','692'};
				load('c:\users\zactus\FeederReductionRepo\StateEstimation\Z_New.mat')
				[c] = FeederReduction_SE(critical_nodes,c,Z);
		end
	case {'clairemont','276'}
        %% clairemont feeder {'clairemont','276'}
        excludeScalingUpIds = [];
        switch lower(configId)
            
            case 'original'
                
                %% cliaremont original
                d = load(which('clairemont.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% clairemont validated
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'clairemont';
                c.transformer(end).Buses = {'sourcebus' '276'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [15000 15000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'clairemont';
                c.regcontrol(end).transformer = 'clairemont';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;
                
                
                % modify the setting of capacitors
                c.capacitor(1).Numsteps = 20;
                c.capacitor(2).Numsteps = 20;
                c.capacitor(3).Numsteps = 20;
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
                
                % Close switches
%                 c.switch(5).Action = 'Close';
%                 c.switch(6).Action = 'Close';
%                 c.switch(11).Action = 'Close';
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % c = addEnergyMeter(c);
                % use constant current load model
                for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                end
                
                
                % assign validated load and fuse to c
                c_validated = dssparse(which('clairemont_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
                c.fuse = c_validated.fuse;
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/clairemont/clairemont_pv_reorder.mat');
                pv = pv.PV_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
             case 'wpv_existing_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
                
           otherwise
                error('Not supported config!');
        end
        
    case {'encinitas','289'}
        %% encinitas feeder {'encinitas','289'}
        excludeScalingUpIds = [];
        switch lower(configId)
            
            case 'original'
                
                %% cliaremont original
                d = load(which('encinitas.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% encinitas validated
                
                % add source bus
                
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                % remove transformer at DP_289_23361
                
                c.transformer(end+1) = c.transformer(end);
                c.transformer(end).Name = 'encinitas';
                c.transformer(end).Buses = {'sourcebus' '289'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [15000 15000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'encinitas';
                c.regcontrol(end).transformer = 'encinitas';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;
                
                
               % modify the setting of capacitors
                c.capacitor(1).Numsteps = 20;
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end


               
                % use constant current load model
                for i = 1:length(c.load);
                
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                end
                
                
                % assign validated load and fuse to c
                c_validated = dssparse(which('encinitas_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
                c.fuse = c_validated.fuse;
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/Encinitas/Encinitas_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
           otherwise
                error('Not supported config!');
        end 
        
    case {'san_marcos_299','299'}
        
        %% san marcos 299 feeder {'san marcos 299','299'}
        excludeScalingUpIds = [];
        switch lower(configId)
            
            case 'original'
                
                %% san marcos 299 original
                d = load(which('san_marcos_299.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% san marcos 299 validated
                                
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'san_marcos_299';
                c.transformer(end).Buses = {'sourcebus' '299'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [15000 15000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'san_marcos_299';
                c.regcontrol(end).transformer = 'san_marcos_299';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;

                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
                % use constant current load model
               for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
               end
               
                % assign validated load and fuse to c
                c_validated = dssparse(which('san_marcos_299_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
                c.fuse = c_validated.fuse;
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/San_Marcos/San_Marcos_299_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
            case 'wpv_existing_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
           otherwise
                error('Not supported config!');
        end 
        
    case {'artesian','1104'}
        
        %% artesian 1104 feeder {'artesian','1104'}
        excludeScalingUpIds = [];
        
        switch lower(configId)
            
            case 'original'
                
                %% artesian original
                d = load(which('artesian.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% artesian validated
                                
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'artesian';
                c.transformer(end).Buses = {'sourcebus' '1104'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [15000 15000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'artesian';
                c.regcontrol(end).transformer = 'artesian';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;
           
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
               
                % use constant current load model
                for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
                end
                
                % assign validated load and fuse to c
                c_validated = dssparse(which('artesian_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
                c.fuse = c_validated.fuse;
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/artesian/artesian_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
                % Since we only mapped 176 pvs in the PV deployment files,
                % but there are actually 180 PV, here we specify PV
                % generation profile Id for the left 4 PVs 
                
                % PV 177 using profile of PV 157
                % PV 178 using profile of PV 157
                % PV 179 using profile of PV 144
                % PV 180 using profile of PV 147
                
                tmp = 1:176;
                opt.fcProfileId = [tmp 157 157 144 147];
                clear tmp;
                
				case 's1b_si'
					[ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                
				% Define Inverter control curve
				c.xycurve(1)=dssxycurve;
				c.xycurve(1).name='vv_curve';
				c.xycurve(1).npts=6;
				c.xycurve(1).Yarray='(1, 1, 0, 0, -1, -1)';
				c.xycurve(1).Xarray='(.5, 0.95, .98, 1.02, 1.05, 1.5)';
				
				%Add inverter class to system
				c.InvControl=dssInvControl; 
				c.InvControl.VoltageChangeTolerance=0.1;
				c.InvControl.VarChangeTolerance=0.1;
				c.InvControl.avgwindowlen='30s';
				c.InvControl.DeltaQ_factor=.1;
				
					
				%Must define both efficiency and PT curve of pvsystems
				%define xycurves 
				%efficiency
				c.xycurve(2)=dssxycurve;
				c.xycurve(2).name='Myeff';
				c.xycurve(2).npts=4;
				c.xycurve(2).Xarray=[.1 .2 .4 1.0];
				c.xycurve(2).yarray=[1 1 1 1];
				
				%P-T
				c.xycurve(3)=dssxycurve;
				c.xycurve(3).name='MyPvsT';
				c.xycurve(3).npts=4;
				c.xycurve(3).Xarray=[0 25 75 100];
				c.xycurve(3).yarray=[1.2 1.0 0.8 0.6];
				
				%Define Loadshape for PV %%% This is for validation
				%purposes only %%%
				c.loadshape=dssloadshape;
				c.loadshape.Name='ValidationLoadShape';
				c.loadshape.npts=1;
				c.loadshape.interval=1;
				c.loadshape.mult=1;
				
							
				%loop over pv and add vaiables
				for i = 1:length(c.pvsystem)
				c.pvsystem(i).EFFCURVE='Myeff';
				c.pvsystem(i).PTCurve='MyPvsT';
				c.pvsystem(i).daily='ValidationLoadShape';
				end
				
				case 's1b_si_50'
					[ c, opt, d, glc] = feederSetup( feederName, 's1b_si');
                load('c:\users\zactus\gridIntegration\SmartInverterSims\ArtesianPvBuses.mat')
				
				for ii=1:length(pvsystPast)
				c.InvControl(ii)=dssInvControl; 
				c.InvControl(ii).Name=['Inv' num2str(ii)];
				c.InvControl(ii).VoltageChangeTolerance=0.1;
				c.InvControl(ii).VarChangeTolerance=0.1;
				c.InvControl(ii).avgwindowlen='30s';
				c.InvControl(ii).DeltaQ_factor=.1;
				c.InvControl(ii).PVSystemList=c.pvsystem(pvsystPast(ii)).name;
				end
				
				otherwise
                error('Not supported config!');
        end   
     
    case {'batiquitos_757','757'}
        
        %% batiquitos 757 feeder {'batiquitos_757','757'}
        excludeScalingUpIds = [];
        
        switch lower(configId)
            
            case 'original'
                
                %% batiquitos 757 original
                d = load(which('batiquitos_757.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% batiquitos 757 validated
                                
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'batiquitos_757';
                c.transformer(end).Buses = {'sourcebus' '757'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [20000 20000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'batiquitos_757';
                c.regcontrol(end).transformer = 'batiquitos_757';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;
           
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
               
                % use constant current load model
                for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
                end
                
                % assign validated load and fuse to c
                c_validated = dssparse(which('batiquitos_757_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
               
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/batiquitos/batiquitos_757_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
           otherwise
                error('Not supported config!');
        end 
        
     case {'batiquitos_1118','1118'}
        
        %% batiquitos 1118 feeder {'batiquitos_1118','1118'}
        
        % pv 122 is about 400kva, which makes it a large PV
        excludeScalingUpIds = 122;
        
        switch lower(configId)
            
            case 'original'
                
                %% batiquitos 1118 original
                d = load(which('batiquitos_1118.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% batiquitos 757 validated
                                
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'batiquitos_1118';
                c.transformer(end).Buses = {'sourcebus' '1118'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [20000 20000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'batiquitos_1118';
                c.regcontrol(end).transformer = 'batiquitos_1118';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;
           
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
               
                % use constant current load model
                for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
                end
                
                % assign validated load and fuse to c
                c_validated = dssparse(which('batiquitos_1118_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
               
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/batiquitos/batiquitos_1118_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
           otherwise
                error('Not supported config!');
        end 
           
     case {'rose_canyon','65'}
        
        %% rose canyon feeder {'rose_canyon','65'}
        excludeScalingUpIds = [];
        
        switch lower(configId)
            
            case 'original'
                
                %% rose canyon original
                d = load(which('rose_canyon.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'rose_canyon';
                c.transformer(end).Buses = {'sourcebus' '65'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [20000 20000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'rose_canyon';
                c.regcontrol(end).transformer = 'rose_canyon';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;
           
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
               
                % use constant current load model
                for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
                end
                
                % assign validated load and fuse to c
                c_validated = dssparse(which('rose_canyon_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
               
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/Rose_Canyon/rose_canyon_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
           otherwise
                error('Not supported config!');
        end   
        
     case {'sampson','80'}
        
        %% sampson feeder {'sampson','80'}
        excludeScalingUpIds = [];
        switch lower(configId)
            
            case 'original'
                
                %% sampson 80 original
                d = load(which('sampson.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% sampson validated
                                
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'sampson';
                c.transformer(end).Buses = {'sourcebus' '80'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [15000 15000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'sampson';
                c.regcontrol(end).transformer = 'san_marcos_299';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;

                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
                % use constant current load model
               for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
               end
               
                % assign validated load and fuse to c
                c_validated = dssparse(which('sampson_imbalanced_loadmove.dss'));
                c.load  = c_validated.load;
                
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/Sampson/sampson_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
            case 'wpv_existing_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
           otherwise
                error('Not supported config!');
        end 
        
     case {'streamview','168'}
        
        %% streamview feeder {'streamview','168'}
        excludeScalingUpIds = [];
        switch lower(configId)
            
            case 'original'
                
                %% streamview 168 original
                d = load(which('streamview.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% sampson validated
                                
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'streamview';
                c.transformer(end).Buses = {'sourcebus' '168'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [15000 15000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'streamview';
                c.regcontrol(end).transformer = 'streamview';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;

                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
                % use constant current load model
               for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
               end
               
                % assign validated load and fuse to c
                c_validated = dssparse('\validation\dss\streamview_imbalanced_loadmove.dss');
                c.load  = c_validated.load;
                
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/Streamview/streamview_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
            case 'wpv_existing_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
           otherwise
                error('Not supported config!');
        end 
        
    case {'olivenhain','1250'}
        
        %% streamview feeder {'streamview','168'}
        excludeScalingUpIds = [];
        switch lower(configId)
            
            case 'original'
                
                %% olivenhain original
                d = load(which('olivenhain.mat')); d = d.obj;
                glc = excel2obj( 'linecode_sdge.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc); ValidatePQ = 0;
                
            case 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                %% olivenhain validated
                                
                % add source bus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and regulator
                c.transformer = dsstransformer;
                c.transformer(end).Name = 'olivenhain';
                c.transformer(end).Buses = {'sourcebus' '1250'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [20000 20000];
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol = dssregcontrol;
                c.regcontrol(end).Name = 'olivenhain';
                c.regcontrol(end).transformer = 'olivenhain';
                
                % Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;

                % modify the setting of capacitors
                for i = 1:length(c.capacitor);
                    c.capacitor(i).Numsteps = 20;
                end
                
                
                for i = 1:length(c.capcontrol);
                    c.capcontrol(i).OFFsetting = 123;
                    c.capcontrol(i).ONsetting  =119;
                end
                % use constant current load model
               for i = 1:length(c.load);
                % c.load(i).model = 5;
                c.load(i).model = 5;
                
                if c.load(i).Phases == 1;
                    c.load(i).kv = 6.9282;
                else c.load(i).kv = 12; end
                
                
               end
               
                % assign validated load and fuse to c
                c_validated = dssparse('\validation\dss\olivenhain_imbalanced_loadmove.dss');
                c.load  = c_validated.load;
                
                
                % remove fuse
                % Typically, load flow is done without fuses
                if isfield(c,'fuse'), c = rmfield(c,'fuse'); end
                
            case 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                
                % add existing PV to the feeder
                pv = load('/dss_new_feeders/mat_files/Olivenhain/olivenhain_pv_reorder.mat');
                pv = pv.pv_re;
                c.pvsystem = pv;

                % reset pvsystem's settings
                c = convertPVModel(c);
                
            case 'wpv_existing_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
           otherwise
                error('Not supported config!');
        end 
    case {'fallbrook','520','avocado','a'}
        %% Fallbrook feeder {'fallbrook','520','avocado','a'}
        excludeScalingUpIds = []; %[44,45];
        switch lower(configId)
            case 'original'
                %% Fallbrook original
                d = excel2obj([conf.gridIntDir '/dssconversion/custdata/520ForENERNEX.xlsx']);
                glc = excel2obj( 'linecode.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc ); validatePQ = 0;
                w=1;
            case 'validated'
                %% Fallbrook Validated
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                
                % add sourcebus
                c.circuit.bus1 = 'sourcebus';
                c.buslist.id = ['sourcebus'; c.buslist.id];
                c.buslist.coord = [c.buslist.coord(1,:); c.buslist.coord];
                c.buslist.coord(1,1) = c.buslist.coord(1,1)-2;
                c.circuit.basekv = 12;
                c.basevoltages = 12;
                
                % add substation transformer and volt regulator
                c.transformer(end+1) = c.transformer(end);
                c.transformer(end).Name = 'Avocado';
                c.transformer(end).Buses = {'sourcebus' '0520'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [28000 28000];% changed from [10000 10000]
                c.transformer(end) = resetTransformer(c.transformer(end));
                
                % add regcontrol
                c.regcontrol(end+1) = c.regcontrol(end);
                c.regcontrol(end).Name = 'Avocado';
                c.regcontrol(end).transformer = 'Avocado';
                
                % %% Changing Transformer Settings
                for i=1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i));
                end
                
                % Changing Regcontrol settings
                for i=1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i));
                    c.regcontrol(i).vreg = 120;
                end
                c.regcontrol(end).vreg = 122;%higher voltage at substation
                % add delays when needed
%                 c.regcontrol(1).delay = 15;
%                 c.regcontrol(2).delay = 45;
%                 c.regcontrol(3).delay = 75;
                
                % modify loads
                loads_1Phase = c.load([c.load.Phases]==1);
                [~, b1] = ismember({loads_1Phase.Name}, {c.load.Name});
                for i=1:length(b1)
                    c.load(b1(i)).Kv = 6.9282;
                    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*2.5;
                end
                loads_3Phase = c.load([c.load.Phases]==3);
                [~, b3] = ismember({loads_3Phase.Name}, {c.load.Name});
                for j=1:length(b3)
                    c.load(b3(j)).Kv = 12;
                    if c.load(b3(j)).Kvar < 189 % the 4 biggest loads stay unchanged
                        c.load(b3(j)).Kvar = c.load(b3(j)).Kvar*2.5;
                    end
                end
                c.load(1731).Kw = 525;
                c.load(1732).Kw = 1313;
                
                % Optimizing capacitor banks to match reactive power
                % Modify existing caps
                c.capacitor(1).Kvar = 1420;% 2900/sqrt(3);
                c.capacitor(1).Numsteps = 20;
                c.capacitor(1).Name = ['cap_' c.capacitor(1).Name];
                c.capcontrol(1).Capacitor = c.capacitor(1).Name;
                c.capacitor(2).Kvar = 1400;%650/sqrt(3);
                c.capacitor(2).Numsteps = 20;
                c.capacitor(2).Name = ['cap_' c.capacitor(2).Name];
                c.capcontrol(2).Capacitor = c.capacitor(2).Name;
                
                % Adding cap banks at appropriate locations
                c.capacitor(3) = c.capacitor(2);
                c.capacitor(3).Bus1 = '05201643A.1.2.3';
                c.capacitor(3).Name = ['cap_' cleanBus(c.capacitor(3).Bus1) 'var'];
                c.capacitor(3).kvar = 1300;
                c.capacitor(3).Numsteps = 25;
                
                c.capacitor(4) = c.capacitor(2);
                c.capacitor(4).Bus1 = '05201947';
                c.capacitor(4).Name = ['cap_' c.capacitor(4).Bus1 'var'];
                c.capacitor(4).kvar = 120;
                
                c.capacitor(5) = c.capacitor(2);
                c.capacitor(5).Bus1 = '05201349';
                c.capacitor(5).Name = ['cap_' c.capacitor(5).Bus1 'var'];
                c.capacitor(5).kvar = 60;
                % Changing Capcontrol settings
                c.capcontrol(1).Capacitor = c.capacitor(1).Name;
                c.capcontrol(1).Vmax = 126;
                
                c.capcontrol(2).Capacitor = c.capacitor(2).Name;
                c.capcontrol(2).Vmax = 126;
                for i=1:length(c.capcontrol)
                    c.capcontrol(i).OFFsetting = 125;
                    c.capcontrol(i).ONsetting = 119;
                    c.capcontrol(i).Vmin = 117;
                    c.capcontrol(i).VoltOverride = 'TRUE';
                    c.capcontrol(i).EventLog = 'yes';
                end
                
                % Changing Generator Settings
                c.generator.Model = 3;
                
                % %% Close all switches
                % seems to be a line to 'nowhere'. when open - voltage is 0V and screws the voltage plots. Does not change
                % the overall results in either mode.
                c.switch(13).Action = 'Close';
                
                %push the voltage regulator to a different bus. it is important for daily
                %simulations. at the end of the feeder the voltage at bus B drops below 93%
                %around 8pm. by pushing the regulator back the voltage in this region can
                %be boosted up to stay above 95%
                c.transformer(3).Buses={'05201400' '05201401'};
                [~, y] = ismember('05201400_05201401', c.line.Name(:));
                n = c.line(y);
                c.line(y) = [];%remove the line to make room for the transformer
                n.bus1 = '05201438';
                n.bus2 = '05201438A';
                n.Name = '05201438_05201438A';
                c.line(end+1) = n;%place a line a the former location of the regulator
                
            case 'wpv_existing'
                %% Fallbrook 'wpv_existing' with existing PV 
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                % add original pv systems
                pv = dssparse(which('f520_pvsystem_pvsystems.dss'));
                c.pvsystem = pv.pvsystem;
                
                % replace the big generator with 2 big pv systems
                c = rmfield(c,'generator');
                PV = c.pvsystem(end);
                PV.Name = 'CG199999_40477_Calle_Roxanne';
                PV.bus1 = '05201644A';
                PV.phases = 3;
                PV.kVA = 999;
                c.pvsystem(end+1) = PV;
                PV.Name = 'CG199999_40177_Calle_Roxanne';
                PV.bus1 = '05201643A';
                PV.kVA = 998.9;
                c.pvsystem(end+1) = PV;
                
                % reset pvsystem's settings
                c = convertPVModel(c);
                for i = 1:length(c.pvsystem)
                    [c.pvsystem(i).bus1, c.pvsystem(i).phases] = findNearestBus(c,c.pvsystem(i).bus1);
                end
                % use constant PQ load model
                for i = 1:length(c.load)
                    c.load(i).model = 1;
                end
                
            case 'wpv_virtual'
                %% Fallbrook 'wpv_virtual' with 432 pv systems
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                % add in the 432 pv systems
                try pv = load(which('scenario2_pvProfiles.mat'));
                catch
                    pv=load([normalizePath('$KLEISSLLAB241') ...
                        'database/gridIntegration/Fallbrook setting/Fallbrook_Scenario2/scenario2_pvProfiles.mat']);
                end
                pv=pv.pv;
                
                % convert it these to dss pv systems
                for i = 1:length(pv)
                    c.pvsystem(i).Name = pv(i).Name;
                    c.pvsystem(i).bus1 = pv(i).bus1;
                    c.pvsystem(i).kVA = pv(i).kVA;
                    [c.pvsystem(i).bus1, c.pvsystem(i).phases] = findNearestBus(c,c.pvsystem(i).bus1);
                end
                c = convertPVModel(c); % convert to model using power as input and output
                
            case {'s1b','balanced_virtualPV_vr_begin'}
                %% Fallbrook {'s1b','balanced_virtualPV_r_begin'}
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_virtual');
                 for i = 1:length(c.regcontrol)
                    c.regcontrol(i).vreg = 1.0*120;
				 end
				 for ii=1:length(c.load)
					 c.load(ii).Model=5;
				 end
				 
			case {'reduced'}
               [ c, opt, d, glc] = feederSetup( feederName, 's1b');
			   critical_nodes={'05201339','05201392','05201643','05202342','05202613','05202098','05202779'};
% 			   [c,~] = FeederReduction(critical_nodes,c);
			   [c,w] = FeederReduction(critical_nodes,c);

				
				case 's1b_si'
                %% s1b with smart inverters!
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
%                 c = to3phases(c);

				%For smart inverters you need to define an xycurve that
				%controls the smart inverter, an xycurve that controls pv
				%efficiency and the xycurve that controls pv P vs. T.
				%Beyond that you need to define the InvControl and the
				%loadshape for PV. have fun.
                
				% Define Inverter control curve
				c.xycurve(1)=dssxycurve;
				c.xycurve(1).name='vv_curve';
				c.xycurve(1).npts=6;
				c.xycurve(1).Yarray='(1, 1, 0, 0, -1, -1)';
				c.xycurve(1).Xarray='(.5, 0.95, .98, 1.02, 1.05, 1.5)';
				
				%Add inverter class to system
				c.InvControl=dssInvControl; 
				c.InvControl.VoltageChangeTolerance=0.1;
				c.InvControl.VarChangeTolerance=0.1;
				c.InvControl.avgwindowlen='30s';
				c.InvControl.DeltaQ_factor=.1;
				
					
				%Must define both efficiency and PT curve of pvsystems
				%define xycurves 
				%efficiency
				c.xycurve(2)=dssxycurve;
				c.xycurve(2).name='Myeff';
				c.xycurve(2).npts=4;
				c.xycurve(2).Xarray=[.1 .2 .4 1.0];
				c.xycurve(2).yarray=[1 1 1 1];
				
				%P-T
				c.xycurve(3)=dssxycurve;
				c.xycurve(3).name='MyPvsT';
				c.xycurve(3).npts=4;
				c.xycurve(3).Xarray=[0 25 75 100];
				c.xycurve(3).yarray=[1.2 1.0 0.8 0.6];
				
				%Define Loadshape for PV %%% This is for validation
				%purposes only %%%
				c.loadshape=dssloadshape;
				c.loadshape.Name='ValidationLoadShape';
				c.loadshape.npts=1;
				c.loadshape.interval=1;
				c.loadshape.mult=1;
				
							
				%loop over pv and add vaiables
				for i = 1:length(c.pvsystem)
				c.pvsystem(i).EFFCURVE='Myeff';
				c.pvsystem(i).PTCurve='MyPvsT';
				c.pvsystem(i).daily='ValidationLoadShape';
				end
            
			case 's1b_si_nolargesystems'	
				%% S1b with smart inverters but removing large system
				[ c, opt, d, glc] = feederSetup( feederName, 's1b_si');
				
				%Remove Inverters from large system
				%Find all PVnames not equal to large systems
				idex = find([c.pvsystem{:}.Pmpp] <800);
				%Need to reformat list to a long string with brackets
				s=[];
				for i=1:length(idex)
				str=c.pvsystem(idex(i)).Name;
				s=[s ' ' str];
				end
				s=['[' s ']'];
				c.InvControl.PVSystemList=s;
				
	
            case 'balanced_virtualPV_vr_middle'
                % modify regcontrol, set voltage level at mid point of each branch to 1.0 p.u
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                mbus = {'05201328','05201366','05201455','05202042','05202424','05202331','05201235'};
                for i = 1:length(c.regcontrol)
                    c.regcontrol(i).vreg = 120;
                    c.regcontrol(i).bus = mbus{i};
                end
%                 c.regcontrol(1).delay = 15;
%                 c.regcontrol(2).delay = 45;
%                 c.regcontrol(3).delay = 75;
                
            case 'balanced_virtualPV_vr_begin_specifyBus'
                %% Fallbrook balanced_vr_begin
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                
                % modify regcontrol, set voltage level at mid point of each branch to 1.0 p.u
                mbus = {'05201317A','05201350A','05201401','05201949A','05202250A','05202315A','0520'};
                for i = 1:length(c.regcontrol)
                    c.regcontrol(i).vreg = 120;
                    c.regcontrol(i).bus = mbus{i};
                    c.regcontrol(i).band = 1.2;
                end
                
            case {'distributed-aggregated','dist-agg'}
                %% Fallbrook distributed PV allocation + aggregated profile
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                
                % aggregated profile
                opt.fcProfileId = 45;
                
                % forming distributed PV allocation by removing 2 big pvsites
                c.pvsystem(44:45) = [];
                %opt.fcProfileId = setdiff(1:432,[44,45]);
                
            case {'distributed-disaggregated','dist-disagg'}
                %% Fallbrook distributed PV allocation + disaggregated profiles
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                
                % forming distributed PV allocation by removing 2 big pvsites
                c.pvsystem(44:45) = [];
                opt.fcProfileId = setdiff(1:432,[44,45]);
                
            case 's0' % 432 PV systems (inluding 2MW site) + single (aggregated) PV profile
                %% Fallbrook S0
                [ c, opt, d, glc] = feederSetup( feederName, 'balanced_virtualPV_vr_begin');
                opt.fcProfileId = 45;
                
            case 's1' % 432 PV systems (inluding 2MW site) + diff (disaggregated) PV profiles
                %% Fallbrook S1
                [ c, opt, d, glc] = feederSetup( feederName, 'balanced_virtualPV_vr_begin');
                
            case 's2' % 45 PV systems (centralized setup) + diff (disaggregated) PV profiles
                %% Fallbrook S2
                [ c, opt, d, glc] = feederSetup( feederName, 'balanced_virtualPV_vr_begin');
                c.pvsystem(46:end) = [];
                
            case 's3' % 430 PV systems (without 2MW site, distributed setup) + diff (disaggregated) PV profiles
                %% Fallbrook S3
                [ c, opt, d, glc] = feederSetup( feederName, 'balanced_virtualPV_vr_begin');
                c.pvsystem(44:45) = [];
                opt.fcProfileId = setdiff(1:432,[44,45]);
                
            case 's0b' % 432 PV systems (inluding 2MW site) + single (aggregated) PV profile
                %% Fallbrook 's0b' S0 with VR set point at beginning of the feeder
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                opt.fcProfileId = 45;
                
            case 's2b' % 45 PV systems (centralized setup) + diff (disaggregated) PV profiles
                %% Fallbrook 's2b' S2 with VR set point at beginning of the feeder
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                c.pvsystem(46:end) = [];
                
            case 's3b' % 430 PV systems (without 2MW site, distributed setup) + diff (disaggregated) PV profiles
                %% Fallbrook 's3b' edS3 with VR set point at beginning of the feeder
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                c.pvsystem(44:45) = [];
                opt.fcProfileId = setdiff(1:432,[44,45]);
                
            case 'distributed_control_45pvs' % with pv systems 1 bus closer to the main line
                %% Fallbrook distributed_control
                % using balanced 3 phase circuit with vr reference point at the beginning of a branch
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                opt.overlapSiteExist = 0;
                % remove the generator
                c = rmfield(c,'generator');
                % reset pvsystems to use simple power dispatch method. Take a look at dsspvsystem doc for more info.
                for i = 1:length(c.pvsystem)
                    c.pvsystem(i).pmpp = c.pvsystem(i).kVA;
                    c.pvsystem(i).irradiance = 1;
                    c.pvsystem(i).pf = 1;
                    c.pvsystem(i).kVA = c.pvsystem(i).kVA;
                end
                % convert to balanced 3 phases
                c = to3phases(c);
                for i = 1:length(c.regcontrol)
                    c.regcontrol(i).vreg = 120*1.05; % set it to the ideal 1.0 pu
                end
                c = mvCompCloser(c,'pvsystem',1); % move pv system one bus closer to the main line
                c = addEnergyStorage(c,c.pvsystem); % add ES to each pvsystem's bus
                c = comp2line(c,c.transformer(1:end-1)); % replace all tranx but the substation one with lines
                % use constant PQ load model
                for i = 1:length(c.load)
                    c.load(i).model = 1;
                end
                
            case {'s1b_nocap'}
                %% Fallbrook 's1b_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
            otherwise
                error('Not supported config!');
        end
    case {'alpine','355','d'}
        %% Alpine feeder {'alpine','355','d'}
        switch lower(configId)
            case 'original'
                %% Alpine original
                d = excel2obj([conf.gridIntDir '/dssconversion/custdata/355ForENERNEX.xlsx']);
                glc = excel2obj( 'linecode.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc );
                w=1;
            case 'validated'
                %% Alpine validated
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                
                % %% Change Circuit settings
                c.circuit.bus1 = 'SourceBus';
                % c.circuit.pu = 1.00;% changed from 1.00
                c.circuit.basekv = 12;
                c.basevoltages = [12 12];
                
                % reset setting of the existing transformer
                id = 1;
                c.transformer(id) = resetTransformer(c.transformer(id));
                c.transformer(id).kvs = [12 12];
                c.transformer(id).Phases = 1;
                c.transformer(id).buses = {'03551328.1','03551328A.1'};
                
                % add 'sourcebus' to buslist
                c.buslist.id{end+1} = 'SourceBus';
                c.buslist.coord(end+1,:) = c.buslist.coord(1,:);
                c.buslist.coord(end,1) = c.buslist.coord(end,1) - 2;
                
                % add substation transformer
                id = 2;
                c.transformer(id).Name = 'ALPINE';
                c.transformer(id).Buses = {'SourceBus' '0355'};
                c.transformer(id).kVs = [12 12];
                c.transformer(id) = resetTransformer(c.transformer(id));
                
                % add regcontrol at substation
                x = dssregcontrol;
                x.Name = 'regSubstation';
                x = resetRegcontrol(x);
                x.vreg = 120*1.01;
                x.transformer = c.transformer(2).Name;
                c.regcontrol(1) = x;
                
                % add another regcontrol for the secondary feeder
                id = 2;
                c.regcontrol(id).Name = 'reg2';
                c.regcontrol(id) = resetRegcontrol(c.regcontrol(id));
                c.regcontrol(id).transformer = c.transformer(1).Name;
                c.regcontrol(id).vreg = 120*1.0;
                
                % change the 2 first line length
                c.line(1).Length = c.line(1).Length*0.3;
                c.line(2).Length = c.line(2).Length*0.55;
                
                c.capacitor(2) = [];
                c.capcontrol(2) =[];
                c.capacitor(1).kvar = 1350;
                c.capacitor.Numsteps = 20;
                
                kVAr_factor = 0.24; % factor used to adjust the kVar consumption of the loads Phase 1 and 3
                kVAr_factor3 = 0.23;
                kW_factor = 0.78;
                
                loads_1Phase = c.load([c.load.Phases]==1);
                [~, b1] = ismember({loads_1Phase.Name}, {c.load.Name});
                for i=1:length(b1)
                    c.load(b1(i)).Kv = 6.9282;
                    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*kVAr_factor;
                    c.load(b1(i)).Kw =c.load(b1(i)).Kw*kW_factor;
                end
                loads_3Phase = c.load([c.load.Phases]==3);
                [~, b3] = ismember({loads_3Phase.Name}, {c.load.Name});
                for j=1:length(b3)
                    c.load(b3(j)).Kv = 12;
                    c.load(b3(j)).Kvar = c.load(b3(j)).Kw*kVAr_factor3;
                    c.load(b3(j)).Kw = c.load(b3(j)).Kw*kW_factor;
                end
                c.load(355).Kvar = c.load(355).Kvar*2.8;
                
                % remove weird connections created by switches
                l = {'SW_0355201_0357901','SW_035560_0357902',...
                    'SW_03552601_1458901','SW_0355302_0356901','TF_03551325_0356902'}; % line names
                lid = find(ismember({c.line.Name},l)); % line id
                x = {};
                for i = 1:length(lid)
                    x = [x; ['Line.' c.line(lid(i)).Name]];
                end
                sid = ismember({c.switch.SwitchedObj},x); %
                c.line(lid) = []; c.switch(sid) = [];
                
                % remove isolated buses
                c = rmBusByName(c,{'0357902','0357901','1458901'});
                
                % reset voltage regulator settings
                for i = 1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i)); 
                end
                for i = 1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i)); 
                end
                w=1;
            case 'wpv_existing'
                %% Alpine with existing PV 
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                CSI_data = excel2obj('data/GIS_CSI_v4.1.xlsx');
                c = addExistingPV(c,feederName,CSI_data);
                % use constant PQ load model
                for i = 1:length(c.load)
                    c.load(i).model = 5;
                end
                w=1;
            case 'wpv_existing_balanced3phase'
                %% Alpine with existing PV in balanced 3-phase config
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                c = to3phases(c);
                w=1;
            case 'wpv_virtual' 
                %% Alpine with all virtual + real pv
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                pv = load([normalizePath('$KLEISSLLAB24-1')  '/database/gridIntegration/PVimpactPaper/virtualPV/f355scenario2_pv.mat']);
                c.pvsystem = pv.pv;
                c = convertPVModel(c);
                for i = 1:length(c.pvsystem)
                    [c.pvsystem(i).bus1, c.pvsystem(i).phases] = findNearestBus(c,c.pvsystem(i).bus1);
                end
                w=1;
            case {'wpv_virtual_balanced_3phase','s1b'}
                %% Alpine {'wpv_virtual_balanced_3phase','s1b'} with all virtual + real PV in balanced 3-phase
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_virtual');
				w=1;
			case {'wpv_virtual_balanced_3phase','s1b_refined'}
                %% Alpine {'wpv_virtual_balanced_3phase','s1b'} with all virtual + real PV in balanced 3-phase
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_virtual');
				for ii=1:length(c.transformer)
					Phases=sum(ismember(char(c.transformer(ii).buses(2)),'\.'));
					if Phases==0
						c.transformer(ii).Phases=3;
					else
						c.transformer(ii).Phases=Phases;
					end
				end
				w=1;
			case {'reduced'}
				[ c, opt, d, glc] = feederSetup( feederName, 's1b_refined');

% 				critical_nodes=c.buslist.id(round((length(c.buslist.id)-2)*rand(kk,1)+1));
% 				critical_nodes=unique(critical_nodes);
% 				while length(critical_nodes)<kk
% 					cnode_tmp=c.buslist.id(round((length(c.buslist.id)-2)*rand(kk,1)+1));
% 					critical_nodes=[critical_nodes;cnode_tmp];
% 					critical_nodes=unique(critical_nodes);
% 				end	
% 				timeStop=tic;
% 				[c] = FeederReduction(critical_nodes,c);
% 				t=toc(timeStop);
% 				c.SimTime=t;
criticalBuses=c.buslist.id([11,13,14]);

for ii=1:length(criticalBuses)
	criticalBuses(ii)=regexprep(criticalBuses(ii),'-','_');
end

[c] = FeederReduction_SE(criticalBuses,c);
% 				critical_nodes=c.buslist.id(round(length(c.buslist.id)*rand(20,1)));
% 				[c] = FeederReduction(critical_nodes,c);

% 				
			case {'reduced_loadmodel_5'}
				[ c, opt, d, glc] = feederSetup( feederName, 'reduced',1);
				for ii=1:length(c.load)
					c.load(ii).Model=5;
				end
			case {'reduced_loadmodel_1'}
				[ c, opt, d, glc] = feederSetup( feederName, 'reduced',1);
				for ii=1:length(c.load)
					c.load(ii).Model=1;
				end
			case {'reduced_loadmodel_2'}
				[ c, opt, d, glc] = feederSetup( feederName, 'reduced',1);
				for ii=1:length(c.load)
					c.load(ii).Model=2;
				end
				
			case 's1b_si'
				%% s1b with smart inverters!
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
				
				% Define Inverter control curve
				c.xycurve(1)=dssxycurve;
				c.xycurve(1).name='vv_curve';
				c.xycurve(1).npts=6;
				c.xycurve(1).Yarray='(1, 1, 0, 0, -1, -1)';
				c.xycurve(1).Xarray='(.5, 0.95, .98, 1.02, 1.05, 1.5)';
				
				%Add inverter class to system
				c.InvControl=dssInvControl; 
				c.InvControl.VoltageChangeTolerance=0.1;
				c.InvControl.VarChangeTolerance=0.1;
				c.InvControl.avgwindowlen='30s';
				c.InvControl.DeltaQ_factor=.1;
				
					
				%Must define both efficiency and PT curve of pvsystems
				%define xycurves 
				%efficiency
				c.xycurve(2)=dssxycurve;
				c.xycurve(2).name='Myeff';
				c.xycurve(2).npts=4;
				c.xycurve(2).Xarray=[.1 .2 .4 1.0];
				c.xycurve(2).yarray=[1 1 1 1];
				
				%P-T
				c.xycurve(3)=dssxycurve;
				c.xycurve(3).name='MyPvsT';
				c.xycurve(3).npts=4;
				c.xycurve(3).Xarray=[0 25 75 100];
				c.xycurve(3).yarray=[1.2 1.0 0.8 0.6];
				
				%Define Loadshape for PV %%% This is for validation
				%purposes only %%%
				c.loadshape=dssloadshape;
				c.loadshape.Name='ValidationLoadShape';
				c.loadshape.npts=1;
				c.loadshape.interval=1;
				c.loadshape.mult=1;
				
							
				%loop over pv and add vaiables
				for i = 1:length(c.pvsystem)
				c.pvsystem(i).EFFCURVE='Myeff';
				c.pvsystem(i).PTCurve='MyPvsT';
				c.pvsystem(i).daily='ValidationLoadShape';
				end
				
				
				case 's1b_si_50'
					[ c, opt, d, glc] = feederSetup( feederName, 's1b_si');
                load('c:\users\zactus\gridIntegration\SmartInverterSims\AlpinePvBuses.mat')
				
				for ii=1:length(pvsystPast)
				c.InvControl(ii)=dssInvControl; 
				c.InvControl(ii).Name=['Inv' num2str(ii)];
				c.InvControl(ii).VoltageChangeTolerance=0.1;
				c.InvControl(ii).VarChangeTolerance=0.1;
				c.InvControl(ii).avgwindowlen='30s';
				c.InvControl(ii).DeltaQ_factor=.1;
				c.InvControl(ii).PVSystemList=c.pvsystem(pvsystPast(ii)).name;
				end
                
            case {'s1b_nocap'}
                %% Alpine 's1b_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
            
            case {'aggregated','agg'}
                %% Alpine aggregated/ single PV profile case
                [ c, opt, d, glc] = feederSetup( feederName, 's1b_nocap');
                
                % aggregated profile
                opt.fcProfileId = 1;
                
            otherwise
                error('Not supported config!');
        end
    case {'pointloma','480','cabrillo','b'}
        %% Point Loma feeder {'pointloma','480','cabrillo','b'}
        switch lower(configId)
            case 'original'
                %% Pointloma original
                d = excel2obj([conf.gridIntDir '/dssconversion/custdata/480ForENERNEX.xlsx']);
                glc = excel2obj( 'linecode.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc );
                validatePQ = 0;
                w=1;
            case 'validated'
                %% Pointloma 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                
                kVAr_factor = 1.5; % factor used to adjust the kVar consumption of the loads Phase 1 and 3
                kW_factor = 1.078;
                pf = 0.95;
                loads_1Phase = c.load([c.load.Phases]==1);
                [~, b1] = ismember({loads_1Phase.Name}, {c.load.Name});
                for i=1:length(b1)
                    c.load(b1(i)).Kv = 6.9282;
                    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*kVAr_factor;
                    c.load(b1(i)).Kw =c.load(b1(i)).Kw*kW_factor;
                    
                end
                loads_3Phase = c.load([c.load.Phases]==3);
                [~, b3] = ismember({loads_3Phase.Name}, {c.load.Name});
                for j=1:length(b3)
                    c.load(b3(j)).Kv = 12;
                    c.load(b3(j)).Kw = c.load(b3(j)).Kw*kW_factor;
                    c.load(b3(j)).Kvar = c.load(b3(j)).Kvar*kVAr_factor;
                end
                
                c.load(210).Kvar = c.load(210).Kw*sqrt(1/(pf*pf)-1);
                c.load(205).Kvar = c.load(205).Kw*sqrt(1/(pf*pf)-1);
                c.load(212).Kvar = c.load(212).Kw*sqrt(1/(pf*pf)-1);
                % Some loads are connected to 1 phase while other connect to all 3 phases
                find([c.load.Phases]==3) % to find 3 phases loads
                
                % %% Change transformer settings
                c.transformer(1).Conns = {'wye' 'wye'};
                c.basevoltages = [12 12 12];
                % %% Delete capcontrolers
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                % %% Optimizing capacitor banks to match reactive power
                % Modify existing caps :
                
                % - 480_70CF
                c.capacitor(1).Kvar = 150;
                c.capacitor(1).Numsteps = 1;
                c.capacitor(1).Name = ['cap_' c.capacitor(1).Name];
                
                c.capacitor(3) = c.capacitor(1);
                c.capacitor(3).Name = 'cap_480_70F_var';
                c.capacitor(3).Kvar = 1120;
                c.capacitor(3).Numsteps = floor(c.capacitor(3).Kvar/40);
                
                % - 480_30CF
                c.capacitor(2).Kvar = 150;
                c.capacitor(2).Numsteps = 1;
                c.capacitor(2).Name = ['cap_' c.capacitor(2).Name];
                
                c.capacitor(4) = c.capacitor(2);
                c.capacitor(4).Name = 'cap_480_30F_var';
                c.capacitor(4).Kvar = 1130;
                c.capacitor(4).Numsteps = floor(c.capacitor(4).Kvar/40);
                
                % Adding cap banks at appropriate locations 
                c.capacitor(5) = c.capacitor(2);
                c.capacitor(5).Bus1 = '048048';
                c.capacitor(5).Name = ['cap_' c.capacitor(5).Bus1];
                c.capacitor(5).Kvar = 60;
                c.capacitor(5).Numsteps = floor(c.capacitor(5).Kvar/10);
                
                c.capacitor(6) = c.capacitor(2);
                c.capacitor(6).Bus1 = '048041';
                c.capacitor(6).Name = ['cap_' c.capacitor(6).Bus1];
                c.capacitor(6).Kvar = 20;
                c.capacitor(6).Numsteps = floor(c.capacitor(6).Kvar/5);
                
                % remove secondary transformer along with the single load it serves.
                % this load is negletible (8 kW) so won't affect the result
                c = rmfield(c,'transformer');
                c = rmElemByName(c,'load','l04804903_Y');
                c = rmElemByName(c,'line','04804902A_04804903');
                c = rmElemByName(c,'line','SW_04808240_00519902');
                c = rmElemByFieldValue(c,'switch','SwitchedObj','Line.SW_04808240_00519902');
                c = rmElemByName(c,'line','SW_048042_00514570');
                c = rmElemByFieldValue(c,'switch','SwitchedObj','Line.SW_048042_00514570');
                c = rmElemByName(c,'line','SW_048037_04819901');
                c = rmElemByFieldValue(c,'switch','SwitchedObj','Line.SW_048037_04819901');
                
                % add substation transformer
                x = dsstransformer;
                x.Name = 'substation';
                x = resetTransformer(x);
                x.kvs = [12 12];
                x.Buses = {'sourcebus', '0480'};
                c.transformer(1) = x;
                
                % add sourcebus 
                id = ismember(c.buslist.id,'0480');
                c.buslist.id{end+1} = 'sourcebus';
                c.buslist.coord(end+1,:) = c.buslist.coord(id,:);
                c.buslist.coord(end,1) = c.buslist.coord(end,1) + 2;
                
                % rewire the substation bus
                c.circuit.bus1 = 'sourcebus';
                
                % create a regcontrol for the new trans
                x = dssregcontrol;
                x.Name = c.transformer(1).Name;
                x = resetRegcontrol(x);
                x.transformer = c.transformer(1).Name;
                x.vreg = .98*120;
                c.regcontrol(1) = x;
                
                % reset voltage regulator settings
                for i = 1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i)); 
                end
                for i = 1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i)); 
                end
                
            case 'wpv_existing'
                %% Pointloma 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                CSI_data = excel2obj('data/GIS_CSI_v4.1.xlsx');
                c = addExistingPV(c,feederName,CSI_data);
                % use constant PQ load model
                for i = 1:length(c.load)
                    c.load(i).model = 1;
                end
                
            case 'wpv_existing_balanced3phase'
                %% Pointloma 'wpv_existing_balanced3phase'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                c = to3phases(c);
                
            case 'wpv_virtual' % with all virtual + real pv
                %% Pointloma 'wpv_virtual'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                pv = load([normalizePath('$KLEISSLLAB24-1')  '/database/gridIntegration/PVimpactPaper/virtualPV/f480scenario2_pv.mat']);
                c.pvsystem = pv.pv;
                c = convertPVModel(c);
                for i = 1:length(c.pvsystem)
                    [c.pvsystem(i).bus1, c.pvsystem(i).phases] = findNearestBus(c,c.pvsystem(i).bus1);
                end
                
            case {'wpv_virtual_balanced_3phase','s1b'}
                %% Pointloma {'wpv_virtual_balanced_3phase','s1b'}
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_virtual');
%                 c = to3phases(c);
                
                % change load model to constant P&Q
                for i = 1:length(c.load)
                    c.load(i).model = 5;
				end
				
			case {'reduced'}
               [ c, opt, d, glc] = feederSetup( feederName, 's1b');
			   critical_nodes={'0480503';'04802323';'04809944';'04801424'};
			   [c,w] = FeederReduction(critical_nodes,c);

				case 's1b_si'
                %% s1b with smart inverters!
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
%                 c = to3phases(c);

				%For smart inverters you need to define an xycurve that
				%controls the smart inverter, an xycurve that controls pv
				%efficiency and the xycurve that controls pv P vs. T.
				%Beyond that you need to define the InvControl and the
				%loadshape for PV. have fun.
                
				% Define Inverter control curve
				c.xycurve(1)=dssxycurve;
				c.xycurve(1).name='vv_curve';
				c.xycurve(1).npts=6;
				c.xycurve(1).Yarray='(1, 1, 0, 0, -1, -1)';
				c.xycurve(1).Xarray='(.5, 0.95, .98, 1.02, 1.05, 1.5)';
				
				%Add inverter class to system
				c.InvControl=dssInvControl; 
				c.InvControl.VoltageChangeTolerance=0.1;
				c.InvControl.VarChangeTolerance=0.1;
				c.InvControl.avgwindowlen='30s';
				c.InvControl.DeltaQ_factor=.1;
				
					
				%Must define both efficiency and PT curve of pvsystems
				%define xycurves 
				%efficiency
				c.xycurve(2)=dssxycurve;
				c.xycurve(2).name='Myeff';
				c.xycurve(2).npts=4;
				c.xycurve(2).Xarray=[.1 .2 .4 1.0];
				c.xycurve(2).yarray=[1 1 1 1];
				
				%P-T
				c.xycurve(3)=dssxycurve;
				c.xycurve(3).name='MyPvsT';
				c.xycurve(3).npts=4;
				c.xycurve(3).Xarray=[0 25 75 100];
				c.xycurve(3).yarray=[1.2 1.0 0.8 0.6];
				
				%Define Loadshape for PV %%% This is for validation
				%purposes only %%%
				c.loadshape=dssloadshape;
				c.loadshape.Name='ValidationLoadShape';
				c.loadshape.npts=1;
				c.loadshape.interval=1;
				c.loadshape.mult=1;
				
							
				%loop over pv and add vaiables
				for i = 1:length(c.pvsystem)
				c.pvsystem(i).EFFCURVE='Myeff';
				c.pvsystem(i).PTCurve='MyPvsT';
				c.pvsystem(i).daily='ValidationLoadShape';
				end
				
				case 's1b_si_50'
					[ c, opt, d, glc] = feederSetup( feederName, 's1b_si');
                load('c:\users\zactus\gridIntegration\SmartInverterSims\PointLomaPvBuses.mat')
				
				for ii=1:length(pvsystPast)
				c.InvControl(ii)=dssInvControl; 
				c.InvControl(ii).Name=['Inv' num2str(ii)];
				c.InvControl(ii).VoltageChangeTolerance=0.1;
				c.InvControl(ii).VarChangeTolerance=0.1;
				c.InvControl(ii).avgwindowlen='30s';
				c.InvControl(ii).DeltaQ_factor=.1;
				c.InvControl(ii).PVSystemList=c.pvsystem(pvsystPast(ii)).name;
				end
                
            case {'s1b_nocap'}
                %% Pointloma s1b_nocap
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
            case {'aggregated','agg'}
                %% Pointloma aggregated/ single PV profile case without capacitors
                [ c, opt, d, glc] = feederSetup( feederName, 's1b_nocap');
                
                % aggregated profile
                opt.fcProfileId = 1;
            otherwise
                error('Not supported config!');
        end
        
    case {'ramona','971','creelman','e'}
        %% Ramona feeder {'ramona','971','creelman','e'}
        switch lower(configId)
            case 'original'
                %% ramona 'original'
                d = excel2obj([conf.gridIntDir '/dssconversion/custdata/971ForENERNEX.xlsx']);
                glc = excel2obj( 'linecode.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc );
                validatePQ = 0;
                w=1;
            case 'validated'
                %% ramona 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                
                % %% Change Circuit settings
                c.circuit.bus1 = 'SourceBus';
                % c.circuit.pu = 1.00;% changed from 1.00
                c.circuit.basekv = 12;
                c.basevoltages = [12 12];
                
                % %% Changing Load Settings
                kVAr_factor1 = 1.294;
                kW_factor1 = 0.98;
                
                loads_1Phase = c.load([c.load.Phases]==1);
                [~, b1] = ismember({loads_1Phase.Name}, {c.load.Name});
                for i=1:length(b1)
                    c.load(b1(i)).Kv = 6.9282;
                    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*kVAr_factor1;
                    c.load(b1(i)).Kw =c.load(b1(i)).Kw*kW_factor1;
                end
                
                for j=1:length(c.switch)
                    c.switch(j).Action = 'Close';
                end
                
                % add substation transformer
                c.transformer(2) = dsstransformer;
                c.transformer(end).Name = c.circuit.Name;
                c.transformer(end).Buses = {'SourceBus' '0971'};
                c.transformer(end).Conns = {'y' 'y'};
                c.transformer(end).kVs = [12 12];
                c.transformer(end).kVAs = [28000 28000];% changed from [10000 10000]
                c.transformer(end).XHL = 1.0871;
                c.transformer(end).sub = 'y';
                c.transformer(end).Rs = [0.103 0.103];
                c.buslist.id(end+1) = {'SourceBus'};
                c.buslist.coord(end+1,:) = c.buslist.coord(1,:);
                
                % add regcontrol to substation transformer
                c.regcontrol(2).transformer = c.transformer(2).Name;
                c.regcontrol(2) = resetRegcontrol(c.regcontrol(2));
                c.regcontrol(2).Name = c.regcontrol(2).transformer ;
                c.regcontrol(2).vreg = 0.99*120;
                
                % change the winding type of existing transformer
                c.transformer(1).Conns = {'y' , 'y'};
                c.transformer(1).Wdg = [];
                c.transformer(1).kVAs = [20000 20000];
                c.transformer(1).Phases = 3;
                
                % modify the secondary tranx's config
                c.regcontrol(1).winding = 2;
                c.regcontrol(1).vreg = 120*.99;
                %c.regcontrol(1).bus = '0971146a';
                %c.regcontrol(1).EventLog = 'true';
                c.regcontrol(1).vlimit = [];
                c.regcontrol(1).ptratio = 57.75;
                c.regcontrol(1).band = 1;
                
                % reset voltage regulator settings
                for i = 1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i)); 
                end
                for i = 1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i)); 
                end
                
            case 'wpv_existing'
                %% ramona 'wpv_existing'
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                CSI_data = excel2obj('data/GIS_CSI_v4.1.xlsx');
                c = addExistingPV(c,feederName,CSI_data);
                % use constant PQ load model
                for i = 1:length(c.load)
                    c.load(i).model = 1;
                end
                
            case 'wpv_existing_balanced3phase'
                %% ramona 'wpv_existing_balanced3phase'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                c = to3phases(c);
                
            case 'wpv_virtual' % with all virtual + real pv
                %% ramona 'wpv_virtual'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                pv = load([normalizePath('$KLEISSLLAB24-1')  '/database/gridIntegration/PVimpactPaper/virtualPV/f971scenario2_pv.mat']);
                c.pvsystem = pv.pv;
                c = convertPVModel(c);
                for i = 1:length(c.pvsystem)
                    [c.pvsystem(i).bus1, c.pvsystem(i).phases] = findNearestBus(c,c.pvsystem(i).bus1);
                end
                
            case {'wpv_virtual_balanced_3phase','s1b'}
                %% ramona {'wpv_virtual_balanced_3phase','s1b'}
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_virtual');
                c = to3phases(c);
                
                % change load model to constant P&Q
                for i = 1:length(c.load)
                    c.load(i).model = 1;
				end
				
				case {'reduced'}
               [ c, opt, d, glc] = feederSetup( feederName, 's1b');
			   critical_nodes=c.buslist.id(round(rand(4,1)*length(c.buslist)));
			   [c,w] = FeederReduction(critical_nodes,c);
			   
				case 's1b_si'
                %% s1b with smart inverters!
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
%                 c = to3phases(c);

				%For smart inverters you need to define an xycurve that
				%controls the smart inverter, an xycurve that controls pv
				%efficiency and the xycurve that controls pv P vs. T.
				%Beyond that you need to define the InvControl and the
				%loadshape for PV. have fun.
                
				% Define Inverter control curve
				c.xycurve(1)=dssxycurve;
				c.xycurve(1).name='vv_curve';
				c.xycurve(1).npts=6;
				c.xycurve(1).Yarray='(1, 1, 0, 0, -1, -1)';
				c.xycurve(1).Xarray='(.5, 0.95, .98, 1.02, 1.05, 1.5)';
				
				%Add inverter class to system
				c.InvControl=dssInvControl; 
				c.InvControl.VoltageChangeTolerance=0.1;
				c.InvControl.VarChangeTolerance=0.1;
				c.InvControl.avgwindowlen='30s';
				c.InvControl.DeltaQ_factor=.1;
				
					
				%Must define both efficiency and PT curve of pvsystems
				%define xycurves 
				%efficiency
				c.xycurve(2)=dssxycurve;
				c.xycurve(2).name='Myeff';
				c.xycurve(2).npts=4;
				c.xycurve(2).Xarray=[.1 .2 .4 1.0];
				c.xycurve(2).yarray=[1 1 1 1];
				
				%P-T
				c.xycurve(3)=dssxycurve;
				c.xycurve(3).name='MyPvsT';
				c.xycurve(3).npts=4;
				c.xycurve(3).Xarray=[0 25 75 100];
				c.xycurve(3).yarray=[1.2 1.0 0.8 0.6];
				
				%Define Loadshape for PV %%% This is for validation
				%purposes only %%%
				c.loadshape=dssloadshape;
				c.loadshape.Name='ValidationLoadShape';
				c.loadshape.npts=1;
				c.loadshape.interval=1;
				c.loadshape.mult=1;
				
						
				
				%loop over pv and add vaiables
				for i = 1:length(c.pvsystem)
				c.pvsystem(i).EFFCURVE='Myeff';
				c.pvsystem(i).PTCurve='MyPvsT';
				c.pvsystem(i).daily='ValidationLoadShape';
				end
				
				case 's1b_si_50'
					[ c, opt, d, glc] = feederSetup( feederName, 's1b_si');
                load('c:\users\zactus\gridIntegration\SmartInverterSims\RamonaPvBuses.mat')
				
				for ii=1:length(pvsystPast)
				c.InvControl(ii)=dssInvControl; 
				c.InvControl(ii).Name=['Inv' num2str(ii)];
				c.InvControl(ii).VoltageChangeTolerance=0.1;
				c.InvControl(ii).VarChangeTolerance=0.1;
				c.InvControl(ii).avgwindowlen='30s';
				c.InvControl(ii).DeltaQ_factor=.1;
				c.InvControl(ii).PVSystemList=c.pvsystem(pvsystPast(ii)).name;
				end
                
            case {'s1b_nocap'}
                %% ramona 's1b_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
            case {'aggregated','agg'}
                %% Ramona aggregated/ single PV profile case without capacitors
                [ c, opt, d, glc] = feederSetup( feederName, 's1b_nocap');
                
                % aggregated profile
                opt.fcProfileId = 1;
            otherwise
                error('Not supported config!');
        end
        
    case {'valleycenter','909','c'}
        %% Valley center feeder {'valleycenter','909','c'}
        switch lower(configId)
            case 'original'
                %% vallaycenter 'original'
                d = excel2obj([conf.gridIntDir '/dssconversion/custdata/909ForENERNEX.xlsx']);
                glc = excel2obj( 'linecode.xlsx' ); glc = glc.LineCode;
                c = dssconversion( d, glc );
                w=1;
            case 'validated'
                %% valleycenter 'validated'
                [ c, opt, d, glc] = feederSetup( feederName, 'original');
                
                % %% Change Circuit settings
                c.circuit.bus1 = 'SourceBus';
                % c.circuit.pu = 1.00;% changed from 1.00
                c.circuit.basekv = 12;
                c.basevoltages = [12 12];
                
                 % add substation transformer
                c.transformer(1) = dsstransformer;
                c.transformer(1).Name = c.circuit.Name;
                c.transformer(1).Buses = {'SourceBus' '0909'};
                c.transformer(1).Conns = {'y' 'y'};
                c.transformer(1).kVs = [12 12];
                c.transformer(1).kVAs = [28000 28000];% changed from [10000 10000]
                c.transformer(1).XHL = 1.0871;
                c.transformer(1).sub = 'y';
                c.transformer(1).Rs = [0.103 0.103];
                c.buslist.id(end+1) = {'SourceBus'};
                id = find(ismember(c.buslist.id,'0909'));
                c.buslist.coord(end+1,:) = c.buslist.coord(id,:);
                
                c.regcontrol = dssregcontrol;
                c.regcontrol = resetRegcontrol(c.regcontrol);
                c.regcontrol.Name = 'reg1';
                c.regcontrol.transformer = c.transformer.Name;
                c.regcontrol.vreg = 1.01*120;
                
                % Change load consumption
                kVAr_factor1 = 0.9; % factor used to adjust the kVar consumption of the loads Phase 1
                kW_factor1 = 1.5;
                kVAr_factor3 = 0.95; % factor used to adjust the kVar consumption of the loads Phase 3
                kW_factor3 = 1.5;
                
                loads_1Phase = c.load([c.load.Phases]==1);
                [~, b1] = ismember({loads_1Phase.Name}, {c.load.Name});
                for i=1:length(b1)
                    c.load(b1(i)).Kv = 6.9282;
                    c.load(b1(i)).Kvar = c.load(b1(i)).Kvar*kVAr_factor1;
                    c.load(b1(i)).Kw =c.load(b1(i)).Kw*kW_factor1;
                end
                loads_3Phase = c.load([c.load.Phases]==3);
                [~, b3] = ismember({loads_3Phase.Name}, {c.load.Name});
                for j=3:length(b3)
                    c.load(b3(j)).Kv = 12;
                    c.load(b3(j)).Kvar = c.load(b3(j)).Kvar*kVAr_factor3;
                    c.load(b3(j)).Kw = c.load(b3(j)).Kw*(kW_factor3);
                end
                c.capacitor.Kvar = 1400;
                c.capacitor.Numsteps = 30;
                kwfac=1.107;
                Kvarfact_end = 1.1;
                c.load(1).Kvar=c.load(1).Kvar * (Kvarfact_end+0.00);
                c.load(2).Kvar=c.load(2).Kvar * (Kvarfact_end+0.00);
                c.load(1).Kw=c.load(1).Kw* (kwfac);
                c.load(2).Kw=c.load(2).Kw* (kwfac);
                
                % remove switches that are not doing anything good
                c = rmElemByName(c,'line','SW_0909100_021637');
                c = rmElemByFieldValue(c,'switch','SwitchedObj','Line.SW_0909100_021637');
                c = rmElemByName(c,'line','SW_090956A_103055AN');
                c = rmElemByFieldValue(c,'switch','SwitchedObj','Line.SW_090956A_103055AN');
                c = rmElemByName(c,'line','SW_090935_103033AN');
                c = rmElemByFieldValue(c,'switch','SwitchedObj','Line.SW_090935_103033AN');
                % remove isolated bus from bus list
                id = ismember(lower(c.buslist.id),{'103033an','103055an'});
                c.buslist.id(id) = [];
                c.buslist.coord(id,:) = [];
                
                % create an extra transformer to validate the voltage profile
                % it can be seen that a trans/ voltage regulator is missing from the voltage comparison plot
                % add a bus next to the first bus of the transformer
                c.buslist.id(end+1) = {'090953a'};
                id = find(ismember(c.buslist.id,'090953'));
                c.buslist.coord(end+1,:) = c.buslist.coord(id,:);
                c.buslist.coord(end,1) = c.buslist.coord(end,1) + 2;
                % create new trans
                x = dsstransformer;
                x.Name = 'tran_2';
                x = resetTransformer(x);
                x.kvs = [12 12];
                x.Buses = {'090953', '090953a'};
                c.transformer(2) = x;
                
                % rewire the line right after it
                [~,id] = findElemByName(c,'line','090953_090954');
                c.line(id).bus1 = '090953a.1.2.3';
                
                % create a regcontrol for the new trans
                x = dssregcontrol;
                x.Name = 'regcon_2';
                x = resetRegcontrol(x);
                x.transformer = c.transformer(2).Name;
                x.vreg = 1.0*120;
                c.regcontrol(2) = x;
                
                % reset voltage regulator settings
                for i = 1:length(c.transformer)
                    c.transformer(i) = resetTransformer(c.transformer(i)); 
                end
                for i = 1:length(c.regcontrol)
                    c.regcontrol(i) = resetRegcontrol(c.regcontrol(i)); 
                end
            case 'wpv_existing'
                %% valleycenter 'wpv_existing' 
                [ c, opt, d, glc] = feederSetup( feederName, 'validated');
                CSI_data = excel2obj('data/GIS_CSI_v4.1.xlsx');
                c = addExistingPV(c,feederName,CSI_data);
                % use constant PQ load model
                for i = 1:length(c.load)
                    c.load(i).model = 1;
                end
                
            case 'wpv_existing_balanced3phase'
                %% valleycenter 'wpv_existing_balanced3phase'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                c = to3phases(c);
                
            case 'wpv_virtual' % with all virtual + real pv
                %% valleycenter 'wpv_virtual'
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_existing');
                pv = load([normalizePath('$KLEISSLLAB24-1')  '/database/gridIntegration/PVimpactPaper/virtualPV/f909scenario2_pv.mat']);
                c.pvsystem = pv.pv;
                c = convertPVModel(c);
                for i = 1:length(c.pvsystem)
                    [c.pvsystem(i).bus1, c.pvsystem(i).phases] = findNearestBus(c,c.pvsystem(i).bus1);
                end
                
            case {'wpv_virtual_balanced_3phase','s1b'}
                %% valleycenter {'wpv_virtual_balanced_3phase','s1b'}
                [ c, opt, d, glc] = feederSetup( feederName, 'wpv_virtual');
                c = to3phases(c);
                
                % change load model to constant P&Q
                for i = 1:length(c.load)
                    c.load(i).model = 1;
				end
				case {'reduced'}
               [ c, opt, d, glc] = feederSetup( feederName, 's1b');
			   critical_nodes=c.buslist.id(round(rand(8,1)*length(c.buslist)));
			   [c,w] = FeederReduction(critical_nodes,c);
				
				case 's1b_si'
                %% s1b with smart inverters!
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
%                 c = to3phases(c);

				%For smart inverters you need to define an xycurve that
				%controls the smart inverter, an xycurve that controls pv
				%efficiency and the xycurve that controls pv P vs. T.
				%Beyond that you need to define the InvControl and the
				%loadshape for PV. have fun.
                
				% Define Inverter control curve
				c.xycurve(1)=dssxycurve;
				c.xycurve(1).name='vv_curve';
				c.xycurve(1).npts=6;
				c.xycurve(1).Yarray='(1, 1, 0, 0, -1, -1)';
				c.xycurve(1).Xarray='(.5, 0.95, .98, 1.02, 1.05, 1.5)';
				
				%Add inverter class to system
				c.InvControl=dssInvControl; 
				c.InvControl.VoltageChangeTolerance=0.1;
				c.InvControl.VarChangeTolerance=0.1;
				c.InvControl.avgwindowlen='30s';
				c.InvControl.DeltaQ_factor=.1;
				c.InvControl.EventLog='No';
					
				%Must define both efficiency and PT curve of pvsystems
				%define xycurves 
				%efficiency
				c.xycurve(2)=dssxycurve;
				c.xycurve(2).name='Myeff';
				c.xycurve(2).npts=4;
				c.xycurve(2).Xarray=[.1 .2 .4 1.0];
				c.xycurve(2).yarray=[1 1 1 1];
				
				%P-T
				c.xycurve(3)=dssxycurve;
				c.xycurve(3).name='MyPvsT';
				c.xycurve(3).npts=4;
				c.xycurve(3).Xarray=[0 25 75 100];
				c.xycurve(3).yarray=[1.2 1.0 0.8 0.6];
				
				%Define Loadshape for PV %%% This is for validation
				%purposes only %%%
				c.loadshape=dssloadshape;
				c.loadshape.Name='ValidationLoadShape';
				c.loadshape.npts=1;
				c.loadshape.interval=1;
				c.loadshape.mult=1;
				
						
				
				%loop over pv and add vaiables
				for i = 1:length(c.pvsystem)
				c.pvsystem(i).EFFCURVE='Myeff';
				c.pvsystem(i).PTCurve='MyPvsT';
				c.pvsystem(i).daily='ValidationLoadShape';
				end
				
				case 's1b_si_50'
					[ c, opt, d, glc] = feederSetup( feederName, 's1b_si');
                load('c:\users\zactus\gridIntegration\SmartInverterSims\ValleyCenterPvBuses.mat')
				
				for ii=1:length(pvsystPast)
				c.InvControl(ii)=dssInvControl; 
				c.InvControl(ii).Name=['Inv' num2str(ii)];
				c.InvControl(ii).VoltageChangeTolerance=0.1;
				c.InvControl(ii).VarChangeTolerance=0.1;
				c.InvControl(ii).avgwindowlen='30s';
				c.InvControl(ii).DeltaQ_factor=.1;
				c.InvControl(ii).PVSystemList=c.pvsystem(pvsystPast(ii)).name;
				end
                
            case {'s1b_nocap'}
                %% valleycenter 's1b_nocap'
                [ c, opt, d, glc] = feederSetup( feederName, 's1b');
                if isfield(c,'capacitor'), c = rmfield(c,'capacitor'); end
                if isfield(c,'capcontrol'), c = rmfield(c,'capcontrol'); end
                
            case {'aggregated','agg'}
                %% ValleyCenter aggregated/ single PV profile case without capacitors
                [ c, opt, d, glc] = feederSetup( feederName, 's1b_nocap');
                
                % aggregated profile
                opt.fcProfileId = 1;
                
            otherwise
                error('Not supported config!');
        end
    otherwise
        error('Not supported feeder!');
end

if validatePQ
    validatepower(feederName,'kvar',0,'',[1 0 0 0 0],c,glc);
    validatepower(feederName,'kw',0,'',[1 0 0 0 0],c,glc);
    validateVoltage(c,feederName);
%     plotFeederProfile(c);
end
% add parameter before saving
opt.excludeScalingUpIds = excludeScalingUpIds;
% save output
save([outDir '/' fn],'c','d','glc','opt');
disp(['Saved file: ' outDir '/' fn]);
end