%% Info: This script is used to run daily/ monthly/ yearly simulation using OpenDSS as the simulation engine
% and forecast and load data as input. The script will run various configurations (e.g. feeders to use, PV penetration levels, etc.)
% which are specified in the local/sim.conf file.

%% initialization
clear; close all;
% add all the files in the repo to running path, assuming you are in the gridIntegration re
def_addpath(1); % force re-add paths  
rmpath([pwd '/tmp']); % remove tmp folder from path
delete('tmp/*.dss');
delete('tmp/*.csv');

% configuration for a simulation is global so keep it global to avoid rerun the getConf function
global conf; conf = getConf();

% configuration for sub-simulation is global in each sub-function so we can keep these global as well
global fName; global fdSetup; global tDay; global pen; global indent; indent = '';

global kk xx stopping
%Overwrite?
ExportToNAS=1; %0=no,1=yes
overwrite=1; %0=no,1=yes
%% simulation main loop
% try

% loop through all the feeders needed to run simulations on
for i = 1:length(conf.feederName)
	fName = conf.feederName{i}; % feederName
	
	%update conf for USI and other things
	[~,~,Config]=xlsread('conf/Configurations.xlsx');
	
	[m,n]=size(Config);
	for ii=2:m
		if strcmp(fName,Config(ii,1))
			for j=3:n
				conf.(char(Config(1,j)))=cell2mat(Config(ii,j));
			end
			continue
		end
	end
	
	conf.fcOutDir = normalizePath(conf.fcOutDir);
	conf.timeDay=cellstr(datestr(datenum(conf.timeStart):1:datenum(conf.timeEnd)));
	%adjust conf from config file
	indent = ''; % start with no indentification
	deploySite = conf.deployment{i}; % get deployment site
	fprintf('%sFeeder name: %s\n',indent,fName);
	
	% loop through different feeder setups or configurations
	for fdid = 1:length(conf.feederSetup)
		fdSetup = conf.feederSetup{fdid}; indent = '   ';
		fprintf('%sFeeder setup: %s\n',indent,fdSetup);
		
		% load feeder's specific config
		fdOpt.fcProfileId=[];
		fdOpt.loadProfId=[];
		fdOpt.excludeScalingUpIds=[];
		[c0, fdOpt] = feederSetup(fName, fdSetup,overwrite);
		% replace the options in conf with customzed ifeeder options if specified
		ffn = fieldnames(fdOpt); cfn = fieldnames(conf); [fdOptId,confId] = ismember(lower(ffn),lower(cfn)); fdOptId = find(fdOptId); confId = confId(fdOptId);
		for k = 1:length(fdOptId), if ~isempty(fdOpt.(ffn{fdOptId(k)})), conf.(cfn{confId(k)}) = fdOpt.(ffn{fdOptId(k)}); end ;end
		
		% loop through days of simulation
		for tId = 1:length(conf.timeDay)
			tDay = conf.timeDay{tId}; indent = '      ';
			fprintf('%sDay: %s\n',indent,tDay);
			
			% all time steps for simulation
			dt = conf.timeStep; % in seconds
			t = datenum(tDay) : dt/24/3600 : (datenum(tDay) + 1 - dt/24/3600); % starting from midnight and end before midnight next day
			
			% load forecast GI profiles (only load the wanted profiles)
			fc = loadForecast( datestr(tDay,'yyyymmdd'), fName, conf.fcProfileId);
			if isempty(fc)
				fprintf('\n\nDay does not exist in forecast!\n\n')
				continue
			else
				% fill in the holes in data with appropriate methods when needed. Look into the fillForecastProfile for more details.
				[fc,emptyProfId] = fillForecastProfile(fc,deploySite,fName,tDay);
				
				% apply forecast scaling factor (to account for seasonal change in solar irradiation) if desired
				fc.profile = fc.profile * conf.fcScaling;
% 				[c0, fdOpt] = feederSetup(fName, fdSetup,overwrite);

				%Determine if it is reduced setup or not
				if isfield(c0,'weights')
					SizeMat=repmat(c0.PvSize',1,length(fc.time))';
					fc_tmp=(SizeMat.*fc.profile)*c0.PvbusMap*c0.weights';
					%Remove zeros so that fc_tmp is time x Pv (not time x nodes)
					fc_tmp(:,find(sum(fc_tmp)==0))=[];
					fc.profile=[];
					SystemSize=SizeMat*c0.PvbusMap*c0.weights';
					SystemSize=SystemSize(1,:);
					SystemSize(find(SystemSize==0))=[];
					for ii=1:length(c0.pvsystem)
						c0.pvsystem(ii).kVA=real(SystemSize(ii));
						c0.pvsystem(ii).pmpp=real(SystemSize(ii));
						c0.pvsystem(ii).kVAr=imag(SystemSize(ii));
						fc.profile(:,ii)=real(fc_tmp(:,ii)./SystemSize(ii));
					end
					% 					fc.profile=real(fc_tmp)./max(max(real(fc_tmp)));
				end
				
				% assign forecast profiles to pv systems
				[c_fc, fc2] = assignProfile(c0, 'pvsystem', fc, conf.fcProfileId, conf.fcTimeZone, t, conf.fcSmoothFactor );
				
				% get load profile data
				[loadProf, rawProf] = getLoadProfile(fName,conf.loadProfileId,tId);
				
				% assign load profiles to loads
				[c_fc, lp2] = assignProfile(c_fc, 'load', loadProf, conf.loadProfileId, conf.loadProfTZone, t, conf.loadProfSmooth);
				
				
				% reset PV levels if needed
				if ~conf.applyPVLevel, conf.PVLevel = getPenLevel(c_fc); end
				
				% loop through the pv penetration levels
				for penId = 1:length(conf.PVLevel)
					pen = conf.PVLevel(penId); indent = '         ';
					fprintf('%sSetup PV Penetration: %d %%',indent,pen);
					
					% check to see if this simulation has been done, if so, reload the saved result and skip this simulation
					fn = sprintf('%s/Res_%s_%s_%s_Pen%03.0f.mat',conf.outputDir,fName,fdSetup,datestr(datenum(tDay),'yyyymmdd'),pen);
					%                     if exist(fn,'file'),
					%                         fprintf(['\n%sSimulation result file exists. Load result and skip this simulation. Move result file if you want to rerun. File: ' fn '\n\n'],indent);
					%                         res = load(fn);
					%                         continue;
					%                     end
					
					% apply penetration level
					if conf.applyPVLevel
						c_pen = applyPenLevel(c_fc, pen, fdOpt.excludeScalingUpIds);
					else
						c_pen = c_fc;
					end
					
					% calculate actual penetration level after
					penLevActual = sum([c_pen.pvsystem.pmpp])/sum([c_pen.load.kw])*100;
					fprintf('... Actual PV Penetration: %3.4f %%\n',penLevActual);
					
					% run simulation using OpenDSS engine
					res = dssSimulation(c_pen, conf.mode, t, tDay, fdSetup);
					
					if stopping==1
						break
					end
					
					% add info to result struct
					res.time = t; res.conf = conf; res.feederName = fName; res.feederSetup = fdSetup; res.feederOption = fdOpt; res.timeDay = tDay; res.penLevel = pen; res.penLevActual = penLevActual;
					
% 					if isfield(c0,'weights')
% 					res.RedTime=c_pen.SimTime;
% 					end
					
					% save result
					save([fn],'-struct','res');
					fprintf(['%sSimulation result saved: ' fn '\n\n',indent]);
					
					% plot sim results
					if conf.doPlot && res.numNotConverged == 0
						plotSimResult(res); close all;
						%plotVoltTSeries(res);
					end
					
					if ExportToNAS==1
						mkdir(['C:\Users\Zactus\gridIntegration\NewRes\'])
						movefile(fn,['C:\Users\Zactus\gridIntegration\NewRes\sim_' num2str(xx) '_' num2str(kk) '.mat']);
					end
				end
			end
		end
	end
end
%  catch err
% %     % notify someone when simulation crashes
% %     ms = sprintf('%s\n%s\n%s\n%s\n%s',err.getReport,'--------------------------------','My software never has bugs. It just develops random features.', 'The definition of an upgrade: Take old bugs out, put new ones in.','The only difference between a bug and a feature is the documentation.');
% %     if conf.sendemail
% %         mailsend(conf.emaillist, sprintf('%s: Simulation crashed!',getenv('USERNAME')),ms);
% %     end
% %     % reset current dir
% %     cd(conf.gridIntDir);
% %     % rethrow the error to help with debugging
% %     rethrow(err);
% end
