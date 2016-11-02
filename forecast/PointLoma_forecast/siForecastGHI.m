%% Solar Resource Assessment
%  Sky Imager Library
%
%  Title: Sky Imager Forecasting, GHI computation
%
%  Author: Bryan Urquhart
%
%
%  Description:
%    This function performs a single advection forecast and shadow mapping. The
%    shadow results are then converted to a GHI estimate
%
%    The processing done by this function has a large potential for improvement.
%    More sophisticated processing techniques can be applied.
%
%
function f_elmt = siForecastGHI( f_elmt , cloudmap , cmPower , ~ , csk , target )
%% Process Input Arguments

% Initialize output
f_elmt.ktavg     = nan ( length( f_elmt.shadow ) , 1 );
f_elmt.isGndShdw = nan ( length( f_elmt.shadow ) , 1 );
f_elmt.shdwFctn  = nan ( length( f_elmt.shadow ) , 2 );
f_elmt.gi		 = nan ( length( f_elmt.shadow ) , 1 );
f_elmt.power	 = nan ( length( f_elmt.shadow ) , 1 );
f_elmt.kt	     = nan ( 3 , 1 );
f_elmt.ktpdf	 = [];
f_elmt.station	 = [];

% Set the minimum fraction of plant coverage considered acceptable
doSensors = isfield(target,'data_type') && strcmpi(target.data_type(1), 'g');
doCombo = isfield(target,'data_type') && strcmpi(target.data_type, 'ghipv');
if doSensors
	minShadowFraction = 0.005;
else
	minShadowFraction = 0.05;
end

% Set a sun-pixel angle limit
limitSPA = 18; %[deg]

% % Integer value corresponding to maximum cloudiness in the map
% thickCld = 2;

%% Generate kt histogram

% We may need to include zenith angle filtering if we change how we process
% images and start including zenith angles higher than 86.

% Extract the kt
if doCombo % Ground sensors + PV
	pdfloops = 3; % Generate ktpdfs for only sensor data, only inverter data, and both combined.
	tmp{1} = vertcat(cmPower.sensor(:));
	tmp{2} = vertcat(cmPower.inverter(:));
	for catidx = 1:size(tmp{1},1)
		tmp{3}(catidx,1).kt = horzcat(tmp{1}(catidx).kt, tmp{2}(catidx).kt);
	end
elseif doSensors % Ground sensors only
	pdfloops = 1;
	% Ground measurement instrument (i.e. from, irradiance sensor). cmPower will have 'sensor' field.
	tmp{1} = vertcat(cmPower.sensor(:));
else % PV only
	pdfloops = 1;
	% Power plant data (i.e. from PV inverters). cmPower will have 'inverter' field.
	tmp{1} = vertcat(cmPower.inverter(:));
end

for pdfidx = 1:pdfloops
	[f_elmt.kt(:,pdfidx), f_elmt.ktpdf{1,pdfidx}] = ktPeaks(tmp{pdfidx}, cloudmap.fraction);
end

%% Determine the average plant kt

% Compute the shadow fraction

if( isfield(target,'data_type') )
	if ( doCombo || doSensors )
		mapname = 'GHI';
	else
		mapname = target.data_type;
	end
	footprint = ~isnan(target.footprint.(mapname));
elseif strcmpi(target.name,'redlands')
	% inverter
	footprint = ~isnan(target.footprint.inverter);
	% combine GHI footprint with inverter footprint
	footprint = footprint || ~isnan(target.footprint.GHI + max(target.footprint.inverter(:)) );
else
	footprint = ~isnan(target.footprint.inverter);
end
npix = sum( footprint(:) );
% for smaller areas, it's faster to use numerical indexing, rather than logical indexing
if(npix < 2e5); footprint = find(footprint); end

if doSensors
	% assume footprint.GHI has stations labeled from 1 to N where N is
	% total number of stations/ PV systems
	station = struct();
	stationId = 1:length(target.footprint.GHInames); % Station IDs of ground sensors
	for ID = stationId
		station(ID).footprint = (target.footprint.(mapname) == ID);
		station(ID).npix = sum(station(ID).footprint(:));
		% for smaller areas, it's faster to use numerical indexing, rather than logical indexing
		if(station(ID).npix < 2e5); station(ID).footprint = find(station(ID).footprint); end
	end
else
	stationId = []; % no stations
end

if doCombo || strcmpi(target.name,'pointloma')
	pv = struct();
	pvId = 1:length(target.footprint.inverterNames);
	for ID = pvId
		pv(ID).footprint = (target.footprint.(mapname) == ( ID + 100 )); % Convention: PV IDs are numbered starting from 101
		pv(ID).npix = sum(pv(ID).footprint(:));
		% for smaller areas, it's faster to use numerical indexing, rather than logical indexing
		if(pv(ID).npix < 2e5); pv(ID).footprint = find(pv(ID).footprint); end
	end
	
	% Stupid hack because of characteristic kt array convention
	if size(f_elmt.kt,2) < 2
		f_elmt.kt(:,2) = f_elmt.kt(:,1);
	end
	
else
	pvId = [];
end

% Compute the fraction of cloud near the sun
mask = (cloudmap.sun.pixel_angle_map <= limitSPA);

% allocate averaging fields for stations and PV arrays
% note that if either stationId or pvId is empty, we won't do anything here, so we don't need to check doSensors/doCombo again
for ID = stationId
	f_elmt.station(ID).ktavg     = nan ( length( f_elmt.shadow ) , 1 );
	f_elmt.station(ID).isGndShdw = true( length( f_elmt.shadow ) , 1 );
	f_elmt.station(ID).shdwFctn  = nan ( length( f_elmt.shadow ) , 2 );
end

for ID = pvId
	f_elmt.pv(ID).ktavg     = nan ( length( f_elmt.shadow ) , 1 );
	f_elmt.pv(ID).isGndShdw = true( length( f_elmt.shadow ) , 1 );
	f_elmt.pv(ID).shdwFctn  = nan ( length( f_elmt.shadow ) , 2 );
end

% Compute kt avg
for sIdx = 1:length( f_elmt.shadow )
	% calculate shadow fractions for individual sensors/pv arrays
	for ID = stationId
		[f_elmt.station(ID).shdwFctn(sIdx,:), f_elmt.station(ID).ktavg(sIdx), f_elmt.station(ID).isGndShdw(sIdx)] = ...
			footprintShadow(f_elmt.shadow{sIdx}, station(ID).footprint, f_elmt.kt(:,1), station(ID).npix, minShadowFraction, mask, f_elmt.advect{sIdx});
	end
	for ID = pvId
		[f_elmt.pv(ID).shdwFctn(sIdx,:), f_elmt.pv(ID).ktavg(sIdx), f_elmt.pv(ID).isGndShdw(sIdx)] = ...
			footprintShadow(f_elmt.shadow{sIdx}, pv(ID).footprint, f_elmt.kt(:,2), pv(ID).npix, minShadowFraction, mask, f_elmt.advect{sIdx});
	end
	% and for the whole plant
	[f_elmt.shdwFctn(sIdx,:), f_elmt.ktavg(sIdx), f_elmt.isGndShdw(sIdx)] = ...
		footprintShadow(f_elmt.shadow{sIdx}, footprint, f_elmt.kt(:,size(f_elmt.kt,2)), npix, minShadowFraction, mask, f_elmt.advect{sIdx});
	
end

%% Determine power output

% Get the time indices
sIdx = zeros(length(f_elmt.ktavg),1);
for idx = 1:length( f_elmt.ktavg )
	sIdx(idx) = find( round(csk.time.*24.*3600) == round(f_elmt.time(idx).*24.*3600), 1); % round(x.*24.*3600) due to precision errors
end
% Compute the gi averages for ground sensors
for ID = stationId
	f_elmt.station(ID).gi = f_elmt.station(ID).ktavg .* csk.gi(sIdx);
end
% PV arrays
for ID = pvId
	csk.pv.gi = ghi2gi( target.footprint.pv(ID).azimuth , target.footprint.pv(ID).tilt , target.footprint.pv(ID).pos , csk.time , csk.ghi );
	clrPwrModel = csk.pv.gi.*target.footprint.pv(ID).scaleFactor;
	clrPwrModel = polyval(target.footprint.pv(ID).eff_coeff, clrPwrModel).*clrPwrModel;
	f_elmt.pv(ID).power = f_elmt.pv(ID).ktavg .* clrPwrModel;
	f_elmt.pv(ID).ghi = f_elmt.pv(ID).ktavg .* csk.ghi;
	f_elmt.pv(ID).gi = f_elmt.pv(ID).ktavg .* csk.pv.gi;
end

% Plant
f_elmt.gi = f_elmt.ktavg .* csk.gi(sIdx);

% Compute the plant power output
f_elmt.power = pwrGIToPower( f_elmt.gi , sum( target.design.power(:) ) );

%% Machine learning output
if isfield(f_elmt, 'brightness')
	% User message
	fprintf( '\n\tGenerating machine learning outputs..........\n\t\t\t\t\t\t\t');
	
	for CircleIdx = 1:length(target.footprint.GHInames)
		% Load individual circular footprints (100m radius, separated so
		% stations do not overlap)
		circlefootprint = load([siNormalizePath('$KLEISSLLAB18-1/database/deployments/UCSD/CircularFootprints/') target.footprint.GHInames{CircleIdx} 'footprint.mat']);
		% Initialize variables
		f_elmt.ravg{CircleIdx} = nan ( length( f_elmt.r ) , 1 );
		f_elmt.gavg{CircleIdx} = nan ( length( f_elmt.g ) , 1 );
		f_elmt.bavg{CircleIdx} = nan ( length( f_elmt.b ) , 1 );
		f_elmt.rbrcslavg{CircleIdx} = nan ( length( f_elmt.rbrcsl ) , 1 );
		f_elmt.brightnessavg{CircleIdx} = nan ( length( f_elmt.brightness ) , 1 );
		f_elmt.MLshdwFctn{CircleIdx} = nan ( length( f_elmt.shadow ) , 2 );
		
		% Average metrics within each footprint
		for AvgIdx = 1:length(f_elmt.ravg{CircleIdx})
			% Convert footprint to binary double
			cfp = double(circlefootprint.GHI > 0);
			f_elmt.rbrcslavg{CircleIdx}(AvgIdx) = nansum(cfp(:).*f_elmt.rbrcsl{AvgIdx}(:)) / sum(cfp(:));
			
			% Convert footprint to binary int16
			cfp = int16(circlefootprint.GHI > 0);
			f_elmt.ravg{CircleIdx}(AvgIdx) = nansum(cfp(:).*f_elmt.r{AvgIdx}(:)) / sum(cfp(:));
			f_elmt.gavg{CircleIdx}(AvgIdx) = nansum(cfp(:).*f_elmt.g{AvgIdx}(:)) / sum(cfp(:));
			f_elmt.bavg{CircleIdx}(AvgIdx) = nansum(cfp(:).*f_elmt.b{AvgIdx}(:)) / sum(cfp(:));
			f_elmt.brightnessavg{CircleIdx}(AvgIdx) = nansum(cfp(:).*f_elmt.brightness{AvgIdx}(:)) / sum(cfp(:));
			
			[f_elmt.MLshdwFctn{CircleIdx}(AvgIdx,:), ~, ~] = footprintShadow(f_elmt.shadow{AvgIdx}, (circlefootprint.GHI == CircleIdx), f_elmt.kt, sum(cfp(:)), minShadowFraction, mask, f_elmt.advect{AvgIdx});
			
		end
		
	end
end

end

function [shdwFctn, ktavg, isGndShdw] = footprintShadow(shadow, footprint, kt, npix, minShadowFraction, mask, advect)
% Mask shadow with footprint
shadow = double(shadow(footprint));

% Compute the shadow fraction
shadow = shadow(shadow >= 0); %remove masked pixels

% Require a minimum level of coverage before we use the ground data
nvalid = length(shadow);
if( nvalid / npix >= minShadowFraction )
	isGndShdw = true;
	shdwFctn = [ sum( shadow == 1 ) sum( shadow == 2 ) ] / numel(shadow);
elseif nvalid==0
	% this case is mainly here to mirror a case from before when we don't have any overlap at all
	isGndShdw = true;
	shdwFctn = nan(1,2);
else
	isGndShdw = false;
	% Compute shadow fraction from SPA masked advect map
	advect = mask .* double(advect+1) - 1;
	advect = advect(advect >= 0); %remove masked pixels
	shdwFctn = [ sum( advect == 1 ) sum( advect == 2 ) ] / numel(advect);
end
% Infer clear fraction
clrfraction = 1 - sum( shdwFctn );

% Compute the average kt
% TODO: make sure this is somehow based on the correct advection if we're not using a ground shadow
ktavg = dot( kt , [ fliplr(shdwFctn) clrfraction ] );

end

function [characteristic_kt, ktpdf] = ktPeaks(powerStruct, cloudFraction)
%% ktPeaks function constants
characteristic_kt = nan(3,1);

% Set the class defaults. Based on 8 months of data.
ktthk = 0.418;
ktthn = 0.700;
ktclr = 1.057;
%ktcle = 1.277;

% Set the number of bins
nbins = 500;

% Set kt min and max values
ktmin = 0.05;
ktmax = 1.4;

% Power validation data limit
datalimit = 2; % In hours
%% Generate kt PDF
% Cap kt data at [datalimit] hours prior to current time step
[tmpsize, ~] = size(powerStruct);
if tmpsize > datalimit*2*60 % 30s image intervals or 2 datapoints for every minute
	powerStruct(1:(tmpsize-datalimit*120)) = []; % Delete first entry/entries if power validation data exceeds [datalimit] hours
end

kt = vertcat(powerStruct(:).kt)';

% Set the kt limits
kt( kt < ktmin ) = nan;
kt( kt > ktmax ) = ktmax;

% Return if kt is all NaN or if only one non-NaN value (minimum 2 points required for binning). Also set characteristic kt values to defaults.
if ( all(isnan(kt(:))) ) || ( sum(~isnan(kt(:))) < 2 ) || ( numel(unique(kt)) < 2 )
	warning('siForecastGHI:Noktvalue','\tNo kt value is reported. Defaults used (0.418, 0.700, 1.057).');
	ktpdf = []; characteristic_kt(1) = ktthk; characteristic_kt(2) = ktthn; characteristic_kt(3) = ktclr;
	return;
end;

% Construct the histogram
ktpdf = slrKtHist( kt , 'NumberOfBins' , nbins );

% Perform normalization (convert to pdf)
ktpdf.binwidth = ktpdf.bins(:,3) - ktpdf.bins(:,1);
ktpdf.integral = sum( ktpdf.binwidth .* ktpdf.binval );

% Filter PDF
ktpdf.nfilt = floor( 0.1/mean(ktpdf.binwidth) );
if(ktpdf.nfilt > nbins); ktpdf.nfilt = nbins; end
ktpdf.binval_filt = lmFilter1( ktpdf.binval , ktpdf.nfilt , 'type ' , 'gaussian' , 'sigma' , 0.1 );

% Find peaks
ktpdf = pwrKtPdfPeaks2( ktpdf );

%% Determine kt for each class

% Pull out the clear and cloud kt. The pwrKtPdfPeaks is hard coded for 3 cloud
% classes but should be written for a more general peak determination. This
% revision should be done in tandem with the capability to detect more than one
% cloud class.

for idx = 1:length(ktpdf.peaks.bins)-1 % Skip cloud enhancement bin
	ktval          = ktpdf.bins( ktpdf.peaks.xmax{idx} , 2 );
	ktpdfval       = ktpdf.binval_filt( ktpdf.peaks.xmax{idx} );
	characteristic_kt(idx) = sum(ktval .* ktpdfval) / sum( ktpdfval );
end

% Use default values if kt pdf produces unreasonable values
if ( characteristic_kt(1) > 0.6 || isnan(characteristic_kt(1)) ), characteristic_kt(1) = ktthk; end % thick cloud should not have kt > .6
if ( characteristic_kt(2) < 0.6 || isnan(characteristic_kt(2)) ), characteristic_kt(2) = ktthn; end % thin cloud should not have kt < .6
if ( characteristic_kt(3) < 0.9 || isnan(characteristic_kt(3)) ), characteristic_kt(3) = ktclr; end % clear sky should not have kt < .9

%% If sky conditions are obviously clear/overcast, use median kt from past minute
if ( cloudFraction < 0.05 ) && ( size(kt, 2) >= 60 ) && ( sum(isnan(kt(end-size(kt, 1)*60+1:end))) < 0.8*size(kt,1)*60 )% Obviously clear
	characteristic_kt(3) = nanmedian(kt(end-size(kt, 1)*60+1:end));
elseif ( cloudFraction > 0.95 ) && ( size(kt, 2) >= 60 ) && ( sum(isnan(kt(end-size(kt, 1)*60+1:end))) < 0.8*size(kt, 1)*60 )% Overcast
	characteristic_kt(1) = nanmedian(kt(end-size(kt, 1)*60+1:end));
end

end