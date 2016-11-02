function plotCeilometerDaily(time,option)
% time: date of interest
% option: 1 for all, 2 for process ceil and METAR and raw heights only (no
% layers)
% NOTE: only works for ceilometer at EBU2 now

%% days of interest
if nargin > 0
	doi = floor(time);
else
	doi = floor(now);
end

%% path to save
p = siNormalizePath('$KleisslLab18-1/infobase/ceilometer_daily_plots');

%% get ceil
c = siCeilometer('EBU2');

%% generate plots
for i = doi
	% UTC time
	t = (i+8/24) : 30/24/3600 : (i+1+8/24);
	% get ceil data
	[hceil, raw] = siGetCloudHeightCeilometer(t);

	%% metar height
	hmet = metGetCloudHeights(c.postion,t);
	
	% save data and plots
	save([p '/data/' datestr(i,'yyyymmdd') '.mat'],'t','hceil','raw','hmet','c');
			
	if ~isempty(raw)
		if ~exist('option','var') || isempty(option) || option == 1
			%% layer height + raw height
			set(0,'defaultaxesfontsize',14);
			figure(10); set(gcf,'position',[150 50 1400 900]);
			plot(raw.time,raw.layerheight,'-','linewidth',2)
			hold on
			plot(t, hceil,'-m','linewidth',2)
			plot(t, hmet,'-g','linewidth',2)
			plot(raw.time,raw.hraw,'.k','linewidth',1), datetick; 
			ylim([0 5000]); xlim([min(t) max(t)]);xlabel('Time (HH:MM)'); ylabel(['Cloud Altitude (Offset: ' sprintf('%.0f',c.altitude) ') (m)']);
			legend('Ceil: Layer 1','Ceil: Layer 2','Ceil: Layer 3','Ceil: Layer 4','Processed Ceil Height','METAR Height','Raw Ceil Heights')
			title(datestr(i));
			hold off

			%%
			fn = [p '/' datestr(i,'yyyymmdd') ''];
			saveas(gcf,[fn '.fig']); saveas(gcf,[fn '.png']);
		elseif option == 2
			set(0,'defaultaxesfontsize',14);
			figure(10); set(gcf,'position',[150 50 1400 900]);
			hold on
			plot(t, hceil,'-m','linewidth',2)
			plot(t, hmet,'-g','linewidth',2)
			plot(raw.time,raw.hraw,'.k','linewidth',1), datetick; 
			ylim([0 5000]); xlim([min(t) max(t)]);xlabel('Time (HH:MM)'); ylabel(['Cloud Altitude (Offset: ' sprintf('%.0f',c.altitude) ') (m)']);
			legend('Processed Ceil Height','METAR Height','Raw Ceil Heights')
			title(datestr(i));
			hold off
		end
	end
end

end