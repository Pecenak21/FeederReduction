function [a samptimes numpts] = interpavg(data, timedata, avgtime, samptimes)
% INTERPAVG(data, timedata, avgtime, samptimes)
% averages timeseries data specified by 'data' and 'time' using 'avgtime' minute averages centered at times 'samptimes'

%% clean up inputs a bit
if( ~exist('samptimes','var') )
	samptimes = timedata;
end
samptimes = samptimes(:);

avgtime = avgtime/60/24/2; % convert from minutes to datenum format, and divide by 2 to get the amount that we add/subtract to do the centering

% want to have the time data be pre-sorted
if( ~issorted(timedata) )
	[timedata idx] = sort(timedata);
	data = data(idx);
end

%% allocate some data
a = zeros(size(samptimes,1),size(data,2)); %zero/zero values will fill in the NaNs when we divide through for regions with no data
numpts = a;

%% compute time bounds
%  this is just a more efficient way to implement (timedata>=(samptimes(i) - avgtime)) & (timedata < (samptimes(i) + avgtime))
ranges = interp1([0 timedata(:)' max(max(timedata(:)),max(samptimes(:)))+2*avgtime], [0 1:length(timedata) length(timedata)+1], [samptimes(:)-avgtime, samptimes(:)+avgtime]);
i_l = ceil(ranges(:,1));
i_h = ceil(ranges(:,2))-1;

%% compute averages
for i=1:length(samptimes)
	m = data(i_l(i):i_h(i),:);
	%m = data((timedata>=(samptimes(i) - avgtime)) & (timedata < (samptimes(i) + avgtime)));
	mm = isnan(m);
	m(mm) = 0;
	numpts(i,:) = sum(~mm);
	a(i) = sum(m)./numpts(i,:);
end

end
