function rr = ramprates( data, time, timestart, timeend, interval )
% timeseries data in 'data' and 'time' (as datenum)
% time bounds in 'timestart' and 'timeend' (as datenums)
% ramp time in 'interval', measured in minutes

%% average the data with a running average
int_days = mode(diff(time));
[d t n] = interpavg(data, time, interval, timestart:int_days:timeend);
s = round(interval/60/24/int_days); % number of points to pass to diff when calculating ramp rates
%% or with a block average
%int_days = interval/60/24;
%[d t n] = interpavg(data, time, interval, timestart-int_days/2:int_days:timeend+int_days/2);
%s = 1;

%% might want to filter out low n times (set to nan or zero?, not empty, otherwise the diff doesn't work properly)

%% calculate the ramp rates

r = d((1+s):end,:)-d(1:end-s,:);
% adjust the time series to match
t = t + int_days*s/2;
t(end-s+1:end) = [];

%% save the data
rr.time = t;
rr.ramp = r;

%% calculate the statistics
xmax = max(abs(r));
[h binx] = hist(r,-xmax:xmax/100:xmax);
% normalize by the number of non-nan values:
h = h/sum(h);
% save the probabilities and the ramp sizes that go with them
rr.pdf = h;
rr.df_x = binx;
% and save a cumulative distribution as well
rr.cdf = cumsum(h); % not sure if this is the desired order, but if not, 1-this will give us the desired order, so...


end
