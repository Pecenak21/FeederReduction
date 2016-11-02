function fc = makeWholeDayProfile(t,fc,dt,smoothFt,type)
global indent;
% convert daily time to indices, starting from 12 am (mid-night)
tid = fc.time;
tid = tid - datenum(datestr(tid(1),'yyyy-mm-dd')); % remove the days
tid = 1+round(tid *24*60*60/dt); % convert to step by every 30 second

% assign profile
prof = zeros(length(t),size(fc.profile,2));
prof(tid,:) = fc.profile;

% fill in the empty holes (nan) values in the middle (missing points) with previous available value
if strcmpi(type,'pvsystem')
    mid = setdiff(tid(1):tid(end),tid); % limit the filling in sun hours (during the day) only for pvsystem
else
    mid = setdiff(1:length(t),tid); % filling in every missing holes
end
fprintf(['%sThere are ' num2str(length(mid)) ' time(s) that forecast data is missing. Fill them with previous available values.\n'],indent);
for j = 1:length(mid)
    if mid(j) == 1, continue; end;
    prof(mid(j),:) = prof(mid(j)-1,:);
end

if smoothFt
    for j = 1:size(prof,2)
        prof(:,j) = smooth(prof(:,j),smoothFt);
    end
end
fc.profile = prof;
fc.time = t;
end