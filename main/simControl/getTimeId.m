function [ tid ] = getTimeId( t, tArray )
[~,tid] = ismember(round(t*3600*24),round(tArray*3600*24));
end

