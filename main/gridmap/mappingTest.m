%% call mapping tool
powerGridMapPath = 'C:\Users\andu\code\tmp\forecast\ucsd microgrid p11.png';
h = mappingTool(powerGridMapPath);
%% predefined buses with lat-lon location
d = get(h.buses,'data');
d{1,2} = 32.878865; d{1,3} = -117.237342; % library woak. chancellor complex
d{2,2} = 32.877252; d{2,3} = -117.233276; % gilman structure
set(h.buses,'data',d);