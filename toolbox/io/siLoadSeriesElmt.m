function [fSeries] = siLoadSeriesElmt(outputDir,datatype,fields)
% load the summary outputs from the most recent forecast run
% output file is indexed from 1 to 999 based on the order of the run
% 
% input:
%			datatype:	(optional) default: 'forecast'. 
%						supported datatype: 'forecast','power'
%
% Note that siLoadSeries tries to load the most recent of each of the series individually, which may not be exactly what you were after

% process input
if ~exist('outputDir','var') || isempty(outputDir)
	conf = readConf(siGetConfPath('forecast.conf'));
	outputDir = conf.outputDir;
end

if ~exist('datatype','var') || isempty(datatype)
	datatype = 'forecast';
end
 
if ~exist('fields','var') || isempty(fields)
	fields = '*'; % load all
end

% Check for previous runs to generate new index
flist = dir( [outputDir '/' datatype '*.mat'] );
fid = max(cellfun(@(s)str2double(s( length(datatype)+2 : length(datatype)+4 )),{flist.name}));

fn = [outputDir '/' datatype '_' sprintf('%03d',fid) '.mat'];
if(exist(fn,'file'))
	if iscell(fields)
		fSeries = load(fn,fields{:});
	else
		fSeries = load(fn,fields);
	end
else
	fSeries = [];
end
	
end